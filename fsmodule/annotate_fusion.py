import pysam

from . import stat_split_reads

def is_overlapped(exon: tuple[str, int, int], candidate_exon: tuple[str, int, int, str]) -> bool:
    """
    check if two exons are overlapped
    Args:
        exon: tuple[str, int, int], chrom, start, end
        candidate_exon: tuple[str, int, int, str], chrom, start, end, symbol
    Returns:
        bool, True if overlapped, False otherwise
    """
    if exon[0] != candidate_exon[0]:
        return False
    return max(0, min(exon[2], candidate_exon[2]) - max(exon[1], candidate_exon[1])) > 0

def split_exon(exon: tuple[str, int, int], all_possible_exons: list[tuple[str, int, int, str]]) -> list[tuple[str, int, int, str]]:
    """
    split the exon into multiple exons
    Args:
        exon: exon to be split, tuple[str, int, int], chrom, start, end
        candidate_exon: all overlapped exons, tuple[str, int, int, str], chrom, start, end, gene name
    Returns:
        tuple[start: int, end: int, symbol: str], start and end of the overlap, symbol of the gene
    """
    split_result = []
    # sort all_possible_exons by start position
    all_possible_exons.sort(key=lambda x: x[1])
    cur = exon[1] + 1 # 0-based to 1-based
    for cand in all_possible_exons:
        s = max(exon[1], cand[1])
        e = min(exon[2], cand[2])
        if cur >= e:
            continue
        if cur < s and s - cur > 5: # hyperpatameter, set here to avoid alignment error
            split_result.append((exon[0], cur, s - 1, "X"))
        split_result.append((exon[0], s, e, cand[3]))
        cur = e + 1
    if cur < exon[2]:
        split_result.append((exon[0], cur, exon[2], "X"))
    return split_result


def annotate_exon_type(exon: tuple[str, int, int], all_candidate_exons: list[tuple[str, int, int, str]]) -> str:
    """
    determine the type of the exon
    Args:
        exon: tuple[str, int, int], chrom, start, end
        all_candidate_exons: list[tuple[str, int, int, str]], chrom, start, end, gene name
    Returns:
        type: str, A, B, X, AX, XB, AXB, ..., A stand for gene1, B stand for gene2, X stand for other sequences
    """
    all_possible_exons = [cand for cand in all_candidate_exons if is_overlapped(exon, cand)]
    if len(all_possible_exons) == 0:
        return 'X'
    split_result = split_exon(exon, all_possible_exons)
    type = ''.join([cand[3] for cand in split_result])
    return type

def stats_exon(reads: list[pysam.AlignedSegment], gene1: str, gene2: str) -> tuple[list[int], list[str], list[tuple[str, int, int]]]:
    """
    stats the exon information for each read and classify fusion type
    Args:
        reads: list[pysam.AlignedSegment]
        gene1: str, gene name
        gene2: str, gene name
    Returns:
        exon_num: list[int], num of gene A, B ,X
        exon_type: list[str], A, B, X, AX, XB, AXB, ...
        exon_coordinates: list[tuple[str, int, int]]
    """
    exon_num = []
    exon_type = []
    exon_coordinates = []
    # sort reads by alignment start position
    reads.sort(key=lambda x: x.reference_start)

    chrom_aligned_range = {} # chrom -> [min_start, max_end]
    for read in reads:
        tmp_exon_coordinates = []
        chrom = read.reference_name
        exon_lengths = stat_split_reads.calculate_exon_lengths(read.cigartuples)
        exon_lengths = [exonlen for exonlen in exon_lengths if exonlen > 0]
        # get exon coordinates
        alignpair = read.get_aligned_pairs() # coordinates: query -> reference
        alignpair = [c for c in alignpair if c[1] != None and c[0] != None]
        accumalignlen = 0
        for exon_len in exon_lengths:
            tmp_exon_coordinates.append((chrom, alignpair[accumalignlen][1], alignpair[accumalignlen + exon_len - 1][1]))
            accumalignlen += exon_len
        if chrom not in chrom_aligned_range:
            chrom_aligned_range[chrom] = [alignpair[0][1], alignpair[-1][1]]
        else:
            chrom_aligned_range[chrom][0] = min(chrom_aligned_range[chrom][0], alignpair[0][1])
            chrom_aligned_range[chrom][1] = max(chrom_aligned_range[chrom][1], alignpair[-1][1])
        if read.is_reverse:
            tmp_exon_coordinates = tmp_exon_coordinates[::-1]
        exon_coordinates.extend(tmp_exon_coordinates)

    # get all possible exons
    all_candidate_exons = []
    for chrom, (min_start, max_end) in chrom_aligned_range.items():
        tmp = stat_split_reads.overlap_with_genes_in_details(chrom, min_start, max_end) # [start, end, gene, strand, [[exon_start, exon_end],...]]
        for gene in tmp:
            if gene[2] not in [gene1, gene2]:
                continue
            symbol = 'A' if gene[2] == gene1 else 'B'
            all_candidate_exons.extend([(chrom, exon[0], exon[1], symbol) for exon in gene[4]])

    # annotate exon type
    exon_type = [annotate_exon_type(exon, all_candidate_exons) for exon in exon_coordinates]

    # calculate exon num from each gene
    exon_num = [len(exon_type)]
    for symbol in ['A', 'B', 'X']:
        tmp = [type for type in exon_type if symbol in type]
        exon_num.append(len(tmp))
    
    return exon_num, exon_type, exon_coordinates

def classify_fusion_type(exon_num: list[int], exon_type: list[str]) -> str:
    """
    classify the fusion type
    Args:
        exon_num: list[int], num of gene A, B ,X
        exon_type: list[str], A, B, X, AX, XB, AXB, ...
    Returns:
        fusionType: str, A, B, X, AX, XB, AXB, ...
    """
    if exon_num[0] == 1:
        type = "single exon"
    elif exon_num[1] == 1 and "A" in exon_type[0]:
        type = "5' fusion"
    elif exon_num[2] == 1 and "B" in exon_type[-1]:
        type = "3' fusion"
    elif exon_num[1] + exon_num[2] > 0.8 * exon_num[0]:
        type = "confident fusion"
    else:
        type = "fusion"

    return type

def annotate_fusion(outpath: str) -> None:
    # stats exon info for each read and classify fusion type
    bam = pysam.AlignmentFile(outpath + 'align_workspace/allalignedseq.bam')

    reads_dict = {} # readName -> list[read object]
    for read in bam.fetch():
        readName = read.query_name
        if readName not in reads_dict:
            reads_dict[readName] = []
        reads_dict[readName].append(read)
    
    # stats reads     
    lines = ""
    # get old breakpoint
    old_breakpoint = {} # ID -> list
    with open(outpath + 'confident_genefusion.txt', 'r') as fi:
        for line in fi:
            if line.startswith('gene1'):
                continue
            ll = line.strip('\n').split('\t') # gene1, gene2, numSupp, chrom1, bp1, chrom2, bp2, ID, readNames
            old_breakpoint[ll[7]] = ll
    with open(outpath + 'confident_genefusion_refined.txt', 'r') as fi:
        for line in fi:
            if line.startswith('ID'):
                lines += f"{line.strip('\n')}\texonNum\texonType\texonCoordinates\tfusionType\n"
                continue
            ll = line.strip('\n').split('\t') # ID, gene1, gene2, numSupp, chrom1, bp1, chrom2, bp2, readNames
            old_ll = old_breakpoint[ll[0]] # gene1, gene2, numSupp, chrom1, bp1, chrom2, bp2, ID, readNames
            readname = f"poa_ctg_{old_ll[7]}_{old_ll[0]}_{old_ll[1]}_{old_ll[3]}_{old_ll[4]}_{old_ll[5]}_{old_ll[6]}"
            if readname not in reads_dict:
                raise ValueError(f"Read {readname} not found in bam file")
            reads = reads_dict[readname]
            # stats exon info for each read
            exon_num, exon_type, exon_coordinates = stats_exon(reads, ll[1], ll[2])
            exon_num_str = '-'.join([str(num) for num in exon_num])
            exon_type_str = '-'.join(exon_type)
            exon_coordinates_str = ','.join([f"{chrom}:{start}-{end}" for chrom, start, end in exon_coordinates])
            # classify fusion type
            fusionType = classify_fusion_type(exon_num, exon_type)
            lines += f"{line.strip('\n')}\t{exon_num_str}\t{exon_type_str}\t{exon_coordinates_str}\t{fusionType}\n"
    
    # replace
    with open(outpath + 'confident_genefusion_refined.txt', 'w') as fo:
        fo.write(lines)