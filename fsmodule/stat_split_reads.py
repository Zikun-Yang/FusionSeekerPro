import pysam
import os

def overlap_with_genes(chrom, start, end, is_single_exon=False):
	"""
	Annotate a genomic segment with overlapping genes using binary search.
	
	Uses binary search to efficiently find genes that overlap with the given
	chromosomal segment. Returns genes that have significant overlap (>50bp or
	>50% of segment length).
	
	Args:
		chrom (str): Chromosome name
		start (int): Start position of the segment
		end (int): End position of the segment
		is_single_exon (bool): Whether the segment is a single exon
	Returns:
		list: List of gene names that overlap with the segment
	"""
	allgenes = geneinfo[chrom] # gene_start, gene_end, gene_name, gene_strand, [[exon_start, exon_end],...]
	pin = int(len(allgenes) / 2)
	pinstart = 0
	pinend = len(allgenes)
	overlapped_candidate_genes = []
	overlapped_genes = []
	notfound = True
	while notfound:
		gene = allgenes[pin]
		ovlplen = min(gene[1], end) - max(gene[0], start)

		if ovlplen >= 0:
			notfound=False
			break
		if start > gene[1]:
			newpin = int((pin + pinend) / 2)
			if newpin == pin:
				break
			else:
				pinstart = pin
				pin = newpin
				continue
		if end < gene[0]:
			newpin = int((pin + pinstart) / 2)
			if newpin == pin:
				break
			else:
				pinend = pin
				pin = newpin
				continue

	for gene in allgenes[max(0, pin - 10) : pin + 10]:
		ovlplen = min(gene[1], end) - max(gene[0], start)
		if ovlplen > 50 or ovlplen > 0.5 * (end - start):
			overlapped_candidate_genes.append(gene) # add gene name

	# remove genes that are included in other genes
	if is_single_exon and len(overlapped_candidate_genes) > 1:
		gene2remove = []
		for idx in range(len(overlapped_candidate_genes)):
			g1 = overlapped_candidate_genes[idx]
			for idx2 in range(idx + 1, len(overlapped_candidate_genes)):
				g2 = overlapped_candidate_genes[idx2]
				if g1[0] <= g2[0] and g1[1] >= g2[1]:
					gene2remove.append(idx2)
				if g2[0] <= g1[0] and g2[1] >= g1[1]:
					gene2remove.append(idx)
		overlapped_candidate_genes = [g for idx, g in enumerate(overlapped_candidate_genes) if idx not in gene2remove]
	
	# check exons
	for gene in overlapped_candidate_genes:
		for exon in gene[4]:
			if min(exon[1], end) - max(exon[0], start) > 0:
				overlapped_genes.append(gene[2])
				break

	return overlapped_genes

def calculate_exon_lengths(cigartuple: list[tuple[int, int]]) -> list[int]:
	"""
	Calculate exon segment lengths from CIGAR tuple.
	
	Args:
		cigartuple (list): List of CIGAR operation tuples (operation, length)
		
	Returns:
		list: List of exon segment lengths (aligned bases)
	"""
	exon_lengths = []
	lastlength = 0
	for pair in cigartuple:
		if pair[0] in [3,7,8]: # match, soft clip, hard clip
			exon_lengths.append(lastlength)
			lastlength =0
		if pair[0] == 0: # insertion
			lastlength += pair[1]
		if pair[0] in [4,5]:
			if lastlength != 0: # soft clip, hard clip
				exon_lengths.append(lastlength)
				lastlength = 0
		if pair[0] in [1,2]: # deletion, intron
			continue

	if lastlength != 0:
		exon_lengths.append(lastlength)
	return exon_lengths


def calculate_soft_clipping(cigartuple: list[tuple[int, int]]) -> list[int]:
	"""
	Extract left and right clipping information from CIGAR tuple.
	
	Identifies soft (4) or hard (5) clipping at the beginning and end
	of the alignment.
	
	Args:
		cigartuple (list): List of CIGAR operation tuples (operation, length)
		
	Returns:
		list: [left_clip_length, -1, right_clip_length], -1 indicates no information is available
	"""
	leftclip = 0
	rightclip = 0
	if cigartuple[0][0] in [4, 5]:
		leftclip = cigartuple[0][1]
	if cigartuple[-1][0] in [4, 5]:
		rightclip = cigartuple[-1][1]
	return [leftclip, -1, rightclip]

def simplify_cigar(cigartuple: list[tuple[int, int]]) -> str:
	"""
	Simplify CIGAR tuple to a string representation with matches and insertions.
	
	Converts CIGAR operations to a simplified string format, combining
	matches (M) and including insertions (I). Skips deletions (D) and
	introns (N).
	
	Args:
		cigartuple (list): List of CIGAR operation tuples (operation, length)
		
	Returns:
		str: Simplified CIGAR string (e.g., "100M,50I,200M")
	"""
	cigarlist = []
	nummatch = 0
	for c in cigartuple:
		if c[0] in [0,7,8]:
			nummatch += c[1]
		if c[0] in [3]:
			cigarlist += [str(nummatch) + 'M']
			nummatch = 0
		if c[0] in [2]:
			continue
		if c[0] in [1]:
			cigarlist += [str(nummatch) + 'M',str(c[1]) + 'I']
			nummatch = 0
		if c[0] in [4,5]:
			if nummatch != 0:
				cigarlist += [str(nummatch) + 'M']
				nummatch = 0
	if nummatch != 0:
		cigarlist += [str(nummatch) + 'M']
		nummatch = 0
	cigarstring = ','.join(cigarlist)
	return cigarstring

def stat_read_info(read, is_record_readseq: bool) -> list:
	"""
	Detect fusion signals within a single read alignment.
	
	Analyzes a read alignment to identify split reads (reads spanning
	multiple exons or genes). Extracts alignment information including
	exon segments, gene annotations, CIGAR information, and sequence data.
	
	Args:
		read: pysam AlignedSegment object representing the read alignment
		is_record_readseq (bool): Whether to record read sequence and quality scores
		
	Returns:
		list: List containing read details [chrom, start, end, read_name, strand,
			             cigar_info, mapq, left_gene, right_gene, exon_segments, sequence,
			             quality, nm_rate, simple_cigar, read_object]
	"""
	alignpair = read.get_aligned_pairs() # get pair of alignment coordinates base by base [(query_base, reference_base),...]
	alignpair = [c for c in alignpair if c[1] != None and c[0] != None]
	chrom = read.reference_name # get chromosome name
	has_sa = read.has_tag('SA')

	# get exon lengths
	cigartuple = read.cigartuples
	exon_lengths = calculate_exon_lengths(cigartuple) # get exon segments
	exon_lengths = [exonlen for exonlen in exon_lengths if exonlen > 0] # remove exon segments with length 0

	# get start and end position
	start = alignpair[0][1]
	end = alignpair[-1][1]

	readinfo = None
	if len(exon_lengths) == 1: # only one exon
		overlapped_genes = overlap_with_genes(chrom, start, end, is_single_exon=True)
		overlapped_segments = [[start, end, exon_lengths[0], overlapped_genes]]
		if len(overlapped_genes) == 0:
			return None
		if len(overlapped_genes) == 1 and not has_sa:
			return None
	else: # multiple exons
		aligned_segments = []
		accumalignlen = 0
		for exon_len in exon_lengths:
			aligned_segments.append([alignpair[accumalignlen][1], alignpair[accumalignlen + exon_len - 1][1], exon_len])
			accumalignlen += exon_len

		overlapped_segments = []
		specifically_aligned_genes = []
		for i in range(len(aligned_segments)):
			overlapped_genes = overlap_with_genes(chrom, aligned_segments[i][0], aligned_segments[i][1])
			overlapped_segments.append(aligned_segments[i] + [overlapped_genes])
			if len(overlapped_genes) == 1:
				specifically_aligned_genes.extend(overlapped_genes)

		num_specifically_aligned_genes = len(list(set(specifically_aligned_genes)))
		if num_specifically_aligned_genes == 0:
			return None
		if num_specifically_aligned_genes == 1 and not has_sa:
			return None


	# get general alignment information
	general_alignment_info = calculate_soft_clipping(cigartuple) # get 5' and 3' clipping length
	general_alignment_info[1] = read.infer_read_length() - general_alignment_info[0] - general_alignment_info[2] # get aligned read length
	# get NM rate
	nmrate = read.get_tag('NM') * 1.0 / read.query_alignment_length # get NM rate
	simplecigar = simplify_cigar(cigartuple) # get simple CIGAR string
	# get read strand
	if read.flag in [0, 2048]:
		readstrand = '+'
	elif read.flag in [16, 2064]:
		readstrand = '-'
	else:
		readstrand = '.' 

	# get read sequence and quality
	if is_record_readseq == False or read.is_secondary or read.is_supplementary:
		readsequence = ''
		readquality = ''
	else:
		readsequence = read.query_sequence
		readqualities = read.query_qualities
		readquality = ''.join([chr(c + 33) for c in readqualities if c is not None])

	readinfo = [chrom,
				start,
				end,
				read.query_name,
				readstrand,
				has_sa,
				general_alignment_info,
				read.mapping_quality,
				overlapped_segments[0][3],
				overlapped_segments[-1][3],
				overlapped_segments,
				readsequence,
				readquality,
				nmrate,
				simplecigar,
				read]
	return readinfo


def sortbyindex(index = 0):
	"""
	Sort multigene read by specified index
	Args:
		index: The index to use as the sort key (default: 0)
	Returns:
		A key function that extracts the specified index from the multigene read
	Example:
		# Sort by first index (default)
		all_multigene_read_info.sort(key=sortmultigeneread())
		
		# Sort by 2nd index
		all_multigene_read_info.sort(key=sortmultigeneread(1))
	"""
	def key_func(a):
		return a[index]
	return key_func

def get_fusion_within_read(readinfo: 
	list[str, int, int, str, str, bool, list[int], int, list[str], list[str], list[list[int, int, int, list[str]]], str, str, float, str, pysam.AlignedSegment]
	) -> list[str, str, str, str, int, int, int, str, list[int], int, int] | None:
	"""
	Detect fusion signals within a single read that spans multiple genes.
	
	For a read that aligns to gene A first, then gene B, this function detects
	the fusion event at the boundary between the two genes.
	
	Args:
		readinfo: Read information list containing exon segments with gene annotations
			Format: [chrom, start, end, read_name, strand, has_sa, cigar_info, mapq,
			         left_gene_list, right_gene_list, overlapped_segments, ...]
			overlapped_segments: [[start, end, length, gene_list], ...]
		
	Returns:
		list: List of fusion candidate information, or empty list if no fusion detected
	"""
	candidate = []
	
	# readinfo structure:
	# [0] chrom, [1] start, [2] end, [3] read_name, [4] strand, [5] has_sa,
	# [6] general_alignment_info, [7] mapq, [8] left_gene_list, [9] right_gene_list,
	# [10] overlapped_segments, [11] readsequence, [12] readquality, [13] nmrate, [14] simplecigar, [15] read object
	
	overlapped_segments = readinfo[10]
	if len(overlapped_segments) < 2:
		return None
	
	# Find the boundary where gene changes
	# Check each adjacent pair of exon segments
	for i in range(len(overlapped_segments) - 1):
		seg1 = overlapped_segments[i]  # [start, end, exon_length, gene_list]
		seg2 = overlapped_segments[i + 1]
		
		genes1 = set(seg1[3]) if seg1[3] else set()
		genes2 = set(seg2[3]) if seg2[3] else set()
		
		# Check if segments align to different genes
		if len(genes1) > 0 and len(genes2) > 0:
			# Find genes that are in one segment but not the other
			genes_only_in_seg1 = genes1 - genes2
			genes_only_in_seg2 = genes2 - genes1
			
			# If there are distinct genes on each side, it's a potential fusion
			if len(genes_only_in_seg1) > 0 and len(genes_only_in_seg2) > 0:
				# Use the gene with longest mapping length for each side
				gene1 = list(genes_only_in_seg1)[0]  # Could be improved to select best gene
				gene2 = list(genes_only_in_seg2)[0]
				
				# Skip if same gene
				if gene1 == gene2:
					continue
				
				# Calculate mapping lengths
				maplen1 = seg1[2]
				maplen2 = seg2[2]
				
				# Filter: each gene should have at least 100bp mapping
				if maplen1 < 50 or maplen2 < 50:
					continue
				
				# Determine breakpoints
				chrom = readinfo[0]
				if readinfo[4] == '+':
					bp1 = seg1[1]  # bp1 < bp2
					bp2 = seg2[0]
				else:
					bp1 = seg1[0]  # bp1 > bp2
					bp2 = seg2[1]
				gapsize = abs(bp2 - bp1)
				
				# Format: [gene1, gene2, 'splitread', chrom1, bp1, chrom2, bp2, read_name, [mapq], maplen1, maplen2, gapsize]
				if bp1 < bp2: # strand +
					candifusion = [gene1, gene2, 'splitread', chrom, bp1, chrom, bp2,
					               readinfo[3], [readinfo[7]], maplen1, maplen2, gapsize]
				else: # strand -
					candifusion = [gene2, gene1, 'splitread', chrom, bp2, chrom, bp1,
					               readinfo[3], [readinfo[7]], maplen2, maplen1, gapsize]
				candidate.append(candifusion)
	return candidate

def get_fusion_from_sameread(sameread: 
							list[list[str, int, int, str, str, bool, list[int], int, list[str], list[str], list[list[int, int, int, list[str]]], str, str, float, str, pysam.AlignedSegment]]
							) -> list[list[str, str, str, str, int, int, int, str, list[int], int, int]]:
	"""
	Get fusion events from a list of same read alignments
	Args:
		sameread: A list of same read alignments
	Returns:
		A list of fusion events
	"""
	candidate=[]
	# First, check for fusions within individual reads
	for read in sameread:
		candidate.append(get_fusion_within_read(read))
	
	# Then, check for fusions between different alignments of the same read
	for i in range(len(sameread)-1):
		read1 = sameread[i]
		for read2 in sameread[i + 1: ]:
			candidate.append(get_fusion_from_readpair(read1,read2))
	return candidate

"""
Entry: get multigene reads signals
"""
def get_split_reads(bampath: str, outpath: str, chrom: str, start: int, end: int, is_record_readseq: bool = True) -> int:
	"""
	Extract split reads from BAM file for a specific chromosome.
	
	Processes all reads in a BAM file for the given chromosome, identifying
	split reads (reads with supplementary alignments). Records detailed
	information about these reads to a file for further fusion detection.
	
	Args:
		bampath (str): Path to input BAM file
		outpath (str): Output directory path
		chrom (str): Chromosome name to process
		start (int): Start position to process
		end (int): End position to process
		is_record_readseq (bool): Whether to record read sequences and quality scores (default: True)

	Returns:
		int: Returns 0 on success
	"""
	bam = pysam.AlignmentFile(bampath, 'rb')
	allread = bam.fetch(chrom, start, end)

	multigene_read_info = []

	for read in allread:
		if read.is_secondary:
			continue
		read_info = stat_read_info(read, is_record_readseq)

		if read_info is not None:
			multigene_read_info.append(read_info)

	lines = ""
	for ll in multigene_read_info:
		exon_info = []
		for mm in ll[10]:
			exon_info.append(f"{mm[0]},{mm[1]},{mm[2]},{':'.join(mm[3])}") # exon_info: start, end, length, gene_name
		exon_info = ';'.join(exon_info)
		lines += f"{ll[0]}\t{ll[1]}\t{ll[2]}\t{ll[3]}\t{ll[4]}\t{ll[5]}\t{','.join(str(mm) for mm in ll[6])}\t{ll[7]}\t{ll[13]}\t{','.join(ll[8])}\t{','.join(ll[9])}\t{exon_info}\t{ll[14]}\t{ll[11]}\t{ll[12]}\n"
	with open(f"{outpath}raw_signal/multigene_read_{chrom}_{start}_{end}.txt", 'w') as fo:
		fo.write('chrom\tstart\tend\tread_name\tstrand\tis_supplementary\tgeneral_alignment_info\tmapq\tnmrate\tleft_gene\tright_gene\texon_info\tcigar\tsequence\tquality\n')
		fo.write(lines)
	
	"""
	polish the fusion events
	"""
	fusion_events = []
	lastname = None
	sameread = []

	for readinfo in multigene_read_info:
		if readinfo[5]: # is supplementary alignment
			pass
			###########################################################
			# TODO: process supplementary alignment
			###########################################################
		else:
			event = get_fusion_within_read(readinfo)
			if event is not None:
				fusion_events.extend(event)

		"""
		if read[3] == lastname:
			sameread.append(read)
			continue
		if lastname is None:
			lastname = read[3]
			sameread = [read]
			continue
		# Process previous group: check for fusions within individual reads and between alignments
		if len(sameread) >= 1:
			splitfusion = get_fusion_from_sameread(sameread)
			fusion_events.append(splitfusion)
		sameread = [read]
		lastname = read[3]
		"""

	"""# Process last group
	if len(sameread) >= 1:
		splitfusion = get_fusion_from_sameread(sameread)
		fusion_events.append(splitfusion)"""

	if not os.path.exists(f"{outpath}rawsignal.txt"):
		with open(f"{outpath}rawsignal.txt", 'w') as fo:
			fo.write('gene1\tgene2\tsplitread\tchrom1\tbp1\tchrom2\tbp2\tread_name\tmapq\tmaplen1\tmaplen2\tgapsize\n')
		
	with open(f"{outpath}rawsignal.txt", 'a') as fo:
		for event in fusion_events:
			#gene1, gene2, 'splitread', chrom1, bp1, chrom2, bp2, read1[3], [leftread[6],rightread[6]], maplen1, maplen2, gapsize
			mapq = [str(c) for c in event[8]]
			fo.write(f"{event[0]}\t{event[1]}\t{event[2]}\t{event[3]}\t{event[4]}\t{event[5]}\t{event[6]}\t{event[7]}\t{','.join(mapq)}\t{event[9]}\t{event[10]}\t{event[11]}\n")
	return 0


##########
##########
##########

def remove_ovlp_exon(oldreadinfo,side,removelength):
	"""
	Remove overlapping exon from a read
	Args:
		oldreadinfo: The old read info
		side: The side to remove the exon
		removelength: The length to remove
	Returns:
		A list of read info
	Example:
		oldreadinfo = [chrom, start, end, read_name, strand, cigar_info, mapq, left_gene, right_gene, exon_info, sequence, quality]
		side = 'right'
		removelength = 10
		readinfo = remove_ovlp_exon(oldreadinfo, side, removelength)
		print(readinfo)
	"""
	readinfo=[]
	exonlist=oldreadinfo[9]
	deletedlen=0

	inspos=oldreadinfo[11]
	if side=='right':
		inspos=inspos[::-1]
	includeins=0
	toremove=removelength
	for c in inspos:
		if toremove<=0:
			break
		if 'M' in c:
			toremove-=int(c[:-1])
			includeins+=int(c[:-1])
		if 'I' in c:
			toremove-=int(c[:-1])
	removelength=includeins
	
	if side=='right':
		while len(exonlist)>0 and removelength>0:
			if removelength>=exonlist[-1][2]-10 :
				removelength-=exonlist[-1][2]
				deletedlen+=exonlist[-1][2]
				exonlist=exonlist[:-1]
				continue
			exonlist[-1]=[exonlist[-1][0],exonlist[-1][1]-removelength,exonlist[-1][2]-removelength,exonlist[-1][3]]
			deletedlen+=removelength
			removelength=0
		if exonlist==[]:
			return []
		readinfo=[oldreadinfo[0],oldreadinfo[1],exonlist[-1][1],oldreadinfo[3],oldreadinfo[4],oldreadinfo[5],oldreadinfo[6],oldreadinfo[7],exonlist[-1][3],exonlist,oldreadinfo[10]]
	else:
		while len(exonlist)>0 and removelength>0:
			if removelength>=exonlist[0][2]-10:
				removelength-=exonlist[0][2]
				deletedlen+=exonlist[0][2]
				exonlist=exonlist[1:]
				continue
			exonlist[0]=[exonlist[0][0]+removelength,exonlist[0][1],exonlist[0][2]-removelength,exonlist[0][3]]
			deletedlen+=removelength
			removelength=0
		if exonlist==[]:
			return []
		readinfo=[oldreadinfo[0],exonlist[0][0],oldreadinfo[2],oldreadinfo[3],oldreadinfo[4],oldreadinfo[5],oldreadinfo[6],exonlist[0][3],oldreadinfo[8],exonlist,oldreadinfo[10]]
	return readinfo

###
###
###
def get_splitgene(readinfo, side):
	"""
	Get the gene name and mapping length for a specific side of a split read.
	
	Args:
		readinfo: Read information list
		side: 'left' or 'right' - which side of the read to analyze
		
	Returns:
		list: [gene_name, mapping_length]
			- gene_name (str): Name of the identified gene (empty string if none)
			- mapping_length (int): Total length of exons mapped to this gene
	"""
	genename = ''
	maplen = 0
	exonlist = readinfo[10]  # overlapped_segments: [[start, end, length, gene_list], ...]
	
	if side == 'right':
		exonlist = exonlist[::-1]  # Reverse for right side
	
	if len(exonlist) == 0:
		return ['', 0]
	
	# Get genes from the first segment (leftmost or rightmost depending on side)
	allcandigene = exonlist[0][3] if exonlist[0][3] else []
	
	if len(allcandigene) == 0:
		return ['', 0]
	
	if len(allcandigene) == 1:
		genename = allcandigene[0]
		if genename == '':
			return [genename, 0]
		# Calculate total mapping length for this gene
		for exon in exonlist:
			if genename in (exon[3] if exon[3] else []):
				maplen += exon[2]
		return [genename, maplen]
	
	# Multiple candidate genes - find the one with longest mapping length
	genemaplen = {}
	for candigene in allcandigene:
		genemaplen[candigene] = 0
	
	for exon in exonlist:
		exon_genes = exon[3] if exon[3] else []
		for candigene in allcandigene:
			if candigene in exon_genes:
				genemaplen[candigene] += exon[2]
	
	# Find gene with longest mapping length
	longestgene = ''
	longestlen = 0
	for candigene in allcandigene:
		if genemaplen[candigene] > longestlen:
			longestgene = candigene
			longestlen = genemaplen[candigene]
	
	genename = longestgene
	maplen = genemaplen[genename]
	return [genename, maplen]

def get_fusion_from_readpair(read1, read2):
	"""
	Detect fusion signals from a read pair
	Args:
		read1: The first read
		read2: The second read
	Returns:
		A list of fusion signals
	Example:
		read1 = [chrom, start, end, read_name, strand, cigar_info, mapq, left_gene, right_gene, exon_info, sequence, quality]
		read2 = [chrom, start, end, read_name, strand, cigar_info, mapq, left_gene, right_gene, exon_info, sequence, quality]
		fusion = get_fusion_readpair(read1, read2)
	"""
	read1=list(read1)
	read2=list(read2)
	read1[5]=list(read1[5])
	read2[5]=list(read2[5])

	candifusion=[]
	if read1[4]=='-':
		temp=read1[5][0]; read1[5][0]=read1[5][2]; read1[5][2]=temp
	if read2[4]=='-':
		temp=read2[5][0];read2[5][0]=read2[5][2];read2[5][2]=temp

	if read1[5][0]<read2[5][0]:
		leftread=read1
		rightread=read2
	else:
		leftread=read2
		rightread=read1

	gapsize=rightread[5][0]-leftread[5][0]-leftread[5][1]

	
	if 0-gapsize > min(100,leftread[5][1]*0.5,rightread[5][1]*0.5):
		return []
	
	if gapsize<-100:
		
		removelength=0-gapsize
		if leftread[10]>rightread[10]:
			if leftread[4]=='+':
				leftread=remove_ovlp_exon(leftread,'right',removelength)
			else:
				leftread=remove_ovlp_exon(leftread,'left',removelength)
		else:
			if rightread[4]=='+':
				rightread=remove_ovlp_exon(rightread,'left',removelength)
			else:
				rightread=remove_ovlp_exon(rightread,'right',removelength)

	if leftread==[] or rightread==[]:
		return []

	if leftread[4]=='+':
		chrom1=leftread[0]
		bp1=leftread[2]
		[gene1,maplen1]=get_splitgene(leftread,'right')
	else:
		chrom1=leftread[0]
		bp1=leftread[1]
		[gene1,maplen1]=get_splitgene(leftread,'left')


	if rightread[4]=='+':
		chrom2=rightread[0]
		bp2=rightread[1]
		[gene2,maplen2]=get_splitgene(rightread,'left')
	else:
		chrom2=rightread[0]
		bp2=rightread[2]
		[gene2,maplen2]=get_splitgene(rightread,'right')

	if gene1 =='' or gene2=='' or gene1==gene2 or maplen1<100 or maplen2<100:
		return []


	if chrom1<chrom2:
		candifusion=[gene1,gene2,'splitread',chrom1,bp1,chrom2,bp2,read1[3],[leftread[6],rightread[6]],maplen1,maplen2,gapsize]
	elif chrom2<chrom1:
		candifusion=[gene2,gene1,'splitread',chrom2,bp2,chrom1,bp1,read1[3],[rightread[6],leftread[6]],maplen2,maplen1,gapsize]
	else:
		if bp1<bp2:
			candifusion=[gene1,gene2,'splitread',chrom1,bp1,chrom2,bp2,read1[3],[leftread[6],rightread[6]],maplen1,maplen2,gapsize]
		else:
			candifusion=[gene2,gene1,'splitread',chrom2,bp2,chrom1,bp1,read1[3],[rightread[6],leftread[6]],maplen2,maplen1,gapsize]


	return [candifusion]

