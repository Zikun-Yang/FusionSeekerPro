from collections import defaultdict
import pybedtools
import pysam
import sys
import os
import re

def extract_gene_id_from_gtf(gtf_interval):
    """Extract gene_id from the attributes field of the GTF interval"""
    if len(gtf_interval.fields) > 8:
        attributes = gtf_interval.fields[8]  # The 9th column (index 8) is the attributes
        keys = ['gene_id', 'gene_name', "name"]
        for key in keys:
            match = re.search(fr'{key}\s+"([^"]+)"', attributes)
            if match:
                return match.group(1)
    return '.'

def extract_multigene_read(input_bam: str, input_gtf: str, output_path: str) -> None:
    """
    Filters a BAM file to retain only those alignments (reads) that overlap
    more than one unique gene feature defined in the GTF file.

    Args:
        input_bam (str): Path to the input Iso-seq BAM file.
        input_gtf (str): Path to the input gene annotation GTF file.
        output_path (str): Path for the output filtered BAM file.
    """

    os.system(f'samtools view -b -F 260 {input_bam} > {output_path}temp.bam') # keep only primary and supplementary alignments

    # --- Step 1: Prepare unique gene intervals from GTF ---
    gene_features = pybedtools.BedTool(input_gtf).filter(lambda x: x[2] == 'gene')
    
    """
    # Merge gene features to ensure each gene corresponds to one merged interval
    # and to properly handle redundancy/overlap within the same gene.
    # The 'd' column (col 4 in BED) is the attribute field in the original GTF,
    # which is kept as a distinct feature after merging.
    try:
        unique_genes_bed = gene_features.merge(s = False)
    except Exception as e:
        sys.exit(f"[ERROR] Error during pybedtools merge/filter: {e}")
    """
    unique_genes_bed = gene_features.each(lambda x: pybedtools.create_interval_from_list([
        x.chrom, 
        str(int(x.start)),
        str(int(x.end)),
        extract_gene_id_from_gtf(x),
        '0',
        x.strand
    ])).saveas()
        
    with open(f'{output_path}log.txt', 'a') as logfile:
        logfile.write(f"[INFO] Unique gene intervals created. Total intervals: {len(unique_genes_bed)}\n")

    # --- Step 2: Intersect reads with unique genes and count overlaps ---
    # split=True: Handle spliced alignments (required for Iso-seq/RNA-seq data)
    intersect_results = pybedtools.BedTool(f"{output_path}temp.bam").intersect(
        unique_genes_bed, wa=True, wb=True, split=True, bed=True
    )

    read_to_genes = defaultdict(set)
    for result in intersect_results:
        read_name = result.name  # The name of the read (field 4)
        gene_key = (result.fields[12], int(result.fields[13]), int(result.fields[14]))  # chr, start, end
        read_to_genes[read_name].add(gene_key)

    # Count the number of unique genes overlapping for each read
    multi_overlap_read_names = set()
    for read_name, gene_set in read_to_genes.items():
        if len(gene_set) >= 2:
            multi_overlap_read_names.add(read_name)
    
    """
    # --- Step 3: Filter for reads overlapping > 1 gene and collect names ---
    multi_overlap_read_names = set()
    
    # The output is a BED file where the last column (field 7) is the count.
    # pybedtools results are 0-indexed lists of fields.
    for read_overlap in intersect_results:
        overlap_count = int(read_overlap.fields[-1]) # Last column is the count
        if overlap_count > 1:
            read_name = read_overlap.name # The name column (field 4)
            multi_overlap_read_names.add(read_name)
    """

    read_count = len(multi_overlap_read_names)
    with open(f'{output_path}log.txt', 'a') as logfile:
        logfile.write(f"[INFO] Found {read_count} unique reads overlapping multiple genes.\n")

    if read_count == 0:
        with open(f'{output_path}log.txt', 'a') as logfile:
            logfile.write("[WARNING] No multi-gene overlapping reads found. Creating an empty BAM file.\n")
        # Proceed to step 4, which will create an empty BAM.
    
    try:
        # Open the input BAM file and the output file
        in_bam = pysam.AlignmentFile(f"{output_path}temp.bam", "rb")
        out_bam = pysam.AlignmentFile(f"{output_path}filtered.bam", "wb", header = in_bam.header)

        # Iterate through the entire input BAM
        for read in in_bam:
            # Check if the read name is in our filtered set
            if read.query_name in multi_overlap_read_names:
                out_bam.write(read)
        
        in_bam.close()
        out_bam.close()
        
    except Exception as e:
        sys.exit(f"Error during pysam BAM handling: {e}")
        
    # make index
    os.system(f'samtools index {output_path}filtered.bam')
    os.system(f"rm {output_path}temp.bam")
