import pybedtools
import pysam
import sys
import os

def extract_multigene_read(input_bam: str, input_gtf: str, output_path: str) -> None:
    """
    Filters a BAM file to retain only those alignments (reads) that overlap
    more than one unique gene feature defined in the GTF file.

    Args:
        input_bam (str): Path to the input Iso-seq BAM file.
        input_gtf (str): Path to the input gene annotation GTF file.
        output_path (str): Path for the output filtered BAM file.
    """

    os.system(f'samtools view -b -F 260 {input_bam} > {output_path}temp.bam')

    # --- Step 1: Prepare unique gene intervals from GTF ---
    gene_features = pybedtools.BedTool(input_gtf).filter(lambda x: x[2] == 'gene')
    
    # Merge gene features to ensure each gene corresponds to one merged interval
    # and to properly handle redundancy/overlap within the same gene.
    # The 'd' column (col 4 in BED) is the attribute field in the original GTF,
    # which is kept as a distinct feature after merging.
    try:
        unique_genes_bed = gene_features.merge(s = False)
    except Exception as e:
        sys.exit(f"Error during pybedtools merge/filter: {e}")
        
    with open(f'{output_path}log.txt', 'a') as logfile:
        logfile.write(f"Unique gene intervals created. Total intervals: {len(unique_genes_bed)}\n")

    # --- Step 2: Intersect reads with unique genes and count overlaps ---

    # Intersect the BAM reads (-a) with the unique gene intervals (-b).
    # c=True: Count the number of overlaps with B for each feature in A.
    # wa=True: Write the original feature A (the read alignment)
    # split=True: Handle spliced alignments (required for Iso-seq/RNA-seq data)
    intersect_results = pybedtools.BedTool(f"{output_path}temp.bam").intersect(
        unique_genes_bed, wa=True, c=True, split=True, bed=True
    )

    # --- Step 3: Filter for reads overlapping > 1 gene and collect names ---
    multi_overlap_read_names = set()
    
    # The output is a BED file where the last column (field 7) is the count.
    # pybedtools results are 0-indexed lists of fields.
    for read_overlap in intersect_results:
        overlap_count = int(read_overlap.fields[-1]) # Last column is the count
        if overlap_count > 1:
            read_name = read_overlap.name # The name column (field 4)
            multi_overlap_read_names.add(read_name)

    read_count = len(multi_overlap_read_names)
    with open(f'{output_path}log.txt', 'a') as logfile:
        logfile.write(f"Found {read_count} unique reads overlapping multiple genes.\n")

    if read_count == 0:
        with open(f'{output_path}log.txt', 'a') as logfile:
            logfile.write("No multi-gene overlapping reads found. Creating an empty BAM file.\n")
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
