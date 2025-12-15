from typing import List, Dict
import re

def sortgeneinfo(index = 0):
	"""
	Sort gene info by specified index
	Args:
		index: The index to use as the sort key (default: 0)
	Returns:
		A key function that extracts the specified index from the gene info
	Example:
		# Sort by first index (default)
		geneinfo[chrom].sort(key=sortgeneinfo())
		
		# Sort by 2nd index
		geneinfo[chrom].sort(key=sortgeneinfo(1))
	"""
	def key_func(a):
		return a[index]
	return key_func

# Pre-compile regex patterns for gene name extraction
GENE_ID_PATTERN = re.compile(r'gene_id\s+"([^"]+)"')
GENE_NAME_PATTERN = re.compile(r'gene_name\s+"([^"]+)"')
GENE_PATTERN = re.compile(r'gene\s+"([^"]+)"')

def _extract_gene_name(attributes: str, use_geneid: bool) -> str:
	"""
	Extract gene name from GTF attributes column using regex.
	
	Args:
		attributes: The 9th column (attributes) from GTF file
		use_geneid: If True, use gene_id; otherwise use gene_name or gene
		
	Returns:
		Gene name string
	"""
	if use_geneid:
		match = GENE_ID_PATTERN.search(attributes)
		if match:
			return match.group(1)
	else:
		match = GENE_NAME_PATTERN.search(attributes)
		if match:
			return match.group(1)
		match = GENE_PATTERN.search(attributes)
		if match:
			return match.group(1)
	return None

def write_geneinfo(geneinfo: dict, outpath: str) -> None:
	"""
	Write gene info to tsv file
	"""
	with open(outpath + 'geneinfo.tsv', 'w') as f:
		for chrom in geneinfo.keys():
			for gene in geneinfo[chrom]:
				f.write(f"{chrom}\t{gene[0]}\t{gene[1]}\t{gene[2]}\t{gene[3]}\t")
				for exon in gene[4]:
					f.write(f"{exon[0]}-{exon[1]};")
				f.write("\n")

def create(gtfinfo: Dict[str, List[str]], 
			chrom_valid: List[str], 
			use_geneid: bool, 
			outpath: str = './fusionseeker_out/',
			writedata: bool = False) -> dict[str, list[list[int, int, str, str, list[list[int, int]]]]]:
	"""
    Construct gene info from GTF file
	"""
	geneinfo = {}
	log_messages = []  # Collect log messages for batch writing
	
	# Pre-filter and parse lines once
	feature_set = {'gene', 'exon'}  # Use set for O(1) lookup
	
	# construct gene info for each chromosome
	for chrom in chrom_valid:
		geneinfo[chrom] = []
		
		# Pre-process lines: split once and filter
		parsed_lines = []
		for line in gtfinfo[chrom]:
			if not line.strip():  # Skip empty lines
				continue
			fields = line.rstrip('\n').split('\t')
			if len(fields) < 9:
				continue
			if fields[0] == chrom and fields[2] in feature_set:
				parsed_lines.append(fields)
		
		lastgene = None
		genepos = []
		exonpos = []
		gene_start = None
		gene_end = None
		gene_strand = None
		
		for fields in parsed_lines:
			feature_type = fields[2]
			start = int(fields[3])
			end = int(fields[4])
			strand = fields[6]
			attributes = fields[8]
			
			# Extract gene name using regex (more efficient)
			genename = _extract_gene_name(attributes, use_geneid)
			
			if genename is None:
				# Fallback to coordinates
				genename = f"{chrom}:{start}-{end}"
				log_messages.append(f"Gene {chrom}:{start}-{end} does not have a {'gene_id' if use_geneid else 'gene_name/gene'}, using coordinates instead\n")
			
			if lastgene is None:
				lastgene = genename
				gene_start = start
				gene_end = end
				gene_strand = strand
			
			if genename == lastgene:
				if feature_type == 'exon':
					exonpos.append([start, end])
				elif feature_type == 'gene':
					# Update gene boundaries
					gene_start = min(gene_start, start)
					gene_end = max(gene_end, end)
			else:
				# Save previous gene
				if lastgene is not None and exonpos:
					genepos = [gene_start, gene_end, lastgene, gene_strand]
					geneinfo[chrom].append(genepos + [exonpos])
				
				# Start new gene
				lastgene = genename
				gene_start = start
				gene_end = end
				gene_strand = strand
				exonpos = []
				if feature_type == 'exon':
					exonpos.append([start, end])
		
		# Save last gene
		if lastgene is not None and exonpos:
			genepos = [gene_start, gene_end, lastgene, gene_strand]
			geneinfo[chrom].append(genepos + [exonpos])

	# sort gene info by position
	for chrom in geneinfo:
		geneinfo[chrom].sort(key=sortgeneinfo())
	
	# Write log messages in batch (more efficient than opening file multiple times)
	if log_messages:
		with open(outpath + 'log.txt', 'a') as f:
			f.writelines(log_messages)
	
	# write gene info to file if set True
	if writedata:
		write_geneinfo(geneinfo, outpath)

	return geneinfo


