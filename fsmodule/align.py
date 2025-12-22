import os
from . import stat_split_reads
from . import cluster

def is_empty(file_path: str) -> bool:
	"""
	Check if a poa fasta file is empty
	Args:
		file_path: str, the path to the file
	"""
	if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
		return True
	lines = open(file_path, 'r').readlines()
	seq = lines[1].strip()
	return seq == ''

def get_fastq_seq_num(fastq_path: str) -> int:
	"""
	Get the number of sequences in a fastq file
	Args:
		fastq_path: str, the path to the fastq file
	"""
	num_seq = 0
	with open(fastq_path, 'r') as fi:
		for line in fi:
			if line.startswith('@'):
				num_seq += 1
	return num_seq

def downsample_fastq(fastq_path: str, new_num_seq: int) -> int:
	"""
	Downsample a fastq file
	Args:
		fastq_path: str, the path to the fastq file
		new_num_seq: int, the number of sequences to downsample to
	"""
	with open(fastq_path, 'r') as fi:
		lines = fi.readlines()
		new_lines = lines[0 : new_num_seq * 4]
	with open(fastq_path, 'w') as fo:
		fo.writelines(new_lines)
	return 0

def poa_all(outpath: str, genomic_regions: list[tuple[str, int, int]]) -> int:
	"""
	Generate aligned consensus sequence for each gene fusion
	Args:
		outpath: str, the path to the output directory
		chrom_valid: list[str], the list of valid chromosomes
	"""
	os.makedirs(outpath + 'align_workspace/', exist_ok=True)
	confident_genefusion = []
	with open(outpath + 'confident_genefusion.txt', 'r') as fi:
		for line in fi:
			if line.startswith('gene1'):
				continue
			confident_genefusion.append(line.strip().split('\t'))

	# get supporting reads' names
	supporting_reads = []
	for event in confident_genefusion:
		supporting_reads.extend(event[8].split(',')) # read names
	
	# get supporting reads' sequences and qualities
	allreadseq = {}
	for chrom, start, end in genomic_regions:
		multigene_read_info = []
		with open(f'{outpath}raw_signal/multigene_read_{chrom}_{start}_{end}.txt', 'r') as fi:
			for line in fi:
				if line.startswith('chrom'):
					continue
				multigene_read_info.append(line.strip('\n').split('\t'))
		for readinfo in multigene_read_info:
			if readinfo[3] in supporting_reads and readinfo[13] != '':
				allreadseq[readinfo[3]] = (readinfo[13], readinfo[14]) # (sequence, quality)
	
	# generate poa sequence for each gene fusion
	with open(f'{outpath}align_workspace/allalignedseq.fa', 'w') as fo:
		for event in confident_genefusion:
			gfinfo = f"{event[7]}_{event[0]}_{event[1]}_{event[3]}_{event[4]}_{event[5]}_{event[6]}" # ID_gene1_gene2_chrom1_bp1_chrom2_bp2
			with open(f'{outpath}align_workspace/suppread_{gfinfo}.fastq','w') as fastq:
				for readname in event[8].split(','): # read names
					if readname not in allreadseq:
						continue
					if allreadseq[readname][1] != '':
						fastq.write(f"@{readname}\n{allreadseq[readname][0]}\n+\n{allreadseq[readname][1]}\n")
					else:
						fastq.write(f"@{readname}\n{allreadseq[readname][0]}\n+\n{'I' * len(allreadseq[readname][0])}\n")
			while not os.path.exists(f'{outpath}align_workspace/poa_{gfinfo}.fa') or is_empty(f'{outpath}align_workspace/poa_{gfinfo}.fa'):
				os.system(f'bsalign poa {outpath}align_workspace/suppread_{gfinfo}.fastq -o {outpath}align_workspace/poa_{gfinfo}.fa > {outpath}align_workspace/poa_{gfinfo}.log 2>&1')
				if is_empty(f'{outpath}align_workspace/poa_{gfinfo}.fa'):
					num_seq = get_fastq_seq_num(f'{outpath}align_workspace/suppread_{gfinfo}.fastq')
					if num_seq > 50: 
						new_num_seq = num_seq - 10
					elif num_seq > 25:
						new_num_seq = num_seq - 5
					else:
						new_num_seq = num_seq - 1
					print(f"Reducing sequence number from {num_seq} to {new_num_seq}")
					if new_num_seq <= 1:
						break
					downsample_fastq(f'{outpath}align_workspace/suppread_{gfinfo}.fastq', new_num_seq)
					with open(f'{outpath}log.txt', 'a') as logfile:
						logfile.write(f"Reducing sequence number from {num_seq} to {new_num_seq}\n, retrying...\n")
			poactg = open(f'{outpath}align_workspace/poa_{gfinfo}.fa','r').read().strip().split('\n') # ['cns_seq', xxxxxxxxx]
			if len(poactg) != 2:
				continue
			fo.write(f">poa_ctg_{gfinfo}\n{poactg[1]}\n")

	return 0

def polish_bp(outpath: str, genomic_regions: list[tuple[str, int, int]], reference: str, datatype: str, polishbp: bool, skipsupplementary: bool) -> int:
	"""
	Polish breakpoints of gene fusion events
	Args:
		outpath: str, the path to the output directory
		chrom_valid: list[str], the list of valid chromosomes
		reference: str, the path to the reference genome
		datatype: str, the type of data (isoseq or long-read)
		polishbp: bool, whether to polish the breakpoints (need reference genome if set to True)
		skipsupplementary: bool, whether to skip the supplementary alignments
	"""
	stat_split_reads.geneinfo = geneinfo

	confident_genefusion = []
	gfinfo = {}
	with open(f'{outpath}confident_genefusion.txt', 'r') as fi:
		for line in fi:
			if line.startswith('gene1'):
				continue
			ll = line.strip('\n').split('\t') # [Gene1, Gene2, NumSupp, Chrom1, Breakpoint1, Chrom2, Breakpoint2, ID, SupportingReads]
			# Convert breakpoint positions to integers for numerical operations in merge_genepair
			ll[4] = int(ll[4])  # Breakpoint1
			ll[6] = int(ll[6])  # Breakpoint2
			confident_genefusion.append(ll)
			gfinfo[f"{ll[7]}_{ll[0]}_{ll[1]}_{ll[3]}_{ll[4]}_{ll[5]}_{ll[6]}"] = ll # ID_gene1_gene2_chrom1_bp1_chrom2_bp2

	if polishbp:
		# align all consensus sequences to reference genome
		if datatype == "isoseq":	
			os.system(f"minimap2 -t 8 -ax splice:hq {reference} --secondary=no {outpath}align_workspace/allalignedseq.fa | samtools sort -o {outpath}align_workspace/allalignedseq.bam")
		else:
			os.system(f"minimap2 -t 8 -ax splice {reference} --secondary=no {outpath}align_workspace/allalignedseq.fa | samtools sort -o {outpath}align_workspace/allalignedseq.bam")
	
		os.system(f"samtools index {outpath}align_workspace/allalignedseq.bam")
		os.system(f"mkdir -p {outpath}align_workspace/raw_signal/")
		# analyse the gene fusion events with aligned consensus sequences
		for chrom, start, end in genomic_regions:
			stat_split_reads.get_split_reads(f"{outpath}align_workspace/allalignedseq.bam",
											f"{outpath}align_workspace/",
											chrom,
											start,
											end,
											is_record_readseq = False)

		# analyze supplementary alignment
		if not skipsupplementary:
			supplementary_read_info = {}
			for chrom, start, end in genomic_regions:
				with open(f"{outpath}align_workspace/raw_signal/multigene_read_{chrom}_{start}_{end}.txt", 'r') as fi:
					for line in fi:
						if line.startswith('chrom'):
							continue
						readinfo = line.strip().split('\t') # chrom, start, end, read_name, strand, is_supplementary, general_alignment_info, mapq, nmrate, left_gene, right_gene, exon_info, cigar, sequence, quality
						readinfo[1], readinfo[2], readinfo[5] = int(readinfo[1]), int(readinfo[2]), bool(readinfo[5])
						if readinfo[5]: # is supplementary alignment
							if readinfo[3] not in supplementary_read_info.keys():
								supplementary_read_info[readinfo[3]] = []
							supplementary_read_info[readinfo[3]].append(readinfo[0:3])
			for read_name, read_pos_list in supplementary_read_info.items():
				stat_split_reads.reanalyze_supplementary_alignment(f"{outpath}align_workspace/allalignedseq.bam",
																	f"{outpath}align_workspace/",
																	read_name,
																	read_pos_list)

		# update
		with open(f"{outpath}align_workspace/rawsignal.txt", "r") as fi:
			for line in fi:
				if line.startswith('gene1'):
					continue
				ll = line.strip('\n').split('\t') # [Gene1, Gene2, splitread, chrom1, bp1, chrom2, bp2, read_name, mapq, maplen1, maplen2, gapsize]
				key = ll[7].replace("poa_ctg_","") # ID_gene1_gene2_chrom1_bp1_chrom2_bp2
				oldinfo = gfinfo[key]
				# Convert breakpoint positions to integers for numerical operations
				newinfo = [ll[0], ll[1], oldinfo[2], ll[3], int(ll[4]), ll[5], int(ll[6]), oldinfo[7], oldinfo[8]] # [Gene1, Gene2, NumSupp, Chrom1, Breakpoint1, Chrom2, Breakpoint2, ID, SupportingReads]
				gfinfo[key] = newinfo

	# merge events with the same gene pair
	merged_events = []
	for event in gfinfo.keys():
		merged_events.append(gfinfo[event])
	merged_events = cluster.merge_genepair(merged_events)

	# write refined gene fusion events
	with open(f"{outpath}confident_genefusion_refined.txt", 'w') as fo:
		fo.write('ID\tgene1\tgene2\tnumSupp\tchrom1\tbp1\tchrom2\tbp2\treadNames\n')
		for event in merged_events:
			fo.write(f"{event[7]}\t{event[0]}\t{event[1]}\t{event[2]}\t{event[3]}\t{event[4]}\t{event[5]}\t{event[6]}\t{event[8]}\n")

	with open(f"{outpath}align_workspace/allalignedseq.fa", 'r') as fi:
		allisoseq = fi.read().split('\n')[:-1]
	
	poaseq={}
	for idx in range(len(allisoseq)):
		seq = allisoseq[idx]
		if seq[0] != '>' or idx + 1 >= len(allisoseq):
			continue
		poaseq[seq.split('_')[2]] = allisoseq[idx + 1] # ID -> sequence

	# write transcript sequences for each gene fusion event
	with open(f"{outpath}confident_genefusion_refined_transcript_sequence.fa", 'w') as fo:	
		for event in merged_events:
			fo.write(f">{event[7]}_{event[0]}_{event[1]}_{event[3]}_{event[4]}_{event[5]}_{event[6]}\n{poaseq[event[7]]}\n")

	return 0