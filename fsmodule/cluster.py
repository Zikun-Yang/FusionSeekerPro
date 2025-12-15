import statistics

def findmax(neighbor: dict[str, list[str]]) -> str:
	"""
	Find the key with the maximum length in the neighbor dictionary.
	Args:
		neighbor: A dictionary of neighbor lists.
	Returns:
		A string of the key with the maximum length.
	"""
	maxlen = 0
	maxkey = list(neighbor.keys())[0]
	for key in neighbor.keys():
		if len(neighbor[key]) > maxlen:
			maxkey = key
			maxlen = len(neighbor[key])
	return maxkey

def sort_by_support_count(a):
	return a[2]

def merge_event_with_same_genepair(event_with_same_genepair: list[list[str, str, int, str, int, str, int, int, str]]
								) -> list[list[str, str, int, str, int, str, int, int, str]]:
	"""
	Merge gene fusion events with the same gene pair.
	Args:
		samepair: A list of gene fusion events with the same gene pair.
	Returns:
		A list of merged gene fusion events.
	"""
	if len(event_with_same_genepair) == 0:
		return []
	if len(event_with_same_genepair) == 1:
		return event_with_same_genepair
	event_with_same_genepair.sort(key = sort_by_support_count, reverse=True)
	candidate = {}
	candidate[1] = event_with_same_genepair[0]
	numcandi=1
	for query in event_with_same_genepair[1:]:
		is_close = False
		for i in range(numcandi):
			reference = candidate[i + 1]
			if abs(query[4]-reference[4]) <= 2000 and abs(query[6] - reference[6]) <= 2000: # bp1 and bp2 are close
				is_close = True
				candidate[i+1] = [reference[0], reference[1], int(reference[2]) + int(query[2]), reference[3], reference[4], reference[5], reference[6], reference[7], reference[8]]
				break
		if not is_close:
			candidate[numcandi + 1] = query
			numcandi += 1
	merged_events = []
	for i in range(numcandi):
		merged_events.append(candidate[i+1])
	return merged_events

def merge_genepair(clustered_candidates: list[list[str, str, int, str, int, str, int, int, str]]
                   	) -> list[list[str, str, int, str, int, str, int, int, str]]:
	"""
	Merge clustered gene fusion events with the same gene pair.
	Args:
		clustered_candidates: A list of clustered gene fusion events.
	Returns:
		A list of merged clustered gene fusion events.
	"""
	candidates = []
	event_with_same_genepair = []
	gene1, gene2 = None, None

	for event in clustered_candidates:
		if event[0] == gene1 and event[1] == gene2:
			event_with_same_genepair.append(event)
			continue
		if event_with_same_genepair:
			candidates += merge_event_with_same_genepair(event_with_same_genepair)
			event_with_same_genepair = []
		event_with_same_genepair.append(event)
		gene1 = event[0]
		gene2 = event[1]
	candidates += merge_event_with_same_genepair(event_with_same_genepair)
	return candidates

def merge_cluster(cluster: list[str], 
		id2signal: dict[str, list[str, str, str, str, int, str, int, str, int, int, int, int]],
		outpath: str) -> list[str, str, int, str, int, str, int, int, str]:
	"""
	Merge gene fusion events in the same cluster.
	Args:
		cluster: A list of gene fusion events in the same cluster.
		id2signal: A dictionary of gene fusion events.
		# id2signal[id] = [gene1, gene2, splitread, chrom1, bp1, chrom2, bp2, read_name, mapq, maplen1, maplen2, gapsize]
		outpath: The path to the output directory.
	Returns:
		A string of merged gene fusion events in the same cluster.
	"""
	signal_list = []
	for c in cluster:
		signal_list.append(id2signal[c])
	s1 = signal_list[0]
	bp1 = [int(c[4]) for c in signal_list]
	bp2 = [int(c[6]) for c in signal_list]
	bp1 = statistics.multimode(bp1)
	bp2 = statistics.multimode(bp2)
	if len(bp1) > 1:
		with open(f'{outpath}log.txt', 'a') as logfile:
			logfile.write(f"[WARNING] Multiple breakpoints found for gene fusion {s1[0]} and {s1[1]}. Using the first most frequent breakpoint.\n")
	if len(bp2) > 1:
		with open(f'{outpath}log.txt', 'a') as logfile:
			logfile.write(f"[WARNING] Multiple breakpoints found for gene fusion {s1[0]} and {s1[1]}. Using the first most frequent breakpoint.\n")
	
	bp1 = bp1[0]
	bp2 = bp2[0]

	quality = [int(q[8]) for q in signal_list]
	quality = int(int(sum(quality)/len(signal_list)/2))
	reads = ','.join([c[7] for c in signal_list])
	numsupp = len(set([c[7] for c in signal_list]))
	merged_signal = [s1[0], s1[1], numsupp, s1[3], bp1, s1[5], bp2, quality, reads] # gene1, gene2, numsupp, chrom1, bp1, chrom2, bp2, quality, reads
	return merged_signal

def cluster_by_dbscan(rawsignals: list[list[str, str, str, str, int, str, int, str, int, int, int, int]], 
					maxdistance: int,
					outpath: str) -> list[list[str, str, int, str, int, str, int, int, str]]:
	"""
	Cluster raw signals into gene fusion candidates using DBSCAN.
	Args:
		rawsignals: A list of raw signal lists.
		maxdistance: The maximum distance between two breakpoints.
		outpath: The path to the output directory.
	Returns:
		A list of clustered gene fusion events.
	"""
	candidates = []
	neighbor = {}
	id2signal = {}
	for i in range(len(rawsignals)):
		id2signal[f'gf{i}'] = rawsignals[i]
		neighbor[f'gf{i}'] = []
	for i in range(len(rawsignals)-1):
		gf1 = id2signal[f'gf{i}']
		for j in range(i + 1, len(rawsignals)):
			gf2 = id2signal[f'gf{j}']
			distance = ( (gf1[4] - gf2[4]) ** 2 + (gf1[6] - gf2[6]) ** 2) ** 0.5
			if distance <= maxdistance:
				neighbor[f'gf{i}'].append(f'gf{j}')
				neighbor[f'gf{j}'].append(f'gf{i}')

	done=[]
	while neighbor != {}:
		startkey = findmax(neighbor)
		cluster = [startkey]
		newadd = neighbor.pop(startkey)
		newadd = [c for c in newadd if c in neighbor]
		while newadd != []:
			if newadd[0] not in neighbor:
				newadd.remove(newadd[0])
				continue
			cluster.append(newadd[0]) # add the new signal to the cluster
			newneighbor = neighbor.pop(newadd[0])
			if len(newneighbor) >= 3: # expand if have > 3 neighbors
				newneighbor = [c for c in newneighbor if c not in cluster + newadd + done] 
				newadd += newneighbor
			newadd.remove(newadd[0]) # remove the first element from the new add list
		candidates.append(merge_cluster(cluster, id2signal, outpath)) # merge the gene fusion events in the same cluster
		done.extend(cluster)
	return candidates

def cluster_bp(outpath: str, maxdistance: int, min_supp: int | None) -> int:
	"""
	Cluster raw signals into gene fusion candidates.
	Args:
		outpath: The path to the output directory.
		maxdistance: The maximum distance between two breakpoints.
		min_supp: The minimum number of supporting reads.
	Returns:
		An integer indicating the success or failure of the clustering.
	"""
	rawsignal = open(outpath + 'rawsignal.txt', 'r').read().split('\n')[:-1]
	rawsignal = rawsignal[1:]  # skip the header
	rawsignal_grouped_by_genepair = {}
	clustered_candidates = []

	if min_supp is None:
		min_supp = 3 + int(len(rawsignal) / 50000)

	# group raw signals by gene pair
	for signal in rawsignal:
		ss = signal.split('\t')
		ss[4] = int(ss[4])
		ss[6] = int(ss[6])
		ss[8] = int(ss[8])
		ss[9] = int(ss[9])
		ss[10] = int(ss[10])
		ss[11] = int(ss[11])

		genepair = f"{ss[0]}\t{ss[1]}"
		if genepair in rawsignal_grouped_by_genepair.keys():
			rawsignal_grouped_by_genepair[genepair].append(ss)
			continue

		genepair = f"{ss[1]}\t{ss[0]}"
		if genepair in rawsignal_grouped_by_genepair.keys():
			ss_reordered = [ss[1], ss[0], ss[2], ss[5], ss[6], ss[3], ss[4], ss[7], ss[8], ss[9], ss[10], ss[11]]
			rawsignal_grouped_by_genepair[genepair].append(ss_reordered)
			continue

		genepair = f"{ss[0]}\t{ss[1]}"
		rawsignal_grouped_by_genepair[genepair] = [ss]

	# cluster raw signals by gene pair
	cluster_id = 1
	for genepair in rawsignal_grouped_by_genepair.keys():
		if len(rawsignal_grouped_by_genepair[genepair]) == 1:
			ss = rawsignal_grouped_by_genepair[genepair][0]
			clustered_candidates.append([ss[0], ss[1], 1, ss[3], ss[4], ss[5], ss[6], ss[8], ss[7]])
			continue
		candidates = cluster_by_dbscan(rawsignal_grouped_by_genepair[genepair], maxdistance, outpath)
		clustered_candidates.extend(candidates)
		cluster_id += 1

	# record all candidate gene fusion events
	clustered_candidates_str = []
	for c in clustered_candidates:
		cc = [str(c) for c in c]
		clustered_candidates_str.append("\t".join(cc))
	with open(outpath + 'clustered_candidate.txt', 'w') as fo:
		fo.write("gene1\tgene2\tnumSupp\tchrom1\tbp1\tchrom2\tbp2\tmapq\treadNames\n")
		fo.write("\n".join(clustered_candidates_str))

	# filter out candidate gene fusion events with less than 3 or min_supp supporting reads
	clustered_candidates = [c for c in clustered_candidates if int(c[2]) >= min_supp]
	
	"""
	# merge gene fusion events with the same gene pair
	clustered_candidates = merge_genepair(clustered_candidates)
	clustered_candidates = [c for c in clustered_candidates if c[2] >= min_supp]
	genepairs = [[c[0], c[1]] for c in clustered_candidates]
	totalgf = len(genepairs)
	genecount = {}
	for c in genepairs:
		if c[0] not in genecount:
			genecount[c[0]] = 1
		else:
			genecount[c[0]] += 1
		if c[1] not in genecount:
			genecount[c[1]] = 1
		else:
			genecount[c[1]] += 1
	repeating = [c for c in genecount if genecount[c] >= max(6, int(totalgf / 20))]
	mergedgf = [c for c in mergedgf if c.split('\t')[0] not in repeating and c.split('\t')[1] not in repeating ]
	"""
	
	with open(outpath + 'confident_genefusion.txt', 'w') as fo:
		# record confident gene fusion events
		fo.write("gene1\tgene2\tnumSupp\tchrom1\tbp1\tchrom2\tbp2\tID\treadNames\n")
		idnum = 1
		for event in clustered_candidates:
			fo.write(f"{event[0]}\t{event[1]}\t{event[2]}\t{event[3]}\t{event[4]}\t{event[5]}\t{event[6]}\tGF{idnum:03d}\t{event[8]}\n")
			idnum += 1

	return 0
