#!/usr/bin/env python3
import argparse
from collections import defaultdict

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fusions', type=str, required=True, action='append', help='The file containing confident gene fusion events')
    parser.add_argument('-o', '--output', type=str, required=True, help='The file to write the merged gene fusion events')
    parser.add_argument('-s', '--sort', type=str, required=False, default="numSupp", help='Sort order: [numSupp, gene]')
    return parser.parse_args()

def mergeGeneFusion(fusion_files: list[str], output_file: str, sort: str = "numSupp") -> None:
    """
    Merge gene fusion events from a file.
    Args:
        fusion_file: The file containing confident gene fusion events
        output_file: The file to write the merged gene fusion events
        sort: Sort order: numSupp, gene
    Returns:
        None
    """
    sortkey = ["numSupp", "gene"]
    if sort not in sortkey:
        raise ValueError(f"Invalid sort order: {sort}. Valid sort orders are: {sortkey}")
    
    genepair2fusion = defaultdict(list)

    # read all fusion files
    fusion = [] # ID, gene1, gene2, numSupp, chrom1, bp1, chrom2, bp2, readNames
    for fusion_file in fusion_files:
        with open(fusion_file, 'r') as f:
            for line in f:
                if line.startswith("ID"):
                    continue
                ll = line.strip().split("\t")
                ll[3] = int(ll[3])
                fusion.append(ll)
                genepair2fusion[(line[1], line[2])].append(ll)

    merged_fusion = []
    # trying to merge fusion events
    for genepair, fusions in genepair2fusion.items():
        if len(fusions) == 1:
            merged_fusion.append(fusions[0])
        else:
            ignored = []
            for i in range(len(fusions)):
                if i in ignored:
                    continue
                f1 = fusions[i]
                for j in range(i + 1, len(fusions)):
                    f2 = fusion[j]
                    if f1[4] == f2[4] and f1[5] == f2[5] and f1[6] == f2[6] and f1[7] == f2[7]:
                        ignored.append(j)
                        f1[8] = ",".join(set(f1[8].split(",") + f2[8].split(",")))
                        f1[3] = len(set(f1[8].split(",")))
                merged_fusion.append(f1)

    # sort merged fusion events
    if sort == "numSupp":
        merged_fusion.sort(key=lambda x: x[3], reverse=True)
    elif sort == "gene":
        merged_fusion.sort(key=lambda x: (x[1], x[2]))

    # write merged fusion events to output file
    cot = 1
    with open(output_file, 'w') as f:
        f.write(f"ID\tgene1\tgene2\tnumSupp\tchrom1\tbp1\tchrom2\tbp2\treadNames\n")
        for fusion in merged_fusion:
            f.write(f"GF{cot:03d}\t{fusion[1]}\t{fusion[2]}\t{fusion[3]}\t{fusion[4]}\t{fusion[5]}\t{fusion[6]}\t{fusion[7]}\t{fusion[8]}\n")
            cot += 1

def main():
    args = parse_args()
    mergeGeneFusion(args.fusions, args.output, args.sort)

if __name__ == '__main__':
    main()