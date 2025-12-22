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

    # check header and determine the file type: clustered, confident, confident_refined
    """
    header:
    clustered: gene1, gene2, numSupp, chrom1, bp1, chrom2, bp2, mapq, readNames
    confident: gene1, gene2, numSupp, chrom1, bp1, chrom2, bp2, ID, readNames
    confident_refined: ID, gene1, gene2, numSupp, chrom1, bp1, chrom2, bp2, readNames, exonNum, exonType, exonCoordinates, fusionType
    """
    file_type = None
    header = None
    for fusion_file in fusion_files:
        with open(fusion_file, 'r') as f:
            line = f.readline()
            if header is None:
                header = line.strip()
                if line.startswith("ID"):
                    file_type = "confident_refined" if "exonNum" in line else "confident"
                else: # start with gene1
                    file_type = "clustered"
            else:
                if header != line.strip():
                    raise ValueError(f"Header mismatch in {fusion_file}, file type should be {file_type}")
    
    # set idx
    idx = {}
    for i, item in enumerate(header.strip().split("\t")):
        idx[item] = i

    # read all fusion files
    print(f"Reading fusion files: {fusion_files}")
    for fusion_file in fusion_files:
        with open(fusion_file, 'r') as f:
            for line in f:
                if line.startswith("ID") or line.startswith("gene1"):
                    continue
                ll = line.strip().split("\t")
                if (ll[idx["gene1"]], ll[idx["gene2"]]) not in genepair2fusion:
                    genepair2fusion[(ll[idx["gene1"]], ll[idx["gene2"]])] = []
                genepair2fusion[(ll[idx["gene1"]], ll[idx["gene2"]])].append(ll)

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
                    f2 = fusions[j]
                    if f1[idx["chrom1"]] == f2[idx["chrom1"]] and f1[idx["bp1"]] == f2[idx["bp1"]] and f1[idx["chrom2"]] == f2[idx["chrom2"]] and f1[idx["bp2"]] == f2[idx["bp2"]]:
                        ignored.append(j)
                        f1[idx["readNames"]] = ",".join(set(f1[idx["readNames"]].split(",") + f2[idx["readNames"]].split(",")))
                        f1[idx["numSupp"]] = str(len(set(f1[idx["readNames"]].split(","))))
                merged_fusion.append(f1)

    # sort merged fusion events
    if sort == "numSupp":
        merged_fusion.sort(key=lambda x: x[idx["numSupp"]], reverse=True)
    elif sort == "gene":
        merged_fusion.sort(key=lambda x: (x[idx["gene1"]], x[idx["gene2"]]))

    # write merged fusion events to output file
    cot = 1
    with open(output_file, 'w') as fo:
        fo.write(header + "\n")
        for fusion in merged_fusion:
            # rename the ID
            if file_type in ["confident_refined", "confident"]:
                fusion[idx["ID"]] = f"GF{cot:03d}"
            fo.write("\t".join(fusion) + "\n")
            cot += 1

def main():
    args = parse_args()
    mergeGeneFusion(args.fusions, args.output, args.sort)

if __name__ == '__main__':
    main()