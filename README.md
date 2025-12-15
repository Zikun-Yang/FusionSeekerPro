# FusionSeeker

A advanced gene fusion caller for long-read single-molecular sequencing data. It is modified from the [FusionSeeker](https://github.com/Maggi-Chen/FusionSeeker) and can better call the fusion events. See [Differences between FusionSeekerPro and FusionSeeker?](#difference) for details.

Authors: [Zikun Yang](https://github.com/Zikun-Yang), [Maggi Chen](https://github.com/Maggi-Chen) (author of FusionSeeker)

## <a name="toc"></a>Table of Contents

- [Quick Start](#quickstart)
- [Introduction](#intro)
- [Differences between FusionSeekerPro and FusionSeeker?](#difference)
- [Installation](#install)
- [Usage](#usage)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Getting Help](#help)
- [Citing FusionSeekerPro](#cite)

## <a name="quickstart"></a> Quick Start
```sh
git clone https://github.com/Zikun-Yang/FusionSeekerPro.git
cd FusionSeekerPro/
./fusionseekerpro -h

# quick gene fusion calling with stored genome
fusionseekerpro --bam merged.sort.bam --datatype isoseq --outpath fusionseekerpro_out/ --human38

# Gene fusion discovery with custom reference (e.g. T2T-CHM13v2.0)
fusionseekerpro --bam merged.sort.bam --datatype nanopore --outpath fusionseekerpro_out/ --ref Human_T2T-CHM13v2.0.fa --gtf Human_T2T-CHM13v2.0.gtf 
```

## <a name="intro"></a> Introduction

FusionSeekerPro is a tool for gene fusion discovery with long-read transcriptome sequencing data. The input should be a sorted BAM (PacBio Iso-Seq, Nanopore, or mixed platform). The output is a list of confident gene fusions and their transcript sequences. By default, FusionSeekerPro uses Human GRCh38 and Ensembl annotation v104 (Homo_sapiens.GRCh38.104.chr.gtf.gz) as reference.<br />
When using custom reference, make sure the chromosome name in BAM, GTF, and reference_genome.fa are identical. By default, FusionSeekerPro only considers gene with valid "gene_name" or "gene" in GTF and skips the remaining genes, unless --geneid is set.<br />

## <a name="difference"></a> Differences between FusionSeekerPro and FusionSeeker?

* reorganized code strcutures, descriptions of funstions and documents
* fixed some known bugs in FusionSeeker
* added bam preprocessing modules that only keep reads aligned to multiple genes. This can shorten the running time
* changed the multiprocessing manner, split the reads by genomic windows rather than chromosomes. This can

## <a name="install"></a> Installation

Dependencies for FusionSeekerPro:

* python3
* setproctitle
* [pysam](https://github.com/pysam-developers/pysam)  (tested with version 0.17.0)
* [minimap2](https://github.com/lh3/minimap2)  (tested with version 2.30)
* [samtools](https://github.com/samtools/samtools)  (tested with version 1.22.1)
* [bedtools](https://github.com/arq5x/bedtools2) (tested with version 2.31.1)
* [pybedtools](https://github.com/daler/pybedtools) (tested with version 0.12.0)
* [bsalign](https://github.com/ruanjue/bsalign)  (tested with version 1.2.1)

```sh
git clone https://github.com/Zikun-Yang/FusionSeekerPro.git
```
Then, please also add this directory to your PATH:
```sh
export PATH=$PWD/FusionSeekerPro/:$PATH
```


 You can use miniforge3 to simplify the environment setup:
```sh
mamba create --name fusions -y
mamba activate fusions
mamba install -c bioconda minimap2=2.24 pysam=0.17 samtools=1.9 setproctitle=1.3.6 bedtools=2.31.1 pybedtools=0.12.0 setproctitle -y
git clone https://github.com/ruanjue/bsalign.git
cd bsalign && make
export PATH=$PWD:$PATH
```

A test dataset is available to verify successful installation:
```
fusionseekerpro --bam testdata/test.bam  -o test_out/ --datatype isoseq --ref testdata/test.fa.gz
```
Output should be identical to confident_genefusion.txt and confident_genefusion_transcript_sequence.fa in the testdata folder, with 1 gene fusion and its transcript sequence. 

## <a name="usage"></a> Usage

```
fusionseekerpro [-h] --bam <sort.bam>

Gene fusion caller for long-read sequencing data

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --bam BAM             Input sorted BAM. index required
  --datatype DATATYPE   Input read type (isoseq, nanopore) [nanopore]
  --gtf GTF             Genome annotation file
  --ref REF             Reference genome. Required for breakpoint polishing
  --geneid              Use Gene ID instead of Gene name [False]
  --human38             Use reference genome and GTF for Human GCRh38 (default)
  --human19             Use reference genome and GTF for Human GCRh37
  -o OUTPATH, --outpath OUTPATH
                        Output directory [./fusionseekerpro_out/]
  -s MINSUPP, --minsupp MINSUPP
                        Minimal reads supporting an event [auto]
  --maxdistance MAXDISTANCE
                        Maximal distance to cluster raw signals [20 for isoseq, 40 for nanopore]
  --keepfile            Keep intermediate files [False]
  --thread THREAD       Number of threads [8]


```

FusionSeekerPro requires a input of read alignment results in BAM format sorted by coordinates. If you start with sequencing reads (Fasta or Fastq format), you may use minimap2 and samtools to map them to a reference genome before you can apply FusionSeekerPro:
```sh
# PacBio Iso-Seq
minimap2 -ax splice:hq reference.fa  isoseq.fastq | samtools sort -o isoseq.bam
samtools index isoseq.bam
# Nanopore
minimap2 -ax splice reference.fa  nanopore.fastq | samtools sort -o nanopore.bam
samtools index nanopore.bam
```

FusionSeekerPro can be applied with built-in Human reference genome (hg38) and annotation (Ensembl v104):
```
fusionseekerpro --bam isoseq.bam --datatype isoseq -o fusionseekerpro_out/
```
Or with custom reference genome and annotation (Make sure the chromosome name in both files are identical to those in BAM file):
```
fusionseekerpro --bam nanopore.bam  --datatype nanopore  -o fusionseekerpro_out/ --gtf annotation.gtf --ref reference.fa
```

By default, FusionSeekerPro uses only gene records with valid "gene_name" in the GTF file. To include all genes in the GTF file, use Gene ID instead:
```
fusionseekerpro --bam isoseq.bam --datatype isoseq -o fusionseekerpro_out/ --geneid 
```


### <a name="parameters"></a> Parameters
#### 1. --minsupp, minimal number of supporting reads
--min_supp is the most important argument for gene fusion candidate filtering of FusionSeekerPro. It is used to remove false-positive signals generated during sequencing or read alignment processes.
By default, FusionSeekerPro estimates the volumn of noise signals from input dataset to assign a resonable --minsupp. If you find number of gene fusions is too few under default settings, you can speficy a lower --minsupp cutoff to allow in more candidates:
```
fusionseekerpro --bam isoseq.bam --datatype isoseq -o test_out/ --minsupp 5 
```
It is not suggested to set a --minsupp below 3, unless the sequencing depth of input dataset is extremely low.

#### 2. --maxdistance, maxminal distance cutoff in density-based spatial clustering of applications with noise
This option adjusts max distance cutoff used for clustering gene fusion raw signals. By default, FusionSeekerPro sets 20 for highly accurate reads (IsoSeq) and 40 for noisy reads (Nanopore).
you can set a larger value of --maxdistance to tolerate more shifts in the breakpoint positions of raw signals:

```
fusionseekerpro --bam isoseq.bam --datatype isoseq -o test_out/ --maxdistance 100
```

#### 3. --ref, reference genome
Input reference genome allows FusionSeekerPro to align transcript sequences and refine breakpoint positions of confident gene fusion calls. Make sure to provide the same reference file used for read alignment.
By default, FusionSeekerPro does NOT refine breakpoint positions when no reference genome is provided. 
(minimap2(>=2.24) is required to map transcript sequences to the reference.)
```
fusionseekerpro --bam isoseq.bam --datatype isoseq -o test_out/  --ref reference.fa
```


## <a name="outputs"></a> Outputs
The output directory includes:
```
# Final results:
confident_genefusion_refined.txt               A list of refined confident gene fusion calls after mapping back to reference genome.
confident_genefusion.txt                       A list of confident gene fusion calls from input BAM file. Includes gene names, breakpoint positions, number and name of fusion-supporting reads.
confident_genefusion_transcript_sequence.fa    Transcript sequences of reported confident gene fusion 

# Intermediate results:
filtered.bam                                   Reads that are overlapped with more than one genes
filtered.bam.bai                               Index of bam file
clustered_candidate.txt                        A full list of gene fusion candidates before applying 
rawsignal.txt                                  A list of all gene fusion raw signals.
log.txt                                        Log file for debug.
(raw_signal/                                   Intermediate files during raw signal detection. Removed by default.)
(align_workspace/                              Intermediate files during transcript sequence          generation with bsalign poa. Removed by default.)
```

## <a name="help"></a>Getting Help

If you have further questions, want to report a bug, or suggest a new feature, please raise an issue at the [issue page](https://github.com/zikun-yang/FusionSeekerPro/issues).

## <a name="cite"></a>Citating FusionSeekerPro

If you use FusionSeekerPro in your work, please cite:
> To be updated

