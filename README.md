# ntsm - Nucleotide Sequence/Sample Matcher
## Summary

We developed PISAD, a tool designed to detect anomalies in cohort samples without requiring reference information. It is primarily divided into two stages. Stage 1: We select low-error data from the cohort and conduct reference-free SNP calling to construct a variant sketch. Stage 2: By comparing the k-mer counts of other cohort data to the variant sketch, we infer the relationships between the sample and other samples to detect the sample swap.

## Dependencies
recommend use conda to install
* GCC (Tested on 8.5.0)
* gperftools(2.10)
* hdf5(1.14.3)
* boost(1.85.0)
## Installation

cloning the PISAD repository to your machine and enter its directory.

```bash
 git clone https://github.com/ZhantianXu/PISAD.git 
 cd pisad/
```
Compiling should be as easy as:

```bash
./configure && make
```

To install in a specified directory:

```bash
./configure --prefix=/PATH && make install
```

## Usage

##### Stage1: SNP callng :

First, we select a low-error-rate sequencing dataset as the target sample for rapid SNP calling.

Example:

```bash
./run.sh 21 /h5 /data/hg002.fastq.gz hg002 /snp 1
```
```bash
    \$1:kmer-size 
    \$2:Address where the dsk generated file is stored
    \$3:*.fq (fq or fq.gz file)
    \$4:Name of the generated file
    \$5:Address for storing the extracted snp file
    \$6:Heterozygosity parameter: if Species heterozygosity <1.5% choose 1,otherwise 2
```

##### Stage1: construct variant sketch:
Next, we convert the called SNPs into a variant sketch.
```bash
./create /snp/hg002_21_2_4_pairex.snp /fa/hg002.fa 21 21
```
```bash
    \$1:snp file 
    \$2:kmer-size
    \$3:Filtering threshold (default: 21)
```

##### Stage2: count the k-mers:
Using this variant sketch we can then count all of these k-mers within a fastq file. Files may be gzipped and multiple threads can be used. Each sample needs a separate run of this command and its own count files.You need to run at least two counts: one for your target sample and one for the sample to be tested.
```bash
./ntsmCount -k 21 -t 2 -s /fa/hg002.fa -n eval/hg002 /data/hg002.fastq.gz
./ntsmCount -k 21 -t 2 -s /fa/hg002.fa -n eval/hg003 /data/hg003.fastq.gz
```
Here, the -s option allows inputting multiple FA files for variant sketching, separated by commas, such as `-s /fa/hg002.fa,/fa/hg001.fa`.
If your input file has a high coverage, you can also add the `-m` parameter to control the reading process and save time, such as `-m 10`.

##### Stage2:Evaluate the samples:
Input the statistics of your target sample and the sample to be tested to calculate their relationship and detect sample swaps.
```bash
./ntsmEval -a /eval/hg002_hg002.txt  /homeb/xuzt/coverage/eval/hg002_hg003.txt > summary.tsv
```



