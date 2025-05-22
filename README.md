# PISAD - Phsaed Intraspecies Sample Anomalies Detection tool

## Summary

We developed PISAD, a tool designed to detect anomalies in cohort samples without requiring reference information. It is primarily divided into two stages. Stage 1: We select low-error data from the cohort and conduct reference-free SNP calling to construct a variant sketch. Stage 2: By comparing the k-mer counts of other cohort data to the variant sketch, we infer the relationships between the sample and other samples to detect the sample swap.

## Dependencies

- GCC (Tested on 8.5.0)
- gperftools(2.10)
- hdf5(1.14.3)
- boost(1.85.0)

## Installation

pisad can be installed with conda using the command:
```bash
conda install bioconda::pisad
```

If cloning directly from the repository run:
```bash
./autogen.sh
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

First, we select a low-error-rate sequencing dataset as the target sample for rapid SNP calling. It supports multi-threaded processing.

Example:

```bash
run.sh -i /data/hg002.fastq.gz -m 0
```

```text
    Required parameters:
      -i:        Input files ( *.fastq or *.fastq.gz files)
      -m:        Heterozygosity parameter (0 for <1.2%, 1 otherwise)
    Optional parameters:
      -k:        kmer-size (default: 21)
      -t:        thread (default: 8)
      -o:        Output prefix (defaults: first input file's prefix)
      -d1:       Directory for dsk files (default: current directory)
      -d2:       Directory for output plot (default: current directory)
      -d3:       Directory for SNP output (default: current directory)
      -h:        Show this help message
    Advanced optional parameters:
      -est:      est_kmercov (default: Estimated by algorithm)
      -cutoff:   cutoff threshold (defaults: 0.95)
      -het:      Initial heterozygosity (defaults: 0/0.12)
      -rho:      Initial rho value (defaults: 0.2)
      -setleft:  Left boundary of the heterozygous region (defaults: Estimated by algorithm)
      -setright: Right boundary of the heterozygous region (defaults: Estimated by algorithm)
```

##### Stage1: construct variant sketch:

Next, we convert the called SNPs into a variant sketch.

```bash
create -i /snp/hg002_21_2_4_pairex.snp
```

```text
    Required parameters:
      -i:        Input files ( .snp file)
    Optional parameters:
      -k:        kmer-size (default: 21)
      -l:        Filtering threshold (default: 21)
      -o:        Output prefix (defaults: current directory)
```

##### Stage2: count the k-mers:

we compare the k-mer counts of other cohort samples to the variant sketch to infer relationships between them. Files may be gzipped and multiple threads can be used.

```bash
pisadCount -s /fa/hg002.fa /data/hg003.fastq.gz
```

```text
    Usage: ./pisadCount -s [FASTA] [OPTION]... [FILES...]
    Required options:
        -s:         variant sketch (one or more)
    Optional options:
        -t:      Number of threads to run (default: Allocate 6 threads for each sequencing file)
        -m:      k-mer coverage threshold for early termination (default: inf)
        -i:      extra debug information
        -k:      k-mer size used (default: 21)
        -o:      Evaluation file path (defaults: current directory)
        -h:      Display this dialog

```

Here, the -s option allows inputting multiple FA files for variant sketching, separated by commas, such as `-s /fa/hg002.fa,/fa/hg001.fa`.
If your input file has a high coverage, you can also add the `-m` parameter to control the reading process and save time, such as `-m 5`.

##### Stage2:Evaluate the samples:

Input the statistics of samples to calculate their relationship and detect sample swaps.

```bash
pisadEval /homeb/xuzt/coverage/eval/hg002_hg003.txt > summary.tsv
```

```text
    Usage: ./pisadEval [OPTION]... [FILES...]
    Optional options:
        -t:      Number of threads to run(default: 1)
        -h:      Display this dialog

```
