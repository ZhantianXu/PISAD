#!/bin/bash
set -e

# less than 6 arguments
if [ $# -lt 6 ]; then
    echo "\$1:kmer-size "
    echo "\$2:Address where the dsk generated file is stored"
    echo "\$3:*.fq (fq or fq.gz file)"
    echo "\$4:Name of the generated file"
    echo "\$5:Address for storing the extracted snp file"
    echo "\$6:Heterozygosity parameter: if Species heterozygosity <1.5% choose 1,otherwise 2"
    exit
fi

# Output the basic parameters for reference
# echo "\$1:kmer-size \$2:tmp \$3:*.fq \$4:prefix \$5:prefix1 \$6:exp"
# echo "$1 $2 ${@:3:$#-5} ${@:(-3):1} ${@:(-2):1} ${@:(-1):1}"

# Set prefix, prefix1, and extra_param using the correct parameter positions
prefix=${@:(-3):1}
prefix1=${@:(-2):1}
exp=${@:(-1):1}

# Concatenate input files from the 3rd argument to the penultimate
input_files=""
for ((i=3; i<=($#-3); i++)); do
    input_files="$input_files ${!i}"
done

# Format the input_files by replacing spaces with commas and removing any trailing commas
input_files=$(echo $input_files | tr ' ' ',' | sed 's/,$//')

# Output the formatted input files
# echo $input_files

if [ ! -f "$2/${prefix}_chr_k$1.h5" ]; then 
    dsk -abundance-min 1 -nb-cores 10 -max-memory 50000 -file $input_files -out-tmp $2 -out-dir $2 -histo 1 -out $2/${prefix}_chr_k$1 -kmer-size $1
fi

if [ ! -f "$2/${prefix}_hete.peak.k$1" ]; then
    python kcov.py $2/${prefix}_chr_k$1.histo $2/${prefix}_hete.peak.k$1 
fi


if [ -f "$2/${prefix}_hete.peak.k$1" ]; then
    n=$(grep -oP '(?<=n:)[0-9]+' "$2/${prefix}_hete.peak.k$1")
    if [ "$n" -eq 2 ]; then
        kcov=$(awk 'NR==2 {print int($1)}' "$2/${prefix}_hete.peak.k$1")
        half_kcov=$(expr $kcov / 2)
        coverage=$(awk 'NR==3 {print int($1)}' "$2/${prefix}_hete.peak.k$1")
    elif [ "$n" -eq 1 ]; then
        kcov=$(awk 'NR==2 {print int($1 / 2)}' "$2/${prefix}_hete.peak.k$1")
        if [ "$exp" -eq 1 ]; then
            #低杂合度
            half_kcov=$(expr $kcov / 2)
            coverage=$(expr $kcov \* 2)
        else
            #高杂合度
            kcov=$(expr $kcov \* 2)
            half_kcov=$(expr $kcov / 2)
            coverage=$(expr $kcov \* 2)
        fi
    else
        echo "Error: n must be 1 or 2. Found: $n"
        exit
    fi

    if [ $kcov -gt 2 ]; then # >2
        left=$(expr $kcov - $half_kcov)
        right=$(expr $kcov + $half_kcov)
    elif [ $kcov -gt 1 ]; then
        left=$(expr $kcov - $half_kcov + 1)
        right=$(expr $kcov + $half_kcov)
    else
        left=1
        right=2
    fi

fi

# Use extra_param (which is the $6 argument)
if [ ! -f "${prefix1}/${prefix}_$1_$left_$right_pairex.snp" ]; then
    command="./kmer2snp --t1 $2/${prefix}_chr_k$1 --c1 $left --c2 $right --k $1 --t 8 --n $prefix1/$prefix"
    echo $command
    ./kmer2snp --t1 $2/${prefix}_chr_k$1 --c1 $left --c2 $right --k $1 --t 8 --n $prefix1/$prefix
fi
