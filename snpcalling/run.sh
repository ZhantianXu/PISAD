#!/bin/bash
set -e

# 初始化变量
input_files=""
exp=""
kmer_size=""
prefix=""
tmp_dir=""
plot_dir=""
snp_dir=""
cutoff=""
rho=""
est=""
het=""
setleft=""
setright=""
thread=""

# 打印使用说明的函数
print_usage() {
    echo "Usage: $0 -i <input_files> -m <heterozygosity> [options...]"
    echo "Required parameters:"
    echo "  -i: Input files (space-separated *.fastq or *.fastq.gz files, no quotes needed)"
    echo "  -m: Heterozygosity parameter (0 for <1.2%, 1 otherwise)"
    echo "Optional parameters:"
    echo "  -k: kmer-size (default: 21)"
    echo "  -t: thread (default: 8)"
    echo "  -o: Output prefix (defaults: first input file's prefix)"
    echo "  -d1: Directory for dsk files (default: current directory)"
    echo "  -d2: Directory for output plot (default: current directory)"
    echo "  -d3: Directory for SNP output (default: current directory)"
    echo "  -h: Show this help message"
    echo "Advanced optional parameters:"
    echo "  -est: est_kmercov (default: Estimated by algorithm)"
    echo "  -cutoff: cutoff threshold (defaults: 0.95)"
    echo "  -het: Initial heterozygosity (defaults: 0/0.12)"
    echo "  -rho: Initial rho value (defaults: 0.2)"
    echo "  -setleft: Left boundary of the heterozygous region (defaults: Estimated by algorithm)"
    echo "  -setright: Right boundary of the heterozygous region (defaults: Estimated by algorithm)"
}

# 检查参数值的函数
check_value() {
    if [ $# -lt 2 ] || [[ "$2" =~ ^- ]] || [ -z "$2" ]; then
        echo "Error: $1 requires a value"
        print_usage
        exit 1
    fi
}

# 手动解析参数
while [ $# -gt 0 ]; do
    case "$1" in
        -m) check_value "$1" "$2"; exp="$2"; shift 2;;
        -i) 
            shift
            input_files=""
            if [ $# -eq 0 ] || [[ "$1" =~ ^- ]]; then
                echo "Error: -i requires at least one input file"
                print_usage
                exit 1
            fi
            while [ $# -gt 0 ] && [[ "$1" != -* ]]; do
                input_files="$input_files $1"
                shift
            done
            ;;
        -k) check_value "$1" "$2"; kmer_size="$2"; shift 2;;
        -t) check_value "$1" "$2"; thread="$2"; shift 2;;
        -o) check_value "$1" "$2"; prefix="$2"; shift 2;;
        -d1) check_value "$1" "$2"; tmp_dir="$2"; shift 2;;
        -d2) check_value "$1" "$2"; plot_dir="$2"; shift 2;;
        -d3) check_value "$1" "$2"; snp_dir="$2"; shift 2;;
        -est) check_value "$1" "$2"; est="$2"; shift 2;;
        -cutoff) check_value "$1" "$2"; cutoff="$2"; shift 2;;
        -het) check_value "$1" "$2"; het="$2"; shift 2;;
        -rho) check_value "$1" "$2"; rho="$2"; shift 2;;
        -setleft) check_value "$1" "$2"; setleft="$2"; shift 2;;
        -setright) check_value "$1" "$2"; setright="$2"; shift 2;;
        -h) print_usage; exit 0;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

# 其余代码保持不变
# 去除 input_files 前后的空格
input_files=$(echo "$input_files" | sed 's/^ *//;s/ *$//')

# 检查必须参数
if [ -z "$input_files" ] || [ -z "$exp" ]; then
    print_usage
    exit 1
fi

# 设置默认值
kmer_size=${kmer_size:-21}
tmp_dir=${tmp_dir:-"."}
snp_dir=${snp_dir:-"."}
thread=${thread:-8}

# 如果 prefix 未通过 -o 指定，则从第一个 input_files 中提取
if [ -z "$prefix" ]; then
    first_file=$(echo "$input_files" | awk '{print $1}')
    prefix=$(basename "$first_file" | sed 's/\(\.fastq\.gz\|\.fastq\|\.fq\.gz\|\.fq\)$//')
fi

# 处理输入文件（将空格分隔的路径转换为逗号分隔）
input_files=$(echo "$input_files" | tr ' ' ',' | sed 's/,$//')
echo "Input files: $input_files"
echo "Prefix: $prefix"

# dsk 处理
if [ ! -f "$tmp_dir/${prefix}.h5" ]; then 
    dsk -nb-cores "$thread" -max-memory 20000 -file "$input_files" -out-tmp "$tmp_dir" -out-dir "$tmp_dir" \
        -histo 1 -out "$tmp_dir/${prefix}" -kmer-size "$kmer_size"
fi

# kcov.py 处理
if ! ls "$tmp_dir/${prefix}"*.png >/dev/null 2>&1; then
    cmd="python kcov.py -i \"$tmp_dir/${prefix}.histo\" -m \"$exp\""
    [ -n "$kmer_size" ] && cmd="$cmd -k \"$kmer_size\""
    [ -n "$est" ] && cmd="$cmd -est \"$est\""
    [ -n "$cutoff" ] && cmd="$cmd -cutoff \"$cutoff\""
    [ -n "$het" ] && cmd="$cmd -het \"$het\""
    [ -n "$rho" ] && cmd="$cmd -rho \"$rho\""
    [ -n "$plot_dir" ] && cmd="$cmd -d2 \"$plot_dir\""
    eval "$cmd"
fi

# PNG 文件处理
png_files=$(ls "$tmp_dir/${prefix}"*.png 2>/dev/null)
if [ -n "$png_files" ]; then
    png_count=$(echo "$png_files" | wc -l)
    if [ "$png_count" -gt 2 ]; then
        echo "Error: More than 2 matching .png files found in $tmp_dir with prefix $prefix:"
        echo "$png_files"
        exit 1
    fi
    first_png=$(echo "$png_files" | head -n 1)
    left=$(basename "$first_png" .png | awk -F'_' '{print $(NF-1)}')
    right=$(basename "$first_png" .png | awk -F'_' '{print $NF}')
    [ -n "$setleft" ] && left="$setleft"
    [ -n "$setright" ] && right="$setleft"
fi


# SNP 处理
if [ ! -f "${snp_dir}/${prefix}_${kmer_size}_${left}_${right}_pairex.snp" ]; then
    command="./kmer2snp -i $tmp_dir/${prefix} -l $left -r $right -k $kmer_size -t $thread -o $snp_dir/"
    echo "$command"
    ./kmer2snp -i "$tmp_dir/${prefix}" -l "$left" -r "$right" -k "$kmer_size" -t "$thread" -o "$snp_dir/"
fi