#!/bin/bash

usage_Help="Usage: ${0##*/} [-W run_hicup] [-Y run_bam2chicago] [-Z run_chicago] [-G run_genomescan]
                 [-P run_plot] [-C run_hicplot] [-I FASTQ_string] [-w scan_window] [-o out_dir] [-p output_prefix] [-t threads]"
run_hicup_Help="* [-W run_hicup]      : \"1\" (default) to run HiCUP pipeline, \"0\" to skip. If skipping, HiCUP_summary_report_*.txt
                         and *R1_2*.hicup.bam need to be in the HiCUP output folder."
run_bam2chicago_Help="* [-Y run_bam2chicago]: \"1\" (default) to run bam2chicago.sh, \"0\" to skip"
run_chicago_Help="* [-Z run_chicago]    : \"1\" (default) to run CHiCAGO pipeline, \"0\" to skip"
run_genomescan_Help="* [-G run_genomescan] : \"1\" (default) to run GenomeScan pipeline, \"0\" to skip"
run_plot_Help="* [-P run_plot]       : \"1\" (default) to generate plots, \"0\" to skip"
run_hicplot_Help="* [-C run_hicplot]    : \"1\" (default) to generate HiC heatmap, \"0\" to skip"
FASTQ_string_Help="* [-I FASTQ_string]: a pair of FASTQ files separated by \",\" (no space is allowed)"
scan_window_Help="* [-w scan_window]    : sliding window size for GenomeScan in base pair. Default: 50000"
out_dir_Help="* [-o out_dir]     : output directory"
output_prefix_Help="* [-p output_prefix]  : output file prefix (filename only, not including the path)"
threads_Help="* [-t threads]        : number of threads to run HiCUP, CHiCAGO, GenomeScan and/or Juicer. Default: 16"
help_Help="* [-h]                : print this help and exit"

printHelpAndExit() {
    echo -e "$usage_Help"
    echo
    echo -e "Required options:"
    echo -e "$FASTQ_string_Help"
    echo -e "$out_dir_Help"
    echo
    echo -e "Optional options:"
    echo -e "$run_hicup_Help"
    echo -e "$run_bam2chicago_Help"
    echo -e "$run_chicago_Help"
    echo -e "$run_genomescan_Help"
    echo -e "$run_plot_Help"
    echo -e "$run_hicplot_Help"
    echo -e "$output_prefix_Help"
    echo -e "$scan_window_Help"
    echo -e "$threads_Help"
    echo -e "$help_Help"
    exit "$1"
}

run_hicup=1 # "1" to run HiCUP pipeline, "0" to skip
run_bam2chicago=1 # "1" to run bam2chicago.sh, "0" to skip
run_chicago=1 # "1" to run CHiCAGO pipeline, "0" to skip
run_genomescan=1 # "1" to run GenomeScan pipeline, "0" to skip
run_plot=1 # "1" to generate plots, "0" to skip
run_hicplot=1 # "1" to generate HiC heatmap, "0" to skip
scan_window=50000 # Sliding window size for GenomeScan in base pair.
threads=16 # Number of threads to run HiCUP, CHiCAGO, GenomeScan and/or Juicer.

while getopts "C:G:hI:o:P:p:t:W:w:Y:Z:" opt; do
    case $opt in
    h) printHelpAndExit 0;;
    W) run_hicup=$OPTARG ;;
    Y) run_bam2chicago=$OPTARG ;;
    Z) run_chicago=$OPTARG ;;
    G) run_genomescan=$OPTARG ;;
    P) run_plot=$OPTARG ;;
    C) run_hicplot=$OPTARG ;;
    I) FASTQ_string=$OPTARG ;;
    w) scan_window=$OPTARG ;;
    o) out_dir=$OPTARG ;;
    p) output_prefix=$OPTARG ;;
    t) threads=$OPTARG ;;
    [?]) printHelpAndExit 1;;
    esac
done

#hash singularity &> /dev/null
command -v singularity &> /dev/null
if [[ $? -ne 0 ]]; then
    echo "Could not find singularity. Please install or include it into the \"PATH\" variable!"
    printHelpAndExit 1
fi

if [ -z "$FASTQ_string" ]; then
    echo "Please provide a pair of FASTQ files separated by \",\" only (no space) (-I)!"
    printHelpAndExit 1
else
    IFS=',' read -a FASTQ <<< "$FASTQ_string"
    for i in "${FASTQ[@]}"; do
        if [ ! -f "$i" ]; then
            echo "$i does not exist (-I)!"
            printHelpAndExit 1
        fi
    done
fi

if [ -z "$output_prefix" ]; then
    output_prefix=$(basename ${FASTQ[0]} | sed 's/^\(.*\)[._]R1.*\.f.*q.*$/\1/')
fi

R1_basename=$(basename ${FASTQ[0]})
R2_basename=$(basename ${FASTQ[1]})
FASTQ_dir=$(dirname ${FASTQ[0]})

if [ -z "$out_dir" ]; then
    echo "Please provide an output directory (-o)!"
    printHelpAndExit 1
fi

# out_dir_abs=`readlink -f $out_dir`
[ -d "$out_dir" ] || mkdir -p $out_dir

cwd=$(dirname $0)

echo "Job is running in the Singularity container ..."
echo "Output: $out_dir"

singularity exec -B $FASTQ_dir:/Oncology/INPUT/ -B $out_dir:/Oncology/OUTPUT/ $cwd/Arima-Oncology-Pipeline-singularity-v0.2.sif bash /Oncology/Arima-Oncology-Pipeline-v0.2.sh -W $run_hicup -Y $run_bam2chicago -Z $run_chicago -G $run_genomescan -P $run_plot -C $run_hicplot -a /usr/bin/bowtie2 -H /HiCUP-0.8.0/ -c /chicagoTools/ -x /Oncology/reference_files/hg38 -d /Oncology/reference_files/Digest_hg38_Arima.txt -s /Oncology/reference_files/hg38.chrom.sizes -e /Oncology/reference_files/hg38_GATC_GANTC.txt -O hg38 -b /Oncology/utils/oncopanel_probes_v1.4.srt.bed -R /Oncology/Arima_files/design/5kb_2Mb/oncopanel_probes_v1.4_hg38_5kb.rmap -B /Oncology/Arima_files/design/5kb_2Mb/oncopanel_probes_v1.4_hg38_5kb.baitmap -D /Oncology/Arima_files/design/5kb_2Mb/ -w $scan_window -g /Oncology/utils/oncopanel_genes_v1.4.srt.bed -t $threads -I /Oncology/INPUT/$R1_basename,/Oncology/INPUT/$R2_basename -o /Oncology/OUTPUT/ -p $output_prefix

echo "Job was finished!"

# Example command using default arguments:
# bash run_Arima-Oncology-Pipeline-singularity-v0.2.sh -I FASTQ_R1.fastq,FASTQ_R2.fastq -o OUT_DIR
