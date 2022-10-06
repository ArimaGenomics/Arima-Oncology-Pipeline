#!/bin/bash

print_usage() {
  printf "Usage: genomescan.sh -i [input.bam] [-g GENE_LIST.bed -w window_size -wf window_file -t threads -p prefix -o output_dir]\n";
}

if [[ "$1" =~ ^((-{1,2})([Hh]$|[Hh][Ee][Ll][Pp])|)$ ]]; then
    print_usage; exit 1
else
    while [[ $# -gt 0 ]]; do
        opt="$1"
        shift;
        current_arg="$1"
        if [[ "$current_arg" =~ ^-{1,2}.* ]]; then
            echo "WARNING: You may have left an argument blank. Double check your command."
        fi
        case "$opt" in
            "-i"|"--input"      ) INPUT_BAM="$1"; shift;;
            "-g"|"--genes"      ) GENE_LIST="$1"; shift;;
            "-w"|"--windowsize" ) WINDOW="$1"; shift;;
            "-wf"|"--windowfile") WINDOW_FILE="$1"; shift;;
            "-t"|"--threads"    ) THREADS="$1"; shift;;
            "-p"|"--prefix"     ) PREFIX="$1"; shift;;
            "-o"|"--out_dir"    ) OUT_DIR="$1"; shift;;
            "-c"|"--cleanup"    ) CLEANUP="$1"; shift;;
            *                   ) echo "ERROR: Invalid option: \""$opt"\"" >&2
                                exit 1;;
                esac
            done
        fi

: "${INPUT_BAM:?"No input file specified."}"
echo "Using input file: "$INPUT_BAM""

: "${GENE_LIST:?"No gene list file specified."}"
echo "Gene regions to analyze from: "$GENE_LIST""

#: "${WINDOW:="50000"}"
echo "Sliding window file to be used: "$WINDOW_FILE""

DEFAULT_PREFIX=`echo $INPUT_BAM | awk -F/ '{print $NF}' | awk -F. '{print $1}'`
: "${PREFIX:=$DEFAULT_PREFIX}"
DEFAULT_OUTPUT=`pwd`
: "${OUT_DIR/${PREFIX}_${WINDOW}/:=${DEFAULT_OUTPUT}/${PREFIX}_${WINDOW}/}"

: "${CLEANUP:=1}"
echo "Intermediate files will be deleted after run. Use '-c 0' to save files."

# : "${WARN:=0}"
# echo "Run will abort in case of file error during scan. Use -a 1 to enable warnings."

# gene_list = bed file (header ok, need line to check for header ^""), check file format (4 column: chr, start, end, gene))
#  -> allow custom bed file for gene gene_list
#  -> default to search Arima oncopanel gene list
#  -> allow gene name lookup against all_gene list (e.g. comma-delim list of all genes for grep)
#  -> allow coord lookup (general, format: chr1:1000-20000) for custom region focus

# input chromosome window files
#   -> preselect from existing sliding window files.
#   -> select custom window file
#   -> (create window file from <genome> with specified characteristics)

bamfilename="$(basename $INPUT_BAM)"
# echo $bamfilename

# Check if sorted file exists, create it not
INPUT_BAM_SRT=$OUT_DIR/$PREFIX/${bamfilename%.bam}.sortn.bam
[ -f ${INPUT_BAM_SRT} ] || samtools sort $INPUT_BAM -@ $THREADS -o $INPUT_BAM_SRT

#check if .bai index file exists, create if not
#[ -f ${INPUT_BAM}.bai ] || samtools index ${INPUT_BAM}
[ -f ${INPUT_BAM_SRT}.bai ] || samtools index ${INPUT_BAM_SRT}

#check if bamtobed file exists, and create if not:
IN_DIR=`echo $INPUT_BAM_SRT | awk '{gsub(/[^\/]*$/,"");print}'`
[ -f ${IN_DIR}/${PREFIX}.bed ] || bedtools bamtobed -i $INPUT_BAM_SRT > ${IN_DIR}/${PREFIX}.bed

# Added by Xiang to reduce memory usage in the grep step
sed s@/[12]@@ ${IN_DIR}/${PREFIX}.bed > ${IN_DIR}/${PREFIX}.modified.bed

#Create output directory tree
mkdir -p $OUT_DIR/$PREFIX/
mkdir -p $OUT_DIR/$PREFIX/tmp/
mkdir -p $OUT_DIR/$PREFIX/gene_results/${WINDOW}/
mkdir -p $OUT_DIR/$PREFIX/gene_results/${WINDOW}/all_genes/

#Add directories for subsequent analysis
mkdir -p $OUT_DIR/$PREFIX/significant_interactions/${WINDOW}/
mkdir -p $OUT_DIR/$PREFIX/significant_interactions/${WINDOW}/chromosome_interactions/
mkdir -p $OUT_DIR/$PREFIX/plots/${WINDOW}/
mkdir -p $OUT_DIR/$PREFIX/plots/${WINDOW}/pdf/
mkdir -p $OUT_DIR/$PREFIX/plots/${WINDOW}/png/

#Add counter variable to monitor progress
counter=0
numberOfGenes=$(cat $GENE_LIST | wc -l | cut -d' ' -f1)
#numberOfGenes=$(sed 1d $GENE_LIST | wc -l | cut -d' ' -f1) # if the file has a header
tenPercentOfGenes=$(($numberOfGenes / 10))

#sed 1d $GENE_LIST | while IFS=$'\t' read -r chrom chromStart chromEnd gene
cat $GENE_LIST | while IFS=$'\t' read -r chrom chromStart chromEnd gene
do
    coord="$chrom:$chromStart-$chromEnd"

    #Get reads aligned to gene region
    samtools view $INPUT_BAM_SRT $coord | cut -f1 > ${OUT_DIR}/${PREFIX}/tmp/${PREFIX}.${gene}.ids.txt

    #Check that id files are not empty
    # if [ -s ${OUT_DIR}/${PREFIX}/tmp/${PREFIX}.${gene}.ids.txt ] then
    #     echo "Gene ids saved for gene $gene."
    # else
    #     echo "No gene ids found. Is the input bam indexed?"

    #     if [[ "$WARN" =~ 0 ]] then
    #         exit 1
    #     fi
    # fi

    #Get read pairs from regions, save to .bed file for intersect
    #LC_ALL=C grep -w -F -f ${OUT_DIR}/${PREFIX}/tmp/${PREFIX}.${gene}.ids.txt < ${IN_DIR}/${PREFIX}.bed > ${OUT_DIR}/${PREFIX}/tmp/${PREFIX}.${gene}.bed

    # Added by Xiang to reduce memory usage
    awk 'NR==FNR{ids[$0]; next} {if($4 in ids) print} f' ${OUT_DIR}/${PREFIX}/tmp/${PREFIX}.${gene}.ids.txt ${IN_DIR}/${PREFIX}.modified.bed > ${OUT_DIR}/${PREFIX}/tmp/${PREFIX}.${gene}.bed

    #Get counts of reads aligned to gene within each sliding window
    bedtools intersect -c -a $WINDOW_FILE -b ${OUT_DIR}/${PREFIX}/tmp/${PREFIX}.${gene}.bed > ${OUT_DIR}/${PREFIX}/gene_results/${WINDOW}/all_genes/${gene}.${WINDOW}.bed

	# Increase gene counter after processing each gene
	counter=$(( $counter + 1 ))
	# Print progress bar after processing every 10% of genes approximately
        if [ $(( $counter % $tenPercentOfGenes )) -eq 0 ] ; then
            echo "Processed $(( ($counter*100) / $numberOfGenes + 1)) % of genes"
        fi

done

# if [[ $CLEANUP =~ 1 ]]; then
#     echo "Deleting intermediate files...";
#     rm -r ${OUT_DIR}/${PREFIX}/tmp/;
#     echo "Deleted ${OUT_DIR}/${PREFIX}/tmp/";
# fi
