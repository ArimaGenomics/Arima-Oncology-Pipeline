#! /bin/bash
########################
# Scheduler directives #
########################
#PBS -q hotel
### Set the name of the job, where jobname is a unique name for your job
#PBS -N GenomeScan_MCF7_S796
#PBS -S /bin/bash
### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=10:00:00
### Inform the scheduler of the number of CPU cores for your job.
#PBS -l nodes=1:ppn=2
### Merge the output and error streams
#PBS -j oe /oasis/tscc/scratch/sofia/GenomeScan/genomescan_MCF7_S740.log
#PBS -M XXXXXX@arimagenomics.com
#PBS -m ae
#################
# Job Execution #
#################

module load R

echo $resolution
echo $outDir
echo $sampleName

#resolution=50kb
#sampleName=MCF7_S740
#padj_method=BH
#outDir=/oasis/tscc/scratch/sofia/GenomeScan
#codeDir=/oasis/tscc/scratch/sofia/GenomeScan
#oncopanelFile=//

echo "Calculate significant interaction events for $sampleName"

# For each gene, identify windows with significant interactions and partner genes
#windowsGenes=$outDir/windows_50kb_genes.bed
Rscript $codeDir/significant_interactions.R  $outDir $currentDir $sampleName $resolution $padj_method $windowsGenes $codeDir

#Rscript $codeDir/significant_interactions.R  $outDir $sampleName $resolution $padj_method $windowsGenes $oncopanelFile $codeDir
