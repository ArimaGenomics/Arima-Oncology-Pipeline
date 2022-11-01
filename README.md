![alt text](https://arimagenomics.com/wp-content/files/2021/08/Arima-Genomics-logo.png "Celebrating Science and Scientist")

# Arima Oncology Pipeline using Arima-HiC<sup>+</sup> Kit

This pipeline is for analyzing capture HiC data with HiCUP and CHiCAGO pipelines for generating Arima Genomics QC metrics and output data files. The pipeline runs the standard HiCUP, CHiCAGO and Juicer "pre" tools from https://www.bioinformatics.babraham.ac.uk/projects/hicup/, https://bioconductor.org/packages/release/bioc/html/Chicago.html and https://github.com/aidenlab/juicer/wiki/Pre with some minor changes to the default parameters. These parameters have been optimized by internal benchmarking and have been found to optimize sensitivity and specificity. It also runs the GenomeScan algorithm developed in house for inter- and intra-chromosomal Structural Variation (SV) detection from Capture HiC data. This pipeline will generate shallow and deep sequencing QC metrics which can be copied into the Arima QC Worksheet for analysis. Additionally, the pipeline automatically generates metaplots for data QC and arc plots of chromatin loops for data visualization.

This pipeline is suited for Capture HiC from any protocol but optimal results are obtained when using the Arima-HiC<sup>+</sup> kit with the Arima-CHiC protocol.

To order Arima-HiC<sup>+</sup> kits, please visit our website: https://arimagenomics.com/

## Getting Started

## <span style="color:red"> Important! </span>

### Although you can install all the tools and dependencies yourself, we strongly recommend using our pre-built Singularity image to run Arima Oncology pipeline!

#### We provide Docker (https://hub.docker.com/repository/docker/arimaxiang/arima_oncology) and Singularity containers (ftp://ftp-arimagenomics.sdsc.edu/pub/ARIMA_Oncology_Pipeline/Arima-Oncology-Pipeline-singularity-v0.3.sif) that allow running Arima Oncology pipeline out of the box. You can mount necessary input and output locations and run Arima Oncology pipeline without dealing with tedious installations or library dependencies.
#### Additionally, a wrapper script (run_Arima-Oncology-Pipeline-singularity-v0.3.sh) is provided in this repository for you to conveniently run the Singularity image by hiding all the mounting steps in a black box so that the only input needed from the users are the FASTQ files and an output location.

## Using the Singularity image with the wrapper script

#### Just download the *.sif file with the corresponding wrapper script.
```
wget ftp://ftp-arimagenomics.sdsc.edu/pub/ARIMA_Oncology_Pipeline/Arima-Oncology-Pipeline-singularity-v0.3.sif
wget https://raw.githubusercontent.com/ArimaGenomics/Arima-Oncology-Pipeline/main/run_Arima-Oncology-Pipeline-singularity-v0.3.sh
```

#### Put them in the same directory and run the command below.
> bash run_Arima-Oncology-Pipeline-singularity-v0.3.sh [-W run_hicup] [-Y run_bam2chicago] [-Z run_chicago] [-G run_genomescan]
                 [-P run_plot] [-C run_hicplot] [-I FASTQ_string] [-w scan_window] [-o out_dir] [-p output_prefix] [-t threads]

```
Required options:
    * [-I FASTQ_string]: a pair of FASTQ files separated by "," (no space is allowed)
    * [-o out_dir]     : output directory

Optional options:
    * [-W run_hicup]      : "1" (default) to run HiCUP pipeline, "0" to skip. If skipping, HiCUP_summary_report_*.txt and *R1_2*.hicup.bam need to be in the HiCUP output folder.
    * [-Y run_bam2chicago]: "1" (default) to run bam2chicago.sh, "0" to skip
    * [-Z run_chicago]    : "1" (default) to run CHiCAGO pipeline, "0" to skip
    * [-G run_genomescan] : "1" (default) to run GenomeScan pipeline, "0" to skip
    * [-P run_plot]       : "1" (default) to generate plots, "0" to skip
    * [-C run_hicplot]    : "1" (default) to generate HiC heatmap, "0" to skip
    * [-p output_prefix]  : output file prefix (filename only, not including the path)
    * [-w scan_window]    : sliding window size for GenomeScan in base pair. Default: 50000
    * [-t threads]        : number of threads to run HiCUP, CHiCAGO, GenomeScan and/or Juicer. Default: 16
    * [-h]                : print this help and exit
```


## Arima Specific Outputs

### Arima Shallow Sequencing QC

#### [output_directory]/[output_prefix]_Arima_QC_shallow.txt
Contents: This file includes QC metrics for assessing the shallow sequencing data for each CHiC library.
- Break down of the number of read pairs
- The target sequencing depth for deep sequencing
- The percentage of long-range cis interactions that overlap the probe regions

### Arima Deep Sequencing QC

#### [output_directory]/[output_prefix]_Arima_QC_deep.txt
Contents: This file includes QC metrics for assessing the deep sequencing data for each CHiC library.
- Break down of the number of read pairs
- The number of loops called
- The percentage of long-range cis interactions that overlap the probe regions

### Loops

#### [output_directory]/chicago/data/[output_prefix].cis_le_2Mb.bedpe
- Significant loops called by CHiCAGO

### Arima Arc Plots

#### [output_directory]/chicago/data/[output_prefix].cis_le_2Mb.arcplot.gz
- bgzipped arcplot file

#### [output_directory]/chicago/data/[output_prefix].cis_le_2Mb.arcplot.gz.tbi
- tabix index of the bgzipped arcplot file

These files can be viewed in the WashU EpiGenome Browser (http://epigenomegateway.wustl.edu/browser/). See the Arima-CHiC Analysis User Guide for more details.

### Metaplots

#### [output_directory]/chicago/data/[output_prefix].bigwig
- bigwig file of all mapped reads. This file can be used to view the sequencing coverage in a genome browser.

#### [output_directory]/chicago/data/[output_prefix].heatmap.pdf
- PDF of metaplot and heatmap of all mapped reads overlapping the probe regions used for loop calling. This file can be used to assess the signal to noise of the CHiC enrichment.

### Arima GenomeScan module with SV Calls and Plots

#### [output_directory]/genomescan/

### Arima HiC heatmap for visualization using Juicebox

#### [output_directory]/plots/[output_prefix]_inter_30.hic


## Test Dataset (hg38)

25 and 200 million Illumina paired-end sequencing datasets and their pipeline output files (based on hg38) can be downloaded from:

ftp://ftp-arimagenomics.sdsc.edu/pub/Oncology/


## Further Reading for Advanced Users

***All the files mentioned in this section are already incorporated into the Docker and Singularity containers!***

On the FTP ( ftp://ftp-arimagenomics.sdsc.edu/pub/ARIMA_Oncology_Pipeline/test_data/design/5kb_2Mb/ ), there are six pre-computed design files for hg38 at 5kb resolution. We found that after the binning, 5kb resolution provides the best replicate reproducibility and sensitivity. This folder contains three pre-computed design files (*.npb, *.poe and *.nbpb), one *.rmap, one *.baitmap, and one *.baitmap_4col.txt file.

The hg38 reference files can also be found on our FTP ( ftp://ftp-arimagenomics.sdsc.edu/pub/ARIMA_Capture_HiC_Settings/hg38/reference/ ), with one reference FASTA file (\*.fa), one HiCUP digest file for Arima's dual-enzyme chemistry and six bowtie2 index files (\*.bt2) in it.

If you would like to run the pipeline with "-C 1" (generate HiC heatmap), these two files are also required: ftp://ftp-arimagenomics.sdsc.edu/pub/Arima_FFPE/Juicer/hg38.chrom.sizes and ftp://ftp-arimagenomics.sdsc.edu/pub/Arima_FFPE/Juicer/hg38_GATC_GANTC.txt.


**Here is a detailed description of each file on our FTP.**

*hg38.fa*
- This is the genome reference file downloaded from UCSC table browser.

*Digest_hg38_Arima.txt*
- This is the HiCUP digest file for Arima's dual-enzyme chemistry, generated using the "hicup_digester" script from HiCUP pipeline.
- Reference: https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html#arima-protocol

*hg38.1.bt2, hg38.2.bt2, hg38.3.bt2, hg38.4.bt2, hg38.rev.1.bt2, hg38.rev.2.bt2*
- These six files are are bowtie2 index files, generated using "bowtie2-build" script from Bowtie 2 pipeline.
- Reference: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer

*.rmap*
- This is a tab-separated file of the format <chr> <start> <end> <numeric ID>, describing the restriction digest (or "virtual digest" if pooled fragments are used). These numeric IDs are referred to as "otherEndID" in CHiCAGO. All read fragments mapping outside of the digest coordinates will be disregarded.
- The file can be created using "create_baitmap_rmap.pl" script from CHiCAGO pipeline (with some pre-processing using custom scripts).
- Reference: https://bitbucket.org/chicagoTeam/chicago/src/master/chicagoTools/

*.baitmap*
- This is a tab-separated file of the format <chr> <start> <end> <numeric ID> <annotation>, listing the coordinates of the baited/captured restriction fragments (should be a subset of the fragments listed in *.rmap file), their numeric IDs (should match those listed in *.rmap file for the corresponding fragments) and their annotations (such as, for example, the names of baited promoters). The numeric IDs are referred to as "baitID" in CHiCAGO.
- The file can be created using "create_baitmap_rmap.pl" script from CHiCAGO pipeline (with some pre-processing using custom scripts).
- Reference: https://bitbucket.org/chicagoTeam/chicago/src/master/chicagoTools/

*.baitmap_4col.txt*
- This is the first 4 columns of the ".baitmap" file.

*.nbpb, .npb, .poe*
- These three files are called design files in CHiCAGO pipeline. They are ASCII files containing the following information:
- NPerBin file (.npb): <baitID> <Total no. valid restriction fragments in distance bin 1> ... <Total no. valid restriction fragments in distance bin N>, where the bins map within the "proximal" distance range from each bait (0; maxLBrownEst] and bin size is defined by the binsize parameter.
- NBaitsPerBin file (.nbpb): <otherEndID> <Total no. valid baits in distance bin 1> ... <Total no. valid baits in distance bin N>, where the bins map within the "proximal" distance range from each other end (0; maxLBrownEst] and bin size is defined by the binsize parameter.
- Proximal Other End (ProxOE) file (.poe): <baitID> <otherEndID> <absolute distance> for all combinations of baits and other ends that map within the "proximal" distance range from each other (0; maxLBrownEst].
- These three files can be created using "makeDesignFiles_py3.py" script from CHiCAGO pipeline (with *rmap and *.baitmap as the input).
- Reference: https://bitbucket.org/chicagoTeam/chicago/src/master/chicagoTools/

*hicup.conf*
- This is the configuration file for HiCUP pipeline. Detailed descriptions on how to tune those parameters can be found at: https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html#running-hicup

*chicago_settings_5kb.txt*
- This is the settings file for CHiCAGO pipeline. Detailed descriptions on how to tune those parameters can be found at: https://rdrr.io/bioc/Chicago/man/defaultSettings.html


## Arima Pipeline Version

0.3

## Support

For Arima customer support, please contact techsupport@arimagenomics.com

## Acknowledgments
#### Authors of HiCUP: https://www.bioinformatics.babraham.ac.uk/projects/hicup/

- Wingett S, et al. (2015) HiCUP: pipeline for mapping and processing Hi-C data F1000Research, 4:1310 (doi: 10.12688/f1000research.7334.1)

#### Authors of CHiCAGO: https://bioconductor.org/packages/release/bioc/html/Chicago.html

- Cairns J, Freire-Pritchett P, Wingett SW, Dimond A, Plagnol V, Zerbino D, Schoenfelder S, Javierre B, Osborne C, Fraser P, Spivakov M (2016). "CHiCAGO: Robust Detection of DNA Looping Interactions in Capture Hi-C data." Genome Biology, 17, 127.

#### Authors of Juicer "Pre": https://github.com/aidenlab/juicer/wiki/Pre

- Neva C. Durand, Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments." Cell Systems 3(1), 2016.
