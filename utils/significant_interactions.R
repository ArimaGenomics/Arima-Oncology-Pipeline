## Description: This script identifies windows with statistically significant interactions ##
## with the gene of interest. The script loops over all genes in the probe ##
## design and outputs a file with information about the window and the corresponding ##
## adjusted p-value. The script currently focuses on inter-chromosomal interactions. ##

## Inputs for script:
## input directory with all results per gene
## output directory
## sample name
## analysis resolution
## method for p value adjustment
## file with windows and corresponding genes
## code directory

require(data.table)
require(ggplot2)
#quartz(width = 10.5, height = 2.5)

# Read inputs
args = commandArgs(trailingOnly = TRUE)
inputDir <- args[1]
outDir <- args[1]
currentDir <- args[2]
sample.id <- args[3]
resolution <- args[4]
padj.method <- args[5]
windows.genes.file <- args[6]
codeDir <- args[7]
analysis_type <- args[8]
distance <- args[9]

# Debugging
#inputDir
#outDir
#sample.id
#resolution
#padj.method
#codeDir
#analysis_type
#distance

#genes[match(c("", ""), genes$gene),]
#all_genes[match(c("", ""), all_genes$gene),]

# Load datasets
load(paste(codeDir,"/oncopanel_genes_v1.4.Rdata", sep = ""))
load(paste(codeDir, "/all_genes.Rdata", sep = ""))

# Process gene dataframe for subsequent analysis
genes$mid_gene <- genes$chrStart + (genes$chrEnd - genes$chrStart)/2
# str(genes)

# Load scripts with additional functions for analysis
source(paste(codeDir, "/cHiC_plots.R", sep =""), chdir = FALSE)
source(paste(codeDir, "/cHiC_fishertest.R", sep =""), chdir = FALSE)
source(paste(codeDir, "/cHiC_analysis.R", sep =""), chdir = FALSE)
#quartz(width = 10.5, height = 2.5)

plot.dir <-  paste(inputDir, "/", sample.id, "/plots/", resolution,  sep = "")

# Process window file with gene information: create file with a comma separated gene list of genes for every window
windows.genes <- read.table(windows.genes.file, header = FALSE, sep = "\t",stringsAsFactors= FALSE)
colnames(windows.genes) <- c("chr","chrStart","chrEnd","gene")
windows_genes_list <- aggregate(gene ~ chr + chrStart + chrEnd, data = windows.genes, paste, collapse = ",")


### Identify significant events for every gene in the panel ###

# Identify all genes in the panel
#all_genes_files <- list.files(paste(inputDir, sample.id, "gene_results", resolution, sep = "/"))

all_genes_files <- list.files(paste(inputDir, sample.id, "gene_results", resolution, currentDir, sep = "/"))
#all_genes_files <- list.files(paste(inputDir, sample.id, "gene_results", resolution, sep = "/"))
#all_genes_files <- list.files(inputDir)
# str(all_genes_files)
# Create dataframe to collectivelty store results from all genes
# The dataframe is initialized with the window coordinates
all_genes_results <- data.frame(chr = character(), chrStart = numeric(), chrEnd = numeric(), abs.pos = numeric(), stringsAsFactors = FALSE)
#colnames(all_genes_results) <- c("chr","chrStart","chrEnd", "abs.pos")

# Create dataframe to store significant results for all genes
all_genes_results_sig <- data.frame(chr = character(), chrStart = numeric(), chrEnd = numeric(), abs.pos = numeric(), stringsAsFactors = FALSE)

# Loop over all genes and perform analysis

# Calculate values that are the same for each gene
# and should be computed only once

# This data frame has one column for each gene and one row for each window.
# The values are the number of reads for each gene and each window.

print("Creating dataframe with reads per gene for every window genome-wide")
df <- getCountMtrx(inputdir = paste(inputDir, sample.id, "gene_results", resolution,"all_genes", sep = "/"))
# str(df)
# length(all_genes_files)

for (i in seq(1,length(all_genes_files),1)){
#	print(paste("I am inside for loop",i, sep = ""))

	focal.gene.file <- all_genes_files[i]
	# The gene name is located before the first "." in the file name
	focal.gene <- strsplit(focal.gene.file, "[.]")[[1]][1]
	# Print gene information
	print(paste("Processing", focal.gene, "reads", sep = " "))

	# Identify gene chromosome - useful for intra-SV analysis
	focal.gene.chr <- genes$chr[genes$gene == focal.gene]
	# focal.gene.chr

	# Calculate p-values for all windows and write results
	if (analysis_type == "inter"){
#		print("I am in the inter padj calculation")
		pvals <- p.adjust(inter.pvals.gene(df, focal.gene), method = padj.method)
	}
	# For intra-chromosomal SVs set parameter to FALSE and provide distance information
	if (analysis_type == "intra") {
		pvals <- p.adjust(inter.pvals.gene(df, focal.gene,FALSE, distance), method = padj.method)
#		print("I am in the intra padj calculation")
	}
#}

	gene_pvals <- cbind(df[,1:4], df[,focal.gene],as.data.frame(pvals))
	colnames(gene_pvals)[5] <- "numberOfReads"
	#colnames(gene_pvals)[5] <- focal.gene
	#colnames(gene_pvals)[6] <- paste("padj", focal.gene, sep = "_")
	colnames(gene_pvals)[6] <- "padj"
#
#	print("all good after padj \n")

	# Ouput file needs to have window coords and adjusted pvalue for each gene
	# Attempt to additionally correct the pvalue for the number of genes in the probe design
#	padj_for_number_of_genes <- 0.05/1404
#	padj_for_number_of_genes
#	gene_pvals_sig <- subset(gene_pvals, gene_pvals$padj <= padj_for_number_of_genes)
	# Previous version
	gene_pvals_sig <- subset(gene_pvals, gene_pvals$padj <= 0.05)

	#gene_pvals_marginal <- subset(gene_pvals, (gene_pvals$padj > 0.05)&&(genes_pvals$padj <= 0.1))
#	print("all good after subsetting \n")
	# Write file with all windows and corresponding p-values
	if (analysis_type == "inter") {
		write.table(gene_pvals, paste(inputDir, "/", sample.id, "/significant_interactions/", resolution,"/all_windows.", focal.gene, ".bed",  sep = ""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep ="\t")
	}else{
		write.table(gene_pvals, paste(inputDir, "/", sample.id, "/significant_interactions/", resolution,"/all_windows.", focal.gene, ".", distance, "_intra.bed",  sep = ""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep ="\t")
	}

	# Write windows with significat pvalues (padj <=0.05) and marginal windows (0.05 <= padj <= 0.1) (TO DO)
	# if there are significant interactions and create genome-wide manhattan plots.
#	print("all good before writing significant file")
#if (FALSE){
	if (nrow(gene_pvals_sig) != 0){
		# Write significant interactions #
		#write.table(gene_pvals_sig, paste(inputDir, "/", sample.id, "/significant_interactions/", resolution,"/significant_events.", focal.gene ,".bed",  sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
		#print("all good after writing significant interactions")

		# Create manhattan plots #
		# Inter-chromosomal
		if (analysis_type == "inter"){
			write.table(gene_pvals_sig, paste(inputDir, "/", sample.id, "/significant_interactions/", resolution,"/significant_events.", focal.gene ,".bed",  sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
			plot_genome <- genomeManhattan(df, gene = focal.gene, pvals,
						plot.title = paste(sample.id, ": ", focal.gene, "\n(", resolution, " sliding windows)", sep = ""),
						calc.pval=TRUE, plot.min = 1, sig.col = "violetred2", annot.col = rgb(65, 144, 243, maxColorValue = 255), pnt.col = "grey20")

			ggsave(plot = plot_genome, filename = paste(sample.id, focal.gene, resolution, "pdf", sep = "."), path = paste(plot.dir, "/pdf/", sep = ""), width = 10.5, height = 2.5, units = "in")
			ggsave(plot = plot_genome, filename = paste(sample.id, focal.gene, resolution, "png", sep = "."), path = paste(plot.dir, "/png/", sep = ""), width = 10.5, height = 2.5, units = "in")
		}
	#	print("all good after plotting\n")

		# Intra-chromosomal
		if (analysis_type == "intra"){
			write.table(gene_pvals_sig, paste(inputDir, "/", sample.id, "/significant_interactions/", resolution,"/significant_events.intra.", focal.gene ,".bed",  sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
			plot_chromosome <- chrManhattan(df, gene = focal.gene, chromosome = focal.gene.chr , distance, pvals,
                                                plot.title = paste(sample.id, ": ", focal.gene, "\n(", resolution, " sliding windows)", sep = ""),
                                                calc.pval=FALSE, plot.min = 1, sig.col = "violetred2", annot.col = rgb(65, 144, 243, maxColorValue = 255), pnt.col = "grey20")
			ggsave(plot = plot_chromosome, filename = paste(sample.id, focal.gene, resolution, distance, "intra","pdf", sep = "."), path = paste(plot.dir, "/pdf/", sep = ""), width = 12.5, height = 2.5, units = "in")
			ggsave(plot = plot_chromosome, filename = paste(sample.id, focal.gene, resolution, distance, "intra", "png", sep = "."), path = paste(plot.dir, "/png/", sep = ""), width = 12.5, height = 2.5, units = "in")
		}

		# Identify interacting partner genes #
		windows_genes_sig <- merge(gene_pvals_sig,windows_genes_list,by = c("chr","chrStart","chrEnd"), all.x = TRUE)
		colnames(windows_genes_sig)[7] <- "partnerGene"
		# Add focal gene coordinates to file
		focal.gene.coords <- genes[genes$gene == focal.gene,c("chr","chrStart","chrEnd")]
		windows_genes_sig$focalGene <- focal.gene
		windows_genes_sig$chrFocalGene <- focal.gene.coords$chr
		windows_genes_sig$startFocalGene <- focal.gene.coords$chrStart
		windows_genes_sig$endFocalGene <- focal.gene.coords$chrEnd
		# Re-order data frame before writing
		windows_genes_sig <- windows_genes_sig[,c("chr","chrStart","chrEnd","chrFocalGene","startFocalGene","endFocalGene","numberOfReads","padj","focalGene","partnerGene", "abs.pos")]
		#windows_genes_sig <- windows_genes_sig[,c("chr","chrStart","chrEnd","chrFocalGene","startFocalGene","endFocalGene","abs.pos","focalGene","partnerGene","numberOfReads","padj")]
		if (analysis_type == "inter") {
			write.table(windows_genes_sig, paste(inputDir, "/", sample.id, "/significant_interactions/", resolution,"/significant_events.", focal.gene, "_annot.bed",  sep = ""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
		}else{
			write.table(windows_genes_sig, paste(inputDir, "/", sample.id, "/significant_interactions/", resolution,"/significant_events.", focal.gene,".",distance, "_intra_annot.bed",  sep = ""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
		}
	#	print("all good after writing annotated results")
	}

	# Create file with p-values for all windows and all genes
 	# and file with significant interactions for all genes together
#	if (i == 1) {
#		all_genes_results <- gene_pvals[,1:4]
#		all_genes_results_sig <- gene_pvals_sig[,1:4]
#	}
#	all_genes_results <- merge(all_genes_results, gene_pvals, by = c("chr","chrStart","chrEnd", "abs.pos"),all = TRUE)
#
#	if (nrow(gene_pvals_sig) != 0){
#		all_genes_results_sig <- merge(all_genes_results_sig, gene_pvals_sig, by = c("chr","chrStart","chrEnd", "abs.pos"),all = TRUE)
#	}

#}
}
