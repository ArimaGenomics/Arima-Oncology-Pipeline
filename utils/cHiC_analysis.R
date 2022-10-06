require(data.table)

getCountMtrx <- function(inputdir, geneSet=1:length(list.files())){
	setwd(inputdir)
	sampleFiles <- list.files()[geneSet]

	if (length(sampleFiles) != 0){
	#Get window positional information
	gene_counts <- fread(sampleFiles[1], header = FALSE, sep = "\t", col.names = c("chr", "chrStart", "chrEnd"), select = 1:3)

	#Get window center position for plotting
	gene_counts$abs.pos <- NA
	gene_counts$abs.pos <- as.data.frame(apply(gene_counts, 1, function(x) mean(c(as.numeric(x[2]), as.numeric(x[3]))) + chr_length[chr_length$chr == x[1],]$chr_start))
	abs.pos <- gene_counts$abs.pos

	#load count data by gene
	#pb <- txtProgressBar(min=0, max=length(sampleFiles), initial = 0, style = 3)
	for(i in 1:length(sampleFiles)){
		gene <- unlist(strsplit(sampleFiles[i], "\\."))[1]

		gene_counts[,gene] <- fread(sampleFiles[i], header = FALSE, sep = "\t", select = 4)

	#	setTxtProgressBar(pb, i)
	}

	return(as.data.frame(gene_counts))
	}
}

getInterchromMtrx <- function(countMatrix){
	interchromCounts <- countMatrix

	for(i in 5:ncol(interchromCounts)){
		intra.chr <- genes$chr[genes$gene == colnames(interchromCounts)[i]]
		interchromCounts[interchromCounts$chr == intra.chr,i] <- NA
	}

	return(interchromCounts)
}



getIntrachromMtrx <- function(countMatrix, distance){
	intrachromCounts <- countMatrix
	intrachromCounts$mid_window <- intrachromCounts$chrStart + (intrachromCounts$chrEnd - intrachromCounts$chrStart)/2
	for(i in 5:ncol(intrachromCounts)-1){
		intra.chr <- genes$chr[genes$gene == colnames(intrachromCounts)[i]]
		# Add distance criteria
		gene_mid_point <- genes$mid_gene[genes$gene == colnames(intrachromCounts)[i]]
		intrachromCounts[(intrachromCounts$chr == intra.chr) & (abs(gene_mid_point - intrachromCounts$mid_window) < distance),i] <- NA
		# Remove inter-chromosomal interactions
		intrachromCounts[intrachromCounts$chr != intra.chr,i] <- NA
	}
	intrachromCounts$mid_window <- NULL
	return(intrachromCounts)
}
