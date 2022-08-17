require(ggplot2)
require(reshape2)

chrHeatPlot <- function(interchrCountMtrx, top_return = 5){
	chr.counts <- aggregate(interchrCountMtrx[,5:ncol(interchrCountMtrx)], by=list(Chromosome=interchrCountMtrx$chr), FUN=sum)
	chr.counts$Chromosome <- factor(chr.counts$Chromosome, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"))
	chr.counts <- cbind(Chromosome = chr.counts$Chromosome, 1000000*chr.counts[,2:ncol(chr.counts)]/chr_length[match(chr.counts$Chromosome, chr_length$chr),2])

	chr.counts.long <- melt(chr.counts)
	colnames(chr.counts.long) <- c("Chromosome", "Gene", "CntPerMb")

	print(noquote("Max value, interchromosomal interactions:"))
	print(head(chr.counts.long[order(chr.counts.long$CntPerMb, decreasing = TRUE),], top_return))

	return(ggplot(chr.counts.long, aes(Chromosome, Gene, fill= CntPerMb)) + geom_tile() + scale_fill_distiller(palette = "Blues"))
}

#################################

genomeManhattan <- function(cntMtrx, gene, pvals = NA, sig.cutoff=0.05, partner.gene = NA, plot.title = NA, calc.pval = FALSE, padj.method = "BH", plot.min = 1, sig.col = "violetred2", annot.col = "darkturquoise", pnt.col = "grey10"){
	if(length(pvals) == 1 & calc.pval == TRUE){
		pvals <- p.adjust(inter.pvals.gene(cntMtrx, gene), method = padj.method)
	}

	plot.df <- cntMtrx[,c("chr", "abs.pos", gene)]
	pvals <- pvals[plot.df[,gene] >= plot.min]				#Removes 0 count rows to reduce drawing time, help with log transforms.
	plot.df <- plot.df[plot.df[,gene] >= plot.min,] 		#Removes 0 count rows to reduce drawing time, help with log transforms.
	max.y = max(plot.df[,gene])
	ylim = c(0.9, 1.25*max.y)

	gplot <- ggplot(plot.df, aes(x = abs.pos, y = plot.df[,3])) +
		geom_rect(aes(xmin=chr_length[1,3], ymin=ylim[1],xmax= chr_length[25,3]+chr_length[25,2], ymax=ylim[2]), fill = "white") +
		geom_rect(aes(xmin=chr_length[2,3], ymin=ylim[1], xmax= chr_length[3,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[4,3], ymin=ylim[1], xmax= chr_length[5,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[6,3], ymin=ylim[1], xmax= chr_length[7,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[8,3], ymin=ylim[1], xmax= chr_length[9,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[10,3], ymin=ylim[1], xmax= chr_length[11,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[12,3], ymin=ylim[1], xmax= chr_length[13,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[14,3], ymin=ylim[1], xmax= chr_length[15,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[16,3], ymin=ylim[1], xmax= chr_length[17,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[18,3], ymin=ylim[1], xmax= chr_length[19,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[20,3], ymin=ylim[1], xmax= chr_length[21,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[22,3], ymin=ylim[1], xmax= chr_length[23,3], ymax=ylim[2]), fill = "gray90") +
		geom_rect(aes(xmin=chr_length[24,3], ymin=ylim[1], xmax= chr_length[25,3], ymax=ylim[2]), fill = "gray90") +
		#scale_x_continuous(name = "Chromosome", breaks = (chr_length$chr_start[1:24]+chr_length$chr_start[2:25])/2, minor_breaks = NULL, labels = c(1:22, "X", "Y"), limits = c(0, chr_length[25,3]+chr_length[25,2]), expand = expansion(mult=0, add=0), position = "bottom") +
		scale_x_continuous(name = "Chromosome", breaks = (chr_length$chr_start[1:24]+chr_length$chr_start[2:25])/2, minor_breaks = NULL, labels = c(1:22, "X", "Y"), limits = c(0, chr_length[25,3]+chr_length[25,2]), position = "bottom") +
		#scale_y_log10(name = "Number of Reads", breaks = c(1, 10, 100, 1000, 10000), minor_breaks = NULL, limits = ylim, expand = expansion(mult=0, add=0)) +
		scale_y_log10(name = "Number of Reads", breaks = c(1, 10, 100, 1000, 10000), minor_breaks = NULL, limits = ylim) +
		labs(title = ifelse(is.na(plot.title), gene, plot.title)) +
		theme(plot.title = element_text(size=8), axis.ticks = element_line(size=0), axis.line = element_line(size=0), axis.text=element_text(size=7), axis.title=element_text(size=8)) +
		geom_vline(xintercept = mean(c(all_genes$absStart[all_genes$gene==gene], all_genes$absEnd[all_genes$gene==gene])), color = annot.col, size=0.3) +
		geom_vline(xintercept = ifelse(is.na(partner.gene), 0, mean(all_genes$absStart[all_genes$gene==partner.gene], all_genes$absEnd[all_genes$gene== partner.gene])), color = annot.col, size = ifelse(is.na(partner.gene), 0, 0.3)) +
		geom_rect(aes(xmin=chr_length[1,3], ymin=ylim[1],xmax= chr_length[25,3]+chr_length[25,2], ymax=ylim[2]), fill=rgb(1, 1, 1, 0), color = "black", size = 0.1) +
		geom_point(col=ifelse(is.na(pvals), pnt.col, ifelse(pvals <= sig.cutoff, sig.col, pnt.col)), size = 0.2) +
		annotate(geom = "text", label= ifelse(is.na(partner.gene), "", partner.gene), x = ifelse(is.na(partner.gene), 0, all_genes$absStart[all_genes$gene==partner.gene]-10000000), y = .9*max.y, hjust = "right", size = 2, color = annot.col) +
		annotate(geom = "text", label= gene, x = all_genes$absEnd[all_genes$gene==gene]+10000000, y = .9*max.y, hjust = "left", size = 2, color = annot.col)

	return(gplot)
}

########################################

chrManhattan <- function(cntMtrx, gene, chromosome, pvals = NA, sig.cutoff=0.05, show.gene = NA, plot.min = 1, plot.title = NA, calc.pval = FALSE, padj.method = "BH", sig.col = "violetred2", annot.col = "darkturquoise", pnt.col = "grey10"){

	if(length(pvals) == 1 & calc.pval == TRUE){
		pvals <- p.adjust(inter.pvals.gene(cntMtrx, gene), method = padj.method)
	}

	plot.df <- cntMtrx[,c("chr", "chrStart", "chrEnd", gene)]
	pvals <- pvals[plot.df$chr == chromosome & plot.df[,gene] >= plot.min]
	plot.df <- plot.df[plot.df$chr == chromosome & plot.df[,gene] >= plot.min,]
	max.y = max(c(plot.df[,gene], sig.cutoff))
	ylim = c(0.9, 1.5*max.y)

	ggplot(plot.df, aes(x=chrStart, y=plot.df[,gene])) +
	scale_x_continuous(name = chromosome, breaks = seq(from=0, to=chr_length$len[chr_length$chr == chromosome], by = 10^7), minor_breaks = NULL, limits = c(0, chr_length$len[chr_length$chr == chromosome]), expand = expansion(mult=0, add=0), position = "bottom") +
	scale_y_log10(name = "Number of Reads", breaks = c(1, 10, 100, 1000, 10000), minor_breaks = NULL, limits = ylim, expand = expansion(mult=0, add=0)) +
	labs(title = ifelse(is.na(plot.title), paste(gene, "interactions with", chromosome, sep = " "), plot.title)) +
	theme(panel.background=element_rect(fill = "white", color = "black", size = 0.1), plot.title = element_text(size=8), axis.ticks.y = element_line(size=0), axis.ticks.x = element_line(size=0.1), axis.line = element_line(size=0), axis.text=element_text(size=7), axis.title=element_text(size=8), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
		geom_rect(aes(xmin= ifelse(all_genes$chr[all_genes$gene==gene]!=chromosome, 0, all_genes$chrStart[all_genes$gene==gene]), ymin=1.1*max.y, xmax= ifelse(all_genes$chr[all_genes$gene==gene]!=chromosome, 0, all_genes$chrEnd[all_genes$gene==gene]), ymax=1.4*max.y), fill = ifelse(all_genes$chr[all_genes$gene==gene]!=chromosome, rgb(1,1,1,0), annot.col)) +
	geom_rect(aes(xmin= ifelse(is.na(show.gene), 0, all_genes$chrStart[all_genes$gene==show.gene]), ymin=1.1*max.y, xmax= ifelse(is.na(show.gene), 0, all_genes$chrEnd[all_genes$gene==show.gene]), ymax=1.4*max.y), fill = ifelse(is.na(show.gene), rgb(1,1,1,0), annot.col)) +
	geom_point(col=ifelse(is.na(pvals), pnt.col, ifelse(pvals <= sig.cutoff, sig.col, pnt.col)), size = 0.2) +
	annotate(geom = "text", label= ifelse(all_genes$chr[all_genes$gene==gene]!=chromosome, "", gene), x = ifelse(all_genes$chr[all_genes$gene==gene]!=chromosome, 0, all_genes$chrEnd[all_genes$gene== gene]+ 500000), y = 1.34*max.y, hjust = "left", vjust = "top", size = 2, color = annot.col) +
	annotate(geom = "text", label= ifelse(is.na(show.gene), "", show.gene), x = ifelse(is.na(show.gene), 0, all_genes$chrStart[all_genes$gene== show.gene]- 500000), y = 1.34*max.y, hjust = "right", vjust = "top", size = 2, color = annot.col)
}

# gene.plot <- function(data, plot.title, x.label, sig.val=10, genes=c("BCR", "ABL1"), chromosome, target.gene.coord, ylim = c(0.9, 1.35*max.y)){

	# plot.df <- data[data$chr == chromosome & data$read.cnt > 1,]
	# max.y = max(plot.df$read.cnt)

	# ggplot(plot.df, aes(x=chr.pos, y=read.cnt)) +
	# scale_x_continuous(name = x.label, minor_breaks = NULL, expand = expansion(mult=0, add=0), limits = c(target.gene.coord[1]-10000, target.gene.coord[2]+10000), position = "bottom") +
	# scale_y_log10(name = "Number of Reads", breaks = c(1, 10, 100, 1000, 10000), minor_breaks = NULL, limits = ylim, expand = expansion(mult=0, add=0)) +
	# labs(title = plot.title) +
	# theme(panel.background=element_rect(fill = "white", color = "black", size = 0.1), plot.title = element_text(size=8), axis.ticks = element_line(size=0), axis.line = element_line(size=0), axis.text=element_text(size=7), axis.title=element_text(size=8), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
	# geom_rect(aes(xmin= target.gene.coord[1], ymin=1.1*max.y, xmax= target.gene.coord[2], ymax=1.3*max.y), fill = "darkturquoise") +
	# geom_point(col=ifelse(plot.df$read.cnt > sig.val & plot.df$chr != "chr22", "violetred2", pnt.col), size = 0.2) +
	# annotate(geom = "text", label= genes[2], x = (target.gene.coord[1]+target.gene.coord[2])/2, y = .91*max.y, hjust = "center", size = 2, color = "darkturquoise") +
	# geom_hline(yintercept= sig.val, color="violetred2", linetype=2, size=0.3)
# }
