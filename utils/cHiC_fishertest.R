inter.pvals.gene <- function(input.df, gene, mk.inter.matrix = TRUE, dist = 0){

	if(mk.inter.matrix == TRUE){
		input.df <- getInterchromMtrx(input.df)
	}

	if(mk.inter.matrix  == FALSE){
		input.df <- getIntrachromMtrx(input.df, dist)
	}

	input.df <- input.df[,5:ncol(input.df)]

	win.sums <- rowSums(input.df, na.rm = TRUE)
	gene.sums <- colSums(input.df, na.rm = TRUE)
	df.total <- sum(input.df, na.rm = TRUE)

#	print(paste("win sums is ", win.sums, sep = ""))
#	print(paste("gene sums is ", gene.sums, sep = ""))
#	print(paste("df.total is ", df.total, sep = ""))

	pval <- numeric(length = nrow(input.df))
	pval[input.df[,gene] == 0] <- 1
	pval[is.na(input.df[,gene])] <- NA

	remaining <- !is.na(input.df[,gene]) & input.df[,gene] != 0
	# test.df <- input.df[remaining, gene]
#	str(remaining)
	pval.remain <- numeric(length = sum(remaining));
#	str(pval.remain)
	if (length(pval.remain) != 0){
	#	pb <- txtProgressBar(min = 0, max = length(pval.remain), style = 3)

		for(i in 1:length(pval.remain)){
			a <- input.df[remaining, gene][i]
			b <- gene.sums[gene] - a
			c <- win.sums[remaining][i] - a
			d <- df.total - b - c + a

			pval.remain[i] <- fisher.test(matrix(c(a, b, c, d), nrow=2))$p.value
	#		setTxtProgressBar(pb, i)
		}
	}
	pval[remaining]<-pval.remain

	return(as.numeric(pval))
}

inter.pvals.all <- function(input.df){
	cont.list <- list()
	win.sums <- rowSums(input.df, na.rm = TRUE)
	gene.sums <- colSums(input.df, na.rm = TRUE)
	df.total <- sum(input.df, na.rm = TRUE)

	#pb<-txtProgressBar(min = 0, max = ncol(input.df), style = 3)

	for(j in 1:ncol(input.df)){
		gene.tbl <- data.frame(a=numeric(length = nrow(input.df)),
							   b=numeric(length = nrow(input.df)),
							   c=numeric(length = nrow(input.df)),
							   d=numeric(length = nrow(input.df)))
		for(i in 1:nrow(input.df)){
			a <- input.df[i, j]
			b <- gene.sums[j] - a
			c <- win.sums[i] - a
			d <- df.total - b - c + a

			gene.tbl[i,] <- c(a, b, c, d)
		}

		cont.list[[j]]<-gene.tbl
		names(cont.list)[j] <- colnames(input.df)[j]
		#setTxtProgressBar(pb, j)
	}

	names(cont.list) <- colnames(input.df)

	#print("Calculating p-values...")

	pval <- lapply(cont.list, FUN = function(y) {apply(y, 1, function(x){ifelse(is.na(x[1]), NA, ifelse(x[1] == 0, 1, fisher.test(matrix(x, nrow=2))$p.value))})})

	return(as.data.frame(pval))
}
