################################################################################
## These are functions that you can use

mor <- function(data){
  #' Median of Ratios RNA Seq normalisation function
  #' 
  #' This function performs a normalisation of RNA-Seq count data using the median of
  #' ratios method.
  #' 
  #' @param data: A data frame of expression counts. rows are genes, columns are samples
  #' 
  #' @return Returns a list of $norm (with the normalised data), $sizefactors (the size factors that were determined), and $pseudoref (the vector that was used as a reference)

  # Define a reference sample (a new sample)
  pseudoref <- apply(data, 1, function(x){exp(mean(log(x[x>0])))}) 
  # Determine ratio of this reference for each gene in each sample
  ratio <- apply(data, 2, function(x){x / pseudoref})
  
  # Size factor per sample is the median ratio across all genes  
  sizefactor <- apply(ratio, 2, function(d){ median(d[d>0])}) 
  sizefactor <- 1 / sizefactor
  
  # Normalize
  mormed <- t(apply(data, 1, function(x){x * sizefactor}))
  mormed <- as.data.frame(mormed)
  
  # Return a list with different elements of the normalization
  list(norm=mormed,
       sizefactors=sizefactor,
       pseudoref=pseudoref)
}

diffex.test.all <- function(form, data, meta){
  #' Differential expression testing function
  #' 
  #' A differential expression test for each gene in a dataframe of normalised gene expression counts.
  #' The test is performed using the negative binomial distribution.
  #' 
  #' @param form: A formula of the 
  #' @param data: A data frame of normalised expression counts. rows are genes, columns are samples.
  #' 
  #' @return Returns a list of $norm (with the normalised data), $sizefactors (the size factors that were determined), and $pseudoref (the vector that was used as a reference)
  
  require(MASS)
  updated.form <- update.formula(form, gene ~ .)
  meta.gene <- meta
  pb <- txtProgressBar(min=0, max=dim(data)[1], initial=0, style=3)
  R <- Reduce(rbind, apply(data, 1, function(expr){
    tryCatch({
      meta.gene$gene <- expr
      fit <- glm.nb(updated.form, data=meta.gene, na.action=na.omit)
      res <- as.data.frame(summary(fit)$coefficients)[2,]
      return(res)
    }, error=function(cond) {
      missing <- as.data.frame(list(NA,NA,NA,NA))
      colnames(missing) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
      return(missing)
    },
    finally={setTxtProgressBar(pb, getTxtProgressBar(pb)+1)})
  }))
  rownames(R) <- rownames(data)
  R$qvalue <- p.adjust(R$`Pr(>|z|)`, "fdr")
  R
}

volcano <- function(diffex.res, q.thresh=0.05, fc.thresh=1){
  diffex.res <- diffex.res[unlist(apply(diffex.res[,c("Estimate","Std. Error", "Pr(>|z|)", 'qvalue')], 1, function(x){all(is.finite(x))})),]
  diffex.res$`Pr(>|z|)t` <- sapply(diffex.res$`Pr(>|z|)`, function(x){max(x,2**-200)})
  p.thresh <- max(diffex.res[diffex.res$qvalue < q.thresh,]$`Pr(>|z|)`)
  significant <- diffex.res[(diffex.res$qvalue < q.thresh) & (abs(diffex.res$Estimate) >= fc.thresh),]
  insignificant <- diffex.res[!((diffex.res$qvalue < q.thresh) & (abs(diffex.res$Estimate) >= fc.thresh)),]
  
  pvals <- diffex.res$`Pr(>|z|)`[is.finite(log2(diffex.res$`Pr(>|z|)`))]
  ylim <- c(0, max(-log2(pvals)))*1.1
  xlim <- c(min(significant$Estimate), max(significant$Estimate))*2
  
  plot(insignificant$Estimate, -log2(insignificant$`Pr(>|z|)`),
       xlim=xlim, ylim=ylim,
       xlab="Log Fold Change",
       ylab="-log p-value")
  points(significant$Estimate, -log2(significant$`Pr(>|z|)`), col=sign(significant$Estimate)+3)
  lines(xlim, -log2(c(p.thresh, p.thresh)), col='black')
  lines(-c(fc.thresh, fc.thresh), ylim, col='black')
  lines(c(fc.thresh, fc.thresh), ylim, col='black')
}


################################################################################
## This is the pre-processing steps

# This is the location of data that we want to download
skin.sun.url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_skin_sun_exposed_lower_leg.gct.gz"
skin.url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_skin_not_sun_exposed_suprapubic.gct.gz"
subject.url <- "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"

# Let's download the files
download.file(skin.url, dest="skin.txt.gz")
download.file(skin.sun.url, dest="skin_sun.txt.gz")
download.file(subject.url, dest="subject_phenotype.txt")

# Open the files and process them a little bit
skin     <- read.table("skin.txt.gz", sep = '\t', skip = 2, header=TRUE, row.names=2)
skin.sun <- read.table("skin_sun.txt.gz", sep = '\t', skip = 2, header=TRUE, row.names=2)
S        <- read.table("subject_phenotype.txt", sep = '\t', header=TRUE)

# Keep the gene names in another list
genenames <- skin[,c('Description')]

#Remove columns that are not gene expression data, and merge the dataframes
skin     <- subset(skin, select=-c(id,Description))
skin.sun <- subset(skin.sun, select=-c(id,Description))

# The data files are a bit large, let's only take 150 samples for each group
n.pergroup <- 150
skin     <- skin[,1:n.pergroup]
skin.sun <- skin.sun[,1:n.pergroup]

E <- cbind(skin,skin.sun)
E <- E[rowSums(E) != 0,] # Remove those with no expression whatsoever

# Make a dataframe with the meta data in it
M <- data.frame(sample.id=c(colnames(skin), colnames(skin.sun)),
                sun=c(rep(FALSE, dim(skin)[2]), rep(TRUE, dim(skin.sun)[2])))
M$SUBJID <- unlist(lapply(M$sample.id, function(x){paste0(c('GTEX-',unlist(strsplit(x, "[.]"))[2]),collapse='')}))
M <- merge(M, S)

rownames(E) <- unlist(lapply(rownames(E), function(x){
  strsplit(x, "[.]")[[1]][1]}))

################################################################################
## From this point onwards, the students must continue.



