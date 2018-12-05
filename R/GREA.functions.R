#' Binning genes into genomic windows.
#' \code{bin.genes.genomic.regions} Adds an extra colum \code{bin.tag} grouping genes to a genomic bins defined by an input \code{window}

#' @param gene.pos Data.frame with "geneid" (ensembl ids), "chr" (only somatic, and as integers (1:22)), "left"(start), "right" (end). For reference please look at data(gene.info).
#' @param window Number of bases to define a genomic loci, 1MB (1000000) is set as defualt
#' @return A dataframe with \code{gene.set}, its genomic position and an extra colum \code{loci} grouping genes to a genomic loci
#' @examples
#' gene.pos.test.bin <- bin.genes.genomic.regions(gene.pos.test, window=1000000)
#'
bin.genes.genomic.regions <- function(gene.pos,
                                      window=500000){
  chrs <- unique(gene.pos$chr)
  if(all(chrs %in% 1:22)==FALSE){
    if("X" %in% chrs){
      gene.pos$chr[which(gene.pos$chr == "X")] <- 23
    }
    if("Y" %in% chrs){
      gene.pos$chr[which(gene.pos$chr == "Y")] <- 24
    }
  }
  gene.pos$chr <- as.integer(gene.pos$chr)
  gene.pos <- gene.pos[order(gene.pos$chr, gene.pos$left),]
  bin.tag.list <- list()

  ## re-define chrs to include new X and Y chromosomes
  chrs <- unique(gene.pos$chr)

  for(i in chrs){
    i.chr.gene.pos <- gene.pos[which(gene.pos$chr == i),]
    min.genes <- min(i.chr.gene.pos$left)
    max.genes <- max(i.chr.gene.pos$right)
    window.bin <- c(min.genes, min.genes+window)
    bin.count <- 1
    chr.bins <- rep(NA, times= nrow(i.chr.gene.pos))
    while(window.bin[1] < max.genes){
      i.bin.name    <- paste0("bin-chr", i, ".", bin.count)
      binned.genes  <- which(
        (i.chr.gene.pos$left >= window.bin[1] & i.chr.gene.pos$right <= window.bin[2]) | # gene within the window.
          (i.chr.gene.pos$left >= window.bin[1] & i.chr.gene.pos$right >= window.bin[2]) #| # gene starts within the window but ends after
        #(i.chr.gene.pos$left <= window.bin[1] & i.chr.gene.pos$right <= window.bin[2])   # gene starts before the window but ends within
      )
      genes.in.bin  <- length(binned.genes)
      chr.bins[binned.genes] <- rep(i.bin.name, times=genes.in.bin)
      bin.count     <- bin.count+1
      window.bin    <- window.bin+window
    }
    bin.tag.list[[i]] <- chr.bins
  }
  bin.tag <- unlist(bin.tag.list)
  gene.pos$bin.tag <- bin.tag
  #gene.pos$bin.tag <- as.factor(gene.pos$bin.tag)


  ### put back chr names
  if(23 %in% gene.pos$chr){
    gene.pos[gene.pos$chr %in% 23, "chr"] <- "X"
  }
  if(24 %in% gene.pos$chr){
    gene.pos[gene.pos$chr %in% 24, "chr"] <- "Y"
  }

  return(gene.pos)
}


#' get.sampled.genes
#' \code{get.sampled.genes} Given a set of \code{input.genes} (which should be in enemsebl) this function wil return a list of \code{n.genes.sets} of genes which represent the same number of genomic regions as the \code{input.genes}

#' @param input.genes Character vector of ensembl ids
#' @param gene.info A data.frame object with the genomic coordinates of the genes used in the study.
#' @param binning.widow Numeric, window size which will define the genomic regions and bining the genes. Default is 125 Kb to each side.
#' @param n.genes.sets Number of random gene sets that are going to be generated.
#' @param parallel Run the sampling iterations in parallel. At the moment this is the only available option.
#' @return A dataframe with \code{gene.set}, its genomic position and an extra colum \code{loci} grouping genes to a genomic loci

#input.genes <- clust.info$ensembl[which(clust.info$cluster == 1)]
get.sampled.genes <- function(input.genes,
                              gene.info.binned,
                              binning.widow=125000,
                              n.excess.genes.limit= 10,
                              n.genes.sets= 1000,
                              parallel = TRUE,
                              n.clusters=3){

  cl <- snow::makeCluster(n.clusters, type = "SOCK")
  doSNOW::registerDoSNOW(cl= cl)
  require(foreach)

  input.genes.gene.info <- gene.info.binned[match(input.genes, as.character(gene.info.binned$gene_id)),]
  ## remove non-somatic genes
  input.genes.gene.info <- input.genes.gene.info[which(input.genes.gene.info$chr %in% c(1:22)),]

  #* bin all these genes
  #input.genes.binned <- bin.genes.genomic.regions(gene.pos = input.genes.gene.info, window = binning.widow)
  n.genomic.regions <- length(unique(input.genes.gene.info$bin.tag))
  n.genes <- length(input.genes)

  if(n.genomic.regions >n.genes ){stop("You cannot have more genomic regions than genes")}
  #permuted.gene.sets <- list()
  #permuted.gene.sets <- vector(mode = "list", length = n.genes.sets)
  avail.bin.tags <- unique(gene.info.binned$bin.tag)

  permuted.gene.sets <- foreach::foreach(j=1:n.genes.sets) %dopar% {

    j.bin <- sample(avail.bin.tags, size = n.genomic.regions, replace = FALSE)
    ## this part is to level up the number of genes with the observed ones.
    j.binned.gene.info <- gene.info.binned[which(gene.info.binned$bin.tag %in% j.bin == TRUE),]
    n.excess.genes <- nrow(j.binned.gene.info)-n.genes
    #remove.genes.counter <- 1
    while(n.excess.genes >= n.excess.genes.limit){
      j.table.bins <- table(j.binned.gene.info$bin.tag)
      bins.with.many.genes <- j.table.bins[names(j.table.bins[which(j.table.bins > 1)])]

      remove.index.per.bin <- unlist(lapply(names(bins.with.many.genes), function(z){
        sample(which(j.binned.gene.info$bin.tag == z), size=1)}))
      j.binned.gene.info <- j.binned.gene.info[-remove.index.per.bin,]
      n.excess.genes <- nrow(j.binned.gene.info)-n.genes
      #remove.genes.counter <- remove.genes.counter+1
    }
    #permuted.gene.sets[[j]] <- j.binned.gene.info[,"ensembl"]
    #if(j %% 500 == 0){cat("[INFO]\t", "Done sampling iteration", j, "\n")}
    # get back only ensembl IDs.
    j.binned.gene.info[,"gene_id"]
  }
}

#' get.genes.within.loci
#' \code{get.genes.within.loci} Given a set of \code{positions}(genomic coordinates), will return a data.frame with the genes within the input \code{positions}

#' @param positions A data.frame with genomic coordinates.
#' @param position.colnames Character vector defining the colnames used for defining the genomic coordinates within the \code{positions} data.frame.
#' @param gene.window Numeric, which will deifine the genomic window that will be used to get neighbouring genes given a certain genomic position defined by \code{positions}.
#' @param gene.info A data.frame with genomic coordinates and genes. For reference please look at data(gene.info)
#' @param pVal.col Character, if provided, then an extra column to the returned data-frame will be added with the pValue of the SNP associated.
#' @return A dataframe

get.genes.within.loci <- function(positions,
                                  position.colnames = c("CHR", "BP", "SNP"),
                                  gene.window=500000,
                                  gene.info,
                                  pVal.col=NULL){

  if(!base::exists(x = "gene.info")){
    stop("Please load gene.info object")
  }
  if(colnames.check(c.names = position.colnames, df = positions) == FALSE){
    stop("Make sure that position.colnames %in%  colnames(df)")
  }
  names(position.colnames) <- c("CHR", "BP", "SNP")
  gene.pos.list <- list( )
  for(i in 1:nrow(positions)){
    i.chr.genes <- gene.info[which(gene.info[,"chr"] == positions[i,position.colnames["CHR"]]),]

    i.chr.gene.index <- c(
      which(i.chr.genes$left > (positions[i, position.colnames["BP"]] - gene.window) &
              i.chr.genes$right < (positions[i,position.colnames["BP"]] + gene.window)), ## all genes within the window

      which(i.chr.genes$left < (positions[i,position.colnames["BP"]] - gene.window) &
              i.chr.genes$right < (positions[i,position.colnames["BP"]] + gene.window) &
              i.chr.genes$right > (positions[i,position.colnames["BP"]] - gene.window)), ## genes with a start before window but ends within

      which(i.chr.genes$left > (positions[i,position.colnames["BP"]] - gene.window) &
              i.chr.genes$right > (positions[i,position.colnames["BP"]] + gene.window)&
              i.chr.genes$left < (positions[i,position.colnames["BP"]] + gene.window)) ## genes with a start within window but ends past window
    )

    if(length(i.chr.gene.index) >= 1){
      i.chr.gene.index <- unique(i.chr.gene.index)
      gene.pos.list[[i]] <- i.chr.genes[i.chr.gene.index,]
    }
  }
  names(gene.pos.list) <- as.character(positions[,position.colnames["SNP"]])
  gene.pos.loci <- do.call(rbind, gene.pos.list)
  gene.pos.loci$SNP <- gsub(rownames(gene.pos.loci), pattern = "\\..*", replacement = "")
  #gene.pos.loci <- reshape2::melt(gene.pos.list, id.vars=c("ensembl", "symbol","biotype","chr","left","right"))
  gene.pos.loci$SNP.pos <- positions[match(gene.pos.loci[,"SNP"],positions[,position.colnames["SNP"]]),position.colnames["BP"]]

  if(!is.null(pVal.col)){
    gene.pos.loci$SNP.pValue <- positions[match(gene.pos.loci[,"SNP"],positions[,position.colnames["SNP"]]),pVal.col]
  }

  return(gene.pos.loci)
}



#' prune.loci.by.proximity
#' \code{prune.loci.by.proximity} Given a set of \code{positions}(genomic coordinates), will return a data.frame with the genes within the input \code{positions}, and if pValue is provided than will trim the original \code{positions} leaving only the SNP with the lowest pValue.

#' @param positions A data.frame with genomic coordinates.
#' @param position.colnames Character vector defining the colnames used for defining the genomic coordinates within the \code{positions} data.frame.
#' @param gene.window Numeric, which will deifine the genomic window that will be used to get neighbouring genes given a certain genomic position defined by \code{positions}.
#' @param gene.info A data.frame with genomic coordinates for genes. For reference please look at data(gene.info)
#' @param pVal.col Character, if provided, then an extra column to the returned data-frame will be added with the pValue of the SNP associated.
#' @return A dataframe
#'
prune.loci.by.proximity <- function(positions,
                                    window="2MB",
                                    position.colnames= c("chr", "pos"),
                                    prune= FALSE, pvalCol= NULL){
  #positions is a dataframe with chr and pos as columns
  if(colnames.check(c.names = position.colnames, df = positions) == FALSE){
    stop("Make sure that position.colnames %in%  colnames(df)")
  }
  window <- as.numeric(gsub(window , pattern = "MB", replacement = ""))* 1000000
  positions$loci <- NA
  lociCounter <- 1
  while(sum(is.na(positions$loci)) != 0){
    nextEvalRow <- which(is.na(positions$loci))[1]
    iWindows <- c(positions[nextEvalRow, position.colnames[2]] - window/2, positions[nextEvalRow , position.colnames[2]] + window/2)
    iChr <- positions[nextEvalRow, position.colnames[1]]
    evalRows <- which(positions[,position.colnames[1]] ==  iChr &
                        positions[,position.colnames[2]] >= iWindows[1] &
                        positions[,position.colnames[2]] <= iWindows[2])
    if(length(evalRows) == 0 || length(evalRows) == 1) {
      positions$loci[nextEvalRow] <- lociCounter
    } else{
      positions$loci[evalRows] <- lociCounter
    }
    lociCounter <- lociCounter+1
  }

  if(prune & !is.null(pvalCol)){
    positions <- do.call(rbind, by(positions, positions$loci, function(x) x[which.min(x[,pvalCol]), ] ))
  } else if (prune & is.null(pvalCol)){
    stop("Please specify the colname")
  }

  return(positions)
}


colnames.check <- function(c.names, df){
  all(c.names %in% colnames(df))
}


#' GREA.test
#' \code{GREA.test} plost the ditribution and observed overlap
#' @param sampled.genes.sets A list where each element is a character vector of ensembl IDs, from the \code{get.sampled.genes} function
#' @param query.gene.set A character vector of ensembl IDs
#' @param genomic.region.gene.sets A list where each element is a data.frame that includes at least SNPs and genes as ensembl IDs.
#' @param plot boolean, if true genrates a plot with the empirical null distribution and the observed overlapped given bye \code{query.gene.set}
#' @return A GREAsummary data.frame
GREA.test <- function(sampled.genes.sets,
                      query.gene.set,
                      genomic.region.gene.sets,
                      plot=TRUE){

  cat("\n[INFO]\t A total of ", length(sampled.genes.sets), "sampled gene sets will be used to generate an empirical null distribution")
  cat("\n[INFO]\t ",length(genomic.region.gene.sets), "gene sets derived from genomic regions are being considerated")

  n.unique.grs <- sapply(genomic.region.gene.sets, function(x){length(unique(x[,"SNP"]))})

  # count number of SNPs have at least 1 neighboring gene is present in query.gene.set
  overlap.per.gr <- sapply(genomic.region.gene.sets,
                              function(x){length(unique(x[,"SNP"][which(x[,"gene_id"] %in% query.gene.set)]))
                              })
  # percentage of count
  perc.per.gr <- overlap.per.gr/n.unique.grs


  # get genes for final report
  overlapped.ensembl.per.gr <- lapply(genomic.region.gene.sets, function(x){query.gene.set[query.gene.set %in% x[,"gene_id"]]})
  overlapped.symbol.per.gr <- lapply(genomic.region.gene.sets, function(x){x$gene_name[x[,"gene_id"] %in% query.gene.set]})

  overlapped.symbol.per.gr.vector <- sapply(overlapped.symbol.per.gr, function(x){paste(x, collapse = "/")})

  ## count number of genes within random sets overlaps with the GWAS "genes"
  null.overlap.per.gr <- lapply(sampled.genes.sets, function(y){
                              sapply(genomic.region.gene.sets, function(x){
                                  length(unique(x[,"SNP"][which(x[,"gene_id"] %in% y )]))
                          })})
  ## get it in percentage
  null.perc.overlap.per.gr <- lapply(null.overlap.per.gr, function(x){x/n.unique.grs})
  null.perc.overlap.per.gr <- do.call(rbind, null.perc.overlap.per.gr)
  rownames(null.perc.overlap.per.gr) <- NULL

  enrich.pvals <- sapply(1:length(genomic.region.gene.sets),
                         function(x){pvalRight(x= perc.per.gr[x], dist = null.perc.overlap.per.gr[,x])})
  #depleted.pvals <- sapply(1:length(full.gwas.genes.list),
  #                        function(x){pvalLeft(x= perc.per.trait[x], dist = null.perc.overlap.per.trait[,x])})

  #Summary report
  GREA.summary <- data.frame(trait= names(perc.per.gr),
                             perc.per.gr= perc.per.gr,
                              enrich.pvals=enrich.pvals,
                              #depleted.pvals=depleted.pvals,
                              mean.permutation.overlap= apply(null.perc.overlap.per.gr, 2, mean),
                              SD.permutation.overlap= apply(null.perc.overlap.per.gr, 2, sd),
                              overlappedGenes=overlapped.symbol.per.gr.vector
                              )

  if(plot){
    null.dist.plot <- null.dist.plotter(null.perc.overlap.per.gr=null.perc.overlap.per.gr,
                                        GREA.summary=GREA.summary)

    return(list(GREA.summary=GREA.summary, null.dist.plot=null.dist.plot))
  } else {
    return(GREA.summary=GREA.summary)
  }

}

#' null.dist.plotter
#' \code{null.dist.plotter} plots the ditribution and observed overlap
#' @param null.perc.overlap.per.gr a list with the overlap between sampled genes and
#' @param GREA.summary a GREA.summary data.frame calculated through the \code{GREA.test} function
#' @return a ggplot object
null.dist.plotter <- function(null.perc.overlap.per.gr,
                              GREA.summary){
  pData <- melt(null.perc.overlap.per.gr)


  vline.data <- GREA.summary[,c("perc.per.gr", "trait", "enrich.pvals")]
  colnames(vline.data) <- c("perc.per.trait", "Var2", "enrich.pvals")
  vline.data$line.color <- "Significant"
  vline.data$line.color[which(enrich.pvals > 0.05)] <- "NonSignificant"

  null.dist.plot <- ggplot(pData, aes(x=value))+
                      geom_bar(stat="count", alpha= 0.7)+
                      facet_wrap(~Var2, scales = "free")+
                      geom_vline(data=vline.data, aes(xintercept= perc.per.trait, color= line.color))+
                      theme_bw()+
                      theme(text= element_text(size=10, family="Helvetica"))

  return(null.dist.plot)
}

#' pvalRight
#' \code{pvalRight} given an empirical distribtion it calculates the quantile in which an observed measurement is located within the distribution
#' @param x numeric, observed value to be compared against a null distribution
#' @param dist numeric vector defining a null distribution
#' @return return the estimated pValue based on the null distrubtion provided
pvalRight <- function(x, dist) {
  sum(dist > x) / length(dist)
}


#' Load and annotate GWAS summary stats
#' \code{load.annotate.genes.from.GWAs.hits} Given common summary stats
#' @param file Character, path to the file were the GWAS summary stats are located
#' @param source Character, either "GWASCatalog" or "Immunobase"
#' @param gene.window Numeric, genomic window
#' @param gene.info A data.frame with genomic coordinates for genes. For reference please look at data(gene.info)
#' @return return ---
load.annotate.genes.from.GWAs.hits <- function(file,
                                               source="Immunobase",
                                               gene.window= 500000,
                                               gene.info){

  if(source== "Immunobase"){

  }

  if(source=="GWASCatalog"){
    #gwas.cat <- fread(file, )
  }

  return()

}

