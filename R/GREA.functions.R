
#' get.sampled.genes
#' \code{get.sampled.genes} Given a set of \code{input.genes} (which should be in enemsebl) this function wil return a list of \code{n.genes.sets} of genes which represent the same number of genomic regions as the \code{input.genes}

#' @param input.genes Character vector of ensembl ids
#' @param gene.info A data.frame object with the genomic coordinates of the genes used in the study.
#' @param binning.widow Window size which will define the genomic regions and bining the genes. Default is 125 Kb to each side.
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


  input.genes.gene.info <- gene.info.binned[match(input.genes, as.character(gene.info.binned$ensembl)),]
  ## remove non-somatic genes
  input.genes.gene.info <- input.genes.gene.info[which(input.genes.gene.info$chr %in% c(1:22)),]

  #* bin all these genes
  #input.genes.binned <- bin.genes.genomic.regions(gene.pos = input.genes.gene.info, window = binning.widow)
  n.genomic.regions <- length(unique(input.genes.gene.info$bin.tag))
  n.genes <- length(input.genes)

  if(n.genomic.regions >n.genes ){stop("You cannot have more genomic regions than genes")}
  #permuted.gene.sets <- list()
  #permuted.gene.sets <- vector(mode = "list", length = n.genes.sets)

  permuted.gene.sets <- foreach::foreach(j=1:n.genes.sets) %dopar% {

    j.bin <- sample(unique(gene.info.binned$bin.tag), size = n.genomic.regions, replace = FALSE)
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
    j.binned.gene.info[,"ensembl"]
  }
}

#' get.genes.within.loci
#' \code{get.genes.within.loci} Given a set of \code{positions}(genomic coordinates), will return a data.frame with the genes within the input \code{positions}

#' @param positions A data.frame with genomic coordinates.
#' @param position.colnames Character vector defining the colnames used for defining the genomic coordinates within the \code{positions} data.frame.
#' @param gene.window Numeric, which will deifine the genomic window that will be used to get neighbouring genes given a certain genomic position defined by \code{positions}.
#' @param gene.info
#' @param pVal.col Character, if provided, then an extra column to the returned data-frame will be added with the pValue of the SNP associated.
#' @return A dataframe

get.genes.within.loci <- function(positions,
                                  position.colnames = c("CHR", "BP", "SNP"),
                                  gene.window=500000,
                                  gene.info, pVal.col=NULL){

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
#' @param gene.info
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


############
GREA.test <- function(sampled.genes.sets, query.gene.set, genomic.region.gene.sets, plot=TRUE){

  cat("\n[INFO]\t A total of ", length(sampled.genes.sets), "sampled gene sets will be used to generate an empirical distribution")
  cat("\n[INFO]\t ", length(genomic.region.gene.sets), "genomic regions")

  n.unique.grs <- sapply(genomic.region.gene.sets, function(x){length(unique(x[,"SNP"]))})

  # count number of SNPs have at least 1 neighboring gene is present in query.gene.set
  overlap.per.gr <- sapply(genomic.region.gene.sets,
                              function(x){length(unique(x[,"SNP"][which(x[,"ensembl"] %in% query.gene.set)]))
                              })
  # percentage of count
  perc.per.gr <- overlap.per.trait/n.unique.grs


  # get genes for final report
  overlapped.ensembl.per.gr <- lapply(genomic.region.gene.sets, function(x){query.gene.set %in% x[,"ensembl"]})
  overlapped.symbol.per.gr <- gene.info$gene_name[match(overlapped.ensembl.per.gr, gene.info$gene_id)]
  overlapped.symbol.per.gr.vector <- sapply(overlapped.symbol.per.trait, function(x){paste(x, collapse = "/")})

  ## count number of genes within random sets overlaps with the GWAS "genes"
  null.overlap.per.gr <- lapply(sampled.genes.sets, function(y){
                              sapply(genomic.region.gene.sets, function(x){
                                  length(unique(x[,"SNP"][which(x[,"ensembl"] %in% y )]))
                          })})
  ## get it in percentage
  null.perc.overlap.per.gr <- lapply(null.overlap.per.gr, function(x){x/n.unique.grs})
  null.perc.overlap.per.gr <- do.call(rbind, null.perc.overlap.per.gr)
  rownames(null.perc.overlap.per.gr) <- NULL

  enrich.pvals <- sapply(1:length(genomic.region.gene.sets),
                         function(x){pvalRight(x= perc.per.trait[x], dist = null.perc.overlap.per.gr[,x])})
  #depleted.pvals <- sapply(1:length(full.gwas.genes.list),
  #                        function(x){pvalLeft(x= perc.per.trait[x], dist = null.perc.overlap.per.trait[,x])})

  #Summary report
  GREA.summary <- data.frame(trait= names(perc.per.trait),
                                          perc.per.trait= perc.per.trait,
                                          enrich.pvals=enrich.pvals,
                                          #depleted.pvals=depleted.pvals,
                                          mean.permutation.overlap= apply(null.perc.overlap.per.gr, 2, mean),
                                          SD.permutation.overlap= apply(null.perc.overlap.per.gr, 2, sd),
                                          overlappedGenes=overlapped.symbol.per.gr.vector
                                          )

  return(GREA.summary)
}




############
load.annotate.genes.from.GWAs.hits <- function(file,  source="Immunobase", gene.window= 500000, gene.info){

  if(source== "Immunobase"){

  }

  if(source=="GWASCatalog"){
    #gwas.cat <- fread(file, )
  }

  return()

}

