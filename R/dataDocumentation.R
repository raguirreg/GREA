#' Gene expression from CD4+ T cells from healthy and Celiac patients from GSE69549
#'
#' A data.frame containing RNASeq based expression levels from GSE69549
#' @docType data
#' @keywords datasets, expression
#' @name count.table
#' @usage data(count.table)
#' @format A data.frame with 23228 genes (by rows) and 74 samples (by column)
"count.table"


#' Sample information from CD4+ T cells from healthy and Celiac patients from GSE69549
#'
#' A data.frame containing sample information from GSE69549
#' @docType data
#' @keywords datasets, sample information
#' @name col.data
#' @usage data(col.data)
#' @format A data.frame with 74 rows and 4 variables
"col.data"


#' List of 7 elements where each is a data.frame with gene sets defined by proximity to GWAS hits from Immunobase studies
#'
#' The list contains gene sets defined by proximity to GWAS hits for these diseases:
#' \itemize{
#' \item "CeD" Celiac Disease
#' \item "CD" Chron's Disease
#' \item "MS" Multiple Sclerosis
#' \item "PBC" Primaire Biliaire Cholangitis
#' \item "RA" Rheumatoid Arthritis
#' \item "T1D" Rheumatoid Arthritis
#' \item "UC" Ulcerative Colitis
#' }
#'
#' Each data.frame of the list has the following information
#' \itemize{
#' \item chr
#' \item left
#' \item right
#' \item width
#' \item SNP
#' \item SNP.pos
#' }
#'
#' @docType data
#' @keywords datasets, gene sets
#' @name ib.gw.traits.genes.list
#' @usage data(ib.gw.traits.genes.list)
#' @format A list with 7
"ib.gw.traits.genes.list"

#' List of 7 elements where each is a data.frame with gene sets defined by proximity to GWAS hits from Immunobase studies (results from immuno-chip platform)
#'
#' The list contains gene sets defined by proximity to GWAS hits for these diseases:
#' \itemize{
#' \item "CeD" Celiac Disease
#' \item "AS" Ankylosing Spondylitis
#' \item "MS" Multiple Sclerosis
#' \item "PBC" Primaire Biliaire Cholangitis
#' \item "RA" Rheumatoid Arthritis
#' \item "T1D" Rheumatoid Arthritis
#' \item "PSO" Psoriasis
#' }
#'
#' Each data.frame of the list has the following information
#' \itemize{
#' \item chr
#' \item left
#' \item right
#' \item width
#' \item SNP
#' \item SNP.pos
#' }
#'
#' @docType data
#' @keywords datasets, gene sets
#' @name ib.traits.genes.list
#' @usage data(ib.gw.traits.genes.list)
#' @format A list with 7
"ib.traits.genes.list"


#' Gene information
#'
#' A data.frame with genomic positions for each gene version GRCh37.75
#' @docType data
#' @keywords datasets, sample information
#' @name gene.info
#' @usage data(gene.info)
#' @format A data.frame with 74 rows and 4 variables
"gene.info"


