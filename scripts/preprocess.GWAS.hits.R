########################################################################
### Load and process Immunobase summary statistics for GREA package.
### Date:
### Author:
########################################################################

library(data.table)
source(file = "./R/GREA.functions.R")
load(file = "./data/gene.info.GRCh37.75.RData")
#colnames(gene.info)[1:3] <- c("chr","left","right")
#save(gene.info, file = "./data/gene.info.GRCh37.75.RData")

gene.binning.window <- 125000

##### Immunobase GWAS hits
# files were downloaded from https://www.immunobase.org/downloads/protected_data/iChip_Data/

ib.files.path <- "/Users/raulaguirre/IELs_celiac/GWAS_enrichment/immunoBase_GWASenrichment/immunoChip/iChipData/"
ib.files <- list.files(ib.files.path, full.names=TRUE, pattern = ".tab$")
#ib.files.names <- sapply(basename(ib.files), function(x){unlist(strsplit(split = "_", x))[4]})
ib.files.names <- c("AS", "CeD", "MS", "PBC", "PSO", "RA", "T1D")

#load data
ib.trait.list <- lapply(ib.files, fread, data.table=FALSE)
names(ib.trait.list) <- ib.files.names

GWAS.threshold.pVal <-  5e-05
# select only SNPs with a given threshold, since we will save these for
ib.trait.list <- lapply(ib.trait.list, function(x){x[which(as.numeric(x[,"PValue"]) < GWAS.threshold.pVal),]})

## remove all SNPs in non-somatic chromosomes
ib.trait.list <- lapply(ib.trait.list, function(x){x[x[,"Chr"] %in% 1:22,]})

### Prune SNPs based on position and pvalue
ib.trait.prunned.list <- list()
for(i.gwas in 1:length(ib.trait.list)){
## make chromosome and position columns as numeric for all traits.
ib.trait.list[[i.gwas]][,"Position"] <- as.numeric(ib.trait.list[[i.gwas]][,"Position"])
ib.trait.list[[i.gwas]][,"Chr"] <- as.numeric(ib.trait.list[[i.gwas]][,"Chr"])

ib.trait.prunned.list[[i.gwas]] <- prune.loci.by.proximity(window = "2.5MB",
                                                            positions = ib.trait.list[[i.gwas]],
                                                            position.colnames = c("Chr", "Position"),
                                                            prune = TRUE,
                                                            pvalCol = "PValue")
}
names(ib.trait.prunned.list) <- names(ib.trait.list)


ib.traits.genes.list <- lapply(ib.trait.prunned.list, function(x){get.genes.within.loci(positions = x,
                                                                                        position.colnames = c("Chr", "Position", "Marker"),
                                                                                        pVal.col = "PValue",
                                                                                        gene.window=gene.binning.window,
                                                                                        gene.info= gene.info)})

###########################################################################
##### Immunobase GWAS hits (genome.wide chip hits)

# files were downloaded from https://www.immunobase.org/downloads/protected_data/iChip_Data/
ib.gw.files.path <- "/Users/raulaguirre/IELs_celiac/GWAS_enrichment/immunoBase_GWASenrichment/immunoChip/GWData/"
ib.gw.files <- list.files(ib.gw.files.path, pattern = ".tab$", full.names = T)
#ib.gw.files.names <- sapply(basename(ib.gw.files),function(x){unlist(strsplit(split = "_", x))[3]})
ib.gw.files.names <- c("CeD", "CD","MS", "PBC", "RA", "T1D", "UC")

#load data
ib.gw.trait.list <- lapply(ib.gw.files, fread, data.table=FALSE)
names(ib.gw.trait.list) <- ib.gw.files.names

GWAS.threshold.pVal <-  5e-05
# select only SNPs with a given threshold, since we will save these for
ib.gw.trait.list <- lapply(ib.gw.trait.list, function(x){x[which(as.numeric(x[,"PValue"]) < GWAS.threshold.pVal),]})

## remove all SNPs in non-somatic chromosomes
ib.gw.trait.list <- lapply(ib.gw.trait.list, function(x){x[x[,"Chr"] %in% 1:22,]})

### Prune SNPs based on position and pvalue
ib.gw.trait.prunned.list <- list()
for(i.gwas in 1:length(ib.gw.trait.list)){
  ## make chromosome and position columns as numeric for all traits.
  ib.gw.trait.list[[i.gwas]][,"Position"] <- as.numeric(ib.gw.trait.list[[i.gwas]][,"Position"])
  ib.gw.trait.list[[i.gwas]][,"Chr"] <- as.numeric(ib.gw.trait.list[[i.gwas]][,"Chr"])

  ib.gw.trait.prunned.list[[i.gwas]] <- prune.loci.by.proximity(window = "2.5MB",
                                                                 ib.gw.trait.list[[i.gwas]],
                                                                position.colnames = c("Chr", "Position"),
                                                                prune = TRUE,
                                                                pvalCol = "PValue")
}
names(ib.gw.trait.prunned.list) <- names(ib.gw.trait.list)


ib.gw.traits.genes.list <- lapply(ib.gw.trait.prunned.list, function(x){get.genes.within.loci(positions = x,
                                                                                                   position.colnames = c("Chr", "Position", "Marker"),
                                                                                                    pVal.col = "PValue",
                                                                                                   gene.window=gene.binning.window,
                                                                                                   gene.info= gene.info)})


##### Save both lists and save them to incorporate them to the package and be used as base enrichments.
save(ib.traits.genes.list, file = "./data/ib.traits.genes.list.RData")
save(ib.gw.traits.genes.list, file = "./data/ib.gw.traits.genes.list.RData")
