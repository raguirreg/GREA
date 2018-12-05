library(data.table)
col.data <- data.frame(fread("/Users/raulaguirre/NewBiologyCourse/RNASeq_1/Data/colData.csv"), row.names = 1)
save(col.data, file = "./data/GSE69549_colData.RData")

count.table <- data.frame(fread("/Users/raulaguirre/NewBiologyCourse/RNASeq_1/Data/countTable.csv"), row.names=1)
save(count.table, file = "./data/GSE69549_countTable.RData")
