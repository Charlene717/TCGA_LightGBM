memory.limit(150000)

TCGA_Meth <- read.csv(
  "D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data/HumanMethylation450.xena",
  sep = "\t",
  header = FALSE
)
#row.names(TCGA_Meth) <- TCGA_Meth[,1]
TCGA_Meth2 <- as.data.frame(t(TCGA_Meth))

TCGA_CNV <- read.csv(
  "D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data/Gistic2_CopyNumber_Gistic2_all_data_by_genes",
  sep = "\t",
  header = FALSE
)
#row.names(TCGA_CNV) <- TCGA_CNV[,1]
TCGA_CNV2 <- as.data.frame(t(TCGA_CNV))

TCGA_SNP <- read.csv(
  "D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data/SNP_and_INDEL.xena",
  sep = "\t",
  header = FALSE
)
#row.names(TCGA_SNP) <- TCGA_SNP[,1]
TCGA_SNP2 <- as.data.frame(t(TCGA_SNP))


TCGA_RNA <- read.csv(
  "D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data/RNASeqV2.geneExp.xena",
  sep = "\t",
  header = FALSE
)
#row.names(TCGA_RNA) <- TCGA_RNA[,1]
TCGA_RNA2 <- as.data.frame(t(TCGA_RNA))

colnames(TCGA_CNV) <- TCGA_CNV[1,]
colnames(TCGA_SNP) <- TCGA_SNP[1,]
colnames(TCGA_Meth) <- TCGA_Meth[1,]
library(dplyr)
TCGA_Com <- left_join(TCGA_CNV, TCGA_SNP, by="Samples")
TCGA_Com <- left_join(TCGA_Com, TCGA_Meth, by="Samples")

