memory.limit(150000)

TCGA_Meth <- read.csv(
  "D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data_LUAD/UCSC_Xena_LUAD_Methlation_450K",
  sep = "\t",
  header = FALSE
)


TCGA_CNV <- read.csv(
  "D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data_LUAD/UCSC_Xena_LUAD_CNV",
  sep = "\t",
  header = FALSE
)


TCGA_SNP <- read.csv(
  "D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data_LUAD/UCSC_Xena_LUAD_SNP_INDEL.txt",
  sep = "\t",
  header = FALSE
)


TCGA_RNA <- read.csv(
  "D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data_LUAD/UCSC_Xena_LUAD_RNA",
  sep = "\t",
  header = FALSE
)



colnames(TCGA_CNV) <- TCGA_CNV[1,]
row.names(TCGA_CNV) <- TCGA_CNV[,1]
TCGA_CNV2 <- as.data.frame(t(TCGA_CNV))
colnames(TCGA_CNV2)[1] <- "sample"

colnames(TCGA_SNP) <- TCGA_SNP[1,]
TCGA_SNP <- TCGA_SNP[!is.na(TCGA_SNP$sample),]
row.names(TCGA_SNP) <- TCGA_SNP[,1]
TCGA_SNP2 <- as.data.frame(t(TCGA_SNP))


colnames(TCGA_Meth) <- TCGA_Meth[1,]
row.names(TCGA_Meth) <- TCGA_Meth[,1]
TCGA_Meth2 <- as.data.frame(t(TCGA_Meth))


colnames(TCGA_RNA) <- TCGA_RNA[1,]
row.names(TCGA_RNA) <- TCGA_RNA[,1]
TCGA_RNA2 <- as.data.frame(t(TCGA_RNA))

##  -----------------------------------------------------------------------------------------  ##

library(dplyr)
TCGA_Com <- inner_join(TCGA_RNA2, TCGA_SNP2, by="sample")
TCGA_Com <- inner_join(TCGA_Com, TCGA_CNV2, by="sample")
TCGA_Com <- inner_join(TCGA_Com, TCGA_Meth2, by="sample")

#
TCGA_Com_TP53 <- inner_join(TCGA_SNP2[,1:1000],TCGA_CNV2[,1:1000], by="sample")
TCGA_Com_TP53 <- inner_join(TCGA_Com_TP53, TCGA_Meth2[,1:1000], by="sample")

TCGA_RNA_TP53 <-　data.frame(sample = TCGA_RNA2[,1],TP53 = TCGA_RNA2[,colnames(TCGA_RNA2) == "TP53"])
TCGA_Com_TP53 <- inner_join(TCGA_Com_TP53, TCGA_RNA_TP53, by="sample")

TCGA_Com_TP53 <- delete.na(TCGA_Com_TP53, 0)
# row.has.na <- apply(final, 1, function(x){any(is.na(x))})

#### Correlation #####
res <- cor.test(my_data$wt, my_data$mpg,  
                 method = "spearman")


TCGA_Com_TP53_cor.lst <- acor(TCGA_Com_TP53[1:(ncol(TCGA_Com_TP53)-1)],
                          TCGA_Com_TP53[ncol(TCGA_Com_TP53)],
                          method = "spearman")


TCGA_Com_TP53_cor_Can.lst <- TCGA_Com_TP53_cor.lst[abs(TCGA_Com_TP53_cor.lst$#???) >= 0.5]

TCGA_Com_TP53_cor <- TCGA_Com_TP53[TCGA_Com_TP53$sample %in% TCGA_Com_TP53_cor_Can.lst,]









                                                         
                                                         
                                                 

