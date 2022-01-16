memory.limit(150000)

## ----------------------------------------
PathName = setwd(getwd())
RVersion = "20211224"
dir.create(paste0(PathName,"/",RVersion))
## ----------------------------------------

##### Load TCGA Data #####
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
##### Combine different type of data #####

library(dplyr)
TCGA_Com <- inner_join(TCGA_RNA2, TCGA_SNP2, by="sample")
TCGA_Com <- inner_join(TCGA_Com, TCGA_CNV2, by="sample")
TCGA_Com <- inner_join(TCGA_Com, TCGA_Meth2, by="sample")

#
TCGA_Com_TP53 <- inner_join(TCGA_SNP2,TCGA_CNV2, by="sample")
TCGA_Com_TP53 <- inner_join(TCGA_Com_TP53, TCGA_Meth2, by="sample")

TCGA_RNA_TP53 <-　data.frame(sample = TCGA_RNA2[,1],TP53 = TCGA_RNA2[,colnames(TCGA_RNA2) == "TP53"])
TCGA_Com_TP53 <- inner_join(TCGA_Com_TP53, TCGA_RNA_TP53, by="sample")

TCGA_Com_TP53 <- delete.na(TCGA_Com_TP53, 0)
# row.has.na <- apply(final, 1, function(x){any(is.na(x))})

#### Correlation #####
# res <- cor.test(my_data$wt, my_data$mpg,  
#                  method = "spearman")
TCGA_Com_TP53_Ori <- TCGA_Com_TP53
row.names(TCGA_Com_TP53) <- TCGA_Com_TP53[,1]
TCGA_Com_TP53 <- TCGA_Com_TP53[,-1]
#NG# TCGA_Com_TP53 <- as.data.frame(as.numeric(unlist(TCGA_Com_TP53)))
TCGA_Com_TP53 <- as.data.frame(lapply(TCGA_Com_TP53, as.numeric))
row.names(TCGA_Com_TP53) <- TCGA_Com_TP53_Ori[,1]
  
TCGA_Com_TP53_cor.lst <- as.data.frame(cor(TCGA_Com_TP53[1:(ncol(TCGA_Com_TP53)-1)],TCGA_Com_TP53[ncol(TCGA_Com_TP53)],method = "spearman"))
TCGA_Com_TP53_cor.lst <- data.frame(Genes= row.names(TCGA_Com_TP53_cor.lst),TP53 = TCGA_Com_TP53_cor.lst)
TCGA_Com_TP53_cor_Can.lst <- TCGA_Com_TP53_cor.lst[abs(TCGA_Com_TP53_cor.lst$TP53) >= 0.2,]  
# TCGA_Com_TP53_cor_Can.lst <- data.frame(TP53=TCGA_Com_TP53_cor.lst[abs(TCGA_Com_TP53_cor.lst$TP53) >= 0.1,])
# TCGA_Com_TP53_cor_Can.lst <- as.data.frame(TCGA_Com_TP53_cor.lst[abs(TCGA_Com_TP53_cor.lst$TP53) >= 0.1,])
TCGA_Com_TP53_cor_Can.lst <- as.data.frame(TCGA_Com_TP53_cor_Can.lst[!is.na(TCGA_Com_TP53_cor_Can.lst$TP53),])
TCGA_Com_TP53_cor <- TCGA_Com_TP53[,colnames(TCGA_Com_TP53) %in% c(TCGA_Com_TP53_cor_Can.lst$Genes,"TP53")]

write.csv(TCGA_Com_TP53_cor ,
           file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor02.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)


set.seed(1) # Fix the seed
TCGA_Com_TP53_cor_Test <- sample_frac(TCGA_Com_TP53_cor,0.1)
TCGA_Com_TP53_cor_Train <- anti_join(TCGA_Com_TP53_cor,TCGA_Com_TP53_cor_Test)

TCGA_Com_TP53_cor_Test2 <- TCGA_Com_TP53_cor_Test
TCGA_Com_TP53_cor_Test2[,ncol(TCGA_Com_TP53_cor_Test2)] <- 0


write.csv(TCGA_Com_TP53_cor_Test2 ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor02_Test.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)
write.csv(TCGA_Com_TP53_cor_Test ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor02_TestAns.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)

write.csv(TCGA_Com_TP53_cor_Train ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor02_Train.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)



# cor 0.1
TCGA_Com_TP53_cor01_Can.lst <- TCGA_Com_TP53_cor.lst[abs(TCGA_Com_TP53_cor.lst$TP53) >= 0.1,]  
TCGA_Com_TP53_cor01_Can.lst <- as.data.frame(TCGA_Com_TP53_cor01_Can.lst[!is.na(TCGA_Com_TP53_cor01_Can.lst$TP53),])


TCGA_Com_TP53_cor01 <- TCGA_Com_TP53[,colnames(TCGA_Com_TP53) %in% c(TCGA_Com_TP53_cor01_Can.lst$Genes,"TP53")]

# cor 0.3
TCGA_Com_TP53_cor03_Can.lst <- TCGA_Com_TP53_cor.lst[abs(TCGA_Com_TP53_cor.lst$TP53) >= 0.3,]  
TCGA_Com_TP53_cor03_Can.lst <- as.data.frame(TCGA_Com_TP53_cor03_Can.lst[!is.na(TCGA_Com_TP53_cor03_Can.lst$TP53),])


TCGA_Com_TP53_cor03 <- TCGA_Com_TP53[,colnames(TCGA_Com_TP53) %in% c(TCGA_Com_TP53_cor03_Can.lst$Genes,"TP53")]




                                                         
#####  cor 0.15  #####
TCGA_Com_TP53_cor015_Can.lst <- TCGA_Com_TP53_cor.lst[abs(TCGA_Com_TP53_cor.lst$TP53) >= 0.15,]  
TCGA_Com_TP53_cor015_Can.lst <- as.data.frame(TCGA_Com_TP53_cor015_Can.lst[!is.na(TCGA_Com_TP53_cor015_Can.lst$TP53),])


TCGA_Com_TP53_cor015 <- TCGA_Com_TP53[,colnames(TCGA_Com_TP53) %in% c(TCGA_Com_TP53_cor015_Can.lst$Genes,"TP53")]

write.csv(TCGA_Com_TP53_cor015 ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor015.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)

set.seed(1) # Fix the seed
TCGA_Com_TP53_cor015_Test <- sample_frac(TCGA_Com_TP53_cor015,0.1)
TCGA_Com_TP53_cor015_Train <- anti_join(TCGA_Com_TP53_cor015,TCGA_Com_TP53_cor015_Test)

TCGA_Com_TP53_cor015_Test2 <- TCGA_Com_TP53_cor015_Test
TCGA_Com_TP53_cor015_Test2[,ncol(TCGA_Com_TP53_cor015_Test2)] <- 0


write.csv(TCGA_Com_TP53_cor015_Test2 ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor015_Test.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)
write.csv(TCGA_Com_TP53_cor015_Test ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor015_TestAns.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)

write.csv(TCGA_Com_TP53_cor015_Train ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor015_Train.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)


## Try
TCGA_RNA_TP53 <-　data.frame(sample = TCGA_RNA2[,1],TP53 = TCGA_RNA2[,colnames(TCGA_RNA2) == "TOP2A"])
TCGA_Com_TP53 <- inner_join(TCGA_Com_TP53, TCGA_RNA_TP53)
TCGA_Com_TP53_cor.lst.pearson <- as.data.frame(cor(TCGA_Com_TP53[1:(ncol(TCGA_Com_TP53)-1)],TCGA_Com_TP53[ncol(TCGA_Com_TP53)],method = "pearson"))
TCGA_Com_TP53_cor.lst.pearson <- data.frame(Genes= row.names(TCGA_Com_TP53_cor.lst.pearson),TP53 = TCGA_Com_TP53_cor.lst.pearson)
TCGA_Com_TP53_cor_Can.lst.pearson <- TCGA_Com_TP53_cor.lst.pearson[abs(TCGA_Com_TP53_cor.lst.pearson$TP53) >= 0.3,] 
TCGA_Com_TP53_cor_Can.lst.pearson <- as.data.frame(TCGA_Com_TP53_cor_Can.lst.pearson[!is.na(TCGA_Com_TP53_cor_Can.lst.pearson$TP53),])


## Try Cor range
TCGA_Com_TP53_cor.lst <- data.frame(TCGA_Com_TP53_cor.lst,Abs = abs(TCGA_Com_TP53_cor.lst$TP53))
TCGA_Com_TP53_cor.lst.sort <- arrange(TCGA_Com_TP53_cor.lst,desc(Abs))

TCGA_Com_TP53_cor.lst.sort2 <- data.frame(Seq = seq(1:nrow(TCGA_Com_TP53_cor.lst.sort)),TCGA_Com_TP53_cor.lst.sort)
TCGA_Com_TP53_cor.lst.TP53 <- TCGA_Com_TP53_cor.lst.sort2[TCGA_Com_TP53_cor.lst.sort2$Genes == "TP53.x",]


TCGA_Com_TP53_cor_Can.lst <- TCGA_Com_TP53_cor.lst[abs(TCGA_Com_TP53_cor.lst$TP53) >= 0.2,]  

##### Candidate 25000 ####
TCGA_Com_TP53_cor.lst <- data.frame(TCGA_Com_TP53_cor.lst,Abs = abs(TCGA_Com_TP53_cor.lst$TP53))
TCGA_Com_TP53_cor.lst.sort <- arrange(TCGA_Com_TP53_cor.lst,desc(Abs))
TCGA_Com_TP53_cor.lst.sort2 <- data.frame(Seq = seq(1:nrow(TCGA_Com_TP53_cor.lst.sort)),TCGA_Com_TP53_cor.lst.sort)

TCGA_Com_TP53_cor_Can_25000.lst <- TCGA_Com_TP53_cor.lst.sort2[1:25000,]                                                   

TCGA_Com_TP53_cor_Can_25000.lst <- as.data.frame(TCGA_Com_TP53_cor_Can_25000.lst[!is.na(TCGA_Com_TP53_cor_Can_25000.lst$TP53),])


TCGA_Com_TP53_cor25000 <- TCGA_Com_TP53[,colnames(TCGA_Com_TP53) %in% c(TCGA_Com_TP53_cor_Can_25000.lst$Genes,"TP53")]

write.csv(TCGA_Com_TP53_cor25000 ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor25000.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)

set.seed(1) # Fix the seed
TCGA_Com_TP53_cor25000_Test <- sample_frac(TCGA_Com_TP53_cor25000,0.1)
TCGA_Com_TP53_cor25000_Train <- anti_join(TCGA_Com_TP53_cor25000,TCGA_Com_TP53_cor25000_Test)

TCGA_Com_TP53_cor25000_Test2 <- TCGA_Com_TP53_cor25000_Test
TCGA_Com_TP53_cor25000_Test2[,ncol(TCGA_Com_TP53_cor25000_Test2)] <- 0


write.csv(TCGA_Com_TP53_cor25000_Test2 ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor25000_Test.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)
write.csv(TCGA_Com_TP53_cor25000_Test ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor25000_TestAns.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)

write.csv(TCGA_Com_TP53_cor25000_Train ,
          file = paste0(PathName,"/",RVersion,"/TCGA_Com_TP53_cor25000_Train.csv"),
          #sep = ",",
          quote = F,
          row.names = F
)
