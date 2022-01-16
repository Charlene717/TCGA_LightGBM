## Create example data.frame with genomic positions
positions <- data.frame(chromosome = c(1,1),
                        start = c(62920, 16863509),
                        end = c(16855942, 16932568))

## load package and set biomaRt dataset
library(biomaRt)
ensembl = useEnsembl(biomart='ensembl', 
                     dataset="hsapiens_gene_ensembl") 

## run query
results <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                 filters = c("chromosome_name", "start", "end"),
                 values = list(positions[,1], positions[,2], positions[,3]),
                 mart = ensembl)

## TCGA_SNP
library(dplyr)
TCGA_SNP_Raw <- read.csv(
  "D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data(Raw)_LUAD/LUAD_mc3.txt",
  sep = "\t",
  header = FALSE)

TCGA_SNP_Raw2 <- TCGA_SNP_Raw
colnames(TCGA_SNP_Raw2) <- TCGA_SNP_Raw2[1,] 
TCGA_SNP_Raw2 <- TCGA_SNP_Raw2[-1,]
TCGA_SNP_Raw3 <- mutate(TCGA_SNP_Raw2,gene_effect = paste0(gene,"_",effect))
TCGA_SNP_Raw3.df <- tidyr::spread(TCGA_SNP_Raw3,gene_effect,DNA_VAF, fill = 0)




TCGA_SNP_Raw3S <- TCGA_SNP_Raw3[1:2000,]
TCGA_SNP_Raw3S.df <- tidyr::spread(t1,gene,effect)

t1 <- TCGA_SNP_Raw3S %>% group_by(sample,gene) %>% slice(1)
TCGA_SNP_Raw3S.df <- tidyr::spread(t1,gene,effect)




load("pollution")
TTT <- tidyr::spread(pollution, size, amount)
unique(TCGA_SNP_Raw3S$effect)
