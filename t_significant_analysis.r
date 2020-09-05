#Rscript
#Author: Yunlong Ma
#2018-06-08
#e-mail: glb-biotech@zju.edu.cn

##read trastuzumab expression data
setwd("C:\\Users\\Administrator\\Desktop\\")
data <- read.table("test_Gene_expression_values.txt",header=TRUE)
head(data)
data1 <- data[,c(-1,-2)]
groupA<- names(data1)[1:51]
groupB<- names(data1)[52:114]
groupC<-names(data1)[115:156]
row.names(data1) <- data[,1] #use the first column (gene ID) to name each row 
head(data1)

#extract all probe in un-treated_trastuzumab(groupA) and treated_trastuzumab(groupB) samples
trastuzumab_expr <- t(data1[, c(groupA, groupB)]) ##extract the column of groupA and groupB, and use t() function to transpose 
group <- c(rep("groupA",51), rep("groupB", 63)) ##generate repeated names
rownames(trastuzumab_expr) <- group
head(trastuzumab_expr)
class(trastuzumab_expr)
class(group)

##perform t test for each probe
Test <- function(x)
{
  T_res <- t.test(x~group)
  stat <- T_res[[1]]
  Pvalue <- T_res[[3]]
  FoldChange <- 2^(T_res[[5]][1]-T_res[[5]][2])
  LogFoldChange <- log(FoldChange,2)
  Test_summary <- c(stat, Pvalue, LogFoldChange,FoldChange)
  return(Test_summary)
}

test_res <- as.data.frame(apply(trastuzumab_expr, 2, Test)) ##apply(x, MARGIN, FUN),MARGIN=1 represents row，MARGIN=2 represents column
rownames(test_res) <- c("t_stat", "Pvalue", "logFoldChange","FoldChange")
test_res[5,] <- p.adjust(test_res[2,], method = "fdr")
test_res[6,] <- p.adjust(test_res[2,], method = "bonferroni")
test_res <- as.data.frame(t(test_res))
names(test_res) <- c("t_stat", "Pvalue", "logFoldChange","FoldChange","FDR","Bonferroni")
ID <- read.csv("Gene_annotation.csv", header=T) ##load gene annotation files
test_res$GeneName <- ID$GeneName[match(rownames(test_res), ID$ProbeID)]##annotate the probe name
write.csv(test_res, file="t-Final-result-33.csv")

#---------END---------------------------------



