# needed packages
install.packages("BiocManager")#one time
install.packages("matrixTests")#one time
BiocManager::install("genefilter")#one time
library(BiocManager)#one time
library(matrixTests)#each time you open your R
library(genefilter)#each time you open your R

#loading the gene expression data for both normal and tumer data
norm = as.matrix(read.csv("lusc-rsem-fpkm-tcga_paired.csv", row.names=1))
tum = as.matrix(read.csv("lusc-rsem-fpkm-tcga-t_paired.csv", row.names=1))

#explore the data, dim function give you the dimension of your data; how 
#many columns(samples) do you have and how many row(genes) do you have
dim(norm)
dim(tum)

#bind the two cases (tumor nd normal) together using column bind function
data = cbind(tum,norm)

#explore the dimension of the whole data
dim(data)

#explore if is there any missing expression value (empty cell)
sum(is.null(data))

#explore the data distribution using the histogram plot
hist(data, col = "orange", main="Histogram")
#scaling the data using log2 transformation to better visulization
# we use (+1) to avoid the infinity character when we log zero valus 
hist(log2(data+1), col = "orange", main="Histogram")

#filter low count genes which have row mean lower than 1
data=data[rowMeans(data) > 1,]

####calculating the fold change####

#caculate the logged mean for each group
tum.mean = apply((log2(data+1))[,1:50], 1, mean)
norm.mean = apply((log2(data+1))[,51:dim(data)[2]], 1, mean)

#calculate the fold change by taking the difference between the two means
#the difference between logged means equl to the fold change insted of using
#division of unlogged data
fold=tum.mean-norm.mean

#view the distribution of the fold change
hist(fold)

####doing the differential expression statstical testing

#creat a phenotype table as its rows contain a phenotype ethier tumor or
#normal corrseponding to the columns in the expression data; so as we know 
#that the first 50 column in the expression data are tumor so the first 50
#row in the phenotype will be labeld tum and the other 50 norm
phenotype = as.data.frame(factor(rep(c("tum","norm"), c(50, 50))))
colnames(phenotype)= "grouping"

#making the hypothesis testing using T test for each row(gene)
t=rowttests(data,phenotype$grouping)

#correct the T test p value using the FDR method
p.adj=p.adjust(t$p.value, method = "fdr")

#save the result in a dataframe contain the fold change and p adjusted valu
result=as.data.frame(cbind(fold , p.adj))

#chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 2
res.deg=result[result$p.adj < 0.05 & abs(result$fold)>2,]

#export the Degs into your current folder for further analysthis
write.csv(as.matrix(res.deg),file="res.degs.csv", quote=F,row.names=T)

#using nonprametric test as wilcoxon test if the data isn't normally distributed
w=row_wilcoxon_twosample(data[,1:50],data[,51:dim(data)[2]])
