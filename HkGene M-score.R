#Housekeeping gene M-score calculation
#include the libraries needed to execute code
library(rJava)
library(xlsx)
library(reshape2)
library(matrixStats)
library(NormqPCR)
library(gplots)
library(ggplot2)
library(Rtsne)
library(limma)
library(plyr)
library(stringr)
library(MASS) 
library(psych)
library(readxl)

options(scipen=999)

# commonly used housekeeping genes.  most studies use four of them.  please check to make sure the list is correct
hklist = c('ACTB','GAPDH', 'PPIA', 'B2M', "RNA18S", "RNA18S5",'PGK1','RPL13A')
###### MODIFY number of housekeeping genes to use - exclude one if it is more variable than others and only use 3
num_hk = 3

#read in well data
assay.dat = read.table(paste0(path, "/", study, "_WellData.txt"), header = TRUE, sep="\t", stringsAsFactors = FALSE)
assay.dat  = assay.dat [, c('Assay','Sample','Ct')]
#tt <- unique(assay.dat$Assay)
assay.dat 
assay.dat
ct = dcast(assay.dat, Assay~Sample, value.var="Ct")
colnames(ct) = gsub(" ","",colnames(ct))
rownames(ct) = ct$Assay
#ct = ct[,-1]

ct = ct[,match(meta$Waf.ID,colnames(ct))]
# dimension to understand number of genes and number of observations
dim(ct)

ct.na = is.na(ct)

#number of genes with any missing ct values	   
na.rows = apply(is.na(ct), 1, any)
sum(na.rows)

#number of genes with all missing ct values
all.na = apply(is.na(ct), 1, all)
sum(all.na)

#remove genes whose ct values are all missing
if(sum(all.na>0)) { 
  ct.nas = ct.na[-which(all.na),]
  na.genes = rownames(ct)[which(all.na)]
  ct = ct[-which(all.na),]
}  

#percentage missing ct
mean(is.na(ct))

#hklist = c('ACTB','GAPDH', 'PPIA', 'B2M', "RNA18S", "RNA18S5",'PGK1','RPL13A')


# most studies have four house keeping genes 
hskgenes = ct[rownames(ct) %in% hklist,]

#hskgenes = ct[rownames(ct) %in% tt,]

##Check the missing in ct values for housekeeping genes. prefer 0 missing
rowSums(is.na(hskgenes))

# replace missing values with row max + 1. please note only use this imputation when the percentage missing should be small
k = which(is.na(ct),arr.ind=TRUE)
ct[k] = (rowMaxs(as.matrix(ct),na.rm=TRUE)+1)[k[,1]]

genes = rownames(ct)

hskgenes = ct[rownames(ct) %in% hklist,]

#M-score calculation
######################################
hk.rank = data.frame(hk = hklist)
# # #get the meanM for each HK genes.
hk.M = data.frame(hk = hklist)
# # 
# # #create qPCR object with expression values
test = get.qPCR(exprs = hskgenes, phenoData = new("AnnotatedDataFrame"), notes = "", verbose = FALSE)
# # #get the ranking of each housekeeping genes  
hk.rank[,"rank"] = match(hklist, selectHKs(test, method = "geNorm", Symbols = featureNames(test), minNrHK = 2, trace = FALSE, log = TRUE)$ranking)
# # #get the average M values of each genes  
mvec = selectHKs(test, method = "geNorm", Symbols = featureNames(test), minNrHK = 2, trace = FALSE, log = TRUE)$meanM
# # #  
hk.M[nas(hk.rank[,'rank']),'M'] = fill(mvec, len=length(hklist))
# # 
hk.M
hk.rank
#########################################
