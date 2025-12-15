#save plink1.9 and all data file in a working file
#from window: cmd;
#from mac: terminal

#go to the working file
cd Data_ST721_xy

plink --merge-list merged.txt --make-bed --out allraw

#missing genotype
plink --bfile allraw --geno 0.05 --make-bed --out cleaned_mydata_geno

#HWD
plink --bfile cleaned_mydata_geno --hwe 0.00000001 --make-bed --out cleaned_mydata_hwe


#Exclude MAF<5%
plink --bfile cleaned_mydata_hwe --maf 0.05 --make-bed --out cleaned_mydata_snp


#individual

plink --bfile cleaned_mydata_snp --check-sex --make-bed --out cleaned_mydata_sex

plink --bfile cleaned_mydata_sex --genome --min 0.05 --make-bed --out cleaned_mydata_relatedness

plink --bfile cleaned_mydata_relatedness --mind 0.1 --make-bed --out all_cleaned



#use R
#set working file directory
setwd("C:/Users/siyaq/Data_ST721")
rm(list=ls())
#Read gene expression data chr19.txt (make sure this txt file is under this directory)
g_ex<-read.table("chr19.txt",header=TRUE)
which(g_ex$GeneID=="ILMN_1729894")

ind_gex<-data.frame(colnames(g_ex)[-1])

#1. get the new fam file change the 6th column
Y<-t(g_ex[which(g_ex$GeneID=="ILMN_1729894"),][-1])

Y<-cbind(ind_gex,Y)



colnames(ind_gex)<-c("IID")
colnames(Y)<-c("IID","ABAC")


SNP_fam<-read.table("all_cleaned.fam",header=FALSE)
colnames(SNP_fam)<-c("FID","IID","PAT","MAT","SEX","PHENOTYPE")

length(unique(SNP_fam$FID))
#139
length(unique(SNP_fam$IID))
#209


#merge with Y
Y_mrg<-merge(SNP_fam, Y, by = "IID", all=TRUE)
all_cleaned.fam<-Y_mrg[,c(2,1,3:5,7)]
colnames(all_cleaned.fam)<-c("FID","IID","PAT","MAT","SEX","PHENOTYPE")

write.table(all_cleaned.fam,"all_cleaned.fam",sep = "\t", col.names=FALSE, row.names = FALSE)

############################
# another phenotype
#1. get the raw fam file, create phenotype
Y<-t(g_ex[c(which(g_ex$GeneID=="ILMN_1729894"),which(g_ex$GeneID=="ILMN_1743205")),])[-1,]

Y_1=rownames(Y)
Y<-cbind(as.numeric(Y[,1]),as.numeric(Y[,2]))
Y<-data.frame(Y_1,Y)

colnames(Y)<-c("IID","ABAC1","ABAC2")


allraw<-read.table("allraw.fam",header=FALSE)
colnames(allraw)<-c("FID","IID","PAT","MAT","SEX","PHENOTYPE")

length(unique(allraw$FID))
#140
length(unique(allraw$IID))
#210


phenotype<-merge(allraw, Y, by = "IID")
phenotype<-phenotype[,c(2,1,7:8)]
sum(is.na(phenotype$ABAC1))
sum(is.na(phenotype$ABAC2))


write.table(phenotype,"phenotype.txt",sep = "\t", col.names=TRUE, row.names = FALSE)

#2. get the FID and IID for the 198 persons

keep.txt<-merge(ind_gex, SNP_fam, by = "IID", all=FALSE)
keep.txt<-keep.txt[,c(2,1)]
write.table(keep.txt,"keep.txt",sep = "\t", row.names = FALSE)

#use plink
#after update the phynotype data
#keep 198 individuals
plink --bfile all_cleaned --keep keep.txt --make-bed --out all_198ind


#QQplot
#Run linear regression to get all p values  
# after all QC
plink --bfile all_198ind --linear 
# the  file produced is: plink.assoc.linear
#use this file to produce qq-plot in R 
# plink.assoc.linear


#use R 
#get qq-plot
#if qqman package is not installed;
#first run  the  commond in R: install.packages("qqman")
library(qqman)
plink.assoc<-read.table("plink.assoc.linear",header=TRUE)

qq(plink.assoc$P)


# run PCA 
plink --bfile all_198ind --pca --out PCA_rslt 



#look at the PCA file in R 
setwd("C:/Users/siyaq/Data_ST721_xy")
rm(list=ls())
PC_vc<-read.table("PCA_rslt.eigenvec",header=FALSE)
colnames(PC_vc)<-c("FID","IID",paste0("PC",c(1:20)))

#find the fid and iid for each individual and cbind with a column as HIS
CEU<-read.table("hapmap_CEU_r23a_filtered.fam",header=FALSE)
CEU<-cbind(CEU[,1:2],c("CEU"))
colnames(CEU)<-c("FID","IID","HIS")

YRI<-read.table("hapmap_YRI_r23a_filtered.fam",header=FALSE)
YRI<-cbind(YRI[,1:2],c("YRI"))
colnames(YRI)<-c("FID","IID","HIS")

JC<-read.table("hapmap_JPT_CHB_r23a_filtered.fam",header=FALSE)
JC<-cbind(JC[,1:2],c("JPT_CHB"))
colnames(JC)<-c("FID","IID","HIS")

HIS<-rbind(CEU,YRI,JC)
PC_vc<-merge(PC_vc,HIS, by=c("FID","IID"), all.x=TRUE)
plot(PC_vc$PC1,PC_vc$PC2)

plot(PC_vc$PC1, PC_vc$PC2, col = c("red", "blue", "green")[as.factor(PC_vc$HIS)])

write.table(PC_vc[,1:22],"pcs.txt",sep = "\t", col.names=TRUE, row.names = FALSE)


##PLINK--run association analysis--chr19
plink --bfile all_198ind --chr 19 --make-bed --out chr19
plink --bfile chr19 --linear --covar pcs.txt


chr19_assoc<-read.table("plink.assoc.linear",header=TRUE)

qq(chr19_assoc$P)
as.data.frame(table(chr19_assoc$CHR))
manhattan(chr19_assoc)

#run association analysis whole genome
plink --bfile all_198ind --linear --covar pcs.txt


chr19_assoc<-read.table("plink.assoc.linear",header=TRUE)

qq(chr19_assoc$P)
as.data.frame(table(chr19_assoc$CHR))
manhattan(chr19_assoc)