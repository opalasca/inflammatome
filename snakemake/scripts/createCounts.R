#!/usr/bin/env Rscript

datadir = "../../data/"

library("DESeq2")
library("tximport")

args = commandArgs(trailingOnly=TRUE)
dataset=args[1]

salmdir = paste0(datadir,"/datasets/",dataset,"/salmon/")

tx2gene <- read.csv(paste(datadir,"ref/tr2gene/Tx2gene.txt",sep=""),sep="\t",header=FALSE)

get_files_salmon<-function(dir,sample_list){
  files <- file.path(dir, sample_list, "quant.sf")
  names(files) <- sample_list
  all(file.exists(files))
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  names(txi) 
  return(txi)
}

samples=list.dirs(path = salmdir, full.names = FALSE, recursive = FALSE)
print(samples)

coldata = as.data.frame(samples)
colnames(coldata)=c("sample")
coldata$condition="ctrl"
#coldata[c(1:n),]$condition="trt"
rownames(coldata)=coldata$sample
coldata$condition=as.factor(coldata$condition)
txi.g<-get_files_salmon(salmdir, samples)
tpm=txi.g$abundance
rownames(coldata)=colnames(txi.g$counts)
coldata<-coldata[colnames(txi.g$counts),]
dds <- DESeqDataSetFromTximport(txi.g, colData = coldata, design = ~ 1)
#dds$condition <- relevel(dds$condition, ref = "ctrl")
cts.raw=counts(dds, normalized=F)
write.table(cts.raw,file=paste(datadir,"/datasets/",dataset,"/counts.tsv",sep=""),sep="\t",quote = F,row.names=T)
