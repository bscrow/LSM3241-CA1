library(oligo)
library(oligoClasses)
library(tidyverse)
library(GEOquery)

gse153922 <- getGEO('GSE153922') 
gse153922 <- gse153922[[1]]

pd <- pData(gse153922)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
gse153922_celdata <- read.celfiles(paste0('C:/Users/Lenovo/Documents/RProjects/GSE153922/',pd$cel_file),phenoData=phenoData(gse153922))

gse153922_eset <- rma(gse153922_celdata)

pData(gse153922_eset)
write.csv(pData(gse153922_eset), paste0(".\\GSE153922_eset_data.csv"))


library(limma)
d<-data.frame(treatment=c("empty","empty","empty","scenescent","scenescent","scenescent", "young", "young", "young", "old", "old", "old"))
design <- model.matrix( ~ d[["treatment"]])
colnames(design)[2] <- "SCGM"
design

fit <- lmFit(gse153922_eset,design)
fitted.ebayes <- eBayes(fit)
topTable(fitted.ebayes)
summary(decideTests(fitted.ebayes[,"SCGM"],lfc=1))
