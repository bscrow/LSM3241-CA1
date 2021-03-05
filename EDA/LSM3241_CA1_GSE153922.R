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

# d<-data.frame(vector=c(1,1,1,1,1,1,0,0,0,0,0,0), age=c(0,0,0,0,0,0,0,0,0,1,1,1), treatment=c(0,0,0,1,1,1,0,0,0,0,0,0))
# design <- model.matrix( ~ d[["treatment"]])
# design
# colnames(design)[2] <- "treatment"
# fit <- lmFit(gse153922_eset,design)
# fitted.ebayes <- eBayes(fit)
# topTable(fitted.ebayes)
# summary(decideTests(fitted.ebayes[,"treatment"],lfc=2))

design <- model.matrix( ~ 0 + gse153922_eset[['source_name_ch1']])
colnames(design) <- levels(as.factor(gse153922_eset[['source_name_ch1']]))
colnames(design) <- c("trf2", "vector", "young", "old")

contrast_matrix1 <- makeContrasts(trf2 - vector, levels=design)
contrast_matrix2 <- makeContrasts(trf2 - young, levels=design)
contrast_matrix3 <- makeContrasts(trf2 - old, levels=design)
contrast_matrix4 <- makeContrasts(vector - young, levels=design)
contrast_matrix5 <- makeContrasts(vector - old, levels=design)
contrast_matrix6 <- makeContrasts(young - old, levels=design)

fit <- lmFit(gse153922_eset,design)
fit2 <- contrasts.fit(fit,contrasts=contrast_matrix6)
#fit2 <- eBayes(fit2)
fitted.ebayes <- eBayes(fit2)
summary(decideTests(fitted.ebayes,lfc=2, p.value = 0.05))
interesting_genes <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=2)

pdf("c6.pdf")
volcanoplot(fitted.ebayes, coef=1, main=sprintf("%d features pass our cutoffs\ncoef:1, lfc:2, p<0.05",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')
dev.off()
