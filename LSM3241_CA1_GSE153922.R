library(GEOquery)
library(hugene20sttranscriptcluster.db)
library(limma)
library(oligo)
library(oligoClasses)
library(tidyverse)

# --------------
# Data retrieval
# --------------

# Path for the raw data files
data_path<-'C:/Users/Lenovo/Documents/RProjects/GSE153922/'

# Retrieve experimental data & RMA
gse153922 <- getGEO('GSE153922') 
gse153922 <- gse153922[[1]]

pd <- pData(gse153922)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
gse153922_celdata <- read.celfiles(paste0(data_path,pd$cel_file),phenoData=phenoData(gse153922))

# No reason to really change default settings
gse153922_eset <- rma(gse153922_celdata)

# Save metadata
write.csv(pData(gse153922_eset), paste0(".\\GSE153922_eset_data.csv"))

# ----------------------------------
# Extraction of most important genes
# ----------------------------------

# One-hot encoding of design
design <- model.matrix( ~ 0 + gse153922_eset[['source_name_ch1']])
colnames(design) <- levels(as.factor(gse153922_eset[['source_name_ch1']]))
colnames(design) <- c("trf2", "vector", "young", "old")
print(design)


# Comparisons of interest:
#   young_trf2: young vs senescence
#   young_old: young vs old
#   old_trf2: old vs senescence
contrasts_matrix <- makeContrasts(
  young_trf2 = young - (trf2 - (vector - young)),
  young_old = young - old,
  old_trf2 = old - (trf2 - (vector - young)),
  levels = design
)
# print(contrasts_matrix)

# Fitting observations to statistical model
fit <- lmFit(gse153922_eset,design)
fit2 <- contrasts.fit(fit,contrasts=contrasts_matrix)
fitted.ebayes <- eBayes(fit2)
summary(decideTests(fitted.ebayes,lfc=1, p.value = 0.05))

# Matching gene symbols to probe ID
probes<-AnnotationDbi::select(hugene20sttranscriptcluster.db,rownames(fitted.ebayes),c("SYMBOL","GENENAME"),keytype="PROBEID")
probes<- probes %>% drop_na()
interesting_genes<-topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=1)
interesting_genes<-merge(probes, interesting_genes, by.x="PROBEID", by.y="row.names")
interesting_genes<-group_by(interesting_genes, PROBEID) %>% slice(1)

# Save comparison data
write.csv(interesting_genes, "./interesting_genes.csv", row.names = FALSE)

# -------------------------------
# Extraction of clusters of genes 
# -------------------------------

# Cluster 1: Increase/decrease expression between both young_old and old_trf2
c1.0<-filter(interesting_genes,young_old>=1 & old_trf2>=1)
c1.1<-filter(interesting_genes,young_old<=-1 & old_trf2<=-1)
c1<-rbind(c1.0, c1.1) # No genes in this category


# Cluster 2: Increase/decrease expression between young_old but decrease/increase expression between old_trf2
c2.0<-filter(interesting_genes,young_old<=-1 & old_trf2>=1)
c2.1<-filter(interesting_genes,young_old>=1 & old_trf2<=-1)
c2<-rbind(c2.0, c2.1)
write.csv(c2, "./cluster2.csv", row.names = FALSE)

# Cluster 3: No significant change between young_old but change in expression between old_trf2
c3.0<-filter(interesting_genes,(young_old>-1 & young_old<1) & old_trf2>=1)
c3.1<-filter(interesting_genes,(young_old>-1 & young_old<1) & old_trf2<=-1)
c3<-rbind(c3.0, c3.1)
write.csv(c3, "./cluster3.csv", row.names = FALSE)

# Cluster 4: Change in expression between young_old but no significant change between old_trf2
c4.0<-filter(interesting_genes,young_old>=1 & (old_trf2>-1 & old_trf2<1))
c4.1<-filter(interesting_genes,young_old<=-1 & (old_trf2>-1 & old_trf2<1))
c4<-rbind(c4.0, c4.1)
write.csv(c4, "./cluster4.csv", row.names = FALSE)

# Cluster 5. No change between young_old and old_trf2 but significant change overall for young_trf2
c5.0<-filter(interesting_genes,(young_old>-1 & young_old<1) & (old_trf2>-1 & old_trf2<1)& young_trf2>=1 )
c5.1<-filter(interesting_genes,(young_old>-1 & young_old<1)& (old_trf2>-1 & old_trf2<1)& young_trf2<=-1)
c5<-rbind(c5.0, c5.1)
write.csv(c5, "./cluster5.csv", row.names = FALSE)

