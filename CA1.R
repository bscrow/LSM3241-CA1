library(GEOquery)
library(oligo)
library(tidyverse)
# it can sometimes give errors, but the resulting elements of the list are easier to work with417
gse1 <- getGEO('GSE153922') # just getting the main data post-rma
#if you download using getGeo,rma(or something else)has already been done


# pull out the ExpressionSet object(3'biased probe set) as the first member of the list 
gse1 <- gse1[[1]]



## In some cases, we may want to get the data itself
## in case we need addictional columns from the GSE records
## or if there is an error in parsing.
## This creates a data structure that resembles the 
## underlying SOFT format.


## Let's look at gse1


# Let's look at gse2

## back to gse1

names(pData(gse1))/#head()
gse1@experimentData
pData(gse1)[,c("culture medium:ch1","supplementary_file")]

### GET OUR RAW DATA


pData(raw_data)







## summarization

# we can't easily run median polish because it's run 56K times
#the resulting value will be after log transformation,
#intensity of the gene different from raw_data after rma
eset <- rma(raw_data)
exprs(eset['1007_s_at',])

pData(eset)
methods(class=class(eset))

pData(eset)[['culture medium:ch1']]
eset[['culture medium:ch1']]

# differential expression

library(limma)
eset[['culture']] <- eset[['culture medium:ch1']]
pData(eset)
eset$culture
~ 0 + eset$culture
~ eset$culture   #intercept is the absolute value of expression
#baseline as MEGM
model <- model.matrix(~ 0 + eset$culture)
colnames(model) <- levels(as.factor(eset$culture))
model
contrasts <- makeContrasts(SCGM - MEGM, levels=model)
# test what's the difference between these two different media conditions.
# now we have a model, and have some contrasts! 

fit <- lmFit(eset,design=model) #t test 54k x 2 estimate of SD rma expression, p value and t test close to 0 it doesn't make sense, the reason that that happens is just because we have so few data points per gene.
methods(class=class(fit))
plotMA(fit)
fitted.contrast <- contrasts.fit(fit,contrasts)
fit2 <- eBayes(fitted.contrast)   #the values that have standard deviation that are actually zero in the original measurement will be raised a little bit.
topTable(fit2,lfc=2)
top_gene <- row.names(topTable(fit2,lfc=2))[1]
exprs(eset)[top_gene,]
exprs(eset)["228335_at",]
dim(topTable(fit2,number=Inf,p.value = 0.05,lfc=2))
# at this point, we've got some genes!

volcanoplot(fit2)#a small P value would be associated with a large fold change either up or down,increase in pvalue y-axis decreas
abline(v=c(-2,2),col='red')
abline(h=1.3,col='red')

interesting_genes <- row.names(topTable(fit2,lfc=1,p.value=0.05,number=Inf))
heatmap(exprs(eset[interesting_genes,]),
        distfun   = function(x) as.dist(1-cor(t(x))))
#generate correlations between the genes and then it's going to take one minus the correlation function and treat that as a distance 
### what about batch effects? 
# This is all the data
#How are they distributed so if you rotate them so that the variation between samples is, if you rotate those six(# samples) dimensions
#otate it such that you only see the "tallest" variation out of the 6 versus the 54K positions
pc <- princomp(exprs(eset))
summary(pc)
?biplot
loadings <- as.data.frame(as.matrix(pc$loadings)[,c(1:2)])
loadings$culture <- eset$culture  
loadings
plot(loadings[,c(1,2)],pch=loadings[,3],
     main="all the data!")

head(exprs(eset))

# select interesting genes and do the same thing!
plot(pc)
#these are eigenvalues of covariance matrix and the eigenvectors are the actual principal components
pc2 <- princomp(exprs(eset[interesting_genes,]))
summary(pc2)
loadings2 <- as.data.frame(as.matrix(pc2$loadings)[,c(1:2)])
loadings2$culture <- eset$culture  
loadings2
plot(loadings2[,c(1,2)],pch=loadings2[,3],main="genes that are differentially expressed")
#variation is actually only about the differences between SCGM and MGM because we've selected for genes that are designed to show that.
#It is the rotation that maximizes all of the variance

