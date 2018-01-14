getwd()
# set file location, according to your computer
setwd("/Users/yanjun/Documents/R data/RNAseq_DESeq2_Stim-TRIP fat tissue/") 
getwd()

## Note that following pacakages need to be installed or upstated according to your computer setting
## try install DESeq2(v3.5) http:// if https:// URLs are not supported

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("IRanges")
library("IRanges")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("XVector")
library("XVector")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
library("GenomicRanges")

install.packages("plyr")
library("plyr")

install.packages("lazyeval")
library("lazyeval")

install.packages("tibble")
library("tibble")

install.packages("htmltools")
library("htmltools")

install.packages("stringi")
library("stringi")

install.packages("data.table")
library("data.table")

install.packages("XML")
library("XML")

install.packages("lattice")
library("lattice") 

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")

install.packages("matrixStats")
library("matrixStats")

source("https://bioconductor.org/biocLite.R")
biocLite("topGO")
library("topGO")

biocLite("org.Dm.eg.db")
library("org.Dm.eg.db")

biocLite("BiocStyle")
library("BiocStyle")

install.packages("rmarkdown")
library("rmarkdown")


library("geneplotter") #Loading required package: lattice
library("ggplot2")

install.packages("LSD")
library("LSD")

library("gplots")
library("RColorBrewer")
library("stringr")

library("genefilter")

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library("biomaRt")

install.packages("bindrcpp")
library("bindrcpp")

install.packages("dplyr")
library("dplyr")

biocLite("Biostrings") #Memory efficient string containers, string matching algorithms, and other utilities, for fast manipulation of large biological sequences or sets of sequences.
library("Biostrings")


biocLite("ShortRead") #This package implements sampling, iteration, and input of FASTQ files. The package includes functions for filtering and trimming reads, and for generating a quality assessment report.
library("ShortRead")

biocLite("EDASeq") #Numerical and graphical summaries of RNA-Seq read data. 
library("EDASeq")


install.packages("fdrtool")
library("fdrtool")

biocLite("vsn") #The method uses a robust variant of the maximum-likelihood estimator for an additive-multiplicative error model and affine calibration.Normally for microarrary data
library("vsn") 

library("Biobase") # for pull out the matrix

install.packages("ggrepel")
library("ggrepel") # for vocalno plot

install.packages("gplots") # for  clustering analysis heatmap
library("gplots")

library("devtools")

devtools::session_info() # recheck your R version etc information

#### 1.loading datasets_the count table
# In DESeq analysis, DESeqDataSet is putforward: The assay(s) (red block) contains the matrix (or matrices) of summarized reads values
#the rowData (blue block) contains information about the genomic ranges: Gene.ID here, or gene name etc
# he colData (purple block) contains information about the samples or experiments: here StimTI On or Off
# http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html is downloaded DESeqDataSet, which include rowData, colData, and count table. 
# here we need construct indicate the count table, and experiment design by ourself. 
# Reading presaved 1551 fat body RNAseq countting table from csv file
getwd()
fbRNAseq<-read.csv("/Users/yanjun/Documents/R data/RNAseq_DESeq2_Stim-TRIP fat tissue/1551_Stim-TRIP fat tissue_gene reads counts.CSV", header=T, sep=",") # as countData, read matrix
head(fbRNAseq)
summary(fbRNAseq) 
class(fbRNAseq)
str(fbRNAseq)
dim(fbRNAseq)
# to pullout the countmatrix from the fbRNAsea


### delete Gene.ID column trial
fbRNAseq1<-fbRNAseq[,2:7] # as count data, without Gene.ID column and 1st row sample information
head(fbRNAseq1)


#### following workflow from http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html
# according to DESeq2 tutorial, fbRNAseq means countData: a matrix of read counts

#### 2. similar as EdgeR, I should also provide a experiment design for colData, change appendix 4 name into "fbRNAseq_design2" accordingly
# A DESeqDataSet object must have an associated design formula, which can be established following the standard procedure
# In my case, we need read the design formular by ourself, which expresses the variables which will be used in modeling such as control and treatment group samples
# the design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model.
fbRNAseq_design=read.table("/Users/yanjun/Documents/R data/RNAseq_DESeq2_Stim-TRIP fat tissue/Stim-TRIP fat tissue_RNAseq experiment design.csv", header = T,sep=",")
fbRNAseq_design # as colData, like experiment design matrix in edgeR, protocol column as protocol

# reading rowname from csv files, trial
fbRNAseq_rowname=fbRNAseq$Gene.ID
head(fbRNAseq_rowname)

#### 3.we can construct a DESeqDataSetFromMatrix based on count matrix(fbRNAseq), sample information_colData(fbRNAseq_design), 
## DESeqDataSetFromMatrix()need to be defined by ourself, similar to 'DESeq2Table' in the online protocol
dds <- DESeqDataSetFromMatrix(countData=fbRNAseq1, colData=fbRNAseq_design, design = ~ Condition )
dds  # One need find a solution to show rownames in DESeqDataSet, which should be Gene.ID names
dds$Condition
colData(dds)
rowData(dds)

#### 4.Important: adding row data to the DESeqDataSetFromMatrix(), this can be used to adding other information, like Gene symbols.
## assign gene.ID as rowname of the DESeqDataSet
# you have to make sure, you have unique gene symbols for every gene_id
# Check NA or duplicated GeneID numbers
sum(is.na(fbRNAseq$Gene.ID))
sum(duplicated(fbRNAseq$Gene.ID))

## original script from DESeqdata set:rownames(dds) <- mcols(dds)$gene # assign Gene.ID from rowData to row names
rownames(dds) <- fbRNAseq$Gene.ID
dds
head(dds)
dim(dds)

####  5.Quality control and Normalization of the count data
### a.Check how many genes with non-zero counts
#fbRNAseq as data.frame, which can not be analysed by counts function
GeneCounts <- counts(dds) # remember to loading correpondent DESeq2 Library
head(GeneCounts)
idx.nz <- apply(GeneCounts, 1, function(x) {all(x>0)}) # only apply to column information of fbRNAseq
sum(idx.nz) # if the non-zero reads are very low, there might be something wrong with RNAseq sequencing or experiments step
head(idx.nz)

### b.random sample from the count matrix, as a snap view
nz.counts <- subset(GeneCounts, idx.nz) 
head(nz.counts)
sam <- sample(dim(nz.counts)[1], 5 )
nz.counts[sam, ]

## c. filtering out low count reads
dim(GeneCounts)
GeneCounts<-GeneCounts[rowMeans(GeneCounts)>5, ]
dim(GeneCounts)
head(GeneCounts)

## filtering on dds directly
dds<-dds[rowMeans(counts(dds))>5,]
dim(dds)

#### 6.Normalization and size factors estimation
####As different libraries will be sequenced to different depths, 
####offsets are built in the statistical model of DESeq2 to ensure that parameters are comparable.
####Generally, the ratios of the size factors should roughly match the ratios of the library sizes.
####Count data sets typically show a trombone shape in average vs meanâdifference plots (MAâplots), reflecting the higher variability of log ratios at lower counts.
####The MAplots are saved to a file in order to avoid cluttering the document.

### a.We first relevel the condition factor to make sure that we have control group samples(Stim-TI Off) as a base level.
dds$Condition
colData(dds)$Condition <-factor(colData(dds)$Condition, levels=c("Stim-TRIP Off","Stim-TRIP On"))
colData(dds)$Condition
## Assigning different sample to two groups: "Treatment" or "control", report errors at moment


### b.We now estimate estimate the size factors for each column  (sample)
dds<-estimateSizeFactors(dds) # &the first sub-function which is called by DESeq
sizeFactors(dds) # quality check need library size correction
dds

### c.we plot the densities of counts for the different samples to see if library size normalisation works or not.
## it is similar to the probability distribution function
#a succesful normalization will lead to overlapping densities.
multidensity( counts(dds, normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 5000)) 

### d.This  same holds true for the empircal cummulative distribution functions (ECDFs) 
###which can be thought of as integrals of the densities and give the probability of observing a certain number of counts equal to x or less given the data.
multiecdf( counts(dds, normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 5000))

### e. pair sample MA plots based on row reads
### To further assess systematic differences between the samples, we can also plot pairwise meanâaverage plots: 
### We plot the average of the un-log transformed counts vs the fold change per gene for each of the sample pairs.
###The combn(x,m) function  helps us to automatically generate all sample pairs(m=2) from x and the function MDPlot from the EDASeq package then generates the the pairwise MA plots.
# MA plot not working due to image showing scripts when use dev.off() function
# The function combn helps us to automatically generate all sample pairs
#the function MDPlot from the EDASeq package then generates the the pairwise Mean-average plots.
#pdf("pairwiseMAs.pdf")
colnames(dds)
colData(dds)$Condition
MA.idx = t(combn(1:6, 2)) # 6 samples, take 2 for all the combinations
MA.idx
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(dds, normalized = T)[idx.nz ,], 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-5,5))
}


#### 7.Sample clustering heatmaps and PCA plot based on sample datasets
###we use it on the regularized log transformed data to calculate alculate the Euclidean distance between samples.
##The aim of the regularized log transform is to stabilize the variance of the data and to make its distribution roughly symmetric 
##methods for clustering and ordination (e.g., principalâcomponent analysis and the like), work best for (at least approximately) homoskedastic data;
##In RNA-Seq data, however, variance grows with the mean. 
##For genes with high counts, the rlog transformation differs not much from an ordinary log2 transformation.
##For genes with lower counts, however, the values are shrunken towards the genesâ averages across all samples.
##Using an empirical Bayesian prior in the form of a ridge penalty, this is done such that the rlogâtransformed data are approximately homoskedastic. 

### a.IMPORTANT: only after rlog, data are applyed to sample similarity analysis

rld <- rlogTransformation(dds, blind=TRUE)  #&& 2nd subfunction of DESeq()
head(rld)
distsRL <- dist(t(assay(rld))) # calculating Elucidean distance valuses
mat <- as.matrix(distsRL)

rownames(mat) <-  colData(rld)$Condition  
colnames(mat) <-  colData(rld)$Sample.ID

### b. Heatmap of sample data clustering analysis
## First we We visualize the Euclidean distance distances in a heatmap, using the function heatmap.2 from the r CRANpkg("gplots") package.
# requires gplot
hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13)) #drawing heatmaps of Eulidean distance
mat

### c. PCA analysis, also use rlog tranformed read counts
# require ggplot2
ntop = 500 # set 500 genes

Pvars <- rowVars(assay(rld)) # calculating variances between rows
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))] # selecting 500 genes with high variance

PCA <- prcomp(t(assay(rld)[select, ]), scale = F) #We then use the function prcomp to compute the principal components. 
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

#The resulting principal component scores will be ploted by the function qplot from ggplot2. 
dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                    sampleNO = colData(rld)$Sample.ID,
                    condition = colData(rld)$Condition)
head(dataGG)


(qplot(PC1, PC2, data = dataGG, color = sampleNO, 
       main = "PC1 vs PC2, top 500 variable genes", size = I(8))
  + labs(x = paste0("PC1, Variance:", round(percentVar[1],4)),
         y = paste0("PC2, Variance:", round(percentVar[2],4)))
  + scale_colour_brewer(type="qual", palette="Paired", direction=1 ) 
  + theme(legend.key = element_blank()) 
  + geom_vline(xintercept=0) + geom_hline(yintercept=0)
  + theme(axis.line.x = element_line(colour="black"),
          axis.line.y = element_line(colour="black"))
  + theme(axis.text.x = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=10,angle=0,hjust=1,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="plain")) # change axis colour and size;
)


#### 8. Differential expression analysis

### a. DESeq2 using negative binomial model to estimate dispersion variances
dds <- estimateDispersions(dds)  # &&&the third sub-function which is called by DESeq
plotDispEsts(dds)
# The black points are the dispersion estimates for each gene as obtained by considering the information from each gene separately.
#the red trend line is fitted, which shows the dispersionsâ dependence on the mean, and 
#then shrink each genE estimate towards the red line to obtain the final estimates (blue points) that are then used in the hypothesis test.
# The blue circles above the main cloud of points are genes which have high gene-wise dispersion estimates which are labelled as dispersion outliers.

### b. Statistical testing of Differential expression
### We can perform the statistical testing for differential expression and extract its results.
### Calling performs the wald test for differential expression, 
### while the call to the function extracts the results of the test and returns adjusted pâvalues according to the BenjaminiâHochberg rule to control the FDR.
### further show the results of logarithmic fold change log2(Stim-TI On/Stim-TI Off)

## a.Statitistical testing with Waldtest 
dds <-  nbinomWaldTest(dds) # &&&& the fourth sub-function which is called by DESeq
head(dds)# important finally, dds were analyzed with waldtest

DESeq2Res <- results(dds, pAdjustMethod = "BH") # FDR calculated by Benjamini-Hochberg Procedure, which is a powerful tool that decreases the false discovery rate
head(DESeq2Res)
head(rownames(DESeq2Res))
dim(DESeq2Res)
# we can summarize the basic parameters of the DESeq2Res
summary(DESeq2Res)
# Check the number of differential expressed genes
table(DESeq2Res$padj < 0.1)

## b.Independent filtering test-optional, required for Stim-TRIP fat tissue RNAseq data 
##  By removing the weakly-expressed genes from the input to the BH-FDR procedure, we can find more genes to be significant among those which we keep, and so improve the power of our test.
## The DESeq2 software-DESeq2 function automatically performs independent filtering which maximizes the number of genes which will have a BH-adjusted p-value less than a critical value (by default, alpha is set to 0.1).
##  this filter is blind to the assignment of samples to the deletion and control group and hence independent under the null hypothesis of equal expression.
##  We can observe how the number of rejections changes for various cutoffs based on mean normalized count.
plot(metadata(DESeq2Res)$filterNumRej, type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")

## c. A histogram of p-values should always be plotted in order to check whether they have been computed correctly. 
##  The null-p-values follow a uniform distribution on the unit interval [0,1] if they are computed using a continuous null distribution. Significant p-values thus become visible as an enrichment of p-values near zero in the histogram.
## Wald test is based on N(0,1) null test 
## Very often, if the assumed variance of the null distribution is too high, we see hill-shaped p-value histogram. If the variance is too low, we get a U-shaped histogram, with peaks at both ends.
hist(DESeq2Res$pvalue, col = "lavender", main = "WT vs Deletion", xlab = "p-values")
# based on L shape p value distribution, we estimate the p value threshold well.
#Thus we don't need further p value independent fitering

##### Attention optional
#if we have e.g. batches or "outlying" samples that are consistently a bit different from others within a group, the dispersion within the experimental group can be different and a single dispersion parameter not be appropriate.
# Fortunately, there is software available to estimate the variance of the null-model from the test statistics. This is commonly referred to as "empirical null modelling".
## Here we use the fdrtool for this using the Wald statistic as input. This packages returns the estimated null variance, as well as estimates of various other FDR-related quantities and the p-values computed using the estimated null model parameters.
DESeq2Res1 <- DESeq2Res[ !is.na(DESeq2Res$padj), ] # remove genes with NA with independent filtering
DESeq2Res1 <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ] # remove dispertion outliers, which has NA p value.
head(DESeq2Res1)
dim(DESeq2Res1)

# We further remove adjusted p values
DESeq2Res1 <- DESeq2Res1[, -which(names(DESeq2Res1) == "padj")]
head(DESeq2Res1)
dim(DESeq2Res1)

#We can now use z-scores returned by DESeq2as input to fdrtool to re-estimate the p-values
dim(dds)
dim(DESeq2Res)
dim(DESeq2Res1)
head(DESeq2Res1)
head(DESeq2Res1$stat)
FDR.DESeq2Res<-fdrtool(DESeq2Res1$stat, statistic= "normal", plot = T)
FDR.DESeq2Res
FDR.DESeq2Res$param[1, "sd"] # specific dispersion
# the estimated null model variance is 0.76, which is higher than 1, the theoretical one in previous Wald test.
# We now add the new BH-adjusted p-values add values to the results data frame.
DESeq2Res1[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
head(DESeq2Res1)
dim(DESeq2Res1)

#We can now re-plot the histogram of the "dispersion factor corrected" p-values.
hist(FDR.DESeq2Res$pval, col = "royalblue4", 
     main = "Stim-TRIP On vs Stim-TRIP Off, correct null model", xlab = "Corrected p-values")
#After independent filtering, the distribution of p value is better 

### b.Exploring differential expressed genes
## we can recheck differential expressed genes
table(DESeq2Res1$padj < 0.1) # check the number of significantly differentially expressed genes
# Since we underestimate dispersion at beginning, the number of DEG is reduced
table(DESeq2Res1$log2FoldChange>0&DESeq2Res1$padj < 0.1) # the number of upregulated genes
table(DESeq2Res1$log2FoldChange<0&DESeq2Res1$padj < 0.1) # the number of downregulated genes

# Extracting differential expressed genes
DGE<-subset(DESeq2Res1,abs(DESeq2Res1$log2FoldChange) >=0.5 & DESeq2Res1$padj <= 0.1 )
dim(DGE)
# Extracting upregulated genes
Upregulated<-subset(DESeq2Res1, DESeq2Res1$log2FoldChange>= 0.5 & DESeq2Res1$padj <= 0.1)
dim(Upregulated)
# Extracting downregulated genes
Downregulated<-subset(DESeq2Res1, DESeq2Res1$log2FoldChange<= -0.5 & DESeq2Res1$padj <= 0.1)
dim(Downregulated)

#Check information about which variables and tests were used by mcols function
head(mcols(DESeq2Res1))
mcols(DESeq2Res1)$description
Readme<-mcols(DESeq2Res1)$description
head(Readme)


# We can order the results table by the smallest adjusted p value
DESeq2ResReordered <- DESeq2Res1[order(DESeq2Res1$padj),]
head(DESeq2ResReordered)

# Export analysed results as csv files
write.csv(as.data.frame(DESeq2Res1), file = "DESeq2Res_Stim-TRIP On-vs-Off gene expression analysis.csv")
write.csv(as.data.frame(DESeq2ResReordered), file = "DESeq2Res_Stim-TRIP gene expression reordered by adjusted p value.csv")
write.csv(as.data.frame(Readme), file = "DESeq2Res_colomn name_Readme.csv")

write.csv(as.data.frame(DGE), file = "DESeq2_Stim-TRIP On-vs-Off differential expressed genes.csv")
write.csv(as.data.frame(Upregulated), file = "DESeq2Res_Stim-TRIP On-vs-Off upregulated genes.csv")
write.csv(as.data.frame(Downregulated), file = "DESeq2Res_Stim-TRIP On-vs-Off downregulated genes.csv")

##try rio package to combine csv files into Excel file
install.packages("rio")
library("rio")

x <- import("DESeq2Res_Stim-TRIP On-vs-Off upregulated genes.csv")
head(x)
y <- import("DESeq2Res_Stim-TRIP On-vs-Off downregulated genes.csv")
head(y)
z <- import("DESeq2_Stim-TRIP On-vs-Off differential expressed genes.csv")
head(z)

export(list(Upregulated = x, Downregulated = y, Differential_regulated=z), "DESeq2_Stim-TRIP On-vs-Off differential expressed genes_2017.xlsx")
##### Attention above scripts are optional, when independent filtering is necessary


#### 9. Further Exploring differential expressed gens

### a. plotMA show the the shrinkage (incorporationn of zero-centered normal prior) of log2fold changes (y axis) over the mean of normalized counts
##It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
## Here in recent DESeq2 versions, we need to use lfcShrink functions
## Points will be colored red if the adjusted p value is smaller than 0.1
## Points fall out of the window will be ploted as open triangles either pointing up or down
resultsNames(dds)
res.shr <- lfcShrink(dds=dds, coef=2, res=DESeq2Res) # res=res is important to keep DESeqdataset format
?lfcShrink()
head(res.shr)
plotMA(res.shr,main="DESeq2_Stim-TRIP On-vs-Off_MA plot", ylim=c(-3,3))

#plotMA(DESeq2Res, main="DESeq2_MA plot", ylim=c(-3, 3))

## one can use identify function to interactively detect the row number of individual genes by clicking on the plot
identify(DESeq2Res$baseMean, DESeq2Res$log2FoldChange)
# one need Esc to escape the running after selecting certain points


##### Optional for extra interests on unshrunken log2 fold change
## one can get the "unshrunken" log2 fold changes with maximum likelihood estimate (MLE)
DESeq2ResMLE <- results(dds, addMLE = TRUE)
head(DESeq2ResMLE)
head(fbRNAseq)
#adding MLE log2FoldChange into DESeq2Res results
write.csv(as.data.frame(DESeq2ResMLE), file = "DESeq2ResMLE_results_DEA.csv")
#### Optional for extra interests



#### 10 based on http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#statistical-testing-of-differential-expression
### Here we can adding Gene Symbols to DESeq2Res
annoRes <- AnnotationDbi::select(org.Dm.eg.db, 
                                 keys=rownames(DESeq2Res), 
                                 columns=c("SYMBOL", "GENENAME"),
                                 keytype="ENSEMBL")
head(annoRes)
colnames(annoRes)[1] <- "Gene.ID"
head(annoRes)
annoRes1<-as.data.frame(annoRes)
head(annoRes1)

head(DESeq2Res)
dim(DESeq2Res)
head(rownames(DESeq2Res))
head(mcols(DESeq2Res))
DESeq2Res[,"Gene.ID"]  <- rownames(DESeq2Res)
head(DESeq2Res)
DESeq2Resdf<-as.data.frame(DESeq2Res)
head(DESeq2Resdf)
dim(DESeq2Resdf)

## merge annoRes1 with DESeq2Resdf
annoDESeq2Res <- merge(DESeq2Resdf, annoRes1, by=c("Gene.ID"))
head(annoDESeq2Res)
dim(annoDESeq2Res)

table(DESeq2Res$padj < 0.1) 
table(annoDESeq2Res$padj < 0.1)

# Reorder dataframe annoDESeq2Res
annoDESeq2ResReorder <- annoDESeq2Res[order(annoDESeq2Res$padj),]
head(annoDESeq2ResReorder, 15)

# ReExtracting differential expressed genes with Gene annotation information
DGE2<-subset(annoDESeq2Res,abs(annoDESeq2Res$log2FoldChange) >=0.5 & annoDESeq2Res$padj <= 0.1 )
DGE2Reorder <- DGE2[order(DGE2$log2FoldChange),]
head(DGE2Reorder)
# ReExtracting upregulated genes with Gene annotation information
Upregulated2<-subset(annoDESeq2Res, annoDESeq2Res$log2FoldChange>= 0.5 & annoDESeq2Res$padj <= 0.1)
head(Upregulated2)
Upregulated2Reorder <- Upregulated2[order(-Upregulated2$log2FoldChange),] # descending
head(Upregulated2Reorder)
# Extracting downregulated genes with Gene annotation information
Downregulated2<-subset(annoDESeq2Res, annoDESeq2Res$log2FoldChange<= -0.5 &annoDESeq2Res$padj <= 0.1)
head(Downregulated2)
Downregulated2Reorder <- Downregulated2[order(Downregulated2$log2FoldChange),] # ascending
head(Downregulated2Reorder)


# Export analysed results as csv files
write.csv(annoDESeq2ResReorder, file = "DESeq2_Stim-TRIP On-vs-Off gene expression reordered by adjusted p value.csv")
write.csv(DGE2Reorder, file = "DESeq2_Stim-TRIP On-vs-Off differential expressed genes.csv")
write.csv(Upregulated2Reorder, file = "DESeq2_Stim-TRIP On-vs-Off upregulated genes.csv")
write.csv(Downregulated2Reorder, file = "DESeq2_Stim-TRIP On-vs-Off Downregulated genes.csv")

##try rio package to combine csv files into Excel file
install.packages("rio")
library("rio")

x <- import("DESeq2_Stim-TRIP On-vs-Off upregulated genes.csv")
head(x)
y <- import("DESeq2_Stim-TRIP On-vs-Off Downregulated genes.csv")
head(y)
z <- import("DESeq2_Stim-TRIP On-vs-Off differential expressed genes.csv")
head(z)

export(list(Upregulated = x, Downregulated = y, Differential_regulated=z), "RNAseq_DESeq2_Stim-TRIP On-vs-Off gene differential expression analysis_2017.xlsx")


### b.It is useful to check the counts of reads for a single gene across the groups by using function plotCounts
# The counts are grouped by the variables in intgroup
# Here, we first select gene with the smallest padjusted value
plotCounts(dds, gene=which.min(DESeq2Res$padj), intgroup = "Condition")
# Genes with largest log2fold changes
plotCounts(dds, gene=which.max(DESeq2Res$log2FoldChange), intgroup = "Condition")
# One can select specific gene with rowname
#(dds, gene=rownames(DESeq2ResReordered[1, ], intgroup = "Condition")

### c.Volcano Plot of DESeq2 analysis
# Using ordered results table by the smallest adjusted p value:DESeq2ResOrdered
## a.judge padjust value based on threshhold 0.1
DESeq2ResReordered1 = as.data.frame(dplyr::mutate(as.data.frame(annoDESeq2ResReorder), sig=ifelse(annoDESeq2ResReorder$padj<0.1, "FDR<0.1", "Not Sig"), row.names=annoDESeq2ResReorder$SYMBOL))
head(DESeq2ResReordered1)
dim(DESeq2ResReordered1)

## c.using ggplot 2 and ggrepel to draw vocalno plot
p = ggplot2::ggplot(DESeq2ResReordered1, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("red", "black")) +
  ggplot2::ggtitle("Volcano Plot of differential expressed genes in Stim-TRIP On fly fat tissue")

p + ggrepel::geom_text_repel(data=DESeq2ResReordered1[1:15, ], ggplot2::aes(label=DESeq2ResReordered1$row.names[1:15])) # labeling 15 most differential regulated genes



###### Backup scripts for heatmap analysis of differental expressed genes
#Plot heatmap (Top 30)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
tiff("heat_map_top30.tiff",width=700,height=700,res=72)
par(mfrow=c(1,3))
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6))
heatmap.2(assay(rld)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
heatmap.2(assay(vsd)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
dev.off()





############################################################################### 
#####Gene Ontology term (biological process) enrichment analysis based on R scripts from Revigo
getwd()
setwd("D:/R data/REVIGO/1551_StimTI_fbRNAseq/")

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0008152","metabolic process",52.035, 3.122, 5.255, 3.772,-2.7124,0.989,0.000),
                     c("GO:0019752","carboxylic acid metabolic process", 3.173,-4.268,-2.731, 2.559,-12.0585,0.283,0.000),
                     c("GO:0030431","sleep", 1.248, 5.565, 2.169, 2.155,-1.4617,0.907,0.054),
                     c("GO:0051156","glucose 6-phosphate metabolic process", 0.141,-4.842, 5.125, 1.230,-2.9489,0.661,0.101),
                     c("GO:1901615","organic hydroxy compound metabolic process", 1.327,-1.837, 6.971, 2.182,-1.9821,0.824,0.102),
                     c("GO:0072524","pyridine-containing compound metabolic process", 0.466, 0.526,-7.982, 1.732,-4.9872,0.694,0.115),
                     c("GO:0005975","carbohydrate metabolic process", 3.683,-0.542, 4.639, 2.623,-4.6253,0.794,0.117),
                     c("GO:0044710","single-organism metabolic process",18.239, 1.639, 1.806, 3.317,-18.9586,0.800,0.125),
                     c("GO:0006732","coenzyme metabolic process", 1.072, 4.695,-1.207, 2.090,-7.3152,0.676,0.126),
                     c("GO:0051186","cofactor metabolic process", 1.327, 4.539,-3.341, 2.182,-5.6819,0.795,0.130),
                     c("GO:0006091","generation of precursor metabolites and energy", 1.696, 3.406,-4.946, 2.288,-4.1163,0.790,0.134),
                     c("GO:0006793","phosphorus metabolic process",10.117, 0.863,-2.973, 3.061,-2.8913,0.748,0.176),
                     c("GO:0006081","cellular aldehyde metabolic process", 0.272,-7.219,-2.729, 1.505,-1.7177,0.597,0.297),
                     c("GO:0015980","energy derivation by oxidation of organic compounds", 1.371,-5.085,-3.209, 2.196,-2.5748,0.494,0.360),
                     c("GO:0044712","single-organism catabolic process", 2.523,-7.103,-1.616, 2.459,-5.1158,0.466,0.390),
                     c("GO:0044255","cellular lipid metabolic process", 2.848,-4.157,-2.048, 2.512,-2.3482,0.465,0.397),
                     c("GO:0044723","single-organism carbohydrate metabolic process", 2.945,-6.437,-0.938, 2.526,-2.6819,0.446,0.399),
                     c("GO:0006629","lipid metabolic process", 4.140,-5.720,-1.064, 2.674,-2.6827,0.499,0.419),
                     c("GO:0019637","organophosphate metabolic process", 4.131,-0.933, 0.969, 2.673,-2.5327,0.624,0.422),
                     c("GO:0072521","purine-containing compound metabolic process", 2.479, 0.456,-6.524, 2.452,-2.4427,0.647,0.436),
                     c("GO:0006575","cellular modified amino acid metabolic process", 0.580,-0.195,-8.026, 1.826,-2.2931,0.698,0.447),
                     c("GO:0019682","glyceraldehyde-3-phosphate metabolic process", 0.097,-5.745, 1.041, 1.079,-2.5286,0.513,0.452),
                     c("GO:0006730","one-carbon metabolic process", 0.211,-6.518,-3.797, 1.398,-1.5368,0.481,0.464),
                     c("GO:0044281","small molecule metabolic process", 8.702,-5.966,-1.835, 2.996,-17.3487,0.464,0.495));


one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
head(one.data);
one.data$description
max(one.data$plot_size)
min(one.data$plot_size)
mediandisp<-median(one.data$dispensability)
mediandisp
meandisp<-mean(one.data$dispensability)
meandisp
# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(1) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, colour = I (alpha ("black", 1) )) + scale_size_area();
p1 <- p1 + scale_size(range=c(5, 25), limits = c(1,  max(one.data$plot_size)) )+ theme_bw(); #  + scale_fill_gradientn(colours = colours = c("blue", "green", "yellow", "red"), limits = c(-300, 0) );
ex <- one.data[one.data$dispensability < meandisp, ]; 
head(ex)
ex$description

p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = ex$description), colour = I(alpha("black", 1)), size = 6 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space Y");
p1 <- p1 + theme(legend.key = element_blank()) ;
# inserted scripts
p1 <- p1 + theme(axis.text.x = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
                 axis.text.y = element_text(colour="black",size=10,angle=0,hjust=1,vjust=0,face="plain"),  
                 axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
                 axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="plain")); # change axis colour and size;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);

# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).
# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");


########################################################################
devtools::session_info() # recheck your R version etc information
