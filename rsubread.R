library(Rsubread)
fastq.files <- list.files( pattern = ".fastq.gz$", full.names = TRUE)
fastq.files
# build index 
buildindex(basename="chr1_mm101",reference="chr1.fa")
# alignement
args(align)
index.files <- list.files( pattern = "chr1_mm101", full.names = TRUE)
index.files
align.stats <- align(index="chr1_mm101",readfile1 = fastq.files)
bam.files <- list.files( pattern = ".BAM$", full.names = TRUE)
bam.files
propmapped(bam.files)
# quality score
qs <- qualityScores(filename="SRR1552450.fastq.gz",nreads=100)
dim(qs)
boxplot(qs)
# Counting
fc <- featureCounts(bam.files, annot.inbuilt="mm10")
# See what slots are stored in fc
names(fc)
## Take a look at the featurecounts stats
fc$stat
## Take a look at the dimensions to see the number of genes
dim(fc$counts)
## Take a look at the first 6 lines
head(fc$counts)
## Take a look at annotation file
head(fc$annotation)
fc.exon <- featureCounts(bam.files, annot.inbuilt="mm10",useMetaFeatures = FALSE)
fc.exon$stat
dim(fc.exon$counts)


# Pre-processing data ----
BiocManager::install("org.Mm.eg.db")
library(Rsubread)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
# Read the sample information into R
sampleinfo <- read.delim("SampleInfo.txt")
sampleinfo
# Read the data into R
seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
head(seqdata)
dim(seqdata)
# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)
colnames(countdata)
# short names of sample
substr("ThisIsAString", start=1, stop=5)
# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata), 1, 7)
# vérification of column
table(colnames(countdata)==sampleinfo$SampleName)

# Filtering to remove lowly expressed genes ----
# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)
# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5,h=10)
# create edgeR object ----
dgeObj <- DGEList(counts.keep)
# have a look at dgeObj
dgeObj
# See what slots are stored in dgeObj
names(dgeObj)
# Library size information is stored in the samples slot
dgeObj$samples
# Quality control
dgeObj$samples$lib.size
# Depth
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)
# Add a title to the plot
title("Barplot of library sizes")
# Normalization and log count ----
# Get log2 counts per million
logcounts <- cpm(dgeObj,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
# MDS Multidimensional scaling plots ----
plotMDS(dgeObj)
# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$CellType)
## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)

# Redo the MDS with cell type colouring
plotMDS(dgeObj,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
# Add a title
title("Cell type")

# Similarly for status
levels(sampleinfo$Status)
col.status <- c("blue","red","dark green")[sampleinfo$Status]
col.status
plotMDS(dgeObj,col=col.status)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")
# There is a sample info corrected file in your data directory
# Old sampleinfo
sampleinfo
# I'm going to write over the sampleinfo object with the corrected sample info
sampleinfo <- read.delim("SampleInfo_Corrected.txt")
sampleinfo
# Redo the MDSplot with corrected information
par(mfrow=c(1,2))
col.cell <- c("purple","orange")[sampleinfo$CellType]
col.status <- c("blue","red","dark green")[sampleinfo$Status]
plotMDS(dgeObj,col=col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
title("Cell type")
plotMDS(dgeObj,col=col.status)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")
# Dimension 3 appears to separate pregnant samples from the rest. Dim4?
plotMDS(dgeObj,dim=c(3,4))
labels <- paste(sampleinfo$SampleName, sampleinfo$CellType, sampleinfo$Status)
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group <- factor(group)
glMDSPlot(dgeObj, labels=labels, groups=group, folder="mds")
# Hierarchical clustering with heatmaps ----
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]

# Plot the heatmap
library(RColorBrewer)
library(gplots)
library(factoextra)
gg <- heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 500 most variable genes across samples",
          ColSideColors=col.cell,scale="row")
gg
# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes\nacross samples",ColSideColors=col.cell,scale="row")
dev.off()


# Normalisation for composition bias
# Apply normalisation to DGEList object
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples
#Before Normalization
par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")
#After normalization
par(mfrow=c(1,2))
plotMD(dgeObj,column = 7)
abline(h=0,col="grey")
plotMD(dgeObj,column = 11)
abline(h=0,col="grey")
# save data for day 2
save(group,dgeObj,sampleinfo,file="preprocessing.Rdata")

# Analysis DE ----

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
load("preprocessing.Rdata")
# Recap of pre-processing ----
## Read the counts from the downloaded data
seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
#
# Remove first two columns from seqdata

countdata <- seqdata[,-(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
countdata
colnames(countdata) <- substr(colnames(countdata), 1, 7)
countdata
## Calculate the Counts Per Million measure
myCPM <- cpm(countdata)
## Identify genes with at least 0.5 cpm in at least 2 samples
thresh <- myCPM > 0.5
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
## Convert to an edgeR object
dgeObj <- DGEList(counts.keep)
## Perform TMM normalisation
dgeObj <- calcNormFactors(dgeObj)
## Obtain corrected sample information
sampleinfo <- read.delim("SampleInfo_Corrected.txt")
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group
# Create the design matrix
# Create the two variables
group <- as.character(group)
type <- sapply(strsplit(group, ".", fixed=T), function(x) x[1])
status <- sapply(strsplit(group, ".", fixed=T), function(x) x[2])
# Specify a design matrix with an intercept term
design <- model.matrix(~ type + status)
design


# Estimating dispersion 
dgeObj <- estimateCommonDisp(dgeObj)
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)
plotBCV(dgeObj)
# testing for DE
# Fit the linear model
fit <- glmFit(dgeObj, design)

names(fit)
head(coef(fit))
# Conduct likelihood ratio tests for luminal vs basal and show the top genes
lrt.BvsL <- glmLRT(fit, coef=2)
topTags(lrt.BvsL)
# Contrast : compare pregnant and virgin
# Suppose we want to find differentially expressed genes between pregnant and virgin. We don’t have a parameter that explicitly will allow us to test that hypothesis. We need to build a contrast:
PvsV <- makeContrasts(statuspregnant-statusvirgin, levels=design)
lrt.pVsV <- glmLRT(fit, contrast=PvsV)
topTags(lrt.pVsV)
# save for annotation
save(lrt.BvsL,dgeObj,group,file="DE.Rdata")

# Annotation and Visualisation of RNA-seq results ----
suppressPackageStartupMessages(library(edgeR))
load("DE.Rdata")
results <- as.data.frame(topTags(lrt.BvsL,n = Inf))
results
# DE and MA plot
dim(results)
summary(de <- decideTestsDGE(lrt.BvsL))
detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.BvsL, de.tags=detags)
# Adding annotation to the edgeR results

library(org.Mm.eg.db)
columns(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
keys(org.Mm.eg.db, keytype="ENTREZID")[1:10]
# verify keys are valid
## Build up the query step-by-step
my.keys <- c("50916", "110308","12293")
my.keys %in% keys(org.Mm.eg.db, keytype="ENTREZID")
all(my.keys %in% keys(org.Mm.eg.db, keytype="ENTREZID"))
## to be filled-in interactively during the class.
ann <- select(org.Mm.eg.db,keys=rownames(results),columns=c("ENTREZID","SYMBOL","GENENAME"))
table(ann$ENTREZID==rownames(results))
results.annotated <- cbind(results, ann)
results.annotated
write.csv(results.annotated,file="B.PregVsLacResults.csv",row.names=FALSE)
# volcano plot
signif <- -log10(results.annotated$FDR)
plot(results.annotated$logFC,signif,pch=16)
points(results.annotated[detags,"logFC"],-log10(results.annotated[detags,"FDR"]),pch=16,col="red")
# Before following up on the DE genes with further lab work, a recommended sanity check is to have a look at the expression levels of the individual samples for the genes of interest. We can quickly look at grouped expression using stripchart. We can use the normalised log expression values in the dgeCounts object (dgeCounts$counts).
library(RColorBrewer)
par(mfrow=c(1,3))
normCounts <- dgeObj$counts
# Let's look at the first gene in the topTable, Krt5, which has a rowname 50916
stripchart(normCounts["110308",]~group)
# This plot is ugly, let's make it better
stripchart(normCounts["110308",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,col=1:6,method="jitter")
# Let's use nicer colours
nice.col <- brewer.pal(6,name="Dark2")
stripchart(normCounts["110308",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="    Krt5")
# Retrieving Genomic Location ----

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
tx <- TxDb.Mmusculus.UCSC.mm10.knownGene
columns(tx)
library(GenomicRanges)
simple.range <-GRanges("1", IRanges(start=1000,end=2000))
simple.range

chrs <- c("chr13", "chr15","chr5")
start <- c(73000000, 101000000, 15000000)
end <- c(74000000,102000000, 16000000)
my.ranges <- GRanges(rep(chrs,3), 
                     IRanges(start=rep(start,each=3),
                             end = rep(end,each=3))
)
my.ranges
keys <- c("50916","110308","12293")
genePos <- select(tx, keys=keys,
                  keytype = "GENEID",
                  columns=c("EXONCHROM","EXONSTART","EXONEND")
)
geneRanges <- GRanges(genePos$EXONCHROM, IRanges(genePos$EXONSTART,genePos$EXONEND), GENEID=genePos$GENEID)
geneRanges
findOverlaps(my.ranges,geneRanges)
seqlevelsStyle(geneRanges)
seqlevelsStyle(simple.range)
# Retrieving Gene Coordinates as GenomicRanges
exo <- exonsBy(tx,"gene")
exo
range(exo[["110308"]])
# Exporting tracks
sigGenes <- results.annotated[detags,]
sigGenes
exoRanges <- unlist(range(exo))
exoRanges
sigRegions <- exoRanges[na.omit(match(sigGenes$ENTREZID, names(exoRanges)))]
sigRegions
mcols(sigRegions) <- sigGenes[match(names(sigRegions), rownames(sigGenes)),]
sigRegions
sigRegions[order(sigRegions$LR,decreasing = TRUE)]
seqlevels(sigRegions)
sigRegions <- keepSeqlevels(sigRegions, paste0("chr", c(1:19,"X","Y")))

Score <- -log10(sigRegions$FDR)
rbPal <-colorRampPalette(c("red", "blue"))
logfc <- pmax(sigRegions$logFC, -5)
logfc <- pmin(logfc , 5)
Col <- rbPal(10)[as.numeric(cut(logfc, breaks = 10))]
mcols(sigRegions)$score <- Score
mcols(sigRegions)$itemRgb <- Col
sigRegions
install.packages("rtracklayer")
library(rtracklayer)
export(sigRegions , con = "topHits.bed")
# extracting reads bam and bai file 
library(GenomicAlignments)
list.files("BAM/")
my.reads <- readGAlignments(file="bam/MCL1.DG.bam",
                            param=ScanBamParam(which=generegion))
my.reads
my.reads <- readGAlignments(file="bam/MCL1.DG.bam",
                            param=ScanBamParam(which=generegion,
                                               what=c("seq","mapq","flag")))
my.reads
# ggplot2 for MA and volcano plot
library(ggplot2)
ggplot(results, aes(x = logCPM, y=logFC)) + geom_point() 
ggplot(results, aes(x = logCPM, y=logFC,col=FDR < 0.05)) + geom_point()
ggplot(results, aes(x = logCPM, y=logFC,col=FDR < 0.05)) + geom_point(alpha=0.4) + scale_colour_manual(values=c("black","red"))
ggplot(results, aes(x = logFC, y=-log10(FDR))) + geom_point()
# introduction to ggbio
install.packages("xfun")
install.packages("htmlTable")
BiocManager::install("biovizBase")
library(htmlTable)
library(ggbio)
library(Gviz)
library(biovizBase)
top200 <- sigRegions[order(sigRegions$LR,decreasing = TRUE)[1:200]]
plotGrandLinear(top200 , aes(y = logFC))
mcols(top200)$UpRegulated <- mcols(top200)$logFC > 0
plotGrandLinear(top200, aes(y = logFC, col = UpRegulated))
autoplot(top200,layout="karyogram",aes(color=UpRegulated,
                                       fill=UpRegulated))
autoplot(tx, which=exo[["110308"]])

# cover read when bam file import
autoplot(bam , stat = "coverage")