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
install.packages("edgeR")
BiocManager::install("org.Mm.eg.db")
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
# Read the sample information into R
sampleinfo <- read.delim("SampleInfo.txt")
View(sampleinfo)
sampleinfo
# Read the data into R
seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
head(seqdata)
View(seqdata)
dim(seqdata)
# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]

View(countdata)
head(countdata)
colnames(countdata)
# short names of sample
substr("ThisIsAString", start=1, stop=5)
# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata), 1, 7)
View(countdata)
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
