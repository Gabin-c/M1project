#Library
library(DESeq2)
library(Biobase)
library(tidyverse)
library(ggrepel)
### import data ----
#Firstly we import the 3 data table which we will use during this analysis. 
#The "counts_table" contains the counting of read mapped on each genes by samples,
#"airway_metadata" contains the information about the samples design 
#and "anno" contains the informations about the genes. 
counts_table <- read.csv("Data/airway_scaledcounts.csv") 
airway_metadata <- read.csv("Data/airway_metadata.csv")
anno <- read.csv("Data/annotables_grch38.csv")
anno <- anno %>% select(ensgene,symbol)

### Create the dds object ---
dds <- DESeqDataSetFromMatrix(counts_table,colData=airway_metadata,design = ~dex,tidy = TRUE)
# Set reference of experience, here "control"
colData(dds)$dex <- relevel(colData(dds)$dex , ref="control")
# To display experiment design.
colData(dds)
# To display column which biological condition is set.
design(dds)

### Explorate data ----
#Transformation of dds counts table in data frame to use ggplot package
count_table_dds <- as.data.frame(counts(dds))
#vizualisation  of count distribution ----
#To facilitate the vizualisation we use the log-freq of each count value "log(count+1)"
for( i in 1:8){
  p <- ggplot(data=count_table_dds, aes(log(count_table_dds[,i]+1))) + geom_histogram(breaks=seq(0,14,1),col="black",fill="grey")+theme_light()+labs(title=colnames(count_table_dds)[i], x="Count value (number of read by genes) in log(count+1)",y="Count frequency") + theme_bw()
  plot(p)
}
#Number for Null for each sample ----
apply(count_table_dds, 2 ,FUN = function(x) sum(x==0))
#vizualisation  of depth for each sample using deph.plot function ----
depth <- colSums(count_table_dds) 
depth <- as.data.frame(depth)
depth$sample <- row.names(depth)
ggplot(depth, aes( x=sample ,y=depth))+ geom_bar(stat="identity",col="black", fill="white")+labs(title = "Depth of each sample", x="Sample", y="Depth")+theme_bw()

#Vizualisation of gene count by sample
count_table_dds[,"name"] <- row.names(count_table_dds)
ggplot(count_table_dds, aes(x=count_table_dds[,"name"], y=count_table_dds[,gene])) + 
  geom_point(size=4,aes(colour=factor(name))) + 
  geom_segment(aes(x=count_table_dds[,"name"], xend=count_table_dds[,"name"], y=0, yend=count_table_dds[,gene]),linetype="dotdash")+ 
  theme(axis.text.x = element_blank() )+ 
  labs(title=paste("Count of",gene,  "for each sample"),x="Samples",y="Counts")+ 
  guides(color= guide_legend(title = "Sample", override.aes = list(size=5))) +
  theme(plot.title = element_text(face = "bold", size= 18)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(legend.text=element_text(size=13)) +
  theme(legend.title=element_blank())

#Differential expression analysis and normalization ----
dds <- DESeq(dds)
#Extraction of the DE result
res <- results(dds,tidy = TRUE)
# The scale factor resulting of normalization are stock in size factor of our design.
colData(dds) 
#Depth of count table after normalization ----
depth_normalize <- colSums(counts(dds, normalized= TRUE))
depth_normalize <- as.data.frame(depth_normalize)
depth_normalize$sample <- row.names(depth_normalize)
ggplot(depth_normalize, aes( x=sample ,y=depth_normalize))+ geom_bar(stat="identity",col="black", fill="white")+labs(title = "Depth of each sample", x="Sample", y="Depth")+theme_minimal()
#vizualisation  of count distribution after normalization ----
count_normalize <- counts(dds, normalized= TRUE)
count_normalize <-as.data.frame(count_normalize)
for(i in 1:8){
  p <- ggplot(data=count_normalize, aes(log(count_normalize[,i]+1))) + geom_histogram(breaks=seq(0,14,1),col="black",fill="grey")+theme_light()+labs(title=colnames(count_normalize)[i], x="Count value (number of read by genes) in log(count+1)",y="Count frequency")
  plot(p)
}


# Dispersion plot ----
# Relationship between dispersion and counts means.
# DESeq2 offer a function that can directly display a plot which describe the relation ship beetween dispersion and count mean.

DESeq2::plotDispEsts(dds, main= "Relationship between dispersion and counts means")
# We obtain a plot that show the final estimate which are obtain after shrunk of genes estimate and we finaly observe the fitted estimate. We can also observed outliers value.

### DE analysis results ----

#Show a summary of DE results at alpha = 0.05

summary(results(dds),0.05)

# Number of genes wich is differential express at 5%
up_regulated <- res %>% filter(padj <= 0.05 & log2FoldChange > 0) %>% nrow()
down_regulated <- res %>% filter(padj <= 0.05 & log2FoldChange < 0) %>% nrow()
tb <- table(res$padj <= 0.05, useNA="always")
tb.DE <- data.frame("No DE" = tb[1], "Down regulated" = down_regulated, "Up regulated" = up_regulated, "NA" = tb[3]  )
row.names(tb.DE) <- ""
tb.DE

# MA plot : relationship between mean count of a gene and it log2 ratio between the two conditions
res <- as.data.frame(res)
res <- res %>% mutate(sig=padj<0.05)

ggplot(res, aes(x = baseMean, y = log2FoldChange, col = sig)) + 
  geom_point() + 
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed",color = "black") + 
  ggtitle("MA plot") + theme_bw()

#Volcano plot at alpha 0.5
#Show ID of the most DE gene 
res_v <- res %>% mutate(sig=padj<0.05) %>%  arrange(padj) %>%
  inner_join(anno,by=c("row"= "ensgene"))


ggplot(res_v, aes(x=log2FoldChange, y=-log10(padj), col=sig)) +
  geom_point() +
  ggtitle("Volcano plot labelling top significant genes") +
  geom_text_repel(data = subset(res_v, (-log10(padj) > 30 | log2FoldChange > 6 | log2FoldChange < -6)),
                  aes(label = symbol),
                  size = 4,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), color = "darkblue") +
  scale_colour_discrete(name="",
                        labels=c("Not significative", "Significative", "NA")) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  geom_vline(xintercept=0,linetype="dashed", color = "red")+
  theme(legend.text=element_text(size=13)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14))

# PCA
# Two PCA with two different transformation
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="dex")

rld <- rlogTransformation(dds, blind = FALSE)
plotPCA(rld, intgroup="dex")

#Distance matrix for sample
library(RColorBrewer)
library(gplots)
library(factoextra)
#assay() used to extracting matrix of normalized values
dists <- get_dist(t(assay(vsdata)),method="pearson")

mat <- as.matrix(dists)
hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat,trace="none",col = rev(hmcol),margin=c(13,13))



#Heatmap of gene expression for 50 better DE gene
library(NMF)
res <- tbl_df(res)
res <- res %>% 
  arrange(padj) %>% 
  inner_join(anno,by=c("row"="ensgene")) %>%
  filter(padj<0.05)
NMF::aheatmap(assay(vsdata)[arrange(res, padj, pvalue)$row[1:25],], 
              labRow=arrange(res, padj, pvalue)$symbol[1:25], 
              scale="row", distfun="pearson", 
              annCol=dplyr::select(airway_metadata, dex, celltype), 
              col=c("green","black","black","red"))
