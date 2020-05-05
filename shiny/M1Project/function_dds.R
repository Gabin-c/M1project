library(DESeq2)
library(Biobase)
library(gplots)
library(RColorBrewer)
library(shiny)
library(tidyverse)
library(NMF)
library(ggrepel)


### Preamble ----
### All of this function can be use after running a DESeq2 workflow.
### You need a DESeq2 dataset to be able to run it.
###   - DESeq2 : count table after run counts() on yout DESeq2 dataset
###      - depth()
###      - count_distribution()
###      - plotcount()
###   - DESeq2 : results object after run DESeq() on your DESeq2 dataset and results() onr DESeq() object 
###      - maplot()
###      - volcanoplot()
###      - number_of_DE 
###   - DESeq2 : after run DESeq() on your DESeq2 datasetn and do vst() or rlogtransformation() on DESeq() object
###      - pca()
###      - heatmap()
###      - clustering_heatmap()
###      - dispersion()

### Depth plot ----
### depth return a barplot of depth sample with ggplot library  
### depth need two argument
###   - dds which is a count table of a RNAseq experience which have in column : sample, and in row :  gene
###   - breaksize which is width of the bar
depth.plot <- function(dds.count,break.width=1){
  depth <- as.data.frame(colSums(dds.count))
  depth$Sample <- row.names(depth)
  
  return(ggplot(depth, aes( x=Sample ,y=depth[,1]))+ 
           geom_bar(stat="identity",fill=brewer.pal(n=length(depth$Sample),name="YlGn"),width = break.width)+
           labs(title = "Depth of each sample", x="Samples", y="Depth")+theme_bw()+
           theme(plot.title = element_text(face = "bold", size= 18)) +
           theme(axis.title.x = element_text(size=14)) +
           theme(axis.title.y = element_text(size=14)))
}


### Count distribution plot ----
### count_distribution return an histogram of count values distribution in the log(count+1) format for one sample
### count_distribution need five arguments
###   - dds which is a count table of an RNAseq experience which have in column : sample, and in row :  gene
###   - breaksize which is width of histogram bar
###   - sample for which we display the count distribution
###   - min which is the min of x axis 
###   - max which is the max of x axis

count.distribution.plot <- function(dds.count, sample,x.min=0,x.max=14,break.width=1){
  
  return(ggplot(data=dds.count, aes(log(dds.count[,sample]+1))) + 
           geom_histogram(breaks=seq(x.min,x.max,break.width),position="identity",alpha=0.5,fill="darkcyan", color="dodgerblue1")+
           theme_classic() +
           labs(title=sample, x="Counts values (number of reads by gene) in log(count+1)",y="Counts frequencies") +
           theme(plot.title = element_text(face = "bold", size= 18)) +
           theme(axis.title.x = element_text(size=14)) +
           theme(axis.title.y = element_text(size=14))
  )
}


### Dispersion plot ----
### dispersion return plot obtained by DESeq2::plotDispEsts() function of DESeq2 package
### dispersion just need one argument : a DESeq() object
dispersion <- function(dds){
  DESeq2::plotDispEsts(dds, main= "Relationship between dispersion and counts means")
}
### nummber of differemtial express gene ----


number.DE.gene <- function(dds.result,p.val = 0.05){
  tb.DE <- as.data.frame(table(dds.result$padj <= p.val ,useNA="always"))
  colnames(tb.DE) = c("DE","Genes")
  return(tb.DE)
}






### Plotcount ----
plotcount <- function(dds,gene){
  dds1 <- dds
  dds1[,"name"] <- row.names(dds1)
  return(
   ggplot(dds1, aes(x=dds1[,"name"], y=dds1[,gene])) + 
     geom_point(size=4,aes(colour=factor(name))) + 
     geom_segment(aes(x=dds1[,"name"], xend=dds1[,"name"], y=0, yend=dds1[,gene]),linetype="dotdash")+ 
     theme(axis.text.x = element_blank() )+ 
     labs(title=paste("Count of",gene,  "for each sample"),x="Samples",y="Counts")+ 
     guides(color= guide_legend(title = "Sample", override.aes = list(size=5))) +
     theme(plot.title = element_text(face = "bold", size= 18)) +
     theme(axis.title.x = element_text(size=14)) +
     theme(axis.title.y = element_text(size=14)) +
     theme(legend.text=element_text(size=13)) +
     theme(legend.title=element_blank())
  )
}


### Maplot ---- 
maplot <- function(dds,padje=0.05){
  res_dif <- dds %>% mutate(sig=padj<padje)
  return(ggplot(res_dif, aes(x = baseMean, y = log2FoldChange, col = sig)) + 
           geom_point() + 
           scale_x_log10() +
           geom_hline(yintercept = 0, linetype = "dashed",color = "black") + 
           theme_bw() +
           scale_colour_discrete(name="",labels=c("Not significative", "Significative", "NA"))+
           guides(color = guide_legend(override.aes = list(size=5))) +
           theme(legend.text=element_text(size=13))+
           theme(axis.title.x = element_text(size=14)) +
           theme(axis.title.y = element_text(size=14)))
}


### VolcanonPlot ---- 
volcanoPlot <-function(dds, annotation=FALSE,anno,padje=0.05,maxlogF=6,minlogF=0,minlogP=30,count){
  res_dif <- dds
  if(annotation == TRUE){
    res_dif <- res_dif %>% mutate(sig=padj<padje) %>%  arrange(padj) %>%
      inner_join(anno,by=c("row"=count[1]))
    return(ggplot(res_dif, aes(x=log2FoldChange, y=-log10(pvalue), col=sig)) +
             geom_point() +
             ggtitle("Volcano plot labelling top significant genes") +
             geom_text_repel(data = subset(res_dif, (-log10(pvalue) > minlogP | log2FoldChange > maxlogF | log2FoldChange < minlogF)),
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
             theme(axis.title.y = element_text(size=14)))
  }else{
    res_dif <- res_dif %>% mutate(sig=padj<padje) %>%  arrange(padj)
    return(ggplot(res_dif, aes(x=log2FoldChange, y=-log10(pvalue), col=sig)) +
             geom_point()+
             scale_colour_discrete(name="",
                                   labels=c("Not significative", "Significative", "NA")) +
             geom_vline(xintercept=0,linetype="dashed", color = "red") +
             guides(colour = guide_legend(override.aes = list(size = 5))) +
             theme(legend.text=element_text(size=13))+
             theme(axis.title.x = element_text(size=14)) +
             theme(axis.title.y = element_text(size=14))) 
  }
}

### PCA ---- 
pca <- function(dds,intgroup){
  return(
    plotPCA(dds, intgroup=intgroup) +
      theme(axis.title.x = element_text(size=14)) +
      theme(axis.title.y = element_text(size=14)) +
      scale_colour_discrete(name="")+
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      theme(legend.text=element_text(size=13))+
      geom_text_repel(
                      aes(label = colnames(dds)),
                      size = 3,
                      box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.3, "lines"), color = "darkblue")
  )
}

### Distance matrix ----
clustering_heatmap <- function(dds){
  dists <- dist(t(assay(dds)))
  mat <- as.matrix(dists)
  hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(100)
  return(heatmap.2(mat,trace="none",col = rev(hmcol),margin=c(13,13)))
}

### Heatmap ----
heatmap <- function(dds, dds2,annotation=FALSE,anno,padje=0.05,metadata,condition,count,min,max){
  res <- dds
  res <- tbl_df(res)
  if(annotation==TRUE){
    res <- res %>% 
      arrange(padj) %>% 
      inner_join(anno,by=c("row"=count[1])) %>%
      filter(padj<padje)
    
    NMF::aheatmap(assay(dds2)[arrange(res, padj, pvalue)$row[min:max],], 
                  labRow=arrange(res, padj, pvalue)$symbol[min:max], 
                  scale="row", distfun="pearson", 
                  annCol=dplyr::select(metadata, condition), 
                  col=c("green","black","black","red"))
    
  }else{
    res <- res %>% 
      arrange(padj) %>% filter(padj<padje)
    
    NMF::aheatmap(assay(dds2)[arrange(res, padj, pvalue)$row[min:max],], 
                  labRow=arrange(res, padj, pvalue)$symbol[min:max], 
                  scale="row", distfun="pearson", 
                  annCol=dplyr::select(metadata, condition), 
                  col=c("green","black","black","red"))
  }
}

