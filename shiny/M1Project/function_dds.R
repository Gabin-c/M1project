library(dashboardthemes)
library(DESeq2)
library(Biobase)
library(gplots)
library(RColorBrewer)
library(shiny)
library(tidyverse)
library(NMF)
library(ggrepel)

depth <- function(dds,breaksize=1){
  depth <- dds
  depth <- as.data.frame(colSums(depth))
  
  depth$Sample <- row.names(depth)
  
  return(ggplot(depth, aes( x=Sample ,y=depth[,1],color=Sample))+ geom_bar(stat="identity",fill=brewer.pal(n=length(depth$Sample),name="YlGn"),width = breaksize)+labs(title = "Depth of each sample", x="Sample", y="Depth")+theme_bw()+scale_fill_brewer(palette="YlGn"))
  
}#ok

count_distribution <- function(dds, sample,min=0,max=14,breaksize=1){
  counts_dds <- dds
  
  return(ggplot(data=counts_dds, aes(log(counts_dds[,sample]+1))) + 
           geom_histogram(breaks=seq(min,max,breaksize),position="identity",alpha=0.5,fill="darkcyan", color="dodgerblue1")+
           theme_classic()+labs(title=sample, x="Count value (number of read by genes) in log(count+1)",y="Count frequency") 
           
         )
}#ok

dispersion <- function(dds){
  DESeq2::plotDispEsts(dds, main= "Relationship between dispersion and counts means")
  
}
number_of_DE <- function(dds,padj = 0.05){
  res_dif <- results(dds, tidy= TRUE)
  return(table(res_dif$padj <= padj, useNA="always"))
}
#Plotcount
plotcount <- function(dds,gene){
  dds1 <- dds
  dds1[,"name"] <- row.names(dds1)
  return(
    ggplot(dds1, aes(x=dds1[,"name"], y=dds1[,gene])) + 
      geom_point(size=4,aes(colour=factor(name))) + geom_segment(aes(x=dds1[,"name"], 
                                                                     xend=dds1[,"name"], 
                                                                     y=0, 
                                                                     yend=dds1[,gene]),linetype="dotdash")+ 
      theme(axis.text.x = element_blank() )+ 
      labs(title=paste("Count of",gene,  "for each sample"), 
           x="Sample",y="Count")+ guides(color= guide_legend(title = "Sample"))
    
  )
}
#Maplot
maplot <- function(dds,padje=0.05){
  
  
  res_dif <- dds %>% mutate(sig=padj<padje)
  return(ggplot(res_dif, aes(x = baseMean, y = log2FoldChange, col = sig)) + 
           geom_point() + 
           scale_x_log10() +
           geom_hline(yintercept = 0, linetype = "dashed",color = "black") + 
           ggtitle("MA plot") + theme_bw() +
           scale_colour_discrete(name="",
                                 labels=c("Not significative", "Significative", "NA")))
}#ok
#VolcanonPlot
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
             geom_vline(xintercept=0,linetype="dashed", color = "red")+theme(legend.text=element_text(size=9)) )
    
    
  }else{
    res_dif <- res_dif %>% mutate(sig=padj<padje) %>%  arrange(padj)
    return(ggplot(res_dif, aes(x=log2FoldChange, y=-log10(pvalue), col=sig)) +
             geom_point()+
             scale_colour_discrete(name="",
                                   labels=c("Not significative", "Significative", "NA")) +
             geom_vline(xintercept=0,linetype="dashed", color = "red")+guides(colour = guide_legend(override.aes = list(size = 5)))+theme(legend.text=element_text(size=7))) 
    
  }
}
#PCA
pca <- function(dds,intgroup){
  
  
  return(plotPCA(dds, intgroup=intgroup))
  
  
  
  
  
}

#heatMap
clustering_heatmap <- function(dds){
  
  dists <- dist(t(assay(dds)))
  mat <- as.matrix(dists)
  hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(100)
  return(heatmap.2(mat,trace="none",col = rev(hmcol),margin=c(13,13)))
  
}

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
  
  

