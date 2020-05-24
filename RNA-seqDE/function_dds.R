library(DESeq2)
library(Biobase)
library(gplots)
library(RColorBrewer)
library(shiny)
library(tidyverse)
library(NMF)
library(ggrepel)
library(factoextra)
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
### count_distribution return an histogram of count values distribution in the log(count+1) (to facilitate vizualisation) format for one sample
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
### number.DE.gene return a table of DE results with number of TRUE,FALSE and NA in function of pvalue
### number.DE.gene need two arguments
###   - dds.result which is a results table obtain by function results() on a DESeq() object
###   - p.val which is pvalue accept un DE analysis, padje is set at 0.05
number.DE.gene <- function(dds.result,p.val = 0.05){
  up_regulated <- dds.result %>% filter(padj <= p.val & log2FoldChange > 0) %>% nrow()
  down_regulated <- dds.result %>% filter(padj <= p.val & log2FoldChange < 0) %>% nrow()
  tb <-table(dds.result$padj <= p.val ,useNA="always")
  tb.DE <- data.frame("No DE" = tb[1], "Down regulated" = down_regulated, "Up regulated" = up_regulated, "NA" = tb[3]  )
  row.names(tb.DE) <- ""
  return(tb.DE)
}

### Plotcount ----
### plotcount return a point plot of read count for one gene in each sample
### plotcount need two argument
###   - dds.count which is a count table of a RNAseq experience which have in column : sample, and in row : gene
###   - gene which is gene which we want to display read count for each sample

gene.count.plot <- function(dds.count,gene){
  dds1 <- dds.count
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
### maplot return a MA plot of DE
### maplot need two arguments
###   - dds.results which is a results table obtain by function results() on a DESeq() object
###   - p.val which is pvalue accept un DE analysis, padje is set at 0.05
ma.plot <- function(dds.results,p.val=0.05){
  dds.res <- dds.results %>% mutate(sig=padj<p.val)
  return(ggplot(dds.res, aes(x = baseMean, y = log2FoldChange, col = sig)) + 
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
### volcano.plot generate a volcano plot of DE
### volcano.plot need 8 arguments
###   - dds.results which is a results table obtain by function results() on a DESeq() object
###   - is.anno if there is an annotation file
###   - anno an annotation 
###   - p.val which is pvalue accept un DE analysis, padje is set at 0.05
###   - count.tb which is a count table of a RNAseq experience which have in column : sample, and in row : gene
###   - maxlogF : max Fold change in x axis which we set gene ID on geom point
###   - minlogF : max Fold change in x axis which we set gene ID on geom point
###   - minlogP : min Log10 in y which we set gene ID on geom point
volcano.plot <-function(dds.results, is.anno=FALSE,anno,p.val=0.05,maxlogF=6,minlogF=0,minlogP=30,count.tb){
  if(is.anno == TRUE){
    dds.res <- dds.results %>% mutate(sig=padj<p.val) %>%  arrange(padj) %>%
      inner_join(anno,by=c("row"=count.tb[1]))
    return(ggplot(dds.res, aes(x=log2FoldChange, y=-log10(padj), col=sig)) +
             geom_point() +
             ggtitle("Volcano plot labelling top significant genes") +
             geom_text_repel(data = subset(dds.res, (-log10(padj) > minlogP | log2FoldChange > maxlogF | log2FoldChange < minlogF)),
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
    dds.res <- dds.results %>% mutate(sig=padj<p.val) %>%  arrange(padj)
    return(ggplot(dds.res, aes(x=log2FoldChange, y=-log10(padj), col=sig)) +
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
### pca.plot generate a pca plot
### pca.plot need 2 arguments
###   - dds.resTransf which is an DESeq2 transformate object with vst() or rLogtransformation()
pca.plot <- function(dds.resTransf,intgroup){
  return(
    plotPCA(dds.resTransf, intgroup=intgroup) +
      theme(axis.title.x = element_text(size=14)) +
      theme(axis.title.y = element_text(size=14)) +
      scale_colour_discrete(name="")+
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      theme(legend.text=element_text(size=13))+
      geom_text_repel(
        aes(label = colnames(dds.resTransf)),
        size = 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"), color = "darkblue")
  )
}

### Distance matrix ----
### distance.matrix.heatmap generate a heatmap of dist between different sample
### distance.matrix.heatmap need just one argument
###   - dds.resTransf which is an DESeq2 transformate object with vst() or rLogtransformation()
distance.matrix.heatmap <- function(dds.resTransf){
  dists <- get_dist(t(assay(dds.resTransf)),method = "pearson")
  mat <- as.matrix(dists)
  hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(100)
  return(heatmap.2(mat,trace="none",col = rev(hmcol),margin=c(13,13)))
}

### Heatmap ----
### gene.expression.heatmap generate a gene expression distance heatmap
### gene.expression.heatmap need 10 argument
###   - dds.results which is a results table obtain by function results() on a DESeq() object
###   - is.anno if there is an annotation file
###   - anno an annotation 
###   - p.val which is pvalue accept un DE analysis, padje is set at 0.05
###   - count.tb which is a count table of a RNAseq experience which have in column : sample, and in row : gene
###   - metadata which a design table of an RNA-seq experience
###   - dds.resTransf which is an DESeq2 transformate object with vst() or rLogtransformation()
###   - condition : design formula of RNA-seq experience
###   - min : num of row in results table starting from the top
###   - max : num of row in results table starting from the top
gene.expression.heatmap <- function(dds.results, dds.resTransf,is.anno=FALSE,anno,p.val=0.05,metadata,condition,count.tb,min,max){
  res <- tbl_df(dds.results)
  if(is.anno==TRUE){
    res <- res %>% 
      arrange(padj) %>% 
      inner_join(anno,by=c("row"=count.tb[1])) %>%
      filter(padj<p.val)
    
    NMF::aheatmap(assay(dds.resTransf)[arrange(res, padj, pvalue)$row[min:max],], 
                  labRow=arrange(res, padj, pvalue)$symbol[min:max], 
                  scale="row", distfun="pearson", 
                  annCol=dplyr::select(metadata, condition), 
                  col=c("green","black","black","red"),treeheight = c(500,100))
    
  }else{
    res <- res %>% 
      arrange(padj) %>% filter(padj<p.val)
    
    NMF::aheatmap(assay(dds.resTransf)[arrange(res, padj, pvalue)$row[min:max],], 
                  labRow=arrange(res, padj, pvalue)$symbol[min:max], 
                  scale="row", distfun="pearson", 
                  annCol=dplyr::select(metadata, condition), 
                  col=c("green","black","black","red"),treeheight = c(500,100))
  }
}

## Theme ----
shinyDashboardThemeDIY <- function(
  appFontFamily, appFontColor, logoBackColor, bodyBackColor, headerButtonBackColor, headerButtonIconColor,
  headerButtonBackColorHover, headerButtonIconColorHover, headerBackColor, headerBoxShadowColor,
  headerBoxShadowSize, sidebarBackColor, sidebarPadding, sidebarShadowRadius, sidebarShadowColor,
  sidebarMenuBackColor, sidebarMenuPadding, sidebarMenuBorderRadius, sidebarUserTextColor, sidebarSearchBackColor,
  sidebarSearchIconColor, sidebarSearchBorderColor,  sidebarTabTextColor, sidebarTabTextSize, sidebarTabBorderStyle,
  sidebarTabBorderColor, sidebarTabBorderWidth, sidebarTabBackColorSelected, sidebarTabTextColorSelected,
  sidebarTabRadiusSelected, sidebarTabTextColorHover, sidebarTabBackColorHover, sidebarTabBorderStyleHover,
  sidebarTabBorderColorHover, sidebarTabBorderWidthHover, sidebarTabRadiusHover, boxBackColor, boxBorderRadius,
  boxShadowSize, boxShadowColor, boxTitleSize, boxDefaultColor, boxPrimaryColor, boxSuccessColor, boxWarningColor,
  boxDangerColor, tabBoxTabColor, tabBoxTabTextSize, tabBoxTabTextColor, tabBoxTabTextColorSelected, tabBoxBackColor,
  tabBoxHighlightColor, tabBoxBorderRadius, buttonBackColor, buttonTextColor, buttonBorderColor, buttonBorderRadius,
  buttonBackColorHover, buttonTextColorHover, buttonBorderColorHover, buttonHeight = 34, buttonPadding = "6px 12px",
  textboxBackColor, textboxBorderColor, textboxBorderRadius, textboxBackColorSelect, textboxBorderColorSelect,
  textboxHeight = 34, textboxPadding = "6px 12px", tableBackColor, tableBorderColor,
  tableBorderTopSize, tableBorderRowSize, primaryFontColor = "auto", successFontColor = "auto",
  warningFontColor = "auto", dangerFontColor = "auto", infoFontColor = "auto", boxInfoColor = "auto", dataFontColor
) {
  
  htmltools::tags$head(
    
    htmltools::tags$style(
      
      htmltools::HTML(
        
        paste0(
          
          '
          /* font: google import [OPTIONAL] */
          /* @import url("//fonts.googleapis.com/css?family=',"Roboto",'"); */
          /* font */
          body, label, input, button, select, box,
          .h1, .h2, .h3, .h4, .h5, h1, h2, h3, h4, h5 {
          font-family: "',appFontFamily,'";
          color: ', appFontColor, ';
          }
          /* font: fix for h6 */
          /* messes up sidebar user section if included above */
          .h6, h6 {
          font-family: "',appFontFamily,'";
          }
          /* sidebar: logo */
          .skin-blue .main-header .logo {
          background: ', logoBackColor, ';
          }
          /* sidebar: logo hover */
          .skin-blue .main-header .logo:hover {
          background: ', logoBackColor, ';
          }
          /* sidebar: collapse button  */
          .skin-blue .main-header .navbar .sidebar-toggle {
          background: ', headerButtonBackColor, ';
          color:', headerButtonIconColor, ';
          }
          /* sidebar: collapse button hover */
          .skin-blue .main-header .navbar .sidebar-toggle:hover {
          background: ', headerButtonBackColorHover, ';
          color:', headerButtonIconColorHover, ';
          }
          /* header */
          .skin-blue .main-header .navbar {
          background: ', headerBackColor, ';
          box-shadow: ', headerBoxShadowSize, ' ', headerBoxShadowColor ,';
          }
          /* sidebar*/
          .skin-blue .main-sidebar {
          background: ', sidebarBackColor, ';
          box-shadow: ', sidebarShadowRadius, " ", sidebarShadowColor, ';
          padding-left: ', sidebarPadding, 'px;
          padding-right: ', sidebarPadding, 'px;
          /* padding-top: 60px; */
          }
          /* sidebar menu */
          .main-sidebar .user-panel, .sidebar-menu, .sidebar-menu>li.header {
          white-space: nowrap;
          background: ', sidebarMenuBackColor, ';
          padding: ', sidebarMenuPadding, 'px;
          border-radius: ', sidebarMenuBorderRadius, 'px;
          }
          /* fix for user panel */
          .user-panel>.info>p, .skin-blue .user-panel>.info {
          color: ', sidebarUserTextColor, ';
          font-size: 12px;
          font-weight: normal;
          }
          section.sidebar .user-panel {
          padding: 10px;
          }
          /* sidebar: tabs */
          .skin-blue .main-sidebar .sidebar .sidebar-menu a {
          color: ', sidebarTabTextColor, ';
          font-size: ', sidebarTabTextSize ,'px;
          border-style: ', sidebarTabBorderStyle, ';
          border-color: ', sidebarTabBorderColor, ';
          border-width: ', sidebarTabBorderWidth, 'px;
          }
          /* sidebar: tab selected */
          .skin-blue .main-sidebar .sidebar .sidebar-menu .active > a {
          color: ', sidebarTabTextColorSelected, ';
          font-size: ', sidebarTabTextSize, 'px;
          border-radius: ', sidebarTabRadiusSelected, ';
          border-style: ', sidebarTabBorderStyleHover, ';
          border-color: ', sidebarTabBorderColorHover, ';
          border-width: ', sidebarTabBorderWidthHover, 'px;
          }
          .skin-blue .sidebar-menu > li:hover > a,
          .skin-blue .sidebar-menu > li.active > a {
          color: ', sidebarTabTextColorSelected, ';
          background: ', sidebarTabBackColorSelected, ';
          border-radius: ', sidebarTabRadiusHover, ';
          }
          /* sidebar: tab hovered */
          .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover {
          background: '
          ,sidebarTabBackColorHover, ';'
          ,'color: ', sidebarTabTextColorHover, ';
          font-size: ', sidebarTabTextSize ,'px;
          border-style: ', sidebarTabBorderStyleHover, ';
          border-color: ', sidebarTabBorderColorHover, ';
          border-width: ', sidebarTabBorderWidthHover, 'px;
          border-radius: ', sidebarTabRadiusHover, ';
          }
          /* sidebar: subtab */
          .skin-blue .sidebar-menu > li > .treeview-menu {
          margin: 0px;
          background: ', sidebarMenuBackColor, ';
          }
          .skin-blue .treeview-menu > li > a {
          background: ', sidebarMenuBackColor, ';
          }
          /* sidebar: subtab selected */
          .skin-blue .treeview-menu > li.active > a,
          .skin-blue .treeview-menu > li > a:hover {
          background: ', sidebarTabBackColorSelected, ';
          }
          /* sidebar: search text area */
          .skin-blue .sidebar-form input[type=text] {
          background: ', sidebarSearchBackColor, ';
          color: ', appFontColor, ';
          border-radius: ', textboxBorderRadius, 'px 0px 0px ', textboxBorderRadius, 'px;
          border-color: ', sidebarSearchBorderColor, ';
          border-style: solid none solid solid;
          }
          /* sidebar: search button */
          .skin-blue .input-group-btn > .btn {
          background: ', sidebarSearchBackColor, ';
          color: ', sidebarSearchIconColor, ';
          border-radius: 0px ', textboxBorderRadius, 'px ', textboxBorderRadius, 'px 0px;
          border-style: solid solid solid none;
          border-color: ', sidebarSearchBorderColor, ';
          }
          /* sidebar form */
          .skin-blue .sidebar-form {
          border-radius: 0px;
          border: 0px none rgb(255,255,255);
          margin: 10px;
          }
          /* body */
          .content-wrapper, .right-side {
          background: ', bodyBackColor, ';
          }
          /* box */
          .box {
          background: ', boxBackColor, ';
          border-radius: ', boxBorderRadius, 'px;
          box-shadow: ', boxShadowSize, ' ', boxShadowColor, ';
          }
          /* box: title */
          .box-header .box-title {
          font-size: ', boxTitleSize, 'px;
          }
          /* tabbox: title */
          .nav-tabs-custom>.nav-tabs>li.header {
          color: ', appFontColor, ';
          font-size: ', boxTitleSize, 'px;
          }
          /* tabbox: tab color */
          .nav-tabs-custom, .nav-tabs-custom .nav-tabs li.active:hover a, .nav-tabs-custom .nav-tabs li.active a {
          background: ', tabBoxTabColor, ';
          color: ', appFontColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          .nav-tabs-custom {
          box-shadow: ', boxShadowSize, ' ', boxShadowColor, ';
          }
          /* tabbox: active tab bg */
          .nav-tabs-custom>.nav-tabs>li.active {
          border-radius: ', tabBoxBorderRadius, 'px;
          border-top-color: ', tabBoxHighlightColor, ';
          # box-shadow: ', boxShadowSize, ' ', boxShadowColor, ';
          }
          /* tabbox: font color */
          .nav-tabs-custom>.nav-tabs>li.active:hover>a, .nav-tabs-custom>.nav-tabs>li.active>a {
          border-bottom-color: ', tabBoxTabColor, ';
          border-top-color: ', tabBoxHighlightColor, ';
          border-right-color: ', tabBoxHighlightColor, ';
          border-left-color: ', tabBoxHighlightColor, ';
          color: ', tabBoxTabTextColorSelected, ';
          font-size: ', tabBoxTabTextSize, 'px;
          border-radius: ', tabBoxBorderRadius, 'px;
          }
          /* tabbox: inactive tabs background */
          .nav-tabs-custom>.nav-tabs>li>a {
          color: ', tabBoxTabTextColor, ';
          font-size: ', tabBoxTabTextSize, 'px;
          }
          /* tabbox: top area back color */
          .nav-tabs-custom, .nav-tabs-custom>.tab-content, .nav-tabs-custom>.nav-tabs {
          border-bottom-color: ', tabBoxHighlightColor, ';
          background: ', tabBoxBackColor, ';
          }
          /* tabbox: top area rounded corners */
          .nav-tabs-custom>.nav-tabs {
          margin: 0;
          border-radius: ', tabBoxBorderRadius, 'px;
          }
          /* infobox */
          .info-box {
          background: ', boxBackColor, ';
          border-radius: ', boxBorderRadius, 'px;
          box-shadow: ', boxShadowSize, ' ', boxShadowColor, ';
          }
          /* valuebox */
          .small-box {
          border-radius: ', boxBorderRadius, 'px;
          box-shadow: ', boxShadowSize, ' ', boxShadowColor, ';
          }
          /* valuebox: main font color */
          .small-box h3, .small-box p {
          color: rgb(255,255,255)
          }
          /* box: default color */
          .box.box-solid>.box-header, .box>.box-header {
          color: ', appFontColor, ';
          }
          .box.box-solid>.box-header {
          border-radius: ', boxBorderRadius, 'px;
          }
          .box.box-solid, .box {
          border-radius: ', boxBorderRadius, 'px;
          border-top-color: ', boxDefaultColor, ';
          }
          /* box: info color */
          .box.box-solid.box-info>.box-header h3, .box.box-info>.box-header h3 {
          color: ', infoFontColor, ';
          }
          .box.box-solid.box-info>.box-header {
          background: ', boxInfoColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          .box.box-solid.box-info, .box.box-info {
          border-color: ', boxInfoColor, ';
          border-left-color: ', boxInfoColor, ';
          border-right-color: ', boxInfoColor, ';
          border-top-color: ', boxInfoColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          /* box: primary color */
          .box.box-solid.box-primary>.box-header h3, .box.box-primary>.box-header h3 {
          color: ', primaryFontColor, ';
          }
          .box.box-solid.box-primary>.box-header {
          background: ', boxPrimaryColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          .box.box-solid.box-primary, .box.box-primary {
          border-color: ', boxPrimaryColor, ';
          border-left-color: ', boxPrimaryColor, ';
          border-right-color: ', boxPrimaryColor, ';
          border-top-color: ', boxPrimaryColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          /* box: success color */
          .box.box-solid.box-success>.box-header h3, .box.box-success>.box-header h3 {
          color: ', successFontColor, ';
          }
          .box.box-solid.box-success>.box-header {
          background: ', boxSuccessColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          .box.box-solid.box-success, .box.box-success {
          border-color: ', boxSuccessColor, ';
          border-left-color: ', boxSuccessColor, ';
          border-right-color: ', boxSuccessColor, ';
          border-top-color: ', boxSuccessColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          /* box: warning color */
          .box.box-solid.box-warning>.box-header h3, .box.box-warning>.box-header h3 {
          color: ', warningFontColor, ';
          }
          .box.box-solid.box-warning>.box-header {
          background: ', boxWarningColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          .box.box-solid.box-warning, .box.box-warning {
          border-color: ', boxWarningColor, ';
          border-left-color: ', boxWarningColor, ';
          border-right-color: ', boxWarningColor, ';
          border-top-color: ', boxWarningColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          /* box: danger color */
          .box.box-solid.box-danger>.box-header h3, .box.box-danger>.box-header h3 {
          color: ', dangerFontColor, ';
          }
          .box.box-solid.box-danger>.box-header {
          background: ', boxDangerColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          .box.box-solid.box-danger, .box.box-danger {
          border-color: ', boxDangerColor, ';
          border-left-color: ', boxDangerColor, ';
          border-right-color: ', boxDangerColor, ';
          border-top-color: ', boxDangerColor, ';
          border-radius: ', boxBorderRadius, 'px;
          }
          /* button */
          .btn-default {
          background: ', buttonBackColor, ';
          color: ', buttonTextColor, ';
          border-color: ', buttonBorderColor, ';
          border-radius: ', buttonBorderRadius, 'px;
          height: ', buttonHeight, 'px;
          padding: ', buttonPadding, ';
          }
          /* button: hover */
          .btn-default:hover {
          background: ', buttonBackColorHover, ';
          color: ', buttonTextColorHover, ';
          border-color: ', buttonBorderColorHover, ';
          }
          /* button: focus */
          .btn-default:focus, .action-button:focus {
          background: ', buttonBackColor, ';
          color: ', buttonTextColor, ';
          border-color: ', buttonBorderColor, ';
          }
          /* button: active */
          .btn-default:active, .action-button:active {
          background: ', buttonBackColor, ';
          color: ', buttonTextColor, ';
          border-color: ', buttonBorderColor, ';
          }
          /* button: visited */
          .btn-default:visited {
          background: ', buttonBackColor, ';
          color: ', buttonTextColor, ';
          border-color: ', buttonBorderColor, ';
          }
          /* textbox */
          .form-control, .selectize-input, .selectize-control.single .selectize-input {
          background: ', textboxBackColor, ';
          color: ', appFontColor, ';
          border-color: ', textboxBorderColor, ';
          border-radius: ', textboxBorderRadius, 'px;
          height: ', textboxHeight, 'px;
          min-height: ', textboxHeight, 'px;
          padding: ', textboxPadding, ';
          }
          /* textbox: selected */
          .form-control:focus, .selectize-input.focus {
          color: ', appFontColor, ';
          background: ', textboxBackColorSelect, ';
          border-color: ', textboxBorderColorSelect, ';
          -webkit-box-shadow: inset 0px 0px 0px, 0px 0px 0px;
          box-shadow: inset 0px 0px 0px, 0px 0px 0px;
          }
          /* multi-row selectize input */
          .selectize-control.multi .selectize-input.has-items {
          height: auto;
          }
          /* verbatim text output */
          .qt pre, .qt code {
          font-family: ', appFontFamily, ' !important;
          }
          pre {
          color: ', appFontColor, ';
          background-color: ', textboxBackColor, ';
          border: 1px solid ', textboxBorderColor, ';
          border-radius: ', textboxBorderRadius, 'px;
          }
          /* drop-down menu */
          .selectize-dropdown, .selectize-dropdown.form-control {
          background: ', textboxBackColor, ';
          border-radius: 4px;
          }
          /* table */
          .table {
          background: ', tableBackColor, ';
          border-radius: ', textboxBorderRadius, 'px;
          }
          /* table: row border color*/
          .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
          border-top: ', tableBorderRowSize, 'px solid ', tableBorderColor, ';
          }
          /* table: top border color*/
          .table>thead>tr>th {
          border-bottom: ', tableBorderTopSize, 'px solid ', tableBorderColor, ';
          }
          /* table: hover row */
          .table-hover>tbody>tr:hover {
          background-color: ', tableBorderColor, ';
          }
          /* table: stripe row */
          .table-striped>tbody>tr:nth-of-type(odd) {
          background-color: ', tableBorderColor, ';
          }
          /* table: body colour */
          table.dataTable tbody tr {
          background-color: ', tableBackColor, ' !important;
          }
          /* table: footer border colour */
          table.dataTable {
          border: 0px !important;
          }
          /* datatable: selected row */
          table.dataTable tr.selected td, table.dataTable td.selected {
          background-color: ', boxSuccessColor, ' !important;
          color: rgb(0,0,0) !important;
          }
          /* datatable: hover row */
          table.dataTable tr.hover td, table.dataTable td.hover {
          background-color: ', tableBorderColor, ' !important;
          }
          table.dataTable.hover tbody tr:hover, table.dataTable.display tbody tr:hover {
          background-color: ', tableBorderColor, ' !important;
          }
          table.dataTable.row-border tbody th, table.dataTable.row-border tbody td,
          table.dataTable.display tbody th, table.dataTable.display tbody td {
          border-top: 1px solid ', tableBorderColor, ' !important;
          }
          /* datatable: stripe row */
          table.dataTable.stripe tbody tr.odd, table.dataTable.display tbody tr.odd {
          background-color: ', tableBorderColor, ' !important;
          }
          /* datatable: page control */
          .dataTables_wrapper .dataTables_paginate .paginate_button {
          background-color: ', tableBorderColor, ' !important;
          color: ', dataFontColor, ' !important;
          }
          /* datatable: table info */
          .dataTables_wrapper .dataTables_paginate .paginate_button.disabled,
          .dataTables_wrapper .dataTables_paginate .paginate_button.disabled:hover,
          .dataTables_wrapper .dataTables_paginate .paginate_button.disabled:active {
          color: ', dataFontColor, ' !important;
          }
          .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter,
          .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,
          .dataTables_wrapper .dataTables_paginate {
          color: ', dataFontColor, ' !important;
          }
          .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter,
          .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,
          .dataTables_wrapper .dataTables_paginate {
          color: ', dataFontColor, ' !important;
          }
          thead {color: ', dataFontColor, '; } tbody {color: ', dataFontColor, '; }
          
          /* datatable search box */
          .dataTables_wrapper .dataTables_filter input {
          background-color: ', textboxBackColor, ';
          border: 1px solid ', textboxBorderColor, ';
          border-radius: ', textboxBorderRadius, 'px;
          }
          thead {color: ', dataFontColor, '; } tbody {color: ', dataFontColor, '; }
          /* notification and progress bar */
          .progress-bar {
          background-color: ', boxSuccessColor, ';
          }
          .shiny-notification {
          height: 80px;
          font-family: ', appFontFamily, ';
          font-size: 15px;
          color: rgb(0,0,0);
          background-color: rgb(225,225,225);
          border-color: rgb(205,205,205);
          border-radius: 10px;
          margin-left: -450px !important;
          }
          /* horizontal divider line */
          hr {
          border-top: 1px solid rgb(215,215,215);
          }
          /* modal */
          .modal-body {
          background-color: ', boxBackColor, ';
          }
          .modal-footer {
          background-color: ', boxBackColor, ';
          }
          .modal-header {
          background-color: ', boxBackColor, ';
          }
          '
        )
      )
    )
  )
}

theme_grey_dark2 <- shinyDashboardThemeDIY(
  ### general
  appFontFamily = "Arial"
  ,appFontColor = "rgb(205,205,205)"
  ,primaryFontColor = "rgb(255,255,255)"
  ,infoFontColor = "rgb(255,255,255)"
  ,successFontColor = "rgb(255,255,255)"
  ,warningFontColor = "rgb(255,255,255)"
  ,dangerFontColor = "rgb(255,255,255)"
  ,bodyBackColor = "rgb(45,55,65)"
  
  ### header
  ,logoBackColor = "rgb(70,80,90)"
  
  ,headerButtonBackColor = "rgb(70,80,90)"
  ,headerButtonIconColor = "rgb(25,35,45)"
  ,headerButtonBackColorHover = "rgb(40,50,60)"
  ,headerButtonIconColorHover = "rgb(0,0,0)"
  
  ,headerBackColor = "rgb(70,80,90)"
  ,headerBoxShadowColor = ""
  ,headerBoxShadowSize = "0px 0px 0px"
  
  ### sidebar
  ,sidebarBackColor = "rgb(52,62,72)"
  ,sidebarPadding = 0
  
  ,sidebarMenuBackColor = "transparent"
  ,sidebarMenuPadding = 0
  ,sidebarMenuBorderRadius = 0
  
  ,sidebarShadowRadius = ""
  ,sidebarShadowColor = "0px 0px 0px"
  
  ,sidebarUserTextColor = "rgb(205,205,205)"
  
  ,sidebarSearchBackColor = "rgb(45,55,65)"
  ,sidebarSearchIconColor = "rgb(153,153,153)"
  ,sidebarSearchBorderColor = "rgb(45,55,65)"
  
  ,sidebarTabTextColor = "rgb(205,205,205)"
  ,sidebarTabTextSize = 14
  ,sidebarTabBorderStyle = "none"
  ,sidebarTabBorderColor = "none"
  ,sidebarTabBorderWidth = 0
  
  ,sidebarTabBackColorSelected = "rgb(70,80,90)"
  ,sidebarTabTextColorSelected = "rgb(255,255,255)"
  ,sidebarTabRadiusSelected = "5px"
  
  ,sidebarTabBackColorHover = "rgb(55,65,75)"
  ,sidebarTabTextColorHover = "rgb(255,255,255)"
  ,sidebarTabBorderStyleHover = "none"
  ,sidebarTabBorderColorHover = "none"
  ,sidebarTabBorderWidthHover = 0
  ,sidebarTabRadiusHover = "5px"
  
  ### boxes
  ,boxBackColor = "rgb(52,62,72)"
  ,boxBorderRadius = 5
  ,boxShadowSize = "0px 0px 0px"
  ,boxShadowColor = ""
  ,boxTitleSize = 16
  ,boxDefaultColor = "rgb(52,62,72)"
  ,boxPrimaryColor = "rgb(200,200,200)"
  ,boxInfoColor = "rgb(80,95,105)"
  ,boxSuccessColor = "rgb(155,240,80)"
  ,boxWarningColor = "rgb(240,80,210)"
  ,boxDangerColor = "rgb(240,80,80)"
  
  ,tabBoxTabColor = "rgb(52,62,72)"
  ,tabBoxTabTextSize = 14
  ,tabBoxTabTextColor = "rgb(205,205,205)"
  ,tabBoxTabTextColorSelected = "rgb(205,205,205)"
  ,tabBoxBackColor = "rgb(52,62,72)"
  ,tabBoxHighlightColor = "rgb(70,80,90)"
  ,tabBoxBorderRadius = 5
  
  ### inputs
  ,buttonBackColor = "rgb(230,230,230)"
  ,buttonTextColor = "rgb(0,0,0)"
  ,buttonBorderColor = "rgb(50,50,50)"
  ,buttonBorderRadius = 5
  
  ,buttonBackColorHover = "rgb(180,180,180)"
  ,buttonTextColorHover = "rgb(50,50,50)"
  ,buttonBorderColorHover = "rgb(50,50,50)"
  
  ,textboxBackColor = "rgb(68,80,90)"
  ,textboxBorderColor = "rgb(76,90,103)"
  ,textboxBorderRadius = 5
  ,textboxBackColorSelect = "rgb(80,90,100)"
  ,textboxBorderColorSelect = "rgb(255,255,255)"
  
  ### tables
  ,tableBackColor = "rgb(52,62,72)"
  ,tableBorderColor = "rgb(70,80,90)"
  ,tableBorderTopSize = 1
  ,tableBorderRowSize = 1
  
  ### datatables
  ,dataFontColor = "rgb(255, 255, 255)"
)

theme_grey_light2 <- shinyDashboardThemeDIY(
  
  ### general
  appFontFamily = "Arial"
  ,appFontColor = "rgb(45,45,45)"
  ,primaryFontColor = "rgb(15,15,15)"
  ,infoFontColor = "rgb(15,15,15)"
  ,successFontColor = "rgb(15,15,15)"
  ,warningFontColor = "rgb(15,15,15)"
  ,dangerFontColor = "rgb(15,15,15)"
  ,bodyBackColor = "rgb(240,240,240)"
  
  ### header
  ,logoBackColor = "rgb(120,120,120)"
  
  ,headerButtonBackColor = "rgb(120,120,120)"
  ,headerButtonIconColor = "rgb(220,220,220)"
  ,headerButtonBackColorHover = "rgb(100,100,100)"
  ,headerButtonIconColorHover = "rgb(60,60,60)"
  
  ,headerBackColor = "rgb(120,120,120)"
  ,headerBoxShadowColor = "#dfdfdf"
  ,headerBoxShadowSize = "3px 5px 5px"
  
  ### sidebar
  ,sidebarBackColor = "rgb(255,255,255)"
  ,sidebarPadding = 0
  
  ,sidebarMenuBackColor = "transparent"
  ,sidebarMenuPadding = 0
  ,sidebarMenuBorderRadius = 0
  
  ,sidebarShadowRadius = "3px 5px 5px"
  ,sidebarShadowColor = "#dfdfdf"
  
  ,sidebarUserTextColor = "rgb(115,115,115)"
  
  ,sidebarSearchBackColor = "rgb(240,240,240)"
  ,sidebarSearchIconColor = "rgb(100,100,100)"
  ,sidebarSearchBorderColor = "rgb(220,220,220)"
  
  ,sidebarTabTextColor = "rgb(100,100,100)"
  ,sidebarTabTextSize = 14
  ,sidebarTabBorderStyle = "none"
  ,sidebarTabBorderColor = "none"
  ,sidebarTabBorderWidth = 0
  
  ,sidebarTabBackColorSelected = "rgb(230,230,230)"
  ,sidebarTabTextColorSelected = "rgb(0,0,0)"
  ,sidebarTabRadiusSelected = "0px"
  
  ,sidebarTabBackColorHover = "rgb(245,245,245)"
  ,sidebarTabTextColorHover = "rgb(0,0,0)"
  ,sidebarTabBorderStyleHover = "none solid none none"
  ,sidebarTabBorderColorHover = "rgb(200,200,200)"
  ,sidebarTabBorderWidthHover = 4
  ,sidebarTabRadiusHover = "0px"
  
  ,boxBackColor = "rgb(248,248,248)"
  ,boxBorderRadius = 5
  ,boxShadowSize = "none"
  ,boxShadowColor = ""
  ,boxTitleSize = 18
  ,boxDefaultColor = "rgb(225,225,225)"
  ,boxPrimaryColor = "rgb(95,155,213)"
  ,boxInfoColor = "rgb(180,180,180)"
  ,boxSuccessColor = "rgb(112,173,71)"
  ,boxWarningColor = "rgb(237,125,49)"
  ,boxDangerColor = "rgb(232,76,34)"
  
  ,tabBoxTabColor = "rgb(248,248,248)"
  ,tabBoxTabTextSize = 14
  ,tabBoxTabTextColor = "rgb(100,100,100)"
  ,tabBoxTabTextColorSelected = "rgb(45,45,45)"
  ,tabBoxBackColor = "rgb(248,248,248)"
  ,tabBoxHighlightColor = "rgb(200,200,200)"
  ,tabBoxBorderRadius = 5
  
  ### inputs
  ,buttonBackColor = "rgb(215,215,215)"
  ,buttonTextColor = "rgb(45,45,45)"
  ,buttonBorderColor = "rgb(150,150,150)"
  ,buttonBorderRadius = 5
  
  ,buttonBackColorHover = "rgb(190,190,190)"
  ,buttonTextColorHover = "rgb(0,0,0)"
  ,buttonBorderColorHover = "rgb(150,150,150)"
  
  ,textboxBackColor = "rgb(255,255,255)"
  ,textboxBorderColor = "rgb(118,118,118)"
  ,textboxBorderRadius = 5
  ,textboxBackColorSelect = "rgb(245,245,245)"
  ,textboxBorderColorSelect = "rgb(108,108,108)"
  
  ### tables
  ,tableBackColor = "rgb(248,248,248)"
  ,tableBorderColor = "rgb(238,238,238)"
  ,tableBorderTopSize = 1
  ,tableBorderRowSize = 1
  
  ### datatables
  ,dataFontColor = "rgb(108,108,108)"
  
)

