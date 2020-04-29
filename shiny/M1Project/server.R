### Server  ----
server <- function(input, output,session) {
  dds <- reactiveValues()
  ### Increase the authorized size for upload ----
  options(shiny.maxRequestSize=30*1024^2)
  
  ### Import the count ----
  count_table <- reactive({
    req(input$file)
    counttable <- read.csv(input$file$datapath, sep = input$sepcount)
  })
  ### Display the count file ----
  output$table <- DT::renderDataTable(count_table(), options = list(pageLength = 20, autoWidth = FALSE,scrollX = TRUE, scrollY = '300px'))
  
  
  
  ### Import the metadata file ---- 
  metadata <- reactive({
    req(input$file2)
    meta_table <- read.csv(input$file2$datapath, sep = input$sepmetadata, row.names=NULL)
  })
  ### Display the metadata file ----
  output$table2 <- DT::renderDataTable(metadata(),options = list(pageLength = 20, autoWidth = FALSE,scrollX = TRUE, scrollY = '300px'))
  
  ### Design condition for DESeq2 ----
  observeEvent(input$file2,{
    updateTextInput(session,"condition", value = paste("~ ",paste(colnames(metadata()), collapse=" + ")))
  })
  
  ### Import annotation file ----
  observeEvent(input$annotation, {
    if(input$annotation== TRUE){
      updateTabsetPanel(session, "params", selected = "annotation")
    }else{
      updateTabsetPanel(session, "params", selected = "nothing")
    }
  })
  anno <- reactive({
    req(input$file3)
    read.csv(input$file3$datapath, sep = input$sepanno)
  })
  output$table3 <- DT::renderDataTable(anno(),options = list(pageLength = 20, autoWidth = FALSE,scrollX = TRUE, scrollY = '300px'))
  
  ### Display parameters for volcano
  observeEvent(input$annotation3, {
    if(input$annotation3== TRUE){
      updateTabsetPanel(session, "param_volc", selected = "Yes")
    }else{
      updateTabsetPanel(session, "param_volc", selected = "No")
    }
  })
  
  
  
  
  ### DDS object  ---- 
  observeEvent(input$deseq2,{
    req(input$deseq2)
    waiter <- waiter::Waiter$new(html = spin_3circles())
    waiter$show()
    on.exit(waiter$hide())
    
    ### Running DESeq1
    dds$dds <- DESeqDataSetFromMatrix(count_table(),colData=metadata(),design=as.formula(input$condition), tidy=TRUE)
    dds$DESeq2 <- DESeq(dds$dds)
    dds$results <- results(dds$DESeq2,tidy=TRUE)
    
    ### Display count table after dds
    output$table4 <- DT::renderDataTable(counts(dds$dds), options = list(pageLength = 5))
    
    ### Choices for count distribution
    updateSelectInput(session,"sample",choices = metadata()[,1])
    
    ### Choices for count by gene
    updateSelectInput(session,"gene",choices = count_table()[,1])
    
    ### Choices for PCA
    updateSelectInput(session,"conditionpca",choices = colnames(metadata()))
    
    ### Choices heatmap
    updateSelectInput(session,"conditionheatmap",choices = colnames(metadata()))
    
    
    ### Counts data frame
    dds$counts_dds <- as.data.frame(counts(dds$DESeq2))
  
  })
  
  ### PCA ----
  observeEvent(input$logaction,{
    if(input$log=="vst"){
      dds$log <- vst(dds$DESeq2, blind=FALSE)
    }else{
      dds$log <- rlogTransformation(dds$DESeq2,blind=FALSE)
    }
  })
  pcaa <- function(){
    pca(dds$log,input$conditionpca)
    
  }
  output$downloadPCA <- downloadHandler(
    filename = "PCA.png",
    content = function(file){
      ggsave(file, plot = pcaa(), device = "png")
    }
  )
  output$pca <- renderPlot({
    withProgress(message = "Running PCA , please wait",{
      
      validate(
        need(dds$log, "Please run DESeq2 and PCA")
      )
      pcaa()
    })})
  
  ### Depth ----
  
  normdepth <- eventReactive(input$normalize1,{
    if(input$normalize1==TRUE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2,normalized=TRUE))
    }
    else if(input$normalize1==FALSE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2))
      
    }
  })
  output$downloadDepth <- downloadHandler(
    filename = "Depth.png",
    content = function(file){
      ggsave(file, plot = depthh(), device = "png")
    }
  )
  depthh <- function(){
    depth(normdepth(),breaksize= input$breaks1)
  }
  output$depth <- renderPlot({
    validate(
      need(dds$counts_dds, "Please run DESeq2")
    )
    depthh()
  })
  
  ### Count distribution ----
  normcount <- eventReactive(input$normalize,{
    if(input$normalize==TRUE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2,normalized=TRUE))
    }
    else if(input$normalize==FALSE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2))
      
    }
  })
  distribution <- function(){count_distribution(normcount(), sample = input$sample,min=input$axis[1],max=input$axis[2],breaksize = input$breaks)
  }
  output$downloadDistribution <- downloadHandler(
    filename = function(){
      paste(input$sample,'.png',sep = '')
    },
    content = function(file){
      ggsave(file, plot = distribution(), device = "png")
    }
  )
  output$count <- renderPlot({
    validate(
      need(dds$DESeq2, "Please run DESeq2")
    )
    distribution()})
  
  ### Count by gene ---
  norm <- eventReactive(input$normalize4,{
    if(input$normalize4==TRUE){
      dds$counts_dds <- as.data.frame(counts(dds$DESeq2,normalized=TRUE))
    }
    else if(input$normalize4==FALSE){
      dds$counts_dds <- as.data.frame(counts(dds$DESeq2))
      
    }
  })
  countg <- function() {
    plotcount(norm(), input$gene)
  }
  
  output$downloadCountgene <- downloadHandler(
    filename = "CountByGene.png",
    content = function(file){
      ggsave(file, plot = countg(), device = "png")
    }
  )
  output$countgene <- renderPlot({
    validate(
      need(dds$DESeq2, "Please run DESeq2")
    )
    countg()})
  
  
  ### MA plot ----
  maplo <- function(){
    maplot(dds$results,padje=input$pvalue)
    
  }
  output$downloadMaplot <- downloadHandler(
    filename = "Maplot.png",
    content = function(file){
      ggsave(file, plot = maplo(), device = "png")
    }
  )
  output$maplot <- renderPlot({
    validate(
      need(dds$results, "Please run DESeq2")
    )
    maplo()
  })
  
  ### Volcano plot ----
  volcan <- function(){
    volcanoPlot(dds$results,annotation = input$annotation3, anno = anno() ,padje=input$pvalue2,minlogF=input$sliderfold[1], maxlogF=input$sliderfold[2], minlogP=input$sliderlog,count=colnames(count_table()))
  }
  output$downloadVulcano <- downloadHandler(
    filename = "Volcanoplot.png",
    content = function(file){
      ggsave(file, plot = volcan(), device = "png")
    }
  )
  
  output$volcano <- renderPlot({
    validate(
      need(dds$results, "Please run DESeq2")
    )
    volcan()
  })
  
  ### Dispersion ----
  dispersion1 <- function(){
    dispersion(dds$DESeq2)
  }
  
  output$downloadDispersion <- downloadHandler(
    filename = "Dispersion.png",
    content = function(file){
      png(file)
      dispersion1()
      dev.off()
    }
  )
  
  output$dispersionPlot <- renderPlot({
    validate(
      need(dds$DESeq2, "Please run DESeq2")
    )
    dispersion1()
    
  })
  
  ### Heat map 1 ----
  observeEvent(input$logaction2,{
    if(input$log1=="vst"){
      dds$log2 <- vst(dds$DESeq2, blind=FALSE)
    }else{
      dds$log2 <- rlogTransformation(dds$DESeq2,blind=FALSE)
    }
  })
  heatmapcluster <- function(){
    clustering_heatmap(dds$log2)
  }
  output$clusteringmap <- renderPlot({
    withProgress(message = "Running heatmap , please wait",{
      validate(
        need(dds$log2, "Please run DESeq2 and Heat map")
      )
      heatmapcluster()
    })})
  
  output$downloadHeatmap1 <- downloadHandler(
    filename = "DistanceMatrix.png",
    content = function(file){
      png(file)
      heatmapcluster()
      dev.off()
    }
  )
  ### Heat map 2 ----
  observeEvent(input$logaction3,{
    if(input$log3=="vst"){
      dds$log3 <- vst(dds$DESeq2, blind=FALSE)
    }else{
      dds$log3 <- rlogTransformation(dds$DESeq2,blind=FALSE)
    }
  })
  
  heatmap2 <- function() {
    heatmap(dds$results,dds$log3,annotation = input$annotation2,metadata=metadata(),condition = input$conditionheatmap,count=colnames(count_table()),min=input$slider2[1],max=input$slider2[2],anno=anno())
  }
  output$clusteringmap2 <- renderPlot({
    validate(
      need(dds$log3, "Please run DESeq2 and Heat map")
    )
    withProgress(message = "Running heatmap , please wait",{
      heatmap2()
    })})
  output$downloadHeatmap2 <- downloadHandler(
    filename = "Heatmap.png",
    content = function(file){
      png(file)
      heatmap2()
      dev.off()
    }
  )
 
}
