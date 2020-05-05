### Server  ----
server <- function(input, output,session) {
  dds <- reactiveValues()
  ### Increase the authorized size for upload ----
  options(shiny.maxRequestSize=30*1024^2)
  
  
  ### Intro : 
  output$countExample <- renderTable({   
    countex <- read.csv("countexample.csv",sep=",")
    countex
  })
  
  output$metadataExample <- renderTable({   
    metaex <- read.csv("metadataexample.csv",sep=",")
    metaex
  }, na = "")
  
  output$annoExample <- renderTable({   
    annoex <- read.csv("annoexample.csv",sep=",")
    annoex
  })
  
  ### Import the count ----
  count_table <- reactive({
    req(input$CountDataTable)
    countTable <- read.csv(input$CountDataTable$datapath, sep = input$separator_Count)
    
  })
  
  
  
  ### Display the count file ----
  output$CountReadTable <- DT::renderDataTable(count_table(), options = list(pageLength = 20, autoWidth = FALSE,scrollX = TRUE, scrollY = '300px'))
  
  
  
  ### Import the metadata file ---- 
  metadata <- reactive({
    req(input$MetadataFile)
    meta_table <- read.csv(input$MetadataFile$datapath, sep = input$separator_Metadata, row.names=NULL)
  })
  ### Display the metadata file ----
  output$MetaTable <- DT::renderDataTable(metadata(),options = list(pageLength = 20, autoWidth = FALSE,scrollX = TRUE, scrollY = '300px'))
  
  ### Design condition for DESeq2 ----
  observeEvent(input$MetadataFile,{
    updateTextInput(session,"DesignDESeq2", value = paste("~ ",paste(colnames(metadata()), collapse=" + ")))
  })
  
  ### Import annotation file ----
  observeEvent(input$CheckAnnotation, {
    if(input$CheckAnnotation == TRUE){
      updateTabsetPanel(session, "paramsAnno", selected = "annotation")
    }else{
      updateTabsetPanel(session, "paramsAnno", selected = "nothing")
    }
  })
  anno <- reactive({
    req(input$AnnotationFile)
    read.csv(input$AnnotationFile$datapath, sep = input$sep_Anno)
  })
  output$AnnoTable <- DT::renderDataTable(anno(),options = list(pageLength = 20, autoWidth = FALSE,scrollX = TRUE, scrollY = '300px'))
  
  ### Display parameters for volcano ---- 
  observeEvent(input$annotationVolcano,{
    if(input$annotationVolcano== TRUE){
      output$SliderFoldVolcano <- renderUI({ 
        sliderInput("sliderfold", "Choose your fold", min=-20, max=20, value=c(-6,6))
      })
      output$SliderLogVolcano <- renderUI({ 
        sliderInput("sliderlog", "Choose your log10", min=0, max=300, value=30)
      })
    }else{
      output$SliderFoldVolcano <- renderUI({})
      output$SliderLogVolcano <- renderUI({})
    }
    
  })
  
  
  
  
  ### DDS object  ---- 
  observeEvent(input$RunDESeq2,{
    req(input$RunDESeq2)
    waiter <- waiter::Waiter$new(html = spin_ball())
    waiter$show()
    
    
    ### Running DESeq1
    dds$dds <- DESeqDataSetFromMatrix(count_table(),colData=metadata(),design=as.formula(input$DesignDESeq2), tidy=TRUE)
    dds$DESeq2 <- DESeq(dds$dds)
    dds$results <- results(dds$DESeq2,tidy=TRUE)
    
    ### Display count table after dds
    output$SuccessMessage <- renderUI({
      box(width = 12, solidHeader = F,
          HTML("<center><h3>DESeq2 workflow successfully completed</h3></center>"))
    })
    
    ### Clustering map
    output$clustering_map <- renderPlot(clustering_heatmap(dds$DESeq2,log="vst"))
    
    ### Choices for count distribution
    updateSelectInput(session,"sample",choices = metadata()[,1])
    
    ### Choices for count by gene
    
    updateSelectizeInput(session,"gene",choices = count_table()[,1], server = TRUE)
    
    
    ### Choices for PCA
    updateSelectInput(session,"conditionpca",choices = colnames(metadata()))
    
    ### Choices heatmap
    updateSelectInput(session,"conditionHeatmap",choices = colnames(metadata()))
    
    
    ### Counts data frame
    dds$counts_dds <-as.data.frame(counts(dds$DESeq2))
    dds$counts_dds_n <-as.data.frame(counts(dds$DESeq2,normalized=TRUE))
    dds$counts_turnup <- as.data.frame(t(dds$counts_dds))
    dds$counts_turnup_n <- as.data.frame(t(dds$counts_dds_n))
    
    
    on.exit(waiter$hide())
    
  })
  
  ### PCA ----
  observeEvent(input$runPCA,{
    if(input$TransformationPCA=="vst"){
      dds$TransformationPCA <- vst(dds$DESeq2, blind=FALSE)
    }else{
      dds$TransformationPCA <- rlogTransformation(dds$DESeq2,blind=FALSE)
    }
  })
  PCAfunction <- function(){
    pca(dds$TransformationPCA,input$conditionpca)
    
  }
  output$downloadPCA <- downloadHandler(
    filename = "PCA.png",
    content = function(file){
      ggsave(file, plot = PCAfunction(), device = "png")
    }
  )
  output$PCAplot <- renderPlot({
    withProgress(message = "Running PCA , please wait",{
      
      validate(
        need(dds$TransformationPCA, "Please run DESeq2 and PCA")
      )
      PCAfunction()
    })})
  
  ### Depth ----
  
  normdepth <- eventReactive(input$normalizeDepth,{
    if(input$normalizeDepth==TRUE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2,normalized=TRUE))
    }
    else if(input$normalizeDepth==FALSE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2))
      
    }
  })
  output$downloadDepth <- downloadHandler(
    filename = "Depth.png",
    content = function(file){
      ggsave(file, plot = depthFunction(), device = "png")
    }
  )
  depthFunction <- function(){
    depth.plot(normdepth(),break.width= input$breaksDepth)
  }
  output$depth <- renderPlot({
    validate(
      need(dds$counts_dds, "Please run DESeq2")
    )
    depthFunction()
  })
  
  ### Count distribution ----
  normcount <- reactive({
    if(input$normalizeDistribution==TRUE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2,normalized=TRUE))
    }
    else if(input$normalizeDistribution==FALSE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2))
      
    }
  })
  distribution <- function(){count.distribution.plot(normcount(), sample = input$sample,x.min=input$axis[1],x.max=input$axis[2],break.width = input$breaksDistribution)
  }
  output$downloadDistribution <- downloadHandler(
    filename = function(){
      paste(input$sample,'.png',sep = '')
    },
    content = function(file){
      ggsave(file, plot = distribution(), device = "png")
    }
  )
  output$CountDistributionPlot <- renderPlot({
    validate(
      need(dds$DESeq2, "Please run DESeq2")
    )
    distribution()})
  
  
  ### Count by gene ----
  
  normCountGene <- eventReactive(input$normalizeCountGene,{
    if(input$normalizeCountGene==TRUE){
      dds$counts_turnup_n
    }
    else if(input$normalizeCountGene==FALSE){
      dds$counts_turnup 
      
    }
  })
  countg <- function() {
    plotcount(normCountGene(), input$gene)
  }
  
  output$downloadCountgene <- downloadHandler(
    filename = "CountByGene.png",
    content = function(file){
      ggsave(file, plot = countg(), device = "png")
    }
  )
  output$CountGenePlot <- renderPlot({
    validate(
      need(dds$DESeq2, "Please run DESeq2")
    )
    countg()})
  
  ### MA plot ----
  MAplotFunction <- function(){
    maplot(dds$results,padje=input$pvalueMAplot)
    
  }
  output$downloadMaplot <- downloadHandler(
    filename = "Maplot.png",
    content = function(file){
      ggsave(file, plot = MAplotFunction(), device = "png")
    }
  )
  output$MAplot <- renderPlot({
    validate(
      need(dds$results, "Please run DESeq2")
    )
    MAplotFunction()
  })
  output$numberDEgenes <- renderTable({
    number_of_DE(dds$results,input$pvalueMAplot)
  })
  
  ### Volcano plot ----
  VolcanoplotFunction <- function(){
    volcanoPlot(dds$results,annotation = input$annotationVolcano, anno = anno() ,padje=input$pvalueVolcano,minlogF=input$sliderfold[1], maxlogF=input$sliderfold[2], minlogP=input$sliderlog,count=colnames(count_table()))
  }
  output$downloadVolcano <- downloadHandler(
    filename = "Volcanoplot.png",
    content = function(file){
      ggsave(file, plot = VolcanoplotFunction(), device = "png")
    }
  )
  
  output$volcanoPlot <- renderPlot({
    validate(
      need(dds$results, "Please run DESeq2")
    )
    VolcanoplotFunction()
  })
  
  ### Dispersion ----
  dispersionFunction <- function(){
    dispersion(dds$DESeq2)
  }
  
  output$downloadDispersion <- downloadHandler(
    filename = "Dispersion.png",
    content = function(file){
      png(file)
      dispersionFunction()
      dev.off()
    }
  )
  
  output$dispersionPlot <- renderPlot({
    validate(
      need(dds$DESeq2, "Please run DESeq2")
    )
    dispersionFunction()
    
  })
  
  ### Heat map 1 ----
  
  observeEvent(input$RunMatrix,{
    if(input$TransformationMatrix=="vst"){
      dds$TransformationMatrix <- vst(dds$DESeq2, blind=FALSE)
    }else{
      dds$TransformationMatrix <- rlogTransformation(dds$DESeq2,blind=FALSE)
    }
  })
  distanceCluster <- function(){
    
    clustering_heatmap(dds$TransformationMatrix)
  }
  output$DistanceMatrixMap <- renderPlot({
    withProgress(message = "Running heatmap , please wait",{
      validate(
        need(dds$TransformationMatrix, "Please run DESeq2 and Heat map")
      )
      distanceCluster()
    })})
  
  output$downloadDistanceMatrix <- downloadHandler(
    filename = "DistanceMatrix.png",
    content = function(file){
      png(file)
      distanceCluster()
      dev.off()
    }
  )
  ### Heat map 2 ----
  observeEvent(input$RunHeatmap,{
    if(input$TransformationHeatmap=="vst"){
      dds$TransformationHeatmap <- vst(dds$DESeq2, blind=FALSE)
    }else{
      dds$TransformationHeatmap <- rlogTransformation(dds$DESeq2,blind=FALSE)
    }
  })
  
  heatmapCluster <- function() {
    input$RunHeatmap
    
    
    heatmap(dds$results,dds$TransformationHeatmap,annotation = input$annotationHeatmap,metadata=metadata(),condition = input$conditionHeatmap,count=colnames(count_table()),min=input$nbGenes[1],max=input$nbGenes[2],anno=anno())
  }
  output$Heatmap <- renderPlot({
    validate(
      need(dds$TransformationHeatmap, "Please run DESeq2 and Heat map")
    )
    
    heatmapCluster()
  })
  output$downloadHeatmap <- downloadHandler(
    filename = "Heatmap.png",
    content = function(file){
      png(file)
      heatmapCluster()
      dev.off()
    }
  )
  
  menuHeatmap <- reactive({
    if(input$RunHeatmap){
      menuSubItem("Heatmap",tabName = "Heatmap", icon = icon("far fa-check-square"))
    }else{
      menuSubItem("Heatmap",tabName = "Heatmap")
      
      
    }
  })
  
  menuDistanceMatrix <- reactive({
    if(input$RunMatrix){
      menuSubItem("Distance matrix",tabName = "DistanceMatrix",icon = icon("far fa-check-square"))
    }else{
      menuSubItem("Distance matrix",tabName = "DistanceMatrix")
      
      
    }
    
  })
  menuPCA <- reactive({
    if(input$runPCA){
      menuSubItem("PCA",tabName = "pca",icon = icon("far fa-check-square"))
    }else{
      menuSubItem("PCA",tabName = "pca")
      
      
    }
    
  })
  
  
  ### Theme ----
  
  observeEvent(input$theme,{
    if(input$theme==TRUE){
      output$themes <- renderUI({
        shinyDashboardThemes("grey_dark")
      })
      
    }else{
      output$themes <-renderUI({
        shinyDashboardThemes("grey_light")
        
      })
    }
  })
  
  
  menu <- reactive({
    if(is.null(input$CountDataTable)==TRUE){
      menuSubItem(text = "1.1 Input count table", tabName = "CountData")
      
    }else{
      menuSubItem(text = "1.1 Input count table", tabName = "CountData", icon = icon("far fa-check-square"))
    }
  })
  output$CountTable <- renderMenu({
    menu()
  })
  
  observeEvent(input$RunDESeq2,{
    output$menuResults <- renderMenu({  menuItem(text = "3 Results", tabName = "deseq2", icon = icon("poll"),startExpanded = TRUE,
                                                 menuSubItem("Count distribution",tabName = "Count_Distribution",icon = icon("far fa-check-square")),
                                                 menuSubItem("Count by gene", tabName = "Count_Gene",icon = icon("far fa-check-square")),
                                                 menuSubItem("Depth of sample",tabName = "Depth",icon = icon("far fa-check-square")),
                                                 menuSubItem("Dispersion",tabName = "Dispersion",icon = icon("far fa-check-square")),
                                                 menuPCA(),
                                                 menuSubItem("MA Plot",tabName = "MAplot",icon = icon("far fa-check-square")),
                                                 menuSubItem("Volcano Plot",tabName = "Volcanoplot",icon = icon("far fa-check-square")),
                                                 menuDistanceMatrix(),
                                                 menuHeatmap()
    )  
    })
  })
  
  menuMetadata <- reactive({
    if(is.null(input$MetadataFile)==TRUE){
      menuSubItem(text = "1.2 Input metadata table", tabName = "Metadata")
      
    }else{
      menuSubItem(text = "1.2 Input metadata table", tabName = "Metadata", icon = icon("far fa-check-square"))
    }
  })
  output$MetadataTable <- renderMenu({
    menuMetadata()
  })
  
  menuAnnotation <- reactive({
    if(is.null(input$AnnotationFile)==TRUE){
      menuSubItem(text = "1.3 Input annotation file", tabName = "Annotation")
      
    }else{
      menuSubItem(text = "1.3 Input annotation file", tabName = "Annotation", icon = icon("far fa-check-square"))
    }
  })
  output$AnnotationTable <- renderMenu({
    menuAnnotation()
  })
  
}
