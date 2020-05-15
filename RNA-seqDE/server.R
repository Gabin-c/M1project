### Server  ----
server <- function(input, output,session) {
  ### Create dds object with no values
  dds <- reactiveValues()
  
  ### Increase the authorized size for upload ----
  options(shiny.maxRequestSize=30*1024^2)
  
  
  ### Display a count table example in the introduction 
  output$countExample <- renderTable({   
    countex <- read.csv("countexample.csv",sep=",")
    countex
  })
  
  ### Display a metadata table example in the introduction
  output$metadataExample <- renderTable({   
    metaex <- read.csv("metadataexample.csv",sep=",")
    metaex
  }, na = "")
  
  ### Display an annotation table example in the introduction
  output$annoExample <- renderTable({   
    annoex <- read.csv("annoexample.csv",sep=",")
    annoex
  })
  
  ### Import the count file ----
  ### csv, tsv or txt files separated by comma, tab or semi-colon
  count_table <- reactive({
    req(input$CountDataTable)
    countTable <- read.csv(input$CountDataTable$datapath, sep = input$separator_Count)
    
  })
  
  ### Display the count file ----
  output$CountReadTable <- DT::renderDataTable(count_table(), options = list(pageLength = 20, autoWidth = FALSE,scrollX = TRUE, scrollY = '300px'))
  
  
  ### Import the metadata file ---- 
  ### csv, tsv or txt files separated by comma, tab or semi-colon
  metadata <- reactive({
    req(input$MetadataFile)
    meta_table <- read.csv(input$MetadataFile$datapath, sep = input$separator_Metadata, row.names=NULL)
  })
  ### Display the metadata file ----
  output$MetaTable <- DT::renderDataTable(metadata(),options = list(pageLength = 20, autoWidth = FALSE,scrollX = TRUE, scrollY = '300px'))
  
  ### Design condition for DESeq2 ----
  ### Corresponds to the columns of the metadata table
  
  ### Verify for each column of metadata() if the values of each row are differents or not
  allValuesDifferent <- reactive({
    apply(metadata(), 2, function(a) length(unique(a))==nrow(metadata()))
  })
  
  ### Select the FALSE of "allValuesDifferent" so the columns with equal values
  notAllValuesDifferent <- reactive({
    metadata()[allValuesDifferent() == FALSE]
  })
  
  ### Verify for each column notAllValuesDifferent if the values of each row are equals
  uniqueValue <- reactive({
    apply(notAllValuesDifferent(), 2, function(a) length(unique(a))==1)
  })
  
  ### Select the false of "uniqueValue" so tje columns without unique values 
  notUniqueValue <- reactive({
    notAllValuesDifferent()[uniqueValue() == FALSE]
  })
  
  ### Update selectinput with columns of notUniqueValue() for the DESeq2 design
  observeEvent(input$MetadataFile,{
    updateSelectInput(session,"DesignDESeq2", choices = paste("~ ",paste(colnames(notUniqueValue()))))
  })
  

  
  
  ### Check if the user has annotation file to upload 
  observeEvent(input$CheckAnnotation, {
    if(input$CheckAnnotation == TRUE){
      updateTabsetPanel(session, "paramsAnno", selected = "annotation")
    }else{
      updateTabsetPanel(session, "paramsAnno", selected = "nothing")
    }
  })
  
  ### Import annotation file ----
  ### csv, tsv or txt files separated by comma, tab or semi-colon
  anno <- reactive({
    req(input$AnnotationFile)
    read.csv(input$AnnotationFile$datapath, sep = input$sep_Anno)
  })
  output$AnnoTable <- DT::renderDataTable(anno(),options = list(pageLength = 20, autoWidth = FALSE,scrollX = TRUE, scrollY = '300px'))
  
  ### Display DESeq2 page after count table uploaded and a message to upload metadata file and choose design
  observeEvent(input$CountDataTable,{
    output$menuDESeq2 <- renderMenu({
      if(input$DesignDESeq2 == ""){
        showNotification("Upload a metadata file and choose a design.")
      }
      menuItem(text = "2 Run DESeq2", tabName = "deseq2", icon = icon("play-circle"))
    })
  })

  ### Running DESeq2 clicking on the button  ---- 
  observeEvent(input$RunDESeq2,{
    if(input$DesignDESeq2 == ""){
      showNotification("Upload a metadata file and choose a design.")
    }
    else{
    req(input$RunDESeq2)
    ### Waiting screen 
    waiter <- waiter::Waiter$new(html = spin_ball())
    waiter$show()
    
    ### DESeq2 process 
    dds$dds <- DESeqDataSetFromMatrix(count_table(),colData=metadata(),design=as.formula(input$DesignDESeq2), tidy=TRUE)
    dds$DESeq2 <- DESeq(dds$dds)
    dds$results <- results(dds$DESeq2,tidy=TRUE)
    
    ### Display success message after running DESeq2
    output$SuccessMessage <- renderUI({
      box(width = 12, solidHeader = F,
          HTML("<center><h3>DESeq2 workflow successfully completed</h3></center>"))
    })
    
    
    ### Samples choices for count distribution
    updateSelectInput(session,"sample",choices = metadata()[,1])
    
    ### Genes choices for count by gene
    updateSelectizeInput(session,"gene",choices = count_table()[,1], server = TRUE)
    
    ### Factor choices for PCA
    updateSelectInput(session,"conditionpca",choices = colnames(metadata()))
    
    ### Factor choices heatmap
    updateSelectInput(session,"conditionHeatmap",choices = colnames(metadata()))
    
    
    ### Counts data frame normalized or not
    dds$counts_dds <-as.data.frame(counts(dds$DESeq2))
    dds$counts_dds_n <-as.data.frame(counts(dds$DESeq2,normalized=TRUE))
    ### Transposed counts data frame normalized or not
    dds$counts_turnup <- as.data.frame(t(dds$counts_dds))
    dds$counts_turnup_n <- as.data.frame(t(dds$counts_dds_n))
    
    
    ### Display "Results" menu when DESeq2 is successfully run
    ### Put a check icon for the menu where the plots are already display
    ### The others (PCA and both heatmap) needs to be run
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
    
    ### End of the waiting screen
    on.exit(waiter$hide())
  }})
  


  
  ### Count distribution ----
  ### Normalize or not the data
  normcount <- reactive({
    if(input$normalizeDistribution==TRUE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2,normalized=TRUE))
    }
    else if(input$normalizeDistribution==FALSE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2))
    }
  })

  ### Display count distribution plot using count.distribution.plot() function (see function_dds.R)
  distribution <- function(){count.distribution.plot(normcount(), sample = input$sample,x.min=input$axis[1],x.max=input$axis[2],break.width = input$breaksDistribution)}

  output$CountDistributionPlot <- renderPlot({
    validate(
      need(dds$DESeq2, "Please run DESeq2")
    )
    distribution()})
  ### Download distribution plot in .png format
  output$downloadDistribution <- downloadHandler(
    filename = function(){
      paste(input$sample,'.png',sep = '')
    },
    content = function(file){
      ggsave(file, plot = distribution(), device = "png")
    }
  )
  
  
  ### Count by gene ----
  ### Normalize or not the data
  normCountGene <- eventReactive(input$normalizeCountGene,{
    if(input$normalizeCountGene==TRUE){
      dds$counts_turnup_n
    }
    else if(input$normalizeCountGene==FALSE){
      dds$counts_turnup 
    }
  })
  ### Display Count by gene plot using gene.count.plot() function (see function_dds.R)
  countg <- function() {
    gene.count.plot(normCountGene(), input$gene)
  }
  output$CountGenePlot <- renderPlot({
    validate(
      need(dds$DESeq2, "Please run DESeq2")
    )
    countg()})
  ### Download Counts by gene plot in .png format
  output$downloadCountgene <- downloadHandler(
    filename = "CountByGene.png",
    content = function(file){
      ggsave(file, plot = countg(), device = "png")
    }
  )
  
  ### Depth ----
  ### Normalize or not the data
  normdepth <- eventReactive(input$normalizeDepth,{
    if(input$normalizeDepth==TRUE){
      dds$counts_dds_n <-as.data.frame(counts(dds$DESeq2,normalized=TRUE))
    }
    else if(input$normalizeDepth==FALSE){
      dds$counts_dds <-as.data.frame(counts(dds$DESeq2))
    }
  })
  ### Display depth plot using depth.plot() function from function_dds.R
  depthFunction <- function(){
    depth.plot(normdepth(),break.width= input$breaksDepth)
  }
  output$depth <- renderPlot({
    validate(
      need(dds$counts_dds, "Please run DESeq2")
    )
    depthFunction()
  })
  ### Download depth file
  output$downloadDepth <- downloadHandler(
    filename = "Depth.png",
    content = function(file){
      ggsave(file, plot = depthFunction(), device = "png")
    }
  )
  
  
  ### Dispersion ----
  ### Display dispersion plot using dispersion() function from function_dds.R
  dispersionFunction <- function(){
    dispersion(dds$DESeq2)
  }
  output$dispersionPlot <- renderPlot({
    validate(
      need(dds$DESeq2, "Please run DESeq2")
    )
    dispersionFunction()
  })
  
  ### Download dispersion plot
  output$downloadDispersion <- downloadHandler(
    filename = "Dispersion.png",
    content = function(file){
      png(file)
      dispersionFunction()
      dev.off()
    }
  )
  
  ### PCA ----
  ### Needs to run PCA 
  ### 2 possibles transformations : vst or rlog
  observeEvent(input$runPCA,{
    if(input$TransformationPCA=="vst"){
      dds$TransformationPCA <- vst(dds$DESeq2, blind=FALSE)
    }else{
      dds$TransformationPCA <- rlogTransformation(dds$DESeq2,blind=FALSE)
    }
  })
  
  ### Display PCA plot usin pca.plot() function from function_dds.R
  PCAfunction <- function(){
    pca.plot(dds$TransformationPCA,input$conditionpca)
    
  }
  
  output$PCAplot <- renderPlot({
    withProgress(message = "Running PCA , please wait",{
      
      validate(
        need(dds$TransformationPCA, "Please run DESeq2 and PCA")
      )
      PCAfunction()
    })})
  ### Possibility to download the PCA plot
  output$downloadPCA <- downloadHandler(
    filename = "PCA.png",
    content = function(file){
      ggsave(file, plot = PCAfunction(), device = "png")
    }
  )
  
  

  ### MA plot ----
  ### Display MA plot using ma.plot() function from function_dds.R
  MAplotFunction <- function(){
    ma.plot(dds$results,p.val=input$pvalueMAplot)
  }
  output$MAplot <- renderPlot({
    validate(
      need(dds$results, "Please run DESeq2")
    )
    MAplotFunction()
  })
  ### Display a table with the number of differential expressed genes
  output$numberDEgenes <- renderTable({
    number.DE.gene(dds$results,input$pvalueMAplot)
  })
  ### Download MAplot in .png format
  output$downloadMaplot <- downloadHandler(
    filename = "Maplot.png",
    content = function(file){
      ggsave(file, plot = MAplotFunction(), device = "png")
    }
  )
  
  ### Volcano plot ----
  ### Display parameters for volcano 
  ### If there is an annotation file, display sliders to choose parameters 
  observeEvent(input$annotationVolcano,{
    if(input$CheckAnnotation== TRUE){
      if(input$annotationVolcano==TRUE){
      output$SliderFoldVolcano <- renderUI({ 
        sliderInput("sliderfold", "Choose your fold", min=-20, max=20, value=c(-6,6))
      })
      output$SliderLogVolcano <- renderUI({ 
        sliderInput("sliderlog", "Choose your log10", min=0, max=300, value=30)
      })
    }}
    if(input$CheckAnnotation== FALSE){
      if(input$annotationVolcano==TRUE){
      output$SliderFoldVolcano <- renderUI({})
      output$SliderLogVolcano <- renderUI({})
      output$AnnoVolcano <- renderUI({
        HTML("<center><h3> You don't have annotation file. </h3></center>")
      })
      }}
    if(input$CheckAnnotation== FALSE){
      if(input$annotationVolcano==FALSE){
        output$SliderFoldVolcano <- renderUI({})
        output$SliderLogVolcano <- renderUI({})
        output$AnnoVolcano <- renderUI({})
      }}
  })
  
  ### Display volcano plot using volcano.plot() function from function_dds.R
  VolcanoplotFunction <- function(){
    volcano.plot(dds$results,is.anno = input$annotationVolcano, anno = anno() ,p.val=input$pvalueVolcano,minlogF=input$sliderfold[1], maxlogF=input$sliderfold[2], minlogP=input$sliderlog,count.tb=colnames(count_table()))
  }
  output$volcanoPlot <- renderPlot({
    validate(
      need(dds$results, "Please run DESeq2")
    )
    VolcanoplotFunction()
  })
  ### Download Volcano plot in .png format
  output$downloadVolcano <- downloadHandler(
    filename = "Volcanoplot.png",
    content = function(file){
      ggsave(file, plot = VolcanoplotFunction(), device = "png")
    }
  )
  

  ### Sample distance matrix heatmap ----
  ### Chose vst or rlog transformation
  observeEvent(input$RunMatrix,{
    if(input$TransformationMatrix=="vst"){
      dds$TransformationMatrix <- vst(dds$DESeq2, blind=FALSE)
    }else{
      dds$TransformationMatrix <- rlogTransformation(dds$DESeq2,blind=FALSE)
    }
  })
  ### Display distance matrix using the fonction distance.matrix.heatmap() from function_dds.R
  distanceCluster <- function(){
    sample.distance.matrix.heatmap(dds$TransformationMatrix)
  }
  output$DistanceMatrixMap <- renderPlot({
    withProgress(message = "Running heatmap , please wait",{
      validate(
        need(dds$TransformationMatrix, "Please run DESeq2 and Heat map")
      )
      distanceCluster()
    })})
  ### Download distance matrix in .png format
  output$downloadDistanceMatrix <- downloadHandler(
    filename = "DistanceMatrix.png",
    content = function(file){
      png(file)
      distanceCluster()
      dev.off()
    }
  )
  
  ### Gene expression heatmap ----
  ### Chose vst or rlog transformation
  observeEvent(input$RunHeatmap,{
    if(input$TransformationHeatmap=="vst"){
      dds$TransformationHeatmap <- vst(dds$DESeq2, blind=FALSE)
    }else{
      dds$TransformationHeatmap <- rlogTransformation(dds$DESeq2,blind=FALSE)
    }
  })
  
  observeEvent(input$annotationHeatmap,{
    if(input$annotationHeatmap==TRUE){
      if(input$CheckAnnotation==FALSE){
        output$Annoheatmap <- renderUI(
          box(width = 12, solidHeader = F,
            HTML("<center><h3> You don't have annotation file. </h3></center>")))
      }}
    else{
      output$Annoheatmap <- renderUI({})
    }
  })
  
  ### Display heatmap using gene.expression.heatmap() function from function_dds.R
  heatmapCluster <- function() {
    input$RunHeatmap
    gene.expression.heatmap(dds$results,dds$TransformationHeatmap,is.anno = input$annotationHeatmap,metadata=metadata(),condition = input$conditionHeatmap,count.tb=colnames(count_table()),min=input$nbGenes[1],max=input$nbGenes[2],anno=anno())
  }
  output$Heatmap <- renderPlot({
    validate(
      need(dds$TransformationHeatmap, "Please run DESeq2 and Heat map")
    )
    heatmapCluster()
  })
  ### Download heatmap in .png format
  output$downloadHeatmap <- downloadHandler(
    filename = "Heatmap.png",
    content = function(file){
      png(file)
      heatmapCluster()
      dev.off()
    }
  )
  
  ### Theme ----
  ### Choice between dark or light theme with a switcher button
  observeEvent(input$theme,{
    if(input$theme==TRUE){
      output$themes <- renderUI({
        theme_grey_dark2
      })
    }else{
      output$themes <-renderUI({
        theme_grey_light2
      })
    }
  })
  
  

  
  ### "Input count table" menu in the sidebar
  ### Change the icon with a check icon when the counts file is imported
  menuCount <- reactive({
    if(is.null(input$CountDataTable)==TRUE){
      menuSubItem(text = "1.1 Input count table", tabName = "CountData")
      
    }else{
      menuSubItem(text = "1.1 Input count table", tabName = "CountData", icon = icon("far fa-check-square"))
    }
  })
  output$CountTable <- renderMenu({
    menuCount()
  })
  
  ### "Input metadata table" menu in the sidebar
  ### Change the icon with a check icon when the metadata file is imported
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
  
  ### "Input annotation file" menu in the sidebar
  ### Change the icon with a check icon when the annotation file is imported
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
  

  
  ### Menu for PCA menu in sidebar
  ### Change the icon with check icon when pca is run successfully
  menuPCA <- reactive({
    if(input$runPCA){
      menuSubItem("PCA",tabName = "pca",icon = icon("far fa-check-square"))
    }else{
      menuSubItem("PCA",tabName = "pca")
    }
  })
  
  ### Menu for Distance matrix in the sidebar
  ### Change the icon with check icon when heatmap is run successfully
  menuDistanceMatrix <- reactive({
    if(input$RunMatrix){
      menuSubItem("Sample distance matrix",tabName = "DistanceMatrix",icon = icon("far fa-check-square"))
    }else{
      menuSubItem("Sample distance matrix",tabName = "DistanceMatrix")
    }
  })
  
  ### Menu for heatmap in the sidebar
  ### Change the icon with check icon when heatmap is run successfully
  menuHeatmap <- reactive({
    if(input$RunHeatmap){
      menuSubItem("Gene expression Heatmap",tabName = "Heatmap", icon = icon("far fa-check-square"))
    }else{
      menuSubItem("Genne expression Heatmap",tabName = "Heatmap")
    }
  })
  
}
