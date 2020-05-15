### Library ----
### here we find all library need for proper functionning of app
library(shiny)
library(shinydashboard)
library(shinyjs)
library(dplyr)
library(tidyverse)
library(vroom)
library(DT)
library(DESeq2)
library(shinyWidgets)
library(shinythemes)
library(waiter)
library(dashboardthemes)
library(shinycssloaders)
library(shinydashboardPlus)
### Files with all the function needed to make plots ----
source("function_dds.R")






### Annotation pannel ----
parameters_Annotation <- tagList(
  tags$style("#paramsAnno { display:none; }"),
### parameter_tabs is a tabsetPanel for input annotation page
### it allow to switch between panel when we use updateTabsetPanel()
### in our case it is when we check right annotation boxes of input annotation page
  tabsetPanel(
    id="paramsAnno",
    tabPanel("nothing"),  
    tabPanel("annotation",
             fluidRow(
               column(width= 6,
                      box(title="Upload annotation file",width = 12, solidHeader = TRUE,collapsible = TRUE,
                          column(width=5,selectInput("sep_Anno", "Separator:", c("Comma" = ",", "Tab" = "\t", "Semi-colon" = ";"))),
                          fileInput("AnnotationFile", "Upload annotation file", accept = c(".csv",".txt",".tsv")))),
               column(width= 6,     
                      box(
                        title = "Accepted files :", width = 12,
                        HTML(
                          "<li> .csv / .tsv / .txt files </li>
                          <li> Separated by tabulation, comma or semi-colon </li>
                          <li> One column with genes symbols named 'symbol'</li>
                          "),
                        height = 160
                        )
                      )
                      )
               )
               )
             )





### User Interface  ----
### UI is a shinydashboard, we use library shinydashboard to set up a dashboard UI
### a dashboard is compound of :
###   - a Header
###   - a sidebar
###   - a body 
ui <- 
  
  tagList(
  
  
  ### Parameters of the dashboard ----
  
  div(
    id = "app",
    dashboardPage(
      
      ### Customize the header ----
      ### header is composate of :
      ###   - a title 
      ###   - a home bouton which on click return user on introduction page
      dashboardHeader(title = "RNA-seq DE analysis", 
                      uiOutput("themes"),

                      ### Home button ----
                      tags$li(a(onclick = "openTab('Intro')",
                                href = NULL,
                                icon("home"),
                                title = "Homepage",
                                style = "cursor: pointer;"),
                              class = "dropdown",
                              tags$script(HTML("
                                               var openTab = function(tabName){
                                               $('a', $('.sidebar')).each(function() {
                                               if(this.getAttribute('data-value') == tabName) {
                                               this.click()
                                               };
                                               });
                                               }")))
      ),
      
      ### Sidebar ----
      ### Sider bar  is composate of a sidebar menu
      ### in this sidebar menu we have menu item which is associate with a tab of body dashboard 
      dashboardSidebar(
        sidebarMenu(id="mysidebar",
                    menuItem(text = "Informations", tabName = "Intro", icon = icon("info-circle")),
                    menuItem(text = "1 Upload data", tabName = "upload", icon = icon("arrow-circle-up"),startExpanded = TRUE,
                             # menuItem which are set up in the server function
                             menuItemOutput("CountTable"),
                             menuItemOutput("MetadataTable"),
                             menuItemOutput("AnnotationTable")),

                    menuItem(text = "2 Run DESeq2", tabName = "deseq2", icon = icon("play-circle")),
                    
                    menuItemOutput("menuResults"),

                    tags$hr(),# a simple line 
                    # this menu generate a switch button to set theme of dashboard
                    # two options are available a light or dark mode
                    menuItem(icon = NULL,
                             materialSwitch(inputId = "theme", label = "Theme", status = "default", value= TRUE)
                    ),tags$hr()
        )
      ),
      ### Dashboard body ----
      ### Organization of the differents pages associate with their menuItem
      
      
      dashboardBody(
        useShinyjs(),
        fluidRow(
          tabItems(
            ### Introduction ----
            ### introduction page associate with menuItem "Informations"
            ### In this page we find all information about application
       
                 ### in particular how it works and different tool used to generate DE analysis
            ### to generate layout html tag provide by shiny library are used 
            ### withSpinner() of shinycssloaders library is use to generate waiting screen during load of img
            tabItem(tabName = "Intro",
                    fluidPage(
                      h2("Introduction"),
                      p("This is an R Shiny web interactive application developed as part of a ", 
                        strong("course project."), "The purpose of this application is to perform a ",
                        strong("differential expression analysis from a counts table"), "in order to help researchers getting interpretable results.",
                        align = "justify"),
                      p("This application uses the package ", 
                        a("DESeq2", href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html"), 
                        "from Bioconductor. It is a package to study differential gene expression analysis based on the negative binomial distribution. It allows a quick visualization of results in order to analyze the counts table data. The results will be given in different forms like graphics, heatmaps, MAplot or even Volcano plot.",
                        align = "justify"),
                      tags$hr(),
                      h3("1. Upload data", style="padding-left: 1em"),
                      p("The input data files accepted for this App are 3 files in '.txt', '.csv' or '.tsv' format separated by comma, tabulation or semi-colon.
                        This App necessarily requires a 'Count Data Table' and a 'Metadata Table'. An optional 'Annotation File' can be added", style="padding-left: 2em", align = "justify"),
                      h4("1.1 Count Data Table", style="padding-left: 3em"),
                      p("The Count Data Table must contain the count for each sample of the experiment for each gene and the first column must be gene ID or gene name as below :",style="padding-left: 5em", align = "justify"),
                      column( 12, style="padding-left: 5em" ,withSpinner(tableOutput("countExample"))),
                      br(),
                      h4("1.2 Metadata Table", style="padding-left: 3em"),
                      p("The Metadata table must contain the information of the experiment with at least 2 columns. The first one corresponds to the samples in the same order as the columns of the Count Table. 
                        The second one is a condition column. You can add as many columns as you have factors in your experiment.",style="padding-left: 5em", align = "justify"),
                      column( 12, style="padding-left: 5em" ,withSpinner(tableOutput("metadataExample"))),
                      h4("1.2  Annotation File", style="padding-left: 3em"),
                      p("The Annotation File contains informations about the genes. If you have one, it must contain a column named 'symbol' in which we can find the symbol of each gene.",style="padding-left: 5em", align = "justify"),
                      column( 12, style="padding-left: 5em" ,withSpinner(tableOutput("annoExample"))),
                      h3("2. Results", style="padding-left: 1em"),
                      p("The results will be display after running DESeq2. You will obtain 9 differents results :", style="padding-left: 2em", align = "justify"),
                      p("- Count distribution",
                        br(), "- Count by gene",
                        br(), "- Depth of sample",
                        br(), "- Dispersion",
                        br(), "- PCA",
                        br(), "- MA plot",
                        br(), "- Volcano plot",
                        br(), "- Sample distance matrix",
                        br(), "- Gene expression Heatmap",style="padding-left: 5em", align = "justify"),
                      p("You can download all the results plots at the bottom of all these pages.",  style="padding-left: 2em", align = "justify")
                      
                      )
                    ),
            ### Upload count table ----


            ### upload page of count table of RNA-seq experience
            ### this page is associate with menuItemOuput("")
            ### On this page we find :
            ###   - a box which countain :
            ###       - a selectInput of separator 
            ###       - a fileInput to upload count table of RNA-seq experience
            ###   - a box which countain :
            ###       - information about file accepted in fileInput 
            ###   - a dataTableOutput of count table input in fileInput
            tabItem(tabName = "CountData",
                    column(width = 6,
                           box(title="Upload count table",width = 12, solidHeader = TRUE,collapsible = TRUE,
                               column(width=5,
                                      selectInput("separator_Count", "Separator:", c("Comma" = ",", "Tab" = "\t", "Semi-colon" = ";"))),
                               fileInput("CountDataTable", "Upload count table", accept = c(".csv",".txt",".tsv"))
                           )),
                    column(width = 6,
                           box(
                             title = "Accepted files :", width = 12,
                             HTML(
                               "<li> .csv / .tsv / .txt files </li>
                               <li> Separated by tabulation, comma or semi-colon </li>
                               <li> First column has to be gene ID or gene name</li>
                               <li> All others columns are count for each sample</li>"),
                             height = 160
                             )),
                    dataTableOutput("CountReadTable")
                    ),
            ### Upload metadata table ----
            ### upload page of metadata file of RNA-seq experience
            ### this page is associate with menuItemOuput("")
            ### On this page we find :
            ###   - a box which countain :
            ###       - a selectInput of separator 
            ###       - a fileInput to upload metadata of RNA-seq experience
            ###   - a box which countain :
            ###       - information about file accepted in fileInput
            ###   - a box which countain :
            ###       - a select input to chose design formula to set for DESeq2 dataset object
            ###   - a dataTableOutput of metadata input in fileInput
            
                  tabItem(tabName = "Metadata",
                    column(width = 6,
                           box(title="Upload metadata table",width = 12, solidHeader = TRUE,collapsible = TRUE,
                               column(width=5,
                                      selectInput("separator_Metadata", "Separator:", c("Comma" = ",", "Tab" = "\t", "Semi-colon" = ";"))),
                               fileInput("MetadataFile", "Upload metadata table", accept = c(".csv",".txt",".tsv"))
                           )),
                    column(width = 6,
                           box(
                             title = "Accepted files :", width = 12,
                             HTML(
                               "<li> .csv / .tsv / .txt files </li>
                               <li> Separated by tabulation, comma or semi-colon </li>
                               <li> At least metadata table contains two columns</li>
                               <li> At least one column has to be factor</li>"),
                             height = 160
                             )),
                    column(width = 12,
                           box(width = 12,
                               selectInput("DesignDESeq2","Choose your design without linear combination", c("")))),
                    dataTableOutput("MetaTable")
                    ),

            ### Upload annotation file ----
            ### upload page of annotation file associate with the RNA-seq experience
            ### this page is associate with menuItemOuput("")
            ### on this page we find :
            ###   - a box with :
            ###       - a checkboxInput :  
            ###            - if box is check on : parameter_tabs is set on annotation
            ###               - On this page we find :
            ###                   - a box which countain :
            ###                       - a selectInput of separator 
            ###                       - a fileInput to upload metadata of RNA-seq experience
            ###                   - a box which countain :
            ###                       - information about file accepted in fileInput
            ###            - if box is check off : parameter_tabs is set on nothing
            ###               - On this page we find : nothing 
              tabItem(tabName = "Annotation",
                    fluidPage(
                      box(width = 12,
                          checkboxInput("CheckAnnotation","Do you have an annotation file ?",value=FALSE)),
                      fluidRow(
                        parameters_Annotation),
                      dataTableOutput("AnnoTable")
                    )
            ),
            ### Run DESeq2 ----
            ### On this page we find little information about how run DESeq2 workflow
            ### a actionButton when press it run DESeq2 workflow 
            ### a uiOutput("") set on server function
            tabItem(tabName = "deseq2",
                    waiter::use_waiter(),
                    fluidPage(
                      box(width = 12, solidHeader = F,
                          HTML(" <center><h3>Here you gonna run DESeq2 workflow.</h3> </pre>
                               <br><p> Check if your design chosen previously is correct.
                                  <br>
                               <br><h7>The process will take a few seconds.</h7></center>")),
                      box(width = 12,

                          actionButton("RunDESeq2","Run DESeq2 Workflow ",icon = icon("fas fa-user-astronaut"), class="btn btn-danger btn-lg btn-block ")),
                      uiOutput("SuccessMessage")
                      
                          )
                    ),
            
            ### Count distribution page ----
            ### On this page we find count distribution plot and it parameters after running DESeq
            ### On this page we find :
            ###     - a box() which countain :
            ###         - selectInput of sample which set sample on count.distribution.plot function
            ###         - sliderInput of break width which set break.width on count.distribution.plot function
            ###         - sliderInput range of x.axis which set x.min and x.max on count.distribution.plot function
            ###     - a box() which countain :
            ###         - return of count.distribution.plot 
            tabItem(tabName = "Count_Distribution",
                    box(title="Count distribution",solidHeader = T, status = "primary",width=12,collapsible = TRUE,
                        column(width = 6,
                               selectInput("sample","Which sample do you want to see ?", choices = c())
                        ),
                        column(width = 6,
                               sliderInput("breaksDistribution","Break size",min=0,max=2,value=1.0,step = 0.25)
                        ),
                        column(width = 6,
                               sliderInput("axis","Axis x",min=0,max=20,value=c(0,14))
                        ),
                        column(width = 6, checkboxInput("normalizeDistribution","Do you want to see distribution after normalisation ?",value=FALSE)
                        )
                    ),
                    box(width=12,status = "primary",withSpinner(plotOutput("CountDistributionPlot"))),
                    column(width= 4,
                           downloadButton("downloadDistribution",'Download plot',class = "btn-warning")
                    )
                    
            ),
            
            ### Count by gene page ----
            ### On this page we find count by gene plot and it parameters after running DESeq
            ### On this page we find :
            ###     - a box() which countain :
            ###         - selectizeInput of gene which set sample on count.gene.plot function
            ###         - a checkBox normalization if checkbox is check ON dds.count use is normalize, else dds.count use is not normalize
            ###     - a box() which countain :
            ###         - return of count.gene.plot 
            ###     - a downloadButton to download  the generate plot
            tabItem(tabName = "Count_Gene",
                    box(title="Count by gene",solidHeader = T, status = "primary",width=12,collapsible = TRUE,
                        column(width = 6,
                               
                               selectizeInput("gene","Which gene do you want to see ?", choices = NULL)
                               
                        ),
                        column(width = 6, checkboxInput("normalizeCountGene","Do you want to see distribution after normalisation ?",value=FALSE)
                        )
                    ),
                    box(width=12,status = "primary",withSpinner(plotOutput("CountGenePlot"))),
                    column(width= 4,
                           downloadButton("downloadCountgene",'Download plot',class = "btn-warning")
                    )
                    
            ),
            ### Depth plot ----
            ### On this page we find depth sample plot and it parameters after running DESeq
            ### On this page we find :
            ###     - a box() which countain :
            ###         - sliderInput of break width which set break.width on depth.plot function
            ###         - a checkBox normalization if checkbox is check ON dds.count use is normalize, else dds.count use is not normalize
            ###     - a box() which countain :
            ###         - return of depth.plot 
            ###     - a downloadButton to download  the generate plot
            tabItem(tabName = "Depth",
                    box(title="Depth of Sample",width = 12,solidHeader = T, status = "primary",collapsible = TRUE,
                        sliderInput("breaksDepth","Bar size",min=0,max=4,value=0.75,step = 0.25),
                        checkboxInput("normalizeDepth","Do you want to see depth after normalisation ?",value=FALSE)
                    ),
                    box(width=12,status = "primary",withSpinner(plotOutput("depth",height = 500))),
                    column(width= 4,
                           downloadButton("downloadDepth",'Download plot',class = "btn-warning")
                    )
            ),
            ### PCA plot ----
            ### On this page we find PCA plot and it parameters after running DESeq
            ### We find :
            ###     - a box() which countain :
            ###         - a selectInput of intgroup for pca.plot function
            ###         - a selectInput to chose transformation for pca.plot 
            ###         - a actionButton to run pca?plot function
            ###     - a box() which countain :
            ###         - return of pca.plot 
            ###     - a downloadButton to download  the generate plot
            tabItem(tabName = "pca",
                    box(width = 12,
                        title = "PCA", solidHeader = T, status = "primary",collapsible = TRUE,
                        selectInput("TransformationPCA",label= "Choose your transformation",choices = c("Variance-stabilizing transformation"="vst","Log transformation"="rld")),
                        selectInput("conditionpca","Choose your intgroup for PCA ?", choices = c()),
                        actionButton("runPCA","Run PCA")
                    ),
                    box(solidHeader = F, status = "primary",width = 12, align = "center",
                        withSpinner(plotOutput("PCAplot",height = 650))
                    ),
                    column(width= 4,
                           downloadButton("downloadPCA",'Download plot',class = "btn-warning")
                    )
            ),
            ### Dispersion plot ----
            ### On this page we find dispersion plot after running DESeq
            ### We find :
            ###     - a box() which countain :
            ###         - return of dispersion.plot 
            ###     - a downloadButton to download  the generate plot
            tabItem(tabName = "Dispersion",
                    box(width = 12,
                        title = "Dispersion", solidHeader = T, status = "primary",collapsible = TRUE,
                        withSpinner(plotOutput("dispersionPlot",height = 650))),
                    column(width= 4,
                           downloadButton("downloadDispersion",'Download plot',class = "btn-warning")
                    )
            ),
            ### MA plot ----
            ### On this page we find MA plot and it parameters after running DESeq
            ### We find :
            ###     - a box() which countain :
            ###         - sliderInput of P.value which set p.val of ma.plot function
            ###         - a tableOutput() of number.DE.gene function
            ###     - a box() which countain :
            ###         - return of ma.plot 
            ###     - a downloadButton to download  the generate plot
            tabItem(tabName = "MAplot",
                    box(width = 12,
                        title = "MA plot", solidHeader = T, status = "primary",collapsible = TRUE,
                        sliderInput("pvalueMAplot", "Choose your pvalue", min=0, max=1, value=0.05),
                        tableOutput("numberDEgenes")
                    ),
                    box(solidHeader = F, status = "primary",width = 12,
                        withSpinner(plotOutput("MAplot",height = 650))),
                    column(width= 4,
                           downloadButton("downloadMaplot",'Download plot',class = "btn-warning"))
            ),
            ### Volcano plot ----
            ### On this page we find MA plot and it parameters after running DESeq
            ### We find :
            ###     - a box() which countain :
            ###         - sliderInput of P.value which set p.val of volcano.plot function
            ###         - a checkbox Yes/No if there an annotation
            ###              - if check Yes uiOutput("sliderFoldVolcano") and uiOutput("SliderLogVolcanon") appear
            ###     - a box() which countain :
            ###         - return of volcano.plot 
            ###     - a downloadButton to download  the generate plot
            tabItem(tabName = "Volcanoplot",
                    fluidPage(
                      box(width = 12,
                          title = "Volcano plot", solidHeader = T, status = "primary",collapsible = TRUE,
                          checkboxInput("annotationVolcano","Do you have an annotation file ?",value=FALSE),
                          sliderInput("pvalueVolcano", "Choose your pvalue", min=0, max=1, value=0.05),
                          uiOutput("AnnoVolcano"),
                          uiOutput("SliderFoldVolcano"),
                          uiOutput("SliderLogVolcano")
                      ),
                      box( solidHeader = F, status = "primary",width = 12,
                           withSpinner(plotOutput("volcanoPlot",height = 650))),
                      column(width= 4,
                             downloadButton("downloadVolcano",'Download plot',class = "btn-warning")))
            ),
            
            
            ### sample distance matrix heat map ----
            ### On this page we find distance matrix heat map and it parameters after running DESeq
            ### We find :
            ###     - a box() which countain :
            ###         - a selectInput to chose transformation for distance.matrix.heatmap  
            ###         - a actionButton to run distance.matrix.heatmap function
            ###     - a box() which countain :
            ###         - return of distance.matrix.heatmap
            ###     - a downloadButton to download  the generate plot
            tabItem(tabName = "DistanceMatrix",
                    
                    box(width = 12,
                        title = "Heat map", solidHeader = T, status = "primary",collapsible = TRUE,
                        selectInput("TransformationMatrix",label= "Choose your transformation",choices = c("Variance-stabilizing transformation"="vst","Log transformation"="rld")),
                        actionButton("RunMatrix","Run Heat map")),
                    box(solidHeader = F, status = "primary",width = 12, align = "center",
                        withSpinner(plotOutput("DistanceMatrixMap",height = 650))
                    ),
                    column(width= 4,
                           downloadButton("downloadDistanceMatrix",'Download plot',class = "btn-warning")
                    )
            ),
            ### Gene expression heatmap ----
            ### On this page we find gene expression heatmap and it parameters after running DESeq
            ### We find :
            ###     - a box() which countain :
            ###         - a selectInput to chose transformation for gene.expression.heatmap  
            ###         - a selectInput to chose condition for gene.expression.heatmap
            ###         - a annotation checkboxInput() for is.Anno of gene.expressions.heatmap function
            ###         - a actionButton to run gene.expression.heatmap function
            ###         - a sliderInput to set the number of genes to display in gene.expression.heatmap
            ###     - a box() which countain :
            ###         - return of distance.matrix.heatmap
            ###     - a downloadButton to download  the generate plot
            
            tabItem(tabName = "Heatmap",
                    waiter::use_waiter(),
                    box(width = 12,
                        title = "Heat map", solidHeader = T, status = "primary",collapsible = TRUE,
                        column(width=6, selectInput("TransformationHeatmap",label= "Choose your transformation",choices = c("Variance-stabilizing transformation"="vst","Log transformation"="rld")),
                               checkboxInput("annotationHeatmap","Do you have a Annotation file ?",value=FALSE)
                        ),
                        column(width=6,
                               selectInput("conditionHeatmap","Choose your condition for Heat map ?", choices = c()),
                               actionButton("RunHeatmap","Run Heat map")),
                        
                        column(width=12,
                               sliderInput("nbGenes",label="Choose the number of genes you want to display", min = 0, 
                                           max = 200, value = c(0, 60)))),
                    uiOutput("Annoheatmap"),
                    box(solidHeader = F, status = "primary",width = 12, align = "center",
                        withSpinner(plotOutput("Heatmap", height = 1000, width = 1000))
                    ),
                    column(width= 4,
                           downloadButton("downloadHeatmap",'Download plot',class = "btn-warning")
                    )
            )
                           )
            )
            )
          )
                    ),
  ### footer settings ----
  tags$footer(
    wellPanel(
      HTML('
           <p align="center" width="4">Developed by <a href="https://www.linkedin.com/in/david-gallien-2096b9193/" target="_blank">David Gallien</a> and <a href="https://www.linkedin.com/in/gabin-coudray-a1941913b/" target="_blank">Gabin Coudray</a>. </p>
           <p align="center" width="4">First year of <a href="http://bioinfo-rennes.fr/" target="_blank">Bioinformatics Master<span>&#39;</span>s degree</a> in Rennes. </p>
           <p align="center" width="4"> <a href="https://www.univ-rennes1.fr/" target="_blank">University of Rennes 1.</a> </p>'
      ), 
      style = 
        "
      position:relative;
      width:100%;
      background-color: rgb(70,80,90);
      color: rgb(205,205,205);"
      ))
  )
