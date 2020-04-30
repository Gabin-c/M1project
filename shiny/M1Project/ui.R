### Library ----
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

### Files with all the function needed to make plots ----
source("function_dds.R")

### Annotation pannel ----
parameter_tabs <- tagList(
  tags$style("#params { display:none; }"),
  tabsetPanel(
    id="params",
    tabPanel("nothing"),  
    tabPanel("annotation",
             fluidRow(
               column(width= 6,
                      box(title="Upload annotation file",width = 12, solidHeader = TRUE,collapsible = TRUE,
                          column(width=5,selectInput("sepanno", "Separator:", c("Comma" = ",", "Tab" = "\t", "Semi-colon" = ";"))),
                          fileInput("file3", "Upload annotation file", accept = c(".csv",".txt",".tsv")))),
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
ui <- tagList(
  
  ### Parameters of the dashboard ----
  div(
    id = "app",
    dashboardPage(
      ### Customize the header ----

      dashboardHeader(title = "RNA-seq DE analysis", 

                      ### Home button----
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
      
      ### Differents menu (page) in the sidebar ----
      dashboardSidebar(
        sidebarMenu(id="mysidebar",
                    menuItem(text = "Informations", tabName = "Intro", icon = icon("info-circle")),
                    menuItem(text = "1 Upload data", tabName = "upload", icon = icon("arrow-circle-up"),startExpanded = TRUE,
                             menuItemOutput("menuCheck"),
                             menuItemOutput("menuCheck1"),
                             menuItemOutput("menuCheck2")),
                    menuItem(text = "2 Run DESeq2", tabName = "deseq2", icon = icon("play-circle")),

                    menuItemOutput("menuResults"),

                    tags$hr(),
                    menuItem(icon = NULL,
                             materialSwitch(inputId = "theme", label = "Theme", status = "default", value= TRUE)
                    ),tags$hr(),
                    helpText("Developed by ", a("David Gallien ", href="https://www.linkedin.com/in/david-gallien-2096b9193/"), "and ", br(), a("Gabin Coudray ", href="https://www.linkedin.com/in/gabin-coudray-a1941913b/"), ", first year of ", br(),
                             a("Bioinformatics Master's degree ", href="http://bioinfo-rennes.fr/"), "in ", br(), "Rennes, ",
                             a("University of Rennes 1 ", href="https://www.univ-rennes1.fr/"), style="padding-left:1em; padding-right:1em;position:absolute; bottom:1em; ", align = "center")
        )
      ),
      
      ### Organization of the differents pages ----
      dashboardBody(
        uiOutput("themes"),
        useShinyjs(),
        fluidRow(
          tabItems(
            ### Introduction ----
            tabItem(tabName = "Intro",
                    fluidPage(
                      h2("Introduction"),
                      p("This is an R Shiny web interactive application developed as part of a ", 
                        strong("course project."), "The purpose of this application is to perform an ",
                        strong("differential expression analysis from a counts table"), "in order to help researchers getting interpretable results.",
                        align = "justify"),
                      p("This application uses the package ", 
                        a("DESeq2", href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html"), 
                        "from Bioconductor. It is a package to study differential gene expression analysis based on the negative binomial distribution. It allows a quick visualization of results in order to analyze the counts table data. The results will be given in different form like graphics, heatmaps, MAplot or even Volcano plot.",
                        align = "justify"),
                      tags$hr(),
                      h3("1. Upload data", style="padding-left: 1em"),
                      p("The input data files accepted for this App are 3 files in '.txt', '.csv' or '.tsv' format separated by comma, tabulation or semi-colon.
                        This App necessarily requires a 'Count Data Table' and a 'Metadata Table'. An optional 'Annotation File' can be added", style="padding-left: 2em", align = "justify"),
                      h4("1.1 Count Data Table", style="padding-left: 3em"),
                      p("The Count Data Table has to contains the count for each sample of the experiment for each gene and the first column has to be gene ID or gene name as below :",style="padding-left: 5em", align = "justify"),
                      column( 12, style="padding-left: 5em" ,withSpinner(tableOutput("countexample"))),
                      br(),
                      h4("1.2 Metadata Table", style="padding-left: 3em"),
                      p("The Metadata table has to contains the informations of the experiment with at least 2 columns. The first one is the samples in the same order as the columns of the Count Table. 
                        The second one is a condition column. You can add as many columns as you have factors in your experiment.",style="padding-left: 5em", align = "justify"),
                      column( 12, style="padding-left: 5em" ,withSpinner(tableOutput("metadataexample"))),
                      h4("1.2  Annotation File", style="padding-left: 3em"),
                      p("The Annotation File contains informations about the genes. If you have one, it must contains a column named 'symbol' in which we can find the symbol of each gene.",style="padding-left: 5em", align = "justify"),
                      column( 12, style="padding-left: 5em" ,withSpinner(tableOutput("annoexample"))),
                      h3("2. Results", style="padding-left: 1em"),
                      p("The results will be display after running DESeq2. You will obtain 9 differents results :", style="padding-left: 2em", align = "justify"),
                      p("- Count distribution",
                        br(), "- Count by gene",
                        br(), "- Depth of sample",
                        br(), "- Dispersion",
                        br(), "- PCA",
                        br(), "- MA plot",
                        br(), "- Volcano plot",
                        br(), "- Distance matrix",
                        br(), "- Heatmap",style="padding-left: 5em", align = "justify"),
                      p("You can download all the results plots at the bottom of all these pages.",  style="padding-left: 2em", align = "justify")
                      
                      )
            ),
            ### Upload count table ----
            tabItem(tabName = "Input",
                    column(width = 6,
                           box(title="Upload count table",width = 12, solidHeader = TRUE,collapsible = TRUE,
                               column(width=5,
                                      selectInput("sepcount", "Separator:", c("Comma" = ",", "Tab" = "\t", "Semi-colon" = ";"))),
                               fileInput("file", "Upload count table", accept = c(".csv",".txt",".tsv"))
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
                    dataTableOutput("table")
                    ),
            ### Upload metadata table ----
            tabItem(tabName = "Input2",
                   column(width = 6,
                           box(title="Upload metadata table",width = 12, solidHeader = TRUE,collapsible = TRUE,
                               column(width=5,
                                      selectInput("sepmetadata", "Separator:", c("Comma" = ",", "Tab" = "\t", "Semi-colon" = ";"))),
                               fileInput("file2", "Upload metadata table", accept = c(".csv",".txt",".tsv"))
                           )),
                    column(width = 6,
                           box(
                             title = "Accepted files :", width = 12,
                             HTML(
                               "<li> .csv / .tsv / .txt files </li>
                               <li> Separated by tabulation, comma or semi-colon </li>
                               <li> At least metadata table contains two column</li>
                               <li> At least one column has to be factor</li>"),
                             height = 160
                             )),
                    column(width = 12,
                           box(width = 12,
                               textInput("condition","Chose your design without linear combination", placeholder = "Conditions"))),
                    dataTableOutput("table2")
                    ),
            ### Uploade annotation file ----
            tabItem(tabName = "Input3",
                    fluidPage(
                      box(width = 12,
                          checkboxInput("annotation","Do you have an annotation file",value=FALSE)),
                      fluidRow(
                        parameter_tabs),
                      dataTableOutput("table3")
                    )
            ),
            ### Run DESeq2 ----
            tabItem(tabName = "deseq2",
                    waiter::use_waiter(),
                    fluidPage(
                      box(width = 12, solidHeader = F,
    
                          
                          
                          HTML(" <center><h3>Here you gonna run DESeq2 workflow.</h3> </pre>
                              
                               <br><h5> Check if  your design chosen previously is good.</h5>
                               <br><h5>If it is not, the application will crash.</h5></center>")),
                      box(width = 12,
                          actionButton("deseq2","Run DESeq2 Workflow ",icon = icon("fas fa-user-astronaut"), class="btn btn-danger btn-lg btn-block ")),
                      uiOutput("table4")
                     
                      )
            ),
            ### Count distribution plot ----
            tabItem(tabName = "count_distribution",
                    box(title="Count distribution",solidHeader = T, status = "primary",width=12,collapsible = TRUE,
                        column(width = 6,
                               selectInput("sample","Which sample do you want to see ?", choices = c())
                        ),
                        column(width = 6,
                               sliderInput("breaks","Break size",min=0,max=2,value=1.0,step = 0.25)
                        ),
                        column(width = 6,
                               sliderInput("axis","Axis x",min=0,max=20,value=c(0,14))
                        ),
                        column(width = 6, checkboxInput("normalize","Do you want to see distribution after normalisation ?",value=FALSE)
                        )
                    ),
                    box(width=12,status = "primary",withSpinner(plotOutput("count"))),
                    column(width= 4,
                           downloadButton("downloadDistribution",'Download plot',class = "btn-warning")
                    )
                    
            ),
            
            ### Count by gene ----
            tabItem(tabName = "count_gene",
                    box(title="Count by gene",solidHeader = T, status = "primary",width=12,collapsible = TRUE,
                        column(width = 6,
                               
                               selectizeInput("gene","Which gene do you want to see ?", choices = NULL)
                               
                        ),
                        column(width = 6, checkboxInput("normalize4","Do you want to see distribution after normalisation ?",value=FALSE)
                        )
                    ),
                    box(width=12,status = "primary",withSpinner(plotOutput("countgene"))),
                    column(width= 4,
                           downloadButton("downloadCountgene",'Download plot',class = "btn-warning")
                    )
                    
            ),
            ### Depth plot ----
            tabItem(tabName = "depth",
                    box(title="Depth of Sample",width = 12,solidHeader = T, status = "primary",collapsible = TRUE,
                        sliderInput("breaks1","Bar size",min=0,max=4,value=0.75,step = 0.25),
                        checkboxInput("normalize1","Do you want to see depth after normalisation ?",value=FALSE)
                    ),
                    box(width=12,status = "primary",withSpinner(plotOutput("depth",height = 500))),
                    column(width= 4,
                           downloadButton("downloadDepth",'Download plot',class = "btn-warning")
                    )
            ),
            ### PCA plot ----
            tabItem(tabName = "pca",
                    box(width = 12,
                        title = "PCA", solidHeader = T, status = "primary",collapsible = TRUE,
                        selectInput("log",label= "Choose your transformation",choices = c("Variance-stabilizing transformation"="vst","Log transformation"="rld")),
                        selectInput("conditionpca","Choose your intgroup for PCA ?", choices = c()),
                        actionButton("logaction","Run PCA")
                    ),
                    box(solidHeader = F, status = "primary",width = 12,
                        withSpinner(plotOutput("pca",height = 650))
                    ),
                    column(width= 4,
                           downloadButton("downloadPCA",'Download plot',class = "btn-warning")
                    )
            ),
            ### Dispersion plot ----
            tabItem(tabName = "dispersion",
                    box(width = 12,
                        title = "Dispersion", solidHeader = T, status = "primary",collapsible = TRUE,
                        withSpinner(plotOutput("dispersionPlot",height = 650))),
                    column(width= 4,
                           downloadButton("downloadDispersion",'Download plot',class = "btn-warning")
                    )
            ),
            ### MA plot ----
            tabItem(tabName = "ma",
                    box(width = 12,
                        title = "MA plot", solidHeader = T, status = "primary",collapsible = TRUE,
                        sliderInput("pvalue", "Chose your pvalue", min=0, max=1, value=0.05),
                        tableOutput("num_DE")
                    ),
                    box(solidHeader = F, status = "primary",width = 12,
                        withSpinner(plotOutput("maplot",height = 650))),
                    column(width= 4,
                           downloadButton("downloadMaplot",'Download plot',class = "btn-warning"))
            ),
            ### Volcano plot ----
            tabItem(tabName = "vulcano",
                    fluidPage(
                      box(width = 12,
                          title = "Volcano plot", solidHeader = T, status = "primary",collapsible = TRUE,
                          checkboxInput("annotation3","Do you have an annotation file",value=FALSE),
                          sliderInput("pvalue2", "Chose your pvalue", min=0, max=1, value=0.05),
                          
                          
                          
                          uiOutput("annotationUi"),
                          uiOutput("annotationUi2")
                          
                      ),
                      
                      box( solidHeader = F, status = "primary",width = 12,
                           withSpinner(plotOutput("volcano",height = 650))),
                      column(width= 4,
                             downloadButton("downloadVulcano",'Download plot',class = "btn-warning")))
            ),
            
            
            ### Heatmap ----
            tabItem(tabName = "heatmap1",
                    
                    box(width = 12,
                        title = "Heat map", solidHeader = T, status = "primary",collapsible = TRUE,
                        selectInput("log1",label= "Choose your transformation",choices = c("Variance-stabilizing transformation"="vst","Log transformation"="rld")),
                        actionButton("logaction2","Run Heat map")),
                    box(solidHeader = F, status = "primary",width = 12,
                        withSpinner(plotOutput("clusteringmap",height = 650))
                    ),
                    column(width= 4,
                           downloadButton("downloadHeatmap1",'Download plot',class = "btn-warning")
                    )
            ),
            ### Heat map 2 ----
            tabItem(tabName = "heatmap2",
                    waiter::use_waiter(),
                    box(width = 12,
                        title = "Heat map", solidHeader = T, status = "primary",collapsible = TRUE,
                        column(width=6, selectInput("log3",label= "Choose your transformation",choices = c("Variance-stabilizing transformation"="vst","Log transformation"="rld")),
                               checkboxInput("annotation2","Do you have a Annotation file ?",value=FALSE)
                        ),
                        column(width=6,
                               selectInput("conditionheatmap","Choose your condition for Heat map ?", choices = c()),
                               actionButton("logaction3","Run Heat map")),

                        column(width=12,
                        sliderInput("slider2",label="Chose the number of genes you want to display", min = 0, 
                                    max = 200, value = c(0, 60)))),
                    box(solidHeader = F, status = "primary",width = 12,
                        withSpinner(plotOutput("clusteringmap2", height = 1000, width = 1000))
                    ),
                    column(width= 4,
                           downloadButton("downloadHeatmap2",'Download plot',class = "btn-warning")
                    )
            )
            )
                    )
                           )
          )
            )
      )

