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
                      box(title="Upload annotation file",width = 12, status = "primary", solidHeader = TRUE,collapsible = TRUE,
                          column(width=5,selectInput("sepanno", "Separator:", c("Comma" = ",", "Tab" = "\t", "Semi-colon" = ";"))),
                          fileInput("file3", "Upload annotation file", accept = c(".csv",".txt",".tsv")))),
               column(width= 6,     
                      box(
                        title = "Accepted files :", width = 12,
                        HTML(
                          "<li> .csv / .tsv / .txt files </li>
                          <li> Separated by tabulation, comma or semi-colon </li>
                          <li> </li>
                          <li> </li>"),
                        height = 160
                        )
                      )
                      )
               )
               )
             )

parameter_volcano <- tagList(
  tags$style("#params { display:none; }"),
  tabsetPanel(
    id="param_volc",
    tabPanel("No"),  
    tabPanel("Yes",
             sliderInput("sliderfold", "Chose your fold", min=-20, max=20, value=c(-6,6)),
             sliderInput("sliderlog", "Chose your log10", min=0, max=300, value=30))
  )
)






### User Interface  ----
ui <- tagList(
  
  ### Parameters of the dashboard ----
  div(
    id = "app",
    dashboardPage(
      ### Customize the header ----
      dashboardHeader(title = "DESEQ DATA COUNT", 
                      ### Home button----
                      tags$li(a(onclick = "openTab('Intro')",
                                href = NULL,
                                icon("home"),
                                title = "Homepage",
                                style = "cursor: pointer;"),
                              class = "dropdown")
                      
      ),
      
      ### Differents menu (page) in the sidebar ----
      dashboardSidebar(
        sidebarMenu(id="mysidebar",
                    menuItem(text = "Informations", tabName = "Intro", icon = icon("info-circle")),
                    menuItem(text = "1 Upload data", tabName = "upload", icon = icon("arrow-circle-up"),startExpanded = FALSE,
                             menuSubItem(text = "1.1 Input count table", tabName = "Input"),
                             menuSubItem(text = "1.2 Input metadata table", tabName = "Input2"),
                             menuSubItem(text = "1.3 Input annotation file", tabName = "Input3")),
                    menuItem(text = "2 Run DESeq2", tabName = "deseq2", icon = icon("play-circle")),
                    menuItem(text = "3 Results", tabName = "deseq2", icon = icon("poll"),startExpanded = FALSE,
                             menuSubItem("Count distribution",tabName = "count_distribution"),
                             menuSubItem("Count by gene", tabName = "count_gene"),
                             menuSubItem("Depth of sample",tabName = "depth"),
                             menuSubItem("Dispersion",tabName = "dispersion"),
                             menuSubItem("PCA",tabName = "pca"),
                             menuSubItem("MA Plot",tabName = "ma"),
                             menuSubItem("Volcano Plot",tabName = "vulcano"),
                             menuSubItem("Distance matrix",tabName = "heatmap1"),
                             menuSubItem("Heatmap",tabName = "heatmap2")
                             
                    ),
                    tags$hr(),
                    menuItem(icon = NULL,
                             materialSwitch(inputId = "theme", label = "Theme", status = "default", value= TRUE)
                    ),tags$hr()
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
                      includeMarkdown("intro.Rmd"))
            ),
            ### Upload count table ----
            tabItem(tabName = "Input",
                    column(width = 6,
                           box(title="Upload count table",width = 12, status = "primary", solidHeader = TRUE,collapsible = TRUE,
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
                           box(title="Upload metadata table",width = 12, status = "primary", solidHeader = TRUE,collapsible = TRUE,
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
                               <li> First column has to be gene ID or gene name</li>
                               <li> All others columns are count for each sample</li>"),
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
                      dataTableOutput("table3"))
                    
            ),
            ### Run DESeq2 ----
            tabItem(tabName = "deseq2",
                    waiter::use_waiter(),
                    fluidPage(
                      box(width = 12, solidHeader = F,
                          HTML("Running DESeq2 will create a dds object .....................
                               <br> Check if the design chosen previously is good.
                               <br>If it is not, the application will crash.")),
                      box(width = 12,
                          actionButton("deseq2","Run DESeq2",icon = icon("fas fa-user-astronaut"), class="btn btn-danger btn-lg btn-block ")),

                      dataTableOutput("table4"))
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
                    box(width=12,status = "primary",plotOutput("count")),
                    column(width= 4,
                           downloadButton("downloadDistribution",'Download plot',class = "btn-warning")
                    )
                    
            ),
            #Count by gene ---
          tabItem(tabName = "count_gene",
                   box(title="Count by gene",solidHeader = T, status = "primary",width=12,collapsible = TRUE,
                       column(width = 6,
                               selectInput("gene","Which gene do you want to see ?", choices = c())
                        ),
                        column(width = 6, checkboxInput("normalize4","Do you want to see distribution after normalisation ?",value=FALSE)
                        )
                    ),
                    box(width=12,status = "primary",plotOutput("countgene")),
                    column(width= 4,
                           downloadButton("downloadCountgene",'Download plot',class = "btn-warning")
                    )
                    
            ),
            ### Depth plot ----
            tabItem(tabName = "depth",
                    box(title="Depth of Sample",width = 12,solidHeader = T, status = "primary",collapsible = TRUE,
                        sliderInput("breaks1","Bar size",min=0,max=4,value=1.0,step = 0.25),
                        checkboxInput("normalize1","Do you want to see depth after normalisation ?",value=FALSE)
                    ),
                    box(width=12,status = "primary",plotOutput("depth")),
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
                        plotOutput("pca")
                    ),
                    column(width= 4,
                           downloadButton("downloadPCA",'Download plot',class = "btn-warning")
                    )
            ),
            ### Dispersion plot ----
            tabItem(tabName = "dispersion",
                    box(width = 12,
                        title = "Dispersion", solidHeader = T, status = "primary",collapsible = TRUE,
                        plotOutput("dispersionPlot")),
                    column(width= 4,
                           downloadButton("downloadDispersion",'Download plot',class = "btn-warning")
                    )
            ),
            ### MA plot ----
            tabItem(tabName = "ma",
                    box(width = 12,
                        title = "MA plot", solidHeader = T, status = "primary",collapsible = TRUE,
                        sliderInput("pvalue", "Chose your pvalue", min=0, max=1, value=0.05)),
                    box(solidHeader = F, status = "primary",width = 12,
                        plotOutput("maplot")),
                    column(width= 4,
                           downloadButton("downloadMaplot",'Download plot',class = "btn-warning"))
            ),
            ### Volcano plot ----
            tabItem(tabName = "vulcano",
                    fluidPage(
                      box(width = 12,
                          title = "Volcano plot", solidHeader = T, status = "primary",collapsible = TRUE,
                          sliderInput("pvalue2", "Chose your pvalue", min=0, max=1, value=0.05)
                      ),
                      box(width = 12,
                          title = "Volcano plot", solidHeader = T, status = "primary",collapsible = TRUE,
                          checkboxInput("annotation3","Do you have an annotation file",value=FALSE),
                          uiOutput("annotationUi"),
                          uiOutput("annotationUi2")
                          
                      ),
                      
                      box( solidHeader = F, status = "primary",width = 12,
                           plotOutput("volcano")),
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
                        plotOutput("clusteringmap")
                    ),
                    column(width= 4,
                           downloadButton("downloadHeatmap1",'Download plot',class = "btn-warning")
                    )
            ),
            ### Heat map 2 ----
            tabItem(tabName = "heatmap2",
                    box(width = 12,
                        title = "Heat map", solidHeader = T, status = "primary",collapsible = TRUE,
                        column(width=6, selectInput("log3",label= "Choose your transformation",choices = c("Variance-stabilizing transformation"="vst","Log transformation"="rld")),
                               checkboxInput("annotation2","Do you have a Annotation file ?",value=FALSE)
                        ),
                        column(width=6,
                               selectInput("conditionheatmap","Choose your condition for Heat map ?", choices = c()),
                               actionButton("logaction3","Run Heat map")),
                        sliderInput("slider2", label = h4("Chose the number of genes you want to display"), min = 0, 
                                    max = 200, value = c(0, 60))),
                    box(solidHeader = F, status = "primary",width = 12,
                        plotOutput("clusteringmap2", height = 1000, width = 1000)
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