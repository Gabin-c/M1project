library(shiny)
runApp("shiny/M1Project")

parameter_volcano <- tagList(
  tags$style("#params { display:none; }"),
  tabsetPanel(
    id="param_volc",
    tabPanel("nothin"),  
    tabPanel("annotat",
             box(width = 12,
                 sliderInput("sliderfold", "Chose your fold", min=-20, max=20, value=c(-6,6)),
                 sliderInput("sliderlog", "Chose your log10", min=0, max=300, value=30)
             )
    )
  )
)
