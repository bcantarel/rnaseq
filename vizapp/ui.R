library(shiny)
library(shinythemes)
shinyUI(list(
  navbarPage(
    title = "RNASeq Analysis",
    id = "nav_shinyapp",
    inverse = TRUE,
    collapsible = TRUE,
    theme = shinytheme("yeti"),
    
    tabPanel("Introduction", uiOutput('ui_intro')),
    tabPanel("Dataset QC", uiOutput('ui_qc')),
    #tabPanel("Sample Sequencing QC", uiOutput('ui_fastqc')),
    tabPanel(
      "Gene Compare",
      id = "panelgc",
      value = "panelgc",
      uiOutput('ui_gc')
    ),
    tabPanel("Gene Alt Splicing", uiOutput('ui_altsplice')),
    tabPanel("DEA", id = "paneldea", value = "paneldea", uiOutput('ui_dea')),
    tabPanel(
      "QuSAGE",
      id = "panelgsea",
      value = "panelgsea",
      uiOutput('ui_gsea')
    )
  ),
  fluidRow(style = "margin-top:40px", tags$footer(class = "footer",
                                                  div(class = "container",
                                                      style = "padding-bottom:0; margin-bottom:0",
                                                        p(icon("envelope", lib = "glyphicon"),
                                                          "brandi.cantarel@utsouthwestern.edu  |  @ UT Southwestern Medical Center",
                                                           style = "margin-bottom:-40px;")
                                                      )
                                                  )
           )
))
