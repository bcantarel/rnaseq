library(shiny)

shinyUI(navbarPage("RNASeq Analysis",
                   id = "nav_shinyapp", inverse = TRUE, collapsible = TRUE,
                   tabPanel("Introduction",uiOutput('ui_intro')),
		   tabPanel("Dataset QC", uiOutput('ui_qc')),
                   #tabPanel("Sample Sequencing QC", uiOutput('ui_fastqc')),
		   tabPanel("Gene Compare", uiOutput('ui_gc')),
                   tabPanel("Gene Alt Splicing", uiOutput('ui_altsplice')),
                   tabPanel("DEA", uiOutput('ui_dea')),
                   tabPanel("QuSAGE", uiOutput('ui_gsea'))
                   ))


