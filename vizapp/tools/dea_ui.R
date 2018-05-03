#######################################
# Shiny interface for data functions
#######################################
# data ui and tabs
output$ui_dea <- renderUI({
  list(fluidPage(
    includeCSS("www/style.css"),
    fluidRow(
      sidebarPanel(
        numericInput("fc.thresh",
                     "LogFold Change Threshold:", 1),
        numericInput("gene.thresh",
                     "Gene Mean Threshold:", 5),
        numericInput("pval.thresh",
                     "P-Value Threshold:", 0.05),
        selectInput(
          "adjust",
          "Choose P-Value Correction",
          choices = c("raw", "FDR", 'BONF'),
          selected = 'FDR'
        ),
        uiOutput("pick.dea"),
        selectInput(
          "heatmap",
          "HeatMap",
          c(
            Top = "top",
            HallmarkGeneSet = "hgeneset",
            CustomGeneSet = "cgeneset"
          )
        ),
        conditionalPanel(condition = "input.heatmap == 'cgeneset'",
                         uiOutput("selectgenes")),
        conditionalPanel(condition = "input.heatmap == 'top'",
                         numericInput("numgenes", "Number Top Genes:", 50)),
        conditionalPanel(condition = "input.heatmap == 'hgeneset'",
                         uiOutput("pick.pathway")),
        radioButtons("cluster", label = "Cluster",
                     choices = list("Display Cluster" = 1, "Hide Cluster" = 2), 
                     selected = 1),
        actionButton("deButton", "Go", class = "btn btn-primary btn-bg centerbtn")
      ),
      column(7,
             tabsetPanel(
               tabPanel(
                 "Differential Gene Set Comparison",
                 br(),br(),
                 dataTableOutput('dge.c', width = '100%'),
                 uiOutput('downloadC')
               ),
               tabPanel(
                 "Heatmap Comparison",
                 #downloadButton('downloadpdf', 'Download Heatmap'),
                 h3(textOutput("hmcomp.top50")),
                 plotOutput("plot.heatmap"),
                 uiOutput("dlheatmap"),
                 h3(textOutput("hmcomp.allgene")),
                 imageOutput("hm.comp", width = "100%", height = "auto"),
                 uiOutput('dlhmcomp'),
                 textOutput("hmcomp.desc")
               )
             ))
    )
  ))
})
