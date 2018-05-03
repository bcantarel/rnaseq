output$ui_qc <- renderUI({
  list(fluidPage(
    includeCSS("www/style.css"),
    fluidRow(
      sidebarPanel(
        uiOutput("dir.qc"),
        p("Click it to load data"),
        actionButton("qcButton", "GO",  class = "btn btn-primary btn-bg centerbtn")
      ),
      column(7,
        tabsetPanel(
          tabPanel(
            "Align Stats",
            br(),br(),
            dataTableOutput('align.stats'),
            br(),
            uiOutput('downloadAlign')
          ),
          
          tabPanel(
            "Gene Type Stats",
            br(),br(),
            dataTableOutput('feature.stats'),
            uiOutput('downloadGT'),
            plotOutput("plot.fs"),
            uiOutput('dlplot'),
            br(),
            textOutput("fs.desc")
          ),
          tabPanel(
            "MDS Analysis",
            br(),br(),
            imageOutput("pca2", width = "100%", height = "auto"),
            uiOutput('dlmds')
          ),
          tabPanel(
            "Sample Distances",
            br(),br(),
            imageOutput("hmap", width = "100%", height = "auto"),
            uiOutput('dlhmap')
          ),
          tabPanel(
            "PCA Analysis",
            br(),br(),
            imageOutput("pca1", width = "100%", height = "auto"),
            uiOutput('dlpca')
          )
        )
      )
    )
    
  ))
})
