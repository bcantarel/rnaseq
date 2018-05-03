output$ui_altsplice <- renderUI({
  list(fluidPage(
    includeCSS("www/style.css"),
    fluidRow(
      sidebarPanel(
        uiOutput("dir.splice"),
        textInput("symsearch", "Search By Gene Symbol", 'IL1B'),
        textInput("enssearch", "Search By ENS ID", ''),
        numericInput("kct", label = "Number of Clusters", value = 2),
        actionButton("altButton", "GO", class = "btn btn-primary btn-bg centerbtn")
      ),
      column(7,
             tabsetPanel(
               tabPanel(
                 "Gene Compare",
                 br(),br(),
                 dataTableOutput('gene.stat'),
                 br(),br(),
                 dataTableOutput('trx.name'),
                 br(),
                 uiOutput('downloadTrx'),
                 plotOutput("plot.cluster"),
                 uiOutput('dlplotcluster'),
                 br(),br(),
                 plotOutput("trx.gene"),
                 uiOutput('dltrxgene'),
                 plotOutput("plot.means", height = 900, width = 'auto'),
                 uiOutput('dlplotmean')
               )
             ))
    )
  ))
  
})

