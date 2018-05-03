output$ui_gc <- renderUI({
  list(fluidPage(
    includeCSS("www/style.css"),
    fluidRow(
      sidebarPanel(
        uiOutput("pick.group"),
        textInput("symsearch", "Search By Gene Symbol", 'IL1B'),
        textInput("enssearch", "Search By ENS ID", ''),
        actionButton("gcButton", "GO", class = "btn btn-primary btn-bg centerbtn")
      ),
      column(7,
             tabsetPanel(
               tabPanel(
                 "Gene Compare",
                 dataTableOutput('gc.stats'),
                 br(),
                 plotOutput("plot.gene"),
                 uiOutput('dlboxplot'),
                 br(),
                 textOutput("bxplot.desc"),
                 br(),
                 plotOutput("violin.gene"),
                 uiOutput('dlviolinplot'),
                 br(),
                 textOutput("violin.desc")
               )
             )),
      tags$style(type='text/css', "#dlboxplot { width:100%;margin-left: 35px;}"),
      tags$style(type='text/css', "#dlviolinplot { width:100%;margin-left: 35px;}")
    )
  ))
  
})
