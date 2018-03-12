output$ui_altsplice <- renderUI({
  fluidPage(
       includeCSS("www/style.css"),
         sidebarLayout(
            sidebarPanel(
                uiOutput("dir.splice"),
                       textInput("symsearch", "Search By Gene Symbol", 'IL1B'),
                       textInput("enssearch", "Search By ENS ID",''),
                       numericInput("kct",label = "Number of Clusters",value = 2),   
                       actionButton("altButton", "GO")
                ),
            mainPanel(
                tabsetPanel(
                       tabPanel("Gene Compare",
                             dataTableOutput('gene.stat'),
                             dataTableOutput('trx.name'),
                             plotOutput("plot.cluster"),
                             plotOutput("trx.gene"),
                             plotOutput("plot.means")
                             )
                    )
                )
            )
        )
    
})

