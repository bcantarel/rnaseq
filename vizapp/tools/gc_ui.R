output$ui_gc <- renderUI({
    fluidPage(
        includeCSS("www/style.css"),
        sidebarLayout(
            sidebarPanel(
                uiOutput("pick.group"),
                textInput("symsearch", "Search By Gene Symbol", 'IL1B'),
                textInput("enssearch", "Search By ENS ID",''),
                actionButton("gcButton", "GO")
                ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Gene Compare",
                             dataTableOutput('gc.stats'),
                             plotOutput("plot.gene"),
                             textOutput("bxplot.desc"),
                             plotOutput("violin.gene"),
                             textOutput("violin.desc")
			     )
                    )
                )
            )
        )
    
})

