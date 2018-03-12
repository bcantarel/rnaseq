output$ui_gsea <- renderUI({
    list(
        includeCSS("www/style.css"),
        sidebarLayout(
            sidebarPanel(
                uiOutput("quresult"),
                uiOutput("gset"),
		actionButton("quButton", "Go")
                ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Gene Set Comparisons",
                             dataTableOutput('gsea.tbl')
                             ),
                    tabPanel("Gene Set Comparison",
                             plotOutput("plot.cis"),
                             plotOutput("plot.genedist")
                             )
                    )
                )
            )
        )
})
