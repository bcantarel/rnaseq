output$ui_gsea <- renderUI({
  list(fluidPage(
    includeCSS("www/style.css"),
    fluidRow(
      sidebarPanel(
        uiOutput("quresult"),
        uiOutput("gset"),
        actionButton("quButton", "Go", class = "btn btn-primary btn-bg centerbtn")
      ),
      column(7,
             tabsetPanel(
               tabPanel(
                 "Gene Set Comparisons",
                 br(),
                 br(),
                 dataTableOutput('gsea.tbl'),
                 uiOutput("downloadGS")
               ),
               tabPanel(
                 "Gene Set Comparison",
                 plotOutput("plot.cis"),
                 uiOutput("dlcis"),
                 plotOutput("plot.genedist"),
                 uiOutput("dlgenedist")
               )
             ))
    )
  ))
})
