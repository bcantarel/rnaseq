output$ui_gc <- renderUI({
  list(fluidPage(
    includeCSS("www/style.css"),
    fluidRow(
      sidebarPanel(
        uiOutput("pick.group"),
        textInput("gc_symsearch", "Search By Gene Symbol", 'IL1B'),
        textInput("gc_enssearch", "Search By ENS ID", ''),
        actionButton("gcButton", "GO", class = "btn btn-primary btn-bg centerbtn")
      ),
      column(7,
             tabsetPanel(
               tabPanel(
                 "Gene Compare",
                 plotOutput("plot.gene"),
                 uiOutput('dlboxplot')
               )
             )),
      tags$style(type='text/css', "#dlboxplot { width:100%;margin-left: 35px;}")
    )
  ))
  
})
