library(markdown)
output$ui_fastqc <- renderUI({
    fluidPage(
    includeCSS("www/style.css"),
        sidebarLayout(
            sidebarPanel(
                uiOutput("dir.fqc"),
                uiOutput("pick.fqc"),
		actionButton("fqcButton", "GO")
                ),
            mainPanel(
	    htmlOutput("inc")
)
)
)
})
