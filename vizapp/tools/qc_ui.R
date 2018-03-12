output$ui_qc <- renderUI({
    fluidPage(
        includeCSS("www/style.css"),
        sidebarLayout(
            sidebarPanel(
                uiOutput("dir.qc"),
                actionButton("qcButton", "GO")
                ),
            mainPanel(
                tabsetPanel(
                     tabPanel("Align Stats",
                             dataTableOutput('align.stats')
			     ),
		tabPanel("Gene Type Stats",
                             dataTableOutput('feature.stats'),
                             plotOutput("plot.fs"),
                             textOutput("fs.desc")
			     ),
		    tabPanel("MDS ANALYSIS",
                             h1("MDS"),
                             imageOutput("pca2",width="1200px",height="1200px")
                             ),
		    tabPanel("Sample Distances",
                             h1("Sample Distances"),
                             imageOutput("hmap",width="1200px",height="1200px")
                             ),
                    tabPanel("PCA ANALYSIS",
                             h1("PCA"),
                             imageOutput("pca1",width="1200px",height="1200px")
                             )
                    )
                )
            )
        )
    
})

