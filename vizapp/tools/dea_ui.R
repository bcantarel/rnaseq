#######################################
# Shiny interface for data functions
#######################################
# data ui and tabs
output$ui_dea <- renderUI({
    list(
        includeCSS("www/style.css"),
        sidebarLayout(
            sidebarPanel(
                numericInput("fc.thresh",
                             "LogFold Change Threshold:", 1),
                numericInput("gene.thresh",
                             "Gene Mean Threshold:", 5),
                numericInput("pval.thresh",
                             "P-Value Threshold:", 0.05),
                selectInput("adjust", "Choose P-Value Correction", choices=c("raw","FDR",'BONF'),selected='FDR'),
                uiOutput("pick.dea"),
	     	selectInput(
			"heatmap", "HeatMap",
      	     		c(Top = "top",
               		HallmarkGeneSet = "hgeneset",
	       		CustomGeneSet= "cgeneset")
	       	),
	     	conditionalPanel(
      	     	condition = "input.heatmap == 'cgeneset'",
	     		  uiOutput("selectgenes")
		),
	     	conditionalPanel(
      	     	condition = "input.heatmap == 'top'",
	      	numericInput("numgenes","Number Top Genes:", 50)
		),
	     	conditionalPanel(
      	     	condition = "input.heatmap == 'hgeneset'",
	     	uiOutput("pick.pathway")
		),
	    	actionButton("deButton", "Go")
	    ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Differential Gene Set Comparison",
                             downloadButton('downloadC', 'Download CSV'),
                             dataTableOutput('dge.c')
                             ),
                    tabPanel("Heatmap",
                             #downloadButton('downloadpdf', 'Download Heatmap'),
			     h1("HeatMap Comparison"),
			     h3("Top 50 User Defined Genes"),
			     plotOutput("plot.heatmap"),
			     h3("All Differentially Expressed Genes"),
                             imageOutput("hm.comp",width="1200px",height="1200px"),
			     textOutput("hmcomp.desc")
                             )
                    )
                )
            )
        )
})
