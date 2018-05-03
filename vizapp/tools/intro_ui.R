output$ui_intro <- renderUI({
  list(
    #includeCSS("www/style.css"),
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
    fluidPage(id = 'intropage',
              fluidRow(
                br(),
                column(
                  8,
                  offset = 2,
                  h1("Welcome to the RNASeq Analysis Portal", align = "center"),
                  br(),
                  br(),
                  div(column(
                    4,
                    div(
                      class = "card text-center",
                      div(class = "card-header",
                          "Gene Compare"),
                      div(
                        class = "card-body",
                        h5(
                          class = "card-title",
                          "To Compare Gene Abundances",
                          align = "center"
                        ),
                        p(
                          class = "card-text",
                          "Here you will be able explore gene abundances by (1) Transmembrane Proteins, (2) Transcription Fractors and (3) Chemocytokines.  Also you can compare the gene abundances of genes by Available Clinicial Data"
                        )
                      ),
                      div(
                        class = "card-footer",
                        # a(class="btn btn-primary btn-sm", 'data-toggle'="tab", 'data-value'="Gene Compare", 'id'="btgc", "Gene Compare")
                        actionButton(
                          class = "btn btn-primary btn-sm",
                          inputId = "actionbtgc",
                          label = "Gene Compare",
                          icon = NULL
                        )
                      )
                    )
                  ),
                  column(
                    4,
                    div(
                      class = "card text-center",
                      div(class = "card-header",
                          "DEA"),
                      div(
                        class = "card-body",
                        h5(class = "card-title", "To Examine Differential Gene Analysis"),
                        p(
                          class = "card-text",
                          "Here you will be able explore Differential Gene Expresss by Group.  You will also be able to compare two comparisons: such as Case1 vs Control and Case2 vs Control."
                        )
                      ),
                      div(
                        class = "card-footer",
                        # a(class="btn btn-primary btn-sm", 'data-toggle'="tab",'data-value'="DEA",id="dea","DEA")
                        actionButton(
                          class = "btn btn-primary btn-sm",
                          inputId = "actionbtdea",
                          label = "DEA",
                          icon = NULL
                        )
                      )
                    )
                  ),
                  column(
                    4,
                    div(
                      class = "card text-center",
                      div(class = "card-header",
                          "QuSAGE"),
                      div(
                        class = "card-body",
                        h5(class = "card-title", "To Examine Gene Set Enrichment Analysis"),
                        p(
                          class = "card-text",
                          "Here you will be able explore Gene Set Enrichment Analysis by Group and Gene Set List."
                        )
                      ),
                      div(
                        class = "card-footer",
                        #a(class="btn btn-primary btn-sm",'data-toggle'="tab", 'data-value'="QuSAGE",id="gsea","GSEA")
                        actionButton(
                          class = "btn btn-primary btn-sm",
                          inputId = "actionbtgsea",
                          label = "QuSAGE"
                        )
                      )
                    )
                  ))
                )
              ),
              fluidRow(column(
                8, offset = 2,
                br(),
                br(),
                div(
                  class = "alert alert-info",
                  a(
                    icon("new-window", lib = "glyphicon"),
                    strong("CITATION LINK"),
                    href = "https://git.biohpc.swmed.edu/BICF/Astrocyte/rnaseq",
                    target = "_blank",
                    style = "color:white"
                  ),
                  p(
                    strong("Note:"),
                    "Error messages can indicate data processing or missing data.  Please wait 30 seconds for the program to catch up to any changes in data loading.",
                    style = "font-size:14px"
                  )
                )
              )))
    
  )
  
})
