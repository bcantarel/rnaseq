output$ui_intro <- renderUI({
    list(
        includeCSS("www/style.css"),
        sidebarLayout(
            sidebarPanel(
                ),
            mainPanel(
                h1("Welcome to the RNASeq Analysis Portal", align = "center"),
                br(),
                br(),
                h4("To Compare Gene Abundances Use Gene Compare", align = "center"),
                p("Here you will be able explore gene abundances by (1) Transmembrane Proteins, (2) Transcription Fractors and (3) Chemocytokines.  Also you can compare the gene abundances of genes by Available Clinicial Data",style = "color:red"),
                h4("To Examine Differential Gene Analysis Use DEA", align = "center"),
                p("Here you will be able explore Differential Gene Expresss by Group.  You will also be able to compare two comparisons: such as Case1 vs Control and Case2 vs Control.",style = "color:red"),
                h4("To Examine Gene Set Enrichment Analysis Use GSEA", align = "center"),
                p("Here you will be able explore Gene Set Enrichment Analysis by Group and Gene Set List.",style = "color:red"),
                br(),
                br(),
                p("Note:Error messages can indicate data processing or missing data.  Please wait 30 seconds for the program to catch up to any changes in data loading."),
                h6("Questions? Contract Brandi Cantarel brandi.cantarel@baylorhealth.edu")
                )
            )
        )
    
})
