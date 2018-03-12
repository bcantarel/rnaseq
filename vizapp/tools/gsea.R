
output$quresult <- renderUI({
    flist <- list.files(data.dir,pattern="*qusage.rda$")
    selectInput('fname', "Choose Comparison", choices=flist)
})

updateData <- function(input) {
    f1 <- paste(data.dir,input$fname,sep='/')
    load(file = f1, envir = .GlobalEnv)
    results <- qsTable(qs.results,number=100,sort.by='fdr')
    return(list(obj=qs.results,tbl=results))
}
values <-  eventReactive(input$quButton,{updateData(input)})

output$gset <- renderUI({
    tbl <- values()$tbl
    choices = setNames(row.names(tbl),tbl$pathway.name)	    
    selectInput("pathname", "Choose Pathway", choices=choices)
})

output$gsea.tbl <- renderDataTable({
    values()$tbl
})

output$plot.cis <- renderPlot({
    qs.results <- values()$obj
    plotCIsGenes(qs.results,path.index=as.numeric(input$pathname))
})

output$plot.genedist <- renderPlot({
    qs.results <- values()$obj
    plotGeneSetDistributions(qs.results,path.index=as.numeric(input$pathname))
})    
