
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
output$downloadGS <- renderUI({
  if(!is.null(values()$tbl)){
    downloadButton('DownloadGS', 'Download CSV')
  }
})
output$DownloadGS <- downloadHandler(
  file <- paste(input$fname,".txt",sep=""),
  content = function(file) {
    write.table(values()$tbl,file,quote=FALSE,row.names=FALSE,sep='\t')
  })
output$plot.cis <- renderPlot({
    qs.results <- values()$obj
    plotCIsGenes(qs.results,path.index=as.numeric(input$pathname), cex.xaxis = 1.3)
})
output$dlcis <- renderUI({
   if (!is.null(values()$obj)) {
  downloadButton('Downloadcis', 'Download PNG')
  }
})
output$Downloadcis <- downloadHandler(
  file = function() {
    paste('cisplot', 'png', sep = ".")
  },
  content <- function(file) {
    png(
      file,
      width = 8 * 120,
      height = 4 * 120,
      res = 120,
      pointsize = 5
    )
    if (is.null(values()$obj)) {
      return()
    }
    qs.results <- values()$obj
    plot <- plotCIsGenes(qs.results,path.index=as.numeric(input$pathname), cex.xaxis = 1.3)
    print(plot)
    dev.off()
  },
  contentType = "image/png"
)
output$plot.genedist <- renderPlot({
    qs.results <- values()$obj
    plotGeneSetDistributions(qs.results,path.index=as.numeric(input$pathname), cex =1.6)
})    

output$dlgenedist <- renderUI({
  if (!is.null(values()$obj)) {
    downloadButton('Downloadgenedist', 'Download PNG')
  }
})
output$Downloadgenedist <- downloadHandler(
  file = function() {
    paste('cisplot', 'png', sep = ".")
  },
  content <- function(file) {
    png(
      file,
      width = 8 * 120,
      height = 4 * 120,
      res = 140,
      pointsize = 5
    )
    if (is.null(values()$obj)) {
      return()
    }
    qs.results <- values()$obj
    plot <- plotGeneSetDistributions(qs.results,path.index=as.numeric(input$pathname), cex =1.6)
    print(plot)
    dev.off()
  },
  contentType = "image/png"
)
