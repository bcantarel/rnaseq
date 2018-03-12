get.qc <- function(var) {
    ctfile <-  paste(data.dir,'countTable.stats.txt',sep='/')
    alfile <-  paste(data.dir,'alignment.summary.txt',sep='/')
    stats <- read.table(file=ctfile,header=TRUE,sep='\t')
    align <- read.table(file=alfile,header=TRUE,sep='\t')
    return(list(featurestat=stats,alignstat=align))
  } 

forstat <- eventReactive(input$qcButton,{get.qc(input)})

output$plot.fs <- renderPlot({
    countTable <- forstat()$featurestat
    par(oma=c(4,4,1,1)) 
    grpnames <- levels(factor(as.character(countTable$Type)))
    boxplot(countTable$ReadPerc ~ countTable$Type,col=rainbow(length(grpnames)),cex.axis=0.7,ylab='Percent of Total Reads',las=2)
  }, height = 400, width = 900)

output$fs.desc <- renderText({ 
      paste("Percent of Total Reads for each gene type.  Boxplots draw to represent the 25th and 75th percentile (the lower and upper quartiles, respectively) as a box with a band in the box representing 50th percentile (the median). The upper whisker is located at the 'smaller' of the maximum x value and Q_3 + 1.5 inner quantile range(IRQ), whereas the lower whisker is located at the 'larger' of the smallest x value and Q_1 - 1.5 IQR",sep=' ')
    })

output$feature.stats <- DT::renderDataTable({
  t1 <- forstat()$featurestat
  t1
},escape=FALSE,filter = 'top',options = list(lengthMenu = c(10, 25, 50, 200, -1)))

output$align.stats <- DT::renderDataTable({
  t1 <- forstat()$alignstat
  t1
},escape=FALSE,filter = 'top',options = list(lengthMenu = c(10, 25, 50, 200, -1)))

output$pca1 <- renderImage({
    list(src=paste(data.dir,'pca.png',sep="/"), alt=paste("Sample Comparison 1"))
},deleteFile=FALSE)
output$pca2 <- renderImage({
    list(src=paste(data.dir,'mds.png',sep="/"), alt=paste("Sample Comparison 2"))
},deleteFile=FALSE)
output$hmap <- renderImage({
    list(src=paste(data.dir,'samples_heatmap.png',sep="/"), alt=paste("Sample Comparison 3"))
},deleteFile=FALSE)
