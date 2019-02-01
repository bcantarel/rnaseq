samtbl <- read_tsv(paste(data.dir,'design.shiny.txt',sep='/'))
wrapper <- function(x, ...) 
 {
   paste(strwrap(x, ...), collapse = "\n")
 }

output$pick.group <- renderUI({
    opts <- names(samtbl)
    selectInput("groupname", "Plot By", choices=opts,selected='SampleGroup')
})

get.ct <- function(var) {
    cttbl <- read_tsv(paste(data.dir,'countTable.logCPM.txt',sep='/'))
    if ("symsearch" %in% names(var)) {
       if (length(as.character(var$symsearch)) > 2) {
       cts <- filter(cttbl,SYMBOL == var$symsearch)
    }
    }
    if ("enssearch" %in% names(var)) {
       if (length(as.character(var$enssearch)) > 2) {
       	  cts <- filter(cttbl,ENSEMBL == var$enssearch)
    	  }
    }
    countTable <- gather(select(cts,4:length(cts)),key="SampleID",value="cts")
    newdf <- inner_join(countTable,samtbl)
    return(list(ctable=newdf))
  }

forct <- eventReactive(input$gcButton,{get.ct(input)})

output$plot.gene <- renderPlot({
    newdf <- forct()$ctable
    p1 <- ggplot(newdf, aes(x=!!as.symbol(input$groupname), y=cts,fill=!!as.symbol(input$groupname))) + geom_boxplot() + theme(legend.position="top",axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + labs(title=paste("Relative Abudance of",input$symsearch),x=input$groupname, y = "Relative Abundance (logCPM)",caption=wrapper("Boxplots draw to represent the 25th and 75th percentile (the lower and upper quartiles, respectively) as a box with a band in the box representing 50th percentile (the median). The upper whisker is located at the 'smaller' of the maximum x value and Q_3 + 1.5 inner quantile range(IRQ), whereas the lower whisker is located at the 'larger' of the smallest x value and Q_1 - 1.5 IQR")) + theme(legend.position = "bottom",plot.margin = margin(15, 15, 15, 15),plot.caption = element_text(size = 10, hjust = 0))
    p2 <- ggplot(newdf, aes(x=!!as.symbol(input$groupname), y=cts,fill=!!as.symbol(input$groupname))) + geom_violin() + theme(legend.position="top",axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + labs(title=paste("Relative Abudance of",input$symsearch),x=input$groupname, y = "Relative Abundance (logCPM)",caption=wrapper("Violin plot shows the kernel probability density of the data at different value. Violin plots include a marker for the median of the data and a box indicating the interquartile range")) + theme(legend.position = "bottom",plot.margin = margin(15, 15, 15, 15),plot.caption = element_text(size = 10, hjust = 0))
    grid.arrange(p1,p2,nrow=2)
}, height = 900, width = 900)

output$dlboxplot <- renderUI({
  if (!is.null(forct()$ctable)) {
    downloadButton('Downloadbp', 'Download PNG')
  }
})

output$Downloadbp <- downloadHandler(
  file = function() {
    paste('boxplot', 'png', sep = ".")
  },
  content <- function(file) {
    png(
      file,
      width = 8 * 120,
      height = 4 * 120,
      res = 120,
      pointsize = 5
    )
    if (is.null(forct()$ctable)) {
      return()
    }
    newdf <- forct()$ctable
    p1 <- ggplot(newdf, aes(x=!!as.symbol(input$groupname), y=cts,fill=!!as.symbol(input$groupname))) + geom_boxplot() + theme(legend.position="top",axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + labs(title=paste("Relative Abudance of",input$symsearch),x=input$groupname, y = "Relative Abundance (logCPM)",caption=wrapper("Boxplots draw to represent the 25th and 75th percentile (the lower and upper quartiles, respectively) as a box with a band in the box representing 50th percentile (the median). The upper whisker is located at the 'smaller' of the maximum x value and Q_3 + 1.5 inner quantile range(IRQ), whereas the lower whisker is located at the 'larger' of the smallest x value and Q_1 - 1.5 IQR")) + theme(legend.position = "bottom",plot.margin = margin(15, 15, 15, 15),plot.caption = element_text(size = 10, hjust = 0))
    p2 <- ggplot(newdf, aes(x=!!as.symbol(input$groupname), y=cts,fill=!!as.symbol(input$groupname))) + geom_violin() + theme(legend.position="top",axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + labs(title=paste("Relative Abudance of",input$symsearch),x=input$groupname, y = "Relative Abundance (logCPM)",caption=wrapper("Violin plot shows the kernel probability density of the data at different value. Violin plots include a marker for the median of the data and a box indicating the interquartile range")) + theme(legend.position = "bottom",plot.margin = margin(15, 15, 15, 15),plot.caption = element_text(size = 10, hjust = 0))
    grid.arrange(p1,p2,nrow=1)
    dev.off()
  },
  contentType = "image/png"
)

