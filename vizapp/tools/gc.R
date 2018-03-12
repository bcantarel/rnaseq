
samfile <-  paste(data.dir,'design.shiny.txt',sep='/')
samtbl <- read.table(samfile,header=TRUE,sep="\t")

output$pick.group <- renderUI({
    opts <- names(samtbl)
    selectInput("groupname", "Plot By", choices=opts,selected='SampleGroup')
})

get.ct <- function(var) {
    ctfile <-  paste(data.dir,'countTable.logCPM.txt',sep='/')
    statfile <- paste(data.dir,'edgeR.results.txt',sep='/')
    wp.df <- data.frame(Column1 = c(0), Column2 = c(0))
    
    if (nchar(as.vector(var$symsearch)) > 2) {
      cts <- read.csv.sql(ctfile, sql = paste("select * from file where symbol ='",var$symsearch,"'",sep=''),sep = "\t")
      if (file.exists(statfile)) {
        wp.df <- read.csv.sql(statfile, sql = paste("select * from file where symbol ='",var$symsearch,"'",sep=''),sep = "\t")
      }
    }
    else {
      if (nchar(as.vector(var$enssearch)) > 2) {
        cts <- read.csv.sql(ctfile, sql = paste("select * from file where ENSEMBL ='",var$enssearch,"'",sep=''),sep = "\t")
        if (file.exists(statfile)) {
          wp.df <- read.csv.sql(statfile, sql = paste("select * from file where ENSEMBL ='",var$enssearch,"'",sep=''),sep = "\t")
        }
      }
    }
    countTable <- cts[,4:length(cts)]
    samples <- colnames(countTable)
    mergetbl <- merge(as.data.frame(samples),samtbl,by.x="samples",by.y="SampleID",all.x=TRUE,sort=FALSE)
    grps <- mergetbl[,var$groupname]
    newdf <- data.frame(ct=t(countTable),grp=as.character(grps))
    names(newdf) = c('cts','grp')
    return(list(ctable=newdf,tbl=cts,design=samtbl,stats=wp.df))
  }

forct <- eventReactive(input$gcButton,{get.ct(input)})

output$plot.gene <- renderPlot({
    countTable <- forct()$ctable
    par(oma=c(4,4,1,1)) 
    grpnames <- levels(factor(as.character(countTable$grp)))
    boxplot(countTable$cts ~ countTable$grp,col=rainbow(length(grpnames)),cex.axis=0.7,ylab='Relative Abundance (logCPM)',las=2,main=input$groupname)
    })

output$bxplot.desc <- renderText({ 
      paste("Relative Abudance of",input$symsearch,"calculated by Log2(Counts Per Million Reads).  Boxplots draw to represent the 25th and 75th percentile (the lower and upper quartiles, respectively) as a box with a band in the box representing 50th percentile (the median). The upper whisker is located at the 'smaller' of the maximum x value and Q_3 + 1.5 inner quantile range(IRQ), whereas the lower whisker is located at the 'larger' of the smallest x value and Q_1 - 1.5 IQR",sep=' ')
    })

output$violin.gene <- renderPlot({
  countTable <- forct()$ctable
  par(oma=c(4,4,1,1))
  p <- ggplot(countTable,aes(x=grp,y=cts)) + geom_violin(trim=FALSE,aes(fill = factor(grp))) + geom_jitter(height = 0) + theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + ylab("Relative Abundance (logCPM)") + xlab("")
  print(p)
}, height = 400, width = 900)

output$violin.desc <- renderText({ 
      paste("Relative Abudance of",input$symsearch,input$enssearch,"calculated by Log2(Counts Per Million Reads).  Violin plot is similar to box plots above, except that it also show the kernel probability density of the data at different value. Violin plots include a marker for the median of the data and a box indicating the interquartile range, as in boxplot above.",sep=' ')
    })
