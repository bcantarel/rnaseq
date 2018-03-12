# UI-elements for DEA

ctfile <-  paste(data.dir,'countTable.logCPM.txt',sep='/')
samfile <-  paste(data.dir,'design.shiny.txt',sep='/')
cts <- read.table(file=ctfile,header=TRUE,sep='\t')
samtbl <- read.table(samfile,header=TRUE,sep="\t")
samples <- colnames(cts[,4:length(cts)])
mergetbl <- merge(as.data.frame(samples),samtbl,by.x="samples",by.y="SampleID",all.x=TRUE,sort=FALSE)
grps <- mergetbl$SampleGroup
grpnames <- levels(factor(grps))
col.blocks <-col.grp(grps,grpnames)
MSIG.geneSets <- read.gmt(paste(data.dir,'geneset.shiny.gmt',sep='/'))

output$pick.dea <- renderUI({
    flist <- list.files(data.dir,pattern="*edgeR.txt$")
    selectInput("file", "Choose Pair", choices=flist)
})
output$pick.pathway <- renderUI({
    pathways <- names(MSIG.geneSets)
    pathchoices = setNames(1:length(pathways),pathways)
    selectInput("deapathname", "Choose Pair", choices=pathchoices)
})

get.data <- function(var) {
    f <- paste(data.dir,var$file,sep='/')
    comp <- read.table(f,header=TRUE,sep='\t')
    comp.filt <- na.omit(comp[abs(comp$logFC) >= var$fc.thresh & comp$rawP <= var$pval.thresh,])
    if (var$adjust == 'FDR') {
        comp.filt <- na.omit(comp[abs(comp$logFC) >= var$fc.thresh & comp$fdr <= var$pval.thresh,])
    }
    if (var$adjust == 'BONF') {
        comp.filt <- na.omit(comp[abs(comp$logFC) >= var$fc.thresh & comp$bonf <= var$pval.thresh,])
    }
    genelist <- as.character(head(comp.filt[order(comp.filt$fdr),]$symbol,n=var$numgenes))
    if (var$heatmap == 'hgeneset') {
        genelist <- unlist(MSIG.geneSets[as.numeric(var$deapathname)])
    }
   if (var$heatmap == 'cgeneset') {
         genelist <- unlist(strsplit(as.character(var$genes), "[;\n]+"))
   }
   return(list(filt=comp.filt,glist=genelist))
}

tbls <-  eventReactive(input$deButton,{get.data(input)})

output$selectgenes <- renderUI({
   symnames <- tbls()$glist
   textAreaInput("genes", "Gene Symbols separated by ';'",value=paste(symnames,collapse=";"), width = '95%',rows=10)
})

output$dge.c <- DT::renderDataTable({
    t1 <- tbls()$filt
    t1$symbol <- paste("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",t1$symbol,'>',t1$symbol,"</a>",sep='')
    t1$ensembl <- paste("<a href=http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",t1$ensembl,'>',t1$ensembl,"</a>",sep='')
    t1
},escape=FALSE,filter = 'top',options = list(lengthMenu = c(10, 25, 50, 200, -1)))

output$downloadC <- downloadHandler(
                                     file <- paste(input$file,".filt.txt",sep=""),
                                     content = function(file) {
                                       write.table(tbls()$filt,file,quote=FALSE,row.names=FALSE,sep='\t')
                                     })

plotHeatmap <- reactive({
	syms <- tbls()$glist
	ct2 <- cts[cts$SYMBOL %in% syms,]
	subset <- ct2[,4:length(ct2)]
	row.names(subset) <- ct2$SYMBOL
	STREE <- hclust(dist(t(subset)))
      	zscores <- scale(t(subset))
      	ngenes <- length(colnames(zscores))
      	textscale <- (1/(ngenes/30))
      	if (textscale > 1) {
           textscale <-1
      	}
      	if (textscale < 0.1) {
           textscale <- 0.1
      	}
	heatmap.2(zscores, col = bluered(100),Rowv = as.dendrogram(STREE), RowSideColors = col.blocks,dendrogram='row', cexCol=textscale,srtRow=45,srtCol=45,trace="none",margins=c(8,16))
	legend("topright",legend=grpnames,col=rainbow(length(grpnames)),pch=20)
})
output$plot.heatmap <- renderPlot({
	plotHeatmap()
})
output$downloadpdf = downloadHandler(
      filename = "output.heatmap.pdf",
      content = function(file) {
      pdf(file = file,paper="letter")
      plotHeatmap()
      dev.off()
      })

output$hm.comp <- renderImage({
  f1 <- paste(data.dir,input$file,sep='/')
  png <- gsub('edgeR.txt','heatmap.edgeR.png',f1)
  list(src=png, alt=paste("HeatMap Comparison"))
},deleteFile=FALSE)

output$hmcomp.desc <- renderText({ 
      paste("Heatmap of all genes with an FDR < 0.05 using EdgeR Results",sep='')
    })
