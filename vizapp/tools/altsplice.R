
 get.bgdata <- function(var) {
   rdafile <- paste(data.dir,'bg.rda',sep='/')
   load(rdafile)
   genetrx <- indexes(bg)$t2g
   genetrx$SYMBOL <- geneNames(bg)

   geneid <- as.character(unique(genetrx[genetrx$SYMBOL %in% var$symsearch,]$g_id))
   if (exists("var$enssearch") && var$enssearch != '') {
       geneid <- as.character(unique(genetrx$g_id[grep(enssearch,genetrx$g_id)]))
     }
   genebg <- subset(bg,paste0("gene_id=='",geneid,"'"))
   rownames(genetrx) <- transcriptNames(bg)
   fpkm <- texpr(bg,meas='FPKM')
   keep <- genetrx[genetrx$g_id == geneid,]$t_id
   cttbl <- fpkm[keep,]
   grps <- pData(bg)$group
   trxnames <-  genetrx[genetrx$g_id == geneid,]
   dtm <- melt(t(cttbl))
   names(dtm) <- c('sample','transcript','value')
   dtm$grps <- rep(grps,length(keep))
   if (length(keep) < 2) {
     dtm$transcript <- rep(keep,length(keep))
   }
   dtm$grptrx <- paste(dtm$grps,dtm$transcript,sep='.')
   test <- stattest(gown=genebg, pData=pData(bg), feature='transcript',covariate='group', libadjust=FALSE,getFC=TRUE)
   if (length(keep) > 1) {
     agg = collapseTranscripts(gene=geneid, gown=bg, k=var$kct, method='kmeans')
     test <- stattest(gowntable=agg$tab, pData=pData(bg), feature='transcript_cluster',covariate='group', libadjust=FALSE,getFC=TRUE)
   }
   return(list(stattbl=test,gid=geneid,obj=genebg,cttbl=dtm,tname=trxnames))
 }

find_sym <- function (sym) {
if (!(toupper(sym) %in% toupper(symsyn$symbol))) {
     syns <- symsyn[grep(input$symsearch,symsyn$synonym,ignore.case=TRUE),]$symbol
     synlist <- paste(as.character(syns),collapse=',')
   if (length(syns) > 1) {
     paste("Please Use Official Gene Symbols",synlist,sep=':')
   }else {"Please Use Official Gene Symbols"}  	 
 }else {NULL}  	 
}	 

getgeneid <- eventReactive(input$altButton,{
 if (input$symsearch == '' & input$enssearch != '') {
    validate (find_sym(input$symsearch))
    }   
get.bgdata(input)
})

output$plot.cluster <- renderPlot({
  par(oma=c(4,4,1,1)) 
  gid <- getgeneid()$gid
  bg <- getgeneid()$obj
  tname <- getgeneid()$tname
  if (nrow(tname) > 1) {
    plotLatentTranscripts(gene=gid, gown=bg, k=input$kct, method='kmeans', returncluster=FALSE)
  }
})

output$gene.stat <- DT::renderDataTable({
  getgeneid()$stattbl
},escape=FALSE)

output$trx.name <- DT::renderDataTable({
  getgeneid()$tname
},escape=FALSE)


output$plot.means <- renderPlot({
  par(oma=c(4,4,1,1))
  par(cex.main=0.75)
  gid <- getgeneid()$gid
  bg <- getgeneid()$obj
  tname <- getgeneid()$tname
  if (nrow(tname) > 1) {
    plotMeans(gid, bg, groupvar='group', meas='cov', colorby='transcript')
  }
}, height = 900, width = 900)

output$trx.gene <- renderPlot({
  countTable <- getgeneid()$cttbl
  par(oma=c(4,4,1,1))
  p <- ggplot(countTable,aes(x=grptrx,y=log2(value+1))) + geom_boxplot(aes(fill = factor(grptrx))) + geom_jitter(height = 0) + theme(legend.position="left",axis.text.x=element_text(angle=45,hjust=1, vjust=1),legend.key.height=unit(0.5,"line"),legend.text=element_text(size=8),legend.title=element_blank()) + ylab("Relative Abundance (FPKM)") + xlab("")
  print(p)
})
