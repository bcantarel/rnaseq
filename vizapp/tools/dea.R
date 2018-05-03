# UI-elements for DEA

ctfile <-  paste(data.dir, 'countTable.logCPM.txt', sep = '/')
samfile <-  paste(data.dir, 'design.shiny.txt', sep = '/')
cts <- read.table(file = ctfile, header = TRUE, sep = '\t')
samtbl <- read.table(samfile, header = TRUE, sep = "\t")
samples <- colnames(cts[, 4:length(cts)])
mergetbl <-
  merge(
    as.data.frame(samples),
    samtbl,
    by.x = "samples",
    by.y = "SampleID",
    all.x = TRUE,
    sort = FALSE
  )
grps <- mergetbl$SampleGroup
grpnames <- levels(factor(grps))
col.blocks <- col.grp(grps, grpnames)
MSIG.geneSets <-
  read.gmt(paste(data.dir, 'geneset.shiny.gmt', sep = '/'))

output$pick.dea <- renderUI({
  flist <- list.files(data.dir, pattern = "*edgeR.txt$")
  selectInput("file", "Choose Pair", choices=flist, width = "100%")
})
output$pick.pathway <- renderUI({
  pathways <- names(MSIG.geneSets)
  pathchoices = setNames(1:length(pathways), pathways)
  selectInput("deapathname", "Choose Pair", choices = pathchoices)
})

get.data <- function(var) {
  f <- paste(data.dir, var$file, sep = '/')
 # f <-paste(f, "edgeR.txt", sep = ".") #This line needs to comment out if input is fullname
  f1 <- paste(data.dir, input$file, sep = '/')
  hmcomp <-paste("Heatmap of all genes with an FDR < 0.05 using EdgeR Results", sep ='')
  top50 <- paste("Top 50 User Defined Genes")
  allgene <- paste("All Differentially Expressed Genes")
  comp <- read.table(f, header = TRUE, sep = '\t')
  comp.filt <-
    na.omit(comp[abs(comp$logFC) >= var$fc.thresh &
                   comp$rawP <= var$pval.thresh,])
  if (var$adjust == 'FDR') {
    comp.filt <-
      na.omit(comp[abs(comp$logFC) >= var$fc.thresh &
                     comp$fdr <= var$pval.thresh,])
  }
  if (var$adjust == 'BONF') {
    comp.filt <-
      na.omit(comp[abs(comp$logFC) >= var$fc.thresh &
                     comp$bonf <= var$pval.thresh,])
  }
  genelist <-
    as.character(head(comp.filt[order(comp.filt$fdr),]$symbol, n = var$numgenes))
  if (var$heatmap == 'hgeneset') {
    genelist <- unlist(MSIG.geneSets[as.numeric(var$deapathname)])
  }
  if (var$heatmap == 'cgeneset') {
    genelist <- unlist(strsplit(as.character(var$genes), "[;\n]+"))
  }
  return(
    list(
      filt = comp.filt,
      glist = genelist,
      f1 = f1,
      hmcomp = hmcomp,
      top50 = top50,
      allgene = allgene
    )
  )
}

tbls <-  eventReactive(input$deButton,
                       {
                         get.data(input)
                       })

output$selectgenes <- renderUI({
  symnames <- tbls()$glist
  textAreaInput(
    "genes",
    "Gene Symbols separated by ';'",
    value = paste(symnames, collapse = ";"),
    width = '95%',
    rows = 10
  )
})

output$dge.c <- DT::renderDataTable({
  t1 <- tbls()$filt
  t1$symbol <-
    paste(
      "<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",
      t1$symbol,
      '>',
      t1$symbol,
      "</a>",
      sep = ''
    )
  t1$ensembl <-
    paste(
      "<a href=http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
      t1$ensembl,
      '>',
      t1$ensembl,
      "</a>",
      sep = ''
    )
  t1
}, escape = FALSE, filter = 'top', options = list(
  lengthMenu = list(c(10, 25, 50, 100, -1), list('10', '25', '50', '100', 'All')),
  autoWidth = TRUE,
  columnDefs = list(list(width = '5%', targets = '0')),
  scrollX = TRUE
))

output$downloadC <- renderUI({
  if (!is.null(tbls()$filt)) {
    downloadButton('OutputFile', 'Download CSV')
  }
})
output$OutputFile <- downloadHandler(
  file <- paste(input$file, ".filt.txt", sep = ""),
  content = function(file) {
    write.table(
      tbls()$filt,
      file,
      quote = FALSE,
      row.names = FALSE,
      sep = '\t'
    )
  }
)
plotHeatmap <- reactive({
  syms <- tbls()$glist
  ct2 <- cts[cts$SYMBOL %in% syms,]
  subset <- ct2[, 4:length(ct2)]
  row.names(subset) <- ct2$SYMBOL
  STREE <- hclust(dist(t(subset)))
  zscores <- scale(t(subset))
  ngenes <- length(colnames(zscores))
  textscale <- (1 / (ngenes/30))
  if (textscale > 1) {
    textscale <- 12
  }
  if (textscale < 0.1) {
    textscale <- 12
  }
  if (input$cluster == 2) {
    heatmap.2(
      zscores,
      col = bluered(100),
      #Rowv = as.dendrogram(STREE),
      Rowv = FALSE,
      Colv = FALSE,
      lwid = c(0.5,3),
      RowSideColors = col.blocks,
      dendrogram = 'row',
      cexCol = 1.1,
      srtRow = 45,
      srtCol = 45,
      trace = "none",
      margins = c(8, 7)
    )
  } else{
    heatmap.2(
      zscores,
      col = bluered(100),
      Rowv = as.dendrogram(STREE),
      RowSideColors = col.blocks,
      dendrogram = 'row',
      cexCol = 1.1,
      srtRow = 45,
      srtCol = 45,
      lwid = c(0.5,3),
      trace = "none",
      margins = c(8, 7)
    )
  }
  legend("top",
         legend = grpnames,
         col = rainbow(length(grpnames)),
         pch = 20)
})
output$plot.heatmap <- renderPlot({
  plotHeatmap()
})

output$dlheatmap <- renderUI({
  if (!is.null(tbls()$glist)) {
    downloadButton('Downloadhp', 'Download PNG')
  }
})

output$Downloadhp <- downloadHandler(
  file = function() {
    paste('mean', 'png', sep = ".")
  },
  content <- function(file) {
    png(
      file,
      width = 8 * 140,
      height = 4 * 100,
      res = 140,
      pointsize = 5
    )
    if (is.null(tbls()$glist)) {
      return()
    }
    syms <- tbls()$glist
    ct2 <- cts[cts$SYMBOL %in% syms,]
    subset <- ct2[, 4:length(ct2)]
    row.names(subset) <- ct2$SYMBOL
    STREE <- hclust(dist(t(subset)))
    zscores <- scale(t(subset))
    ngenes <- length(colnames(zscores))
    textscale <- (1 / (ngenes/30))
    if (textscale > 1) {
      textscale <- 12
    }
    if (textscale < 0.1) {
      textscale <- 12
    }
    if (input$cluster == 2) {
      plot <- heatmap.2(
        zscores,
        col = bluered(100),
        Rowv = FALSE,
        Colv = FALSE,
        lwid = c(0.5,3),
        RowSideColors = col.blocks,
        dendrogram = 'row',
        cexCol = 1.1,
        srtRow = 45,
        srtCol = 45,
        trace = "none",
        margins = c(8, 7)
      )
    } else{
     plot <- heatmap.2(
        zscores,
        col = bluered(100),
        Rowv = as.dendrogram(STREE),
        RowSideColors = col.blocks,
        dendrogram = 'row',
        cexCol = 1.1,
        srtRow = 45,
        srtCol = 45,
        lwid = c(0.5,3),
        trace = "none",
        margins = c(8, 7)
      )
    }
    legend("top",
           legend = grpnames,
           col = rainbow(length(grpnames)),
           pch = 20)
    print(plot)
    dev.off()
  },
  contentType = "image/png"
)

output$downloadpdf = downloadHandler(
  filename = "output.heatmap.pdf",
  content = function(file) {
    pdf(file = file, paper = "letter")
    plotHeatmap()
    dev.off()
  }
)

output$hm.comp <- renderImage({
  f1 <- paste(tbls()$f1)
  png <- gsub('edgeR.txt', 'heatmap.edgeR.png', f1)
  list(
    src = png,
    alt = paste("HeatMap Comparison"),
    width = "100%",
    height = "auto"
  )
}, deleteFile = FALSE)

output$dlhmcomp <- renderUI({
  if (!is.null(tbls()$f1)) {
    downloadButton('Downloadhm', 'Download PNG')
  }
})
output$Downloadhm <- downloadHandler(
  file = function() {
    paste('heatmap_comp', 'png', sep = ".")
  },
  content <- function(file) {
    f1 <- paste(tbls()$f1)
    png <- gsub('edgeR.txt', 'heatmap.edgeR.png', f1)
    file.copy(png, file)
  }
)

output$hmcomp.allgene <- renderText({
  paste(tbls()$allgene)
})
output$hmcomp.top50 <- renderText({
  paste(tbls()$top50)
})
output$hmcomp.desc <- renderText({
  paste(tbls()$hmcomp)
})
