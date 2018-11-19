get.qc <- function(var) {
  ctfile <-  paste(data.dir, 'countTable.stats.txt', sep = '/')
  alfile <-  paste(data.dir, 'alignment.summary.txt', sep = '/')
  stats <- read.table(file = ctfile,header = TRUE,sep = '\t')
  align <- read.table(file = alfile,header = TRUE,sep = '\t')
  fsdes <- paste(
      "Percent of Total Reads for each gene type.  Boxplots draw to represent the 25th and 75th percentile (the lower and upper quartiles, respectively) as a box with a band in the box representing 50th percentile (the median). The upper whisker is located at the 'smaller' of the maximum x value and Q_3 + 1.5 inner quantile range(IRQ), whereas the lower whisker is located at the 'larger' of the smallest x value and Q_1 - 1.5 IQR",
      sep = ' '
    )
  pca <- paste(data.dir, 'pca.png', sep = "/")
  mds <- paste(data.dir, 'mds.png', sep = "/")
  sheat <- paste(data.dir, 'samples_heatmap.png', sep = "/")
  return(
    list(
      featurestat = stats,
      alignstat = align,
      fsdes = fsdes,
      pca = pca,
      mds = mds,
      sheat = sheat
    )
  )
}

forstat <- eventReactive(input$qcButton, {get.qc(input)})

output$plot.fs <- renderPlot({
  countTable <- forstat()$featurestat
  grpnames <- levels(factor(as.character(countTable$Type)))
  boxplot(
    countTable$ReadPerc ~ countTable$Type,
    col = rainbow(length(grpnames)),
    cex.axis = 0.7,
    ylab = 'Percent of Total Reads',
    las = 2,
    par(mar = c(12, 5, 2, 2))
  )
}, height = "auto", width = 'auto')

#Download Button for plot.fs
output$dlplot <- renderUI({
  if (!is.null(forstat()$featurestat)) {
    downloadButton('Downloadplot', 'Download PNG')
  }
})
output$Downloadplot <- downloadHandler(
  file = function() {
    paste('plot', 'png', sep = ".")
  },
  content <- function(file) {
    png(
      file,
      width = 8 * 120,
      height = 5 * 120,
      res = 200,
      pointsize = 6
    )
    if (is.null(forstat()$featurestat)) {
      return()
    }
    countTable <- forstat()$featurestat
    par(oma = c(4, 4, 1, 1))
    grpnames <- levels(factor(as.character(countTable$Type)))
    plot <- boxplot(
      countTable$ReadPerc ~ countTable$Type,
      col = rainbow(length(grpnames)),
      cex.axis = 0.7,
      ylab = 'Percent of Total Reads',
      las = 2,
      par(mar = c(12, 5, 2, 2))
    )
    print(plot)
    dev.off()
  },
  contentType = "image/png"
)

#display description for the plot
output$fs.desc <- renderText({
  paste(forstat()$fsdes)
})

output$feature.stats <- DT::renderDataTable({
  t1 <- forstat()$featurestat
  DT::datatable(
    t1,
    escape = FALSE,
    filter = 'top',
    options = list(lengthMenu = list(
      c(10, 25, 50, 100,-1), list('10', '25', '50', '100', 'All')
    ))
  )
})

output$align.stats <- DT::renderDataTable({
  t1 <- forstat()$alignstat
  DT::datatable(
    t1,
    escape = FALSE,
    filter = 'top',
    options = list(lengthMenu = list(
      c(10, 25, 50, 100,-1), list('10', '25', '50', '100', 'All')
    ))
  )
})
output$downloadAlign <- renderUI({
  if (!is.null(forstat()$alignstat)) {
    downloadButton('DownloadAlign', 'Download CSV')
  }
})
output$DownloadAlign <- downloadHandler(
  file <- paste('alignment.summary.txt'),
  content = function(file) {
    write.table(
      forstat()$alignstat,
      file,
      quote = FALSE,
      row.names = FALSE,
      sep = '\t'
    )
  }
)

output$downloadGT <- renderUI({
  if (!is.null(forstat()$featurestat)) {
    downloadButton('DownloadGT', 'Download CSV')
  }
})
output$DownloadGT <- downloadHandler(
  file <- paste('countTable.stats.txt'),
  content = function(file) {
    write.table(
      forstat()$featurestat,
      file,
      quote = FALSE,
      row.names = FALSE,
      sep = '\t'
    )
  }
)
output$pca1 <- renderImage({
  list(
    src = forstat()$pca,
    alt = paste("Sample Comparison 1"),
    width = "100%",
    height = "100%"
  )
}, deleteFile = FALSE)

output$dlpca <- renderUI({
  if (!is.null(forstat()$pca)) {
    downloadButton('Downloadpca', 'Download PNG')
  }
})
output$Downloadpca <- downloadHandler(
  file = function() {
    paste('pca_analysis', 'png', sep = ".")
  },
  content <- function(file) {
    file.copy(forstat()$pca, file)
  }
)

output$pca2 <- renderImage({
  list(
    src = forstat()$mds,
    alt = paste("Sample Comparison 2"),
    width = "100%",
    height = "auto"
  )
}, deleteFile = FALSE)

output$dlmds <- renderUI({
  if (!is.null(forstat()$mds)) {
    downloadButton('Downloadmds', 'Download PNG')
  }
})
output$Downloadmds <- downloadHandler(
  file = function() {
    paste('mds', 'png', sep = ".")
  },
  content <- function(file) {
    file.copy(forstat()$mds, file)
  }
)
output$hmap <- renderImage({
  list(
    src = forstat()$sheat,
    alt = paste("Sample Comparison 3"),
    width = "100%",
    height = "auto"
  )
}, deleteFile = FALSE)

output$dlhmap <- renderUI({
  if (!is.null(forstat()$sheat)) {
    downloadButton('Downloadhmap', 'Download PNG')
  }
})
output$Downloadhmap <- downloadHandler(
  file = function() {
    paste('sample_distance', 'png', sep = ".")
  },
  content <- function(file) {
    file.copy(forstat()$sheat, file)
  }
)