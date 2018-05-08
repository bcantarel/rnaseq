library(shiny)
library(Vennerable)
library(qusage)
library(DT)
library(ggplot2)
library(ballgown)
library(sqldf)
library(reshape2)
library("gplots")

shinyServer(function(input, output, session) {

    data.dir <- Sys.getenv('outputDir')
    rda.dir <- data.dir
    symsyn <- read.table(file='symbol2synonym.txt',header=FALSE)
    names(symsyn) <- c('symbol','synonyms')
    source('functions.R', local = TRUE)
    source('tools/intro.R', local = TRUE)
    source('tools/intro_ui.R', local = TRUE)
    source('tools/gc.R', local = TRUE)
    source('tools/qc.R', local = TRUE)
    source('tools/altsplice.R', local = TRUE)
    source('tools/dea.R', local = TRUE)
    source('tools/gsea.R', local = TRUE)
    #source('tools/fastqc.R', local = TRUE)
    source('tools/gc_ui.R', local = TRUE)
    source('tools/qc_ui.R', local = TRUE)
    source('tools/dea_ui.R', local = TRUE)
    source('tools/gsea_ui.R', local = TRUE)
    #source('tools/fastqc_ui.R', local = TRUE)
    source('tools/altsplice_ui.R', local = TRUE)

})

