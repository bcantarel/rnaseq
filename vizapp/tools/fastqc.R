output$dir.fqc <- renderUI({
    datadir <- dir(data.dir)
    selectInput("dset.fqc", "Choose Dataset", choices=datadir,selected='coarse_cellpops')
})
output$pick.fqc <- renderUI({
    flist <- list.files(paste(data.dir,input$dset.fqc,"fastqc",sep='/'),pattern="*html$")
    flist <- gsub(".sort_fastqc.html",'',flist)
    selectInput("samp", "Choose Sample", choices=flist)
})

getpage <- function(var) {
	return(includeHTML(paste(data.dir,"/",var$dset.fqc,"/fastqc/",var$samp,".sort_fastqc.html",sep='')))
}
output$inc <- renderUI({getpage(input)})