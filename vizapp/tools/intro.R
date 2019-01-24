# x <- reactiveVal(1)
observeEvent(input$actionbtgc, {
  #x(x() + 1)
  # if (input$actionbtgc > 0){
  # updateNavbarPage(session, "rnaseqnav", selected = "panelgc")
  #   output$Text <-renderText({"lailailailailaillai"})
  #}
  updateNavbarPage(session, 'nav_shinyapp', 'panelgc')
})

observeEvent(input$actionbtdea, {
  updateNavbarPage(session, 'nav_shinyapp', 'paneldea')
})

observeEvent(input$actionbtgsea, {
  updateNavbarPage(session, 'nav_shinyapp', 'panelgsea')
})
# output$tb1 <- renderPrint({x()})
#  output$tb1 <- renderText({
#   updateNavbarPage(session, '')
# })
# output$ui_intro <-renderUI({
#  source('tools/gc.R', local = TRUE)
#})
