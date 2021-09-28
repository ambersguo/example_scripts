library(ggplot2)
library(dplyr)
library(tidyverse)
library(plotly)
library(shiny)

set.seed(100)

mydata<-read.table("Pathak_HeLa_CPSF6_IS_larger_than_5_and_RNAseq_nonzero_binned.txt",head=T)
p<-mydata %>% tidyr::gather(sample, IS_count, 1:6) %>% 
  ggplot(., aes(HeLa_RNAseq_log10, IS_count))+
  geom_point(aes(color=sample),alpha=0.2)

fig <- ggplotly(p)

ui <- fluidPage(
  plotlyOutput("plot"),
  verbatimTextOutput("hover"),
  verbatimTextOutput("click"),
  verbatimTextOutput("brush"),
  verbatimTextOutput("zoom")
  
)

server <- function(input, output, session) {
  
  output$plot <- renderPlotly({
    p<-mydata %>% tidyr::gather(sample, IS_count, 1:6) %>% 
      ggplot(., aes(HeLa_RNAseq_log10, IS_count))+
      geom_point(aes(color=sample),alpha=0.2)
    
    fig <- ggplotly(p)
    fig
  })
  
  output$hover <- renderPrint({
    d <- event_data("plotly_hover")
    if (is.null(d)) "Hover events appear here (unhover to clear)" else d
  })
  
  output$click <- renderPrint({
    d <- event_data("plotly_click")
    if (is.null(d)) "Click events appear here (double-click to clear)" else d
  })
  
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
  })
  
  output$zoom <- renderPrint({
    d <- event_data("plotly_relayout")
    if (is.null(d)) "Relayout (i.e., zoom) events appear here" else d
  })
  
}

shinyApp(ui, server)
