library(shiny)
library(datasets)
library(ggplot2)
library(stringr)
library(plyr)
library(tidyverse)
library(reshape2)
library(ggrepel)
library(periscope)
library(shinybusy)
options(shiny.maxRequestSize=100*1024^2) 
#setwd(gsub("shiny.R", "",rstudioapi::getActiveDocumentContext()$path))

ui <-shinyUI(fluidPage(# 1
  titlePanel("Multiple phenotypes Manhattan plots"),
  # First screen
  tabsetPanel(#2
    tabPanel("Two-way Manhatan", #3
             titlePanel("Two-way Manhatan"),
             sidebarLayout(  
               sidebarPanel( 
                 fileInput('file1', 'Choose CSV File',
                           accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')), #5
                 tags$br(),
                 checkboxInput('header', 'Header', TRUE),
                 radioButtons('sep', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
                 radioButtons('quote', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"), '"'),
                 selectInput('marker_ID', 'Marker_ID', ""),
                 selectInput('posi', 'Posi', ""),
                 selectInput('pvalue', 'P value', "", selected = ""),
                 selectInput('chromosome', 'Chromosome', "", selected = ""),
                 checkboxInput('highlight', 'Highlight', TRUE),
                 selectInput('trait1', 'Trait 1', "", selected = ""),
                 selectInput('trait2', 'Trait 2', "", selected = ""),
                 numericInput("obs", "Threshold:", 1),
                 numericInput("ylim", "ylim:", 1),
                 numericInput("point1", "Point size:", 0.1),
                 sliderInput("aspect.ratio1", ("aspect.ratio"),
                             min = 0.1, max = 2, value = .5),
                 sliderInput("scale1", ("Scale"),
                             min = 0.5, max = 3, value = 1),
                 textInput("xlab","Xlab"," "),
                 textInput("ylab","Ylab","p value"),
                 downloadButton("downloadData", "Download example data")
                 
               ),
               mainPanel(plotOutput('MyPlot'), add_busy_spinner(spin = "fading-circle"),
                         downloadButton('downloadPlot','Download Plot'))
             )),
    # Second screen
    tabPanel("PheWAS plot",
             pageWithSidebar(
               headerPanel('PheWAS plot'),
               sidebarPanel(
                 fileInput('file2', 'Choose CSV File',
                           accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')), #5
                 tags$br(),
                 checkboxInput('header', 'Header', TRUE),
                 radioButtons('sep', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
                 radioButtons('quote', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"), '"'),
                 selectInput('marker_ID2', 'Marker_ID', ""),
                 selectInput('pheno_cor', 'Type', ""),
                 selectInput('pheno_name', 'Pheno', ""),
                 selectInput('pvalue2', 'P valxue', "", selected = ""),
                 checkboxInput('highlight2', 'Highlight', TRUE),
                 numericInput("obs2", "Threshold:", 1),
                 numericInput("ylim2", "ylim:", 1),
                 numericInput("point2", "Point size:", 0.1),
                 sliderInput("aspect.ratio2", ("aspect.ratio"),
                             min = 0.1, max = 2, value = .5),
                 sliderInput("scale2", ("Scale"),
                             min = 0.5, max = 3, value = 1),
                 numericInput("ncols", ("Ncols"),5),
                 textInput("xlab2","Xlab"," "),
                 textInput("ylab2","Ylab","p value"),
                 downloadButton("downloadData2", "Download example Phewas data")
                 
               ),
               mainPanel(plotOutput('PheWAS_plot'),add_busy_spinner(spin = "fading-circle"),
                         downloadButton('downloadPlot2','Download Plot'))
             
             ))
  )#3
)#2
)#1



server <- shinyServer(function(input, output, session) {
  data <- reactive({ 
    req(input$file1) 
    inFile <- input$file1 
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    
    updateSelectInput(session, inputId = 'posi', label = 'Marker position',
                      choices = names(df), selected = names(df)[8])
    updateSelectInput(session, inputId = 'pvalue', label = 'p value',
                      choices = names(df), selected = names(df)[5])
    updateSelectInput(session, inputId = 'chromosome', label = 'Chromosome',
                      choices = names(df), selected = names(df)[7])
    updateSelectInput(session, inputId = 'trait1', label = 'Trait 1',
                      choices = names(df), selected = names(df)[1])
    updateSelectInput(session, inputId = 'trait2', label = 'Trait 2',
                      choices = names(df), selected = names(df)[6])
    updateSelectInput(session, inputId = 'marker_ID', label = 'Marker_ID',
                      choices = names(df), selected = names(df)[2])
    return(df)})
  data2 <- reactive({ 
    req(input$file2) 
    inFile2 <- input$file2 
    df2 <- read.csv(inFile2$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    updateSelectInput(session, inputId = 'pheno_cor', label = 'Phenotype group',
                      choices = names(df2), selected = names(df2)[9])
    updateSelectInput(session, inputId = 'pheno_name', label = 'Phenotype ID',
                      choices = names(df2), selected = names(df2)[1])
    updateSelectInput(session, inputId = 'pvalue2', label = 'p value',
                      choices = names(df2), selected = names(df2)[5])
    updateSelectInput(session, inputId = 'marker_ID2', label = 'Marker_ID',
                      choices = names(df2), selected = names(df2)[2])
    return(df2)})
# Tab 1
  output$MyPlot <- renderPlot({
    dados <- data()[, c(input$posi, input$pvalue, input$chromosome,input$trait1, input$trait2, input$marker_ID)]
    colnames(dados)<-c("posi", "pvalue", "chromosome", "trait1", "trait2", "Marker_ID")
    chromossos<-unique(dados$chromosome)
    chromossos_size=c()
    for (i in 1:length(chromossos)){chromossos_size[i]=max(dados$posi[dados$chromosome==chromossos[i]])}
    for (i in 2:length(chromossos)){chromossos_size[i]<-chromossos_size[i]+chromossos_size[i-1]}
    chrom_max_size<-c()
    for (i in 1:length(chromossos)){
      if (i == 1){chrom_max_size[i]<-max(dados$posi[dados$chromosome==chromossos[i]])/2 } else {
        chrom_max_size[i] <- chromossos_size[i-1] + max(dados$posi[dados$chromosome==chromossos[i]])/2
      }}
    
    for (i in 2:length(chromossos)){
      dados$posi[dados$chromosome==chromossos[i]]<-dados$posi[dados$chrom==chromossos[i]]+chromossos_size[i-1]
    }
    data_1 <- reactiveValues()
    
    a <-ggplot(dados, aes(x = posi, y=pvalue, colour = chromosome)) +
      geom_point(size = input$point1)+facet_grid(as.formula(paste("trait1~ trait2")))+
      theme_bw()  + ylab(input$ylab)+ xlab(input$xlab)+ ylim(0, input$ylim)+
      scale_x_continuous(breaks=chrom_max_size,labels=chromossos) + 
      theme(legend.position = "none",axis.text.x = element_text(angle = 90, size= 7, hjust = 0), aspect.ratio = input$aspect.ratio1,
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if(input$highlight == TRUE){data_1$plot <- a + geom_hline(yintercept = input$obs, color="red", alpha = 0.75, size = 0.5, linetype = 2)+ 
                                geom_text_repel(aes(label=ifelse(pvalue>input$obs,as.character(Marker_ID),'')), size = 1.5, max.overlaps = 40) # alternative 2
      }
    else{
      data_1$plot <- a
    }
    output$downloadPlot <- downloadHandler(
          filename = function(){paste("Manhattan Plot",'.pdf',sep='')},
          content = function(file){
    ggsave(file,plot= data_1$plot, device = "pdf", width = 7*input$scale1,height = 7*input$scale1,bg = "white", dpi = 300)
   
    output$downloadData <- downloadHandler(
      filename = function() {paste("example_manhatan.csv")},
      content = function(file) {
        write.csv(datasetInput(), file, row.names = FALSE)
      }
    )
  
    
      })
    
    data_1$plot
      })
  
  output$downloadData <- downloadHandler(
    filename = function() {paste("example_manhatan.csv")},
    content = function(file) {write.csv(read.csv("teste_Manhatan.csv"), file, row.names = FALSE)})
  
# Tab 2
  output$PheWAS_plot <- renderPlot({
    dados2 <- data2()[, c(input$pheno_cor, input$pheno_name, input$pvalue2,input$marker_ID2)]
    colnames(dados2)<-c("pheno_cor", "pheno_name", "pvalue2", "marker_ID2")
    data_2 <- reactiveValues()
    dados2$Index<-paste(dados2$pheno_cor,dados2$pheno_name, sep = "_")
    dados2<-dados2[order(dados2$Index),]
    dados2$Index <- as.numeric(as.factor(dados2$Index))
    pheno_cor <- unique(dados2$pheno_cor) ; breaks <- c()
    for (i in 1:length(pheno_cor)){breaks[i]<- median(dados2[dados2$pheno_cor==pheno_cor[i], "Index"])}
    temp_graph<- ggplot(dados2, aes(x = Index, y=pvalue2, colour = pheno_cor)) +
        geom_point(size = input$point2)+ theme_bw() +facet_wrap(~marker_ID2, ncol =  as.integer(input$ncols))+
        ylab(input$ylab2)+ xlab(input$xlab2)+
        scale_x_continuous(breaks=breaks,labels=unique(dados2$pheno_cor)) +
        theme(aspect.ratio = input$aspect.ratio2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position ="none" )
    
    if(input$highlight2){ data_2$plot <- temp_graph + geom_hline(yintercept = input$obs2, color="red",alpha = 0.75, size = 0.5, linetype = 2)+ 
      geom_text_repel(aes(label=ifelse(pvalue2>input$obs2,as.character(pheno_name),'')), size = 1.5, max.overlaps = 40)

      }
    else{data_2$plot <- temp_graph}
    
    
    output$downloadPlot2 <- downloadHandler(
      filename = function(){paste("PHEWAS_Plot",'.pdf',sep='')},
      content = function(file){ggsave(file,plot= data_2$plot, device = "pdf", width = 7*input$scale2,height = 7*input$scale2, dpi = 300)}) 
    data_2$plot
          })

  output$downloadData2 <- downloadHandler(
    filename = function() {paste("PheWAS_test.csv")},
    content = function(file) {write.csv(read.csv("PheWAS_test.csv"), file, row.names = FALSE)})
  
})

shinyApp(ui, server)