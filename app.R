#### Install packages not yet installed ####

cran_packages <- c("rlang", "shiny","shinyjs", "DT", "ggplot2", "RColorBrewer", "heatmaply", "plotly", "reshape", "tidyverse","bslib", "scales", "ggrepel")
cran_installed_packages <- cran_packages %in% rownames(installed.packages())
bioconductor_packages <- c("limma", "edgeR")
bioconductor_installed_packages <- bioconductor_packages %in% rownames(installed.packages())
if (any(cran_installed_packages == FALSE)) {
  install.packages(cran_packages[!cran_installed_packages])
}

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") 
  }

if (any(bioconductor_installed_packages == FALSE)) {
  
  suppressMessages({
    BiocManager::install(bioconductor_packages[!bioconductor_installed_packages], quiet = TRUE)
  })
}


#### Packages loading ####
invisible(lapply(cran_packages, library, character.only = TRUE))
invisible(lapply(bioconductor_packages, library, character.only = TRUE))

source("functions.R")

feature1_style='border-left: 1.25px solid #f0f0f0;z-index:1;'
feature2_style='border-left: 1.25px solid #f0f0f0;z-index:2;'
title_style='height:60px;border-left: 1.25px solid #f0f0f0;border-bottom: 1.25px solid #f0f0f0;z-index:2;'
space1_style='height:40px;z-index:2;'
space2_style='height:20px;border-left: 1.25px solid #f0f0f0;z-index:2;'
space3_style='height:20px;border-left: 1.25px solid #f0f0f0;border-top: 1.25px solid #f0f0f0;z-index:2;'

ui <- fluidPage(
  theme = bs_theme(bootswatch = "minty", base_font = font_google("Roboto")),
  tags$head(
    tags$style(
      HTML("
      @import url('https://fonts.googleapis.com/css2?family=Roboto');
        body {font-family: 'Roboto';}
        h1 {font-family: 'Roboto';}
        h2 {font-family: 'Roboto';}
        h5 {font-family: 'Roboto';}
        strong {font-family: 'Roboto';}",
        "#inTabset {
          position:fixed;
          width: 100%;
          height: 4.7%;
          background-color: white;
          margin-top: -30px;
          z-index:3;
          }", 
        ".tab-content {
        margin-top: 25px;
        z-index:1;
        }",
        
        )
      )
    ),
  
  div(
    titlePanel(
      h2(strong("RNA interference and CRISPR/Cas9 screen results"))
      ), style = "background-color:#ffffff;position:fixed; z-index: 3; top: -10px;width:100%;height:80px;"),
  br(), 
  br(),
  br(),
  sidebarLayout(
    sidebarPanel(
      style = "position:fixed;width:32%;background-color:#dfede7;",
      fileInput('files', strong('Upload workflow outputs'), accept = c(".rds"), multiple = T),
      selectInput("contrast", "Contrast", choices = ""),
      numericInput("FDR", "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("logFC", "LogFC threshold", value = 1.5, min = 0, max = 15, step = 0.5), 
      selectInput("export_gene", strong("Select gene to generate report on:"), choices = ""),
      downloadButton("report", "Generate report")
    ),
    mainPanel(
      tabsetPanel(id = "inTabset",
                  type = "tabs",
                  tabPanel(icon("house-user"),
                           fluidRow(column(12, style = 'height:20px')),
                           h5(strong("Welcome!")),
                           h6("This site allows the exploration and visualization of the snakemake workflow results from RNA interference and CRISPR/Cas9 screens. 
                              More information on how to use this workflow can be found here: https://github.com/zifornd/shrnaseq"),
                           h6("Once you have successfully run the snakemake workflow on your data, upload the workflow output where indicated on the left panel of this site. 
                             This workflow output can be found in the snakemake results directory, named “shiny.rds”."),
                           h6("You can set a threshold for false discovery rate (FDR) and log fold-change (FC) where indicated on the left, then explore the findings through the tabs above.
                              The plots throughout this site are interactive, hover over features of the plots for more specific details."),
                           h6("You can also generate a report of results on a gene of interest in your data.")
                           ),
                  tabPanel("Quality control", style = "position:relative;",
                           fluidRow(column(12, style = 'height:20px')),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('indexcounts'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Guide RNA counts per sample. Hover over each bar to get specific counts per sample"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('shcounts'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Density of guide RNA counts across all samples"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('bcv'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Variability in the screen as a functon of guide RNA abundance"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(6, plotlyOutput('mds'), style = feature1_style),
                                    column(6, plotlyOutput('correctedmds'))),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Multidimensional plot to visualise the relationships between samples, A) uncorrected and B) corrected for batch"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(6, plotlyOutput('pca'), style = feature1_style),
                                    column(6, plotlyOutput('correctedpca'))),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("First two principal components to visualise the relationships between samples, A) uncorrected and B) corrected for batch"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(6, plotlyOutput('sampledist'), style = feature1_style),
                                    column(6, plotlyOutput('correctedsampledist'))),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Distances between samples, A) uncorrected and B) corrected for batch"), style = title_style)),
                           fluidRow(column(12, style = space1_style))
                           ),
                  tabPanel("Differential expression",
                           fluidRow(column(12, style = 'height:20px')),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('hist'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Histogram of guide RNA p values for the selected contrast, hover over bars to get further details"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                        
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(7, plotlyOutput('plotsmear'), style = feature1_style),
                                    column(5, plotlyOutput('volcano'))),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(7, strong("Log fold-change (FC) against log counts-per-million (CPM) of the guide RNAs for the selected contrast. 
                                                 Hover over points to get specific information on each guide RNA")
                                           , style = title_style),
                                    column(5, strong("Volcano plot of guide RNA for the selected contrast 
                                                     Hover over points to get specific information on each guide RNA")
                                           , style = 'height:60px;border-bottom: 1.25px solid #f0f0f0;')),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(4, numericInput("topguides", "Select number of top guides", value = 5, min = 2, max = 30), style = feature2_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(6, plotlyOutput('de'), style = feature1_style),
                                    column(6, plotlyOutput('correctedde'))),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Differential expression (logCPM) across samples for the top guide RNAs, A) uncorrected and B) corrected for batch"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, selectInput("corrected", "Select uncorrected or corrected to view in the below table:", 
                                                           choices=c("Uncorrected", "Corrected")), style = feature2_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(12, DTOutput('detable'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Table of differential expression (logCPM) across samples for each guide RNA."), style = title_style)),
                           fluidRow(column(12, style = space1_style))
                           ),
                  tabPanel("Gene level",
                           fluidRow(column(12, style = 'height:20px')),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('camerarank'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Rank plot of competitive gene set  findings for the selected contrast"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput('camera'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Table of competitive gene set test"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, plotlyOutput('generank'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Gene level mean logFC rank plot across genes for the selected contrast. Change the logFC threshold on the left-hand side panel to change highlighted genes on plot"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput('genelevel'), style = feature1_style)),
                           fluidRow(column(12, strong("Table of gene level statistics for the selected contrast"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(6, selectInput("barcode_gene", "Select Gene", choices = ""), style = feature2_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(12, plotOutput('barcode'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Enrichment of guide RNAs across selected gene for the selected contrast"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           ),
                  tabPanel("Enrichment",
                           fluidRow(column(12, style = 'height:20px')),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(4, numericInput("topgos", "Select number of top terms to view"
                                                           , value = 10, min = 1, max =30), style = feature2_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(6, plotlyOutput("upgoplot"), style = 'height:600px;border-left: 1.25px solid #f0f0f0;'),
                                    column(6, plotlyOutput("downgoplot"), style = 'height:600px;')),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(6, strong("Top up-regulated gene onotology terms for the selected contrast"), style = title_style),
                                    column(6, strong("Top down-regulated gene onotology terms for the selected contrast"), style = 'height:60px;border-bottom: 1.25px solid #f0f0f0;')),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput('gotable'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong('Table of gene ontology (GO) results for the selected contrast. 
                                                      The "Ont" refers to the ontology the GO term belongs to; Biological Process (BP), Cellular Component (CC) or Molecular Function (MF)'), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(4, numericInput("topkeggs", "Select number of top pathways to view"
                                                           , value = 10, min = 1, max = 30), style = feature2_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(6, plotlyOutput("upkeggplot"), style = 'height:600px;border-left: 1.25px solid #f0f0f0;'),
                                    column(6, plotlyOutput("downkeggplot"), style = 'height:600px;')),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(6, strong("Top up-regulated KEGG pathways for the selected contrast"), style = title_style),
                                    column(6, strong("Top down-regulated KEGG pathways for the selected contrast"), style = 'height:60px;border-bottom: 1.25px solid #f0f0f0;')),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput("keggtable"), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Table of KEGG pathway results for selected contrast"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           ),
                  tabPanel("Comparing contrasts",
                           fluidRow(column(12, style = 'height:20px')),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(4, selectInput("contrast_gene", "Select Gene", choices=""), style = feature2_style)),
                           fluidRow(column(12, style = space2_style)),
                           fluidRow(column(12, plotlyOutput('guidecontrast'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Plot comparing guide RNA logFC between contrasts"), style = title_style)),
                           fluidRow(column(12, style = space1_style)),
                           
                           fluidRow(column(12, style = space3_style)),
                           fluidRow(column(12, DTOutput('FCcontrast'), style = feature1_style)),
                           fluidRow(column(12, style= space2_style)),
                           fluidRow(column(12, strong("Table of mean log fold-change of genes across contrasts"), style = title_style)),
                           fluidRow(column(12, style = space1_style))
                           )
      )
    )
  )
)

server <- function(input, output, session) {
  
##### Reactive input #####
  
  data <- reactive({
    file <- input$files
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "rds", label = "Please upload rds file"))
    readRDS(file$datapath)
  })

  contrasts <- reactive({
    name <- names(data())[which(str_detect(names(data()), "_lrt"))]
    split <- str_split(name, "_")
    as.vector(sapply(split,"[[",1))
  })
  observe({
    updateSelectInput(session, "contrast", choices = contrasts())
  })
  inputcontrast <- reactive({
    validate(need(input$contrast != "", "Please select contrast in side panel"))
    paste(input$contrast)
    })
  FCthres <- reactive({
    validate(need(input$logFC != "", "Please specify logFC threshold in side panel"),
             need(input$logFC >= 0, "Please specify a positive logFC value"))
    as.numeric(input$logFC)})
  FDRthres <- reactive({
    validate(need(input$FDR != "", "Please specify FDR threshold in side panel"),
             need((input$FDR <= 1 & input$FDR >=0), "Please specifiy a FDR threshold value between 0 and 1"))
    as.numeric(input$FDR)})
  topguides <- reactive({
    validate(need(input$topguides != "", "Please specify number of top guides to view"),
             need((input$topguides >= 2 & input$topguides <= 30), "Please specify a number between 2 and 30"))
    as.numeric(input$topguides)})
  topgos <- reactive({
    validate(need(input$topgos != "", "Please specify number of top GO terms to view"),
             need((input$topgos >= 1 & input$topgos <= 30), "Please specify a number between 1 and 30"))
    as.numeric(input$topgos)})
  topkeggs <- reactive({
    validate(need(input$topkeggs != "", "Please specify number of top KEGG pathways to view"),
             need((input$topkeggs >= 1 & input$topkeggs <= 30), "Please specify a number between 1 and 30"))
    as.numeric(input$topkeggs)})
  
##### Quality control tab #####
  
  output$indexcounts <- renderPlotly({
    indexcounts(data())
  })
  output$shcounts <- renderPlotly({
    guidecounts(data())
  })
  output$bcv <- renderPlotly({
    bcv(data())
  })
  output$mds <- renderPlotly({
    mds(data())
  })
  output$correctedmds <- renderPlotly({
    cor_mds(data())
  })
  output$pca <- renderPlotly({
    pca(data())
  })
  output$correctedpca <- renderPlotly({
    cor_pca(data())
  })
  output$sampledist <- renderPlotly({
    sampledist(data())
  })
  output$correctedsampledist <- renderPlotly({
    cor_sampledist(data())
  })
  
##### Differential expression tab #####
  
  output$hist <- renderPlotly({
    hist(data(), inputcontrast())
  })
  output$plotsmear <- renderPlotly({
    plotsmear(data(), inputcontrast(), FCthres(), FDRthres())
  })
  output$volcano <- renderPlotly({
    volcano(data(), inputcontrast(), FDRthres())
  })
  output$de <- renderPlotly({
    de(data(), inputcontrast(), topguides())
  })
  output$correctedde <- renderPlotly({
    cor_de(data(), inputcontrast(), topguides())
  })
  output$detable <- renderDT({
    datatable(detable(data(), inputcontrast(), input$corrected))
  })
  
##### Gene level tab #####
  
  output$camerarank <- renderPlotly({
    s <- input$camera_rows_selected
    camerarank(data(), inputcontrast(),FDRthres(), s)
  })
  output$camera <- renderDT({
    datatable(camera(data(), inputcontrast()))
  })
  output$generank <- renderPlotly({
    s <- input$genelevel_rows_selected
    generank(data(), inputcontrast(), FCthres(), s)
  })
  output$genelevel <- renderDT({
    datatable(genelevel(data(), inputcontrast()))
  })
  observe({
    updateSelectInput(session,
                      "barcode_gene",
                      choices = c(unique(data()$x$genes$Gene)))
  })
  barcodegene <- reactive({
    validate(need(input$barcode_gene != "", "Please select gene to view"))
    paste(input$barcode_gene)
  })
  output$barcode<- renderPlot({
    barcode(data(), inputcontrast(), barcodegene())
  })
  
##### Enrichment tab #####
  
  output$upgoplot <- renderPlotly({
    upgoplot(data(), inputcontrast(), topgos())
  })
  output$downgoplot <- renderPlotly({
    downgoplot(data(), inputcontrast(), topgos())
  })
  output$gotable <- renderDT({
    datatable(gotable(data(), inputcontrast()))
  })
  output$upkeggplot <- renderPlotly({
    upkeggplot(data(), inputcontrast(), topkeggs())
  })
  output$downkeggplot <- renderPlotly({
    downkeggplot(data(), inputcontrast(), topkeggs())
  })
  output$keggtable <- renderDT({
    datatable(keggtable(data(), inputcontrast()))
  })
  
#### Comparing contrasts tab #####
  
  observe({
    updateSelectInput(session, "contrast_gene", choices = c(unique(data()$x$genes$Gene)))
  })
  contrastgene <- reactive({
    validate(need(input$contrast_gene != "", "Please select gene to view"))
    paste(input$contrast_gene)
  })
  output$guidecontrast <- renderPlotly({
    guidecontrast(data(), contrastgene())
  })
  output$FCcontrast <- renderDT({
    datatable(FCcontrast(data(), inputcontrast()),
              selection = "single")
  })
  observe({
    updateSelectInput(session, "export_gene", choices = c(unique(data()$x$genes$Gene)))
  })
  exportgene <- reactive({
    paste(input$export_gene)
  })
  output$report <- downloadHandler(
    filename <- function() {
      paste0(exportgene(), ".html")
      },
    content <- function(file) {
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      params <- base::list(gene = exportgene(),
                           data=data())
      rmarkdown::render(tempReport, 
                        output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

shinyApp(ui,server)
