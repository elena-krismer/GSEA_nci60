#setwd("~/Documents/GitHub/nci60_gsea/")

library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(shiny)

source("code/functions.R")

# angiogensesis GOBP_BLOOD_VESSL_MORPHOGENESIS
# https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=GOBP_BLOOD_VESSEL_MORPHOGENESIS



#------------------------------------------------------------------------------
# RUN GSEA
#------------------------------------------------------------------------------

run_gsea <-  function(gene_set, FC.vec, scoreType){
  gsea_results <- fgseaSimple(pathways = gene_set,
                     stats = FC.vec,
                     scoreType = scoreType,
                     nperm = 1000) # permutations depending on the p-value
  return(gsea_results)
}


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Shiny
# ----------------------------------------------------------------------------

# input variables 
active_pathway_options <-read.csv("data/clustering_pathwayenrichment.txt", 
                                  sep = "\t") %>%
  select(Description) # %>%
 # as.list()
  

gene_set_database_options <- c("Hallmark", "C5_gobp", "C5_gomf", "C5_gocc_genset")


Hallmark = msigdbr(species = "Homo sapiens", category = "H")
C5_gobp = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
C5_gomf = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
C5_gocc = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")


ui <- fluidPage(
  titlePanel("Gene set enrichment Analysis"),
  sidebarLayout(
    # --------------------------------------------------------------------------
    # Side bar Panel
    # --------------------------------------------------------------------------
    sidebarPanel(
      tags$form(
        selectInput("selected_pathway", "Select active pathway", 
                    active_pathway_options),
        selectInput("selected_geneset", "Select geneset database", 
                    gene_set_database_options)
      )
    ),
    # --------------------------------------------------------------------------
    # Main Panel
    # --------------------------------------------------------------------------
    mainPanel(
      # subsetting main panel into tabs
      tabsetPanel(
        type = "tabs",
        # plots
        tabPanel("Significant Gene set", dataTableOutput("gsea_table")),
        # reaction description
        tabPanel("GSEA Plot", plotOutput("gsea_plot")),
        actionButton("button", strong("Submit"))
        )
      )
    )
  )

server <- function(input, output) {
  
  selected_pathway <- eventReactive(input$button,{
    return(input$selected_pathway)
  })
  gene_set <- eventReactive(input$button, {
    return(input$selected_geneset)
  })
    
  FC.vec <- geneexpression_data_preprocessing(selected_pathway)
  
  # set score type
  scoreType <- "std"
  # if both negative score type is negative
  if(min(FC.vec) < 0 && max(FC.vec) < 0){
    scoreType <- "pos"
  }
  if(min(FC.vec) > 0 && max(FC.vec) > 0){
    scoreType <- "neg"
  }
  
  gsea_results <- run_gsea(gene_set, FC.vec, scoreType)
  
  # gsea table
  output$gsea_table <- renderDataTable(
    gsea_signficiant(gsea_results)
  )
  # gsea plot
  output$gsea_plot <- renderPlot({
    gsea_results <- run_gsea(gene_set, FC.vec, scoreType)
    gsea_plot_func(gsea_results)
  })
}


shinyApp(ui, server)

