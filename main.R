# setwd("~/Documents/GitHub/nci60_gsea/")

library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(shiny)
library(dplyr)
library(fgsea)
library(plyr)
library(data.table)
library(vroom)
library(qdapTools)
library(corpcor)
library(qvalue)
library(splines)

source("code/functions.R")
source("code/lct_functions.R")
source("code/correlation.R")
source("code/gene_sets.R")

# angiogensesis GOBP_BLOOD_VESSL_MORPHOGENESIS
# https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=GOBP_BLOOD_VESSEL_MORPHOGENESIS

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Shiny
# ----------------------------------------------------------------------------

#

path_task_analysis_results <- "data/flux_table.txt"
# input variables
active_pathway_options <- vroom(path_task_analysis_results,
  delim = "\t"
) %>%
  mutate(
    sum =
      dplyr::select(., is_numeric) %>%
        rowSums(na.rm = TRUE)
  ) %>%
  dplyr::filter(sum != 61000) %>%
  dplyr::filter(sum != 0) %>%
  dplyr::select(Description) %>%
  as.list() %>%
  unlist()

names(active_pathway_options) <- NULL

gene_set_database_options <- c("Hallmark", "C5_gobp", "C5_gomf", "C5_gocc",
                               "C2_cgp", "combined_gs")


Hallmark <- msigdbr(species = "Homo sapiens", category = "H")
C5_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
C5_gomf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
C5_gocc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
C2_cgp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
combined_gs <- get_combined_gs()


ui <- fluidPage(
  titlePanel("Gene Set Enrichment Analysis"),
  sidebarLayout(
    # --------------------------------------------------------------------------
    # Side bar Panel
    # --------------------------------------------------------------------------
    sidebarPanel(
      tags$form(
        selectInput(
          inputId = "selected_pathway", label = "Select active pathway",
          choices = active_pathway_options, selected = "CMP-N-acetylneuraminate synthesis"
        ),
        selectInput(
          inputId = "selected_geneset", label = "Select geneset database",
          choices = gene_set_database_options, selected = "Hallmark"
        ) # ,
        # actionButton("button", strong("Submit"))
      ),
      tags$form(
        selectInput(
          inputId = "selected_pathway_lct", label = "Select active pathway",
          choices = active_pathway_options, selected = "CMP-N-acetylneuraminate synthesis"
        ),
        selectInput(
          inputId = "selected_geneset_lct", label = "Select geneset database",
          choices = gene_set_database_options, selected = "Hallmark"
        ),  actionButton("button", "Submit")# ,
        # actionButton("button", strong("Submit"))
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
        tabPanel("Significant Gene set GSEA", dataTableOutput("gsea_table")),
        # reaction description
        tabPanel("GSEA Plot", plotOutput("gsea_plot")),
        tabPanel("Linear Combination Test", dataTableOutput("lct_table")),
        tabPanel("Correlation", dataTableOutput("correlation_df"))
      )
    )
  )
)

server <- function(input, output) {
  
  path_task_analysis_results <- "data/flux_results_models_publication.txt"
  
  #------------------------------------------------------------------------------
  # run gsea
  #-----------------------------------------------------------------------------

  run_gsea <- function() {
    pathway <- input$selected_pathway
    gene_set_1 <- input$selected_geneset

    # find better way
    if (gene_set_1 == "Hallmark") {
      gene_set <- Hallmark
    } else if (gene_set_1 == "C5_gobp") {
      gene_set <- C5_gobp
    } else if (gene_set_1 == "C5_gomf") {
      gene_set <- C5_gomf
    } else if (gene_set_1 == "combined_gs") {
      gene_set <- combined_gs
    } else if (gene_set_1 == "C2_cgp") {
      gene_set <- C2_cgp
    } else {
      gene_set <- C5_gocc
    }

    FC.vec <- geneexpression_data_preprocessing(pathway)

    # set score type
    scoreType <- "std"
    # if both negative score type is negative
    if (min(FC.vec) < 0 && max(FC.vec) < 0) {
      scoreType <- "pos"
    }
    if (min(FC.vec) > 0 && max(FC.vec) > 0) {
      scoreType <- "neg"
    }
    gene_set <- gene_set %>%
      geneset_database_formating()
    gsea_results <- fgseaSimple(
      pathways = gene_set,
      stats = FC.vec,
      scoreType = scoreType,
      nperm = 1000
    ) # permutations depending on the p-value
    return(gsea_results)
  }

  #------------------------------------------------------------------------------
  # output gsea
  #-----------------------------------------------------------------------------
  # gsea table
  output$gsea_table <- renderDataTable({
    gsea_results <- run_gsea()
    gsea_signficiant(gsea_results)
  })
  # gsea plot
  output$gsea_plot <- renderPlot({
    gsea_results <- run_gsea()
    gsea_plot_func(gsea_results)
  })

  
  #-----------------------------------------------------------------------------
  # run lct
  #-----------------------------------------------------------------------------
  
  pathway_activity <- read.csv("data/flux_results_models_publication.txt",
                               sep = "\t",
                               colClasses = c(rep("character", 3), rep("double", 59))) %>%
    dplyr::select(!c("System", "Subsystem")) %>%
    column_to_rownames("Description") %>%
    mutate_all(function(x) as.numeric(as.character(x))) %>%
    apply(1, function(x)(x-min(x))/(max(x)-min(x))) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Description")
  
  geneexpression <- read.csv("data/nci_60_expression.txt", sep = "\t") %>%
    # remove ME MDA N as no expression data is available
    dplyr::select(-"ME_MDA_N", -starts_with("..")) %>%
    # filter ensembl ids
    dplyr::filter(grepl("ENSG", Identifier)) %>%
    plyr::ddply("Identifier", numcolwise(mean))
  
  geneset_lct_preprocessing <- function() {
    # converts geneset dataframe into matrix
    # rows - genes
    # columns - geneset
    # if gene in geneset == 1, else 0
    gene_set_1 <- input$selected_geneset_lct
    
    if (gene_set_1 == "Hallmark") {
      gene_set <- Hallmark
    } else if (gene_set_1 == "C5_gobp") {
      gene_set <- C5_gobp
    } else if (gene_set_1 == "C5_gomf") {
      gene_set <- C5_gomf
    } else if (gene_set_1 == "combined_gs") {
      gene_set <- combined_gs
    } else if (gene_set_1 == "C2_cgp") {
      gene_set <- C2_cgp
    } else {
      gene_set <- C5_gocc
    }
    
    gene_set_lct <- gene_set %>% select(gs_name, ensembl_gene)
    l <- split(gene_set_lct, gene_set_lct$gs_name, drop = TRUE)
    nestedlist <- sapply(l, function(x) select(x, ensembl_gene) %>% as.list())
    # create binary matrix out of lists
    geneset_matrix <- mtabulate(nestedlist) != 0
    # rows - genes; columns  - gene_sets
    geneset_matrix <- geneset_matrix %>% t()
    # convert TRUE FALSE to 1 and 0
    geneset_matrix <- 1 * geneset_matrix
    
    return(geneset_matrix)
  }
  
  phenotype_data_preprocessing <- function(phenotype) {
    pathway <- input$selected_pathway_lct
    
    phenotype_data <- phenotype %>%
      filter(Description == pathway) %>%
      pivot_longer(!Description, names_to = "ID", values_to = "phenotype")
    return(phenotype_data)
  }
  
  geneexpression_lct <- geneexpression
  rownames(geneexpression_lct) <- geneexpression_lct$Identifier
  geneexpression_lct$Identifier <- NULL
  
  
  #------------------------------------------------------------------------------
  # output lct
  #-----------------------------------------------------------------------------
  observeEvent(input$button, {
  output$lct_table <- renderDataTable({
    geneset_matrix <- geneset_lct_preprocessing()
    phenotype_data <- phenotype_data_preprocessing(phenotype = pathway_activity)
    lct_df <- linear_combination_test(
      expressiondata = geneexpression_lct,
      geneset = geneset_matrix,
      phenotype_data = phenotype_data,
      nbPermutations = 1000
    )
  })
  })
  #------------------------------------------------------------------------------
  # output correlation
  #-----------------------------------------------------------------------------
  output$correlation_df <- renderDataTable({
    correlation_df()
  })
}


shinyApp(ui, server)
