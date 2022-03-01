# RUN LCT
setwd("~/Documents/GitHub/nci60_gsea/GSEA_nci60/")
source("code/lct_functions.R")

pathway_activity <- read.csv("data/flux_table.txt", sep= "\t", colClasses = c(rep("character",3),rep("double",59))) %>%
  select(!c("System", "Subsystem"))

geneexpression <- read.csv("data/nci_60_expression.txt", sep = "\t") %>%
  # remove ME MDA N as no expression data is available
  dplyr::select(-"ME_MDA_N", -starts_with("..")) %>%
  # filter ensembl ids
  dplyr::filter(grepl("ENSG", Identifier)) %>%
  plyr::ddply("Identifier", numcolwise(mean))

geneset_lct_preprocessing <- function(){
  # converts geneset dataframe into matrix
  # rows - genes
  # columns - geneset
  # if gene in geneset == 1, else 0
  #gene_set_1 <- input$selected_geneset_lct
  gene_set_1 <- "Hallmark"
  Hallmark <- msigdbr(species = "Homo sapiens", category = "H")
  if (gene_set_1 == "Hallmark") {
    gene_set <- Hallmark
  } else if (gene_set_1 == "C5_gobp") {
    gene_set <- C5_gobp
  } else if (gene_set_1 == "C5_gomf") {
    gene_set <- C5_gomf
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
  geneset_matrix <- 1*geneset_matrix
  
  return(geneset_matrix)
}

phenotype_data_preprocessing <- function(phenotype){
  
  #pathway <- input$selected_pathway_lct
  pathway <- "CMP-N-acetylneuraminate synthesis"
  
  phenotype_data <- phenotype %>%
    filter(Description == pathway) %>%
    pivot_longer(!Description,names_to = "ID", values_to = "phenotype") 
  return(phenotype_data)
}


rownames(geneexpression) <- geneexpression$Identifier
geneexpression$Identifier <- NULL
geneset_matrix <- geneset_lct_preprocessing()
phenotype_data <- phenotype_data_preprocessing(phenotype = pathway_activity)
linear_combination_test(expressiondata = geneexpression, 
                        geneset = geneset_matrix,
                        phenotype_data = phenotype_data,
                        nbPermutations = 1000)