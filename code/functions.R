

#------------------------------------------------------------------------------
# foramt for gsea
#-----------------------------------------------------------------------------

# gsea database
# get all genes into list for one hallmark
# needs to be in list for GSEA

geneset_database_formating <- function(gene_set_selected) {
  gene_set <- gene_set_selected %>%
    select(gs_name, ensembl_gene) %>%
    group_by(gs_name) %>%
    summarise(all.genes = list(ensembl_gene)) %>%
    deframe()
  return(gene_set)
}


#-----------------------------------------------------------------------------
# get gene expression data
# distinguish them by pathway activity INACTIVE/ACTIVE
# gives back vector with fold changes for gsea
#-----------------------------------------------------------------------------
geneexpression_data_preprocessing <- function(selected_pathway, 
                                              path_task_analysis = path_task_analysis_results) {
  geneexpression <- vroom("data/nci_60_expression.txt", delim = "\t") %>%
    # remove ME MDA N as no expression data is available
    select(-ME_MDA_N, -starts_with("..")) %>%
    # filter ensembl ids
    dplyr::filter(grepl("ENSG", Identifier)) %>%
    pivot_longer(!Identifier, names_to = "cellline", values_to = "expression") %>%
    dplyr::filter(!grepl("X", cellline)) # remove obsolet x rows
  
  pathway_activity <- read.csv("data/flux_results_models_publication.txt",
                               sep = "\t",
                               colClasses = c(rep("character", 3), rep("double", 59))) %>%
    dplyr::select(!c("System", "Subsystem")) %>%
    column_to_rownames("Description") %>%
    mutate_all(function(x) as.numeric(as.character(x))) %>%
    apply(1, function(x)(x-min(x))/(max(x)-min(x))) %>%
    as.data.frame() %>%
    rownames_to_column("Description")

  # get activity from model analysis
  inactive <- pathway_activity %>%
    dplyr::filter(Description == selected_pathway) %>%
    dplyr::select(Description, celllines) %>%
    pivot_longer(!Description, names_to = "cellline", values_to = "activity") %>%
    dplyr::select(cellline, activity) %>%
    group_by(activity) %>%
    dplyr::filter(activity == 0)


  inactive <- inactive$cellline
  active <- setdiff(celllines, inactive)

  # add column for two groups
  geneexpression <- geneexpression %>%
    mutate(group = case_when(
      cellline %in% inactive ~ "INACTIVE",
      cellline %in% active ~ "ACTIVE",
    )) %>%
    na.omit() %>%
    dplyr::select(-cellline) %>%
    pivot_wider(names_from = group, values_from = expression)

  # get mean of rows
  geneexpression$ACTIVE_mean <- sapply(
    strsplit(gsub(
      "[c()]", "",
      geneexpression$ACTIVE
    ), ","),
    function(x) mean(as.numeric(x), na.rm = TRUE)
  )
  geneexpression$INACTIVE_mean <- sapply(
    strsplit(gsub(
      "[c()]", "",
      geneexpression$INACTIVE
    ), ","),
    function(x) mean(as.numeric(x), na.rm = TRUE)
  )

  geneexpression$ACTIVE <- NULL
  geneexpression$INACTIVE <- NULL


  # calculate Fold changes
  # extract expression data
  FC <- geneexpression %>%
    # Calculate inactive minus active fold change
    # beause exprssion values are log2, substratction is the same as divison
    # https://www.biostars.org/p/209790/
    # mutate(delta = INACTIVE-ACTIVE) %>%
    # Calculate mean fold change per gene
    group_by(Identifier) %>%
    summarise(mean.delta = ACTIVE_mean - INACTIVE_mean) %>%
    # arrange by descending fold changes
    arrange(desc(mean.delta))

  # get fold changes into vector
  FC.vec <- FC$mean.delta
  names(FC.vec) <- FC$Identifier

  return(FC.vec)
}


#------------------------------------------------------------------------------
# perform gsea
#-----------------------------------------------------------------------------

# filter gsea table by signficant adjusted p value
gsea_signficiant <- function(gsea) {
  signficant_df <- gsea %>%
    filter(padj <= 0.05) %>%
    select(!leadingEdge)
  return(signficant_df)
}

gsea_plot_func <- function(gsea) {
  plot <- gsea %>%
    filter(padj <= 0.05) %>%
    ggplot(aes(
      x = reorder(pathway, NES),
      y = NES
    )) +
    geom_col() +
    theme_classic() +
    # force equal max min
    lims(y = c(-3.2, 3.2)) +
    labs(
      y = "Normalized enrichment score (NES)",
      x = "Gene set",
      title = "Hallmark GSEA (FDR < 0.05"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  return(plot)
}

#------------------------------------------------------------------------------
# perform LCT
#-----------------------------------------------------------------------------

lct_data_preprocessing <- function(path_task_analysis = path_task_analysis_results){
  
  # geneset: from dataframe to binary matrix 
  # pathway_activity_ranking: order according to expression data
  

  
  geneexpression <- read.csv("data/nci_60_expression.txt", sep = "\t") %>%
    # remove ME MDA N as no expression data is available
    select(-ME_MDA_N) %>%
    # filter ensembl ids
    dplyr::filter(grepl("ENSG", Identifier)) %>%
    ddply("Identifier", numcolwise(mean))
  
  row.names(geneexpression) <- geneexpression$Identifier
  geneexpression$Identifier <- NULL
  cellline_order <- colnames(geneexpression) %>% as.vector()
  
  pathway_activity_ranking <- pathway_activity <- read.csv(path_task_analysis, sep = "\t") %>%
    select(!c("System", "Subsystem")) %>%
    filter(Description == "CMP-N-acetylneuraminate synthesis")
  row.names(pathway_activity_ranking) <- pathway_activity_ranking$Description
  pathway_activity_ranking$Description <- NULL
  # order df according to expression data
  pathway_activity_ranking <-  pathway_activity_ranking %>%
    setcolorder(cellline_order)
  
}
