#------------------------------------------------------------------------------
# Correlation matrix
#------------------------------------------------------------------------------

`%notin%` <- Negate(`%in%`)

correlation_df <- function(path_task_analysis = path_task_analysis_results){
  
  Hallmark <- msigdbr(species = "Homo sapiens", category = "H")
  
  geneexpression <- read.csv("data/nci_60_expression.txt", sep= "\t") %>%
    # remove ME MDA N as no expression data is available
    dplyr::select(-c("ME_MDA_N"), -starts_with("X"))%>%
    # filter ensembl ids
    dplyr::filter(grepl("ENSG", Identifier)) %>%
    tidyr::pivot_longer(!Identifier, names_to = "cellline", values_to = "expression", 
                        values_ptypes = list(val = 'character')) %>%
    mutate(ensembl_gene = Identifier) %>%
    inner_join(Hallmark, by = "ensembl_gene") %>%
    select(cellline, expression, gs_name) %>%
    mutate(gs_name = gsub("HALLMARK_", "", gs_name)) %>%
    transform(expression = as.numeric(expression)) %>%
    pivot_wider(names_from = gs_name, values_from = expression) %>%
    filter(across(everything(), ~ !grepl('NA', .)))

  cellline <- geneexpression$cellline
  geneexpression$cellline <- NULL

  gene_sets <- colnames(geneexpression)
  for (set in gene_sets) {
    (
      geneexpression[[set]] <- sapply(
        strsplit(gsub(
         "[c()]", "",
          geneexpression[[set]]
        ), ","),
        function(x) mean(as.numeric(x), na.rm = TRUE)
      )
    )
  }
  rownames(geneexpression) <- cellline

  pathway_activity <- vroom(path_task_analysis, delim = "\t", 
                            col_types = cols(.default = "c")) %>%
    select(!c("System", "Subsystem")) %>%
    pivot_longer(cols = c(-Description), names_to = "X") %>%
    pivot_wider(names_from = c(Description)) %>%
    as.data.frame()

  rownames(pathway_activity) <- pathway_activity$X
  pathway_activity$X <- NULL
  fluxomics <- colnames(pathway_activity)
  rownames(pathway_activity) <- rownames(geneexpression)

  transcriptomics_fluxomics <- merge(geneexpression, pathway_activity, by = "row.names") %>%
    select(!Row.names) %>%
    na.omit()  %>%
    mutate_all(function(x) as.numeric(as.character(x))) 
  rownames(transcriptomics_fluxomics) <- rownames(pathway_activity)
  transcriptomics_fluxomics$Row.names <- NULL
  correlation <- cor(transcriptomics_fluxomics)
  
  correlation$pathway <- rownames(correlation)

  df <- correlation %>%
    as.data.frame() %>%
  select(!fluxomics) %>%
  filter(row.names(correlation) %in% fluxomics) %>%
  na.omit() %>%
  as.data.frame() %>%
  select(pathway, everything())
  #as.matrix() %>%
  #heatmap()
  return(df)
}

