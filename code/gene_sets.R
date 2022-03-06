#------------------------------------------------------------------------------
# Create new combined Genesets for Angiogenesis
#-----------------------------------------------------------------------------

# GOBP_BLOOD_VESSEL_Morphognesis
# GOBP_PLACENTA_BLOOD_VESSEL_DEVELOPMENT C5
# GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_DIFFERENTIATION
# GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION
# GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS
# GOBP_BLOOD_VESSEL_MATURATION
# GOBP_BLOOD_VESSEL_MORPHOGENESIS
# GOBP_BLOOD_VESSEL_REMODELING
# GOBP_BRANCHING_INVOLVED_IN_BLOOD_VESSEL_MORPHOGENESIS
# GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS
# GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION
# GOBP_RETINAL_BLOOD_VESSEL_MORPHOGENESIS
# GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT
# HALLMARK_ANGIOGENESIS
# LU_TUMOR_ANGIOGENESIS_UP C2 CGP
# LU_TUMOR_VASCULATURE_UP C2 CGP
# RANKIN_ANGIOGENIC_TARGETS_OF_VHL_HIF2A_UP C2 CGP

# WP_PLATELETMEDIATED_INTERACTIONS_WITH_VASCULAR_AND_CIRCULATING_CELLS C2 cp CP:WIKIPATHWAYS

# group GOBP togther + subgroups

get_combined_gs <- function() {
  Hallmark <- msigdbr(species = "Homo sapiens", category = "H")
  C5_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  C5_gomf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
  C5_gocc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
  C2_cgp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")



  blood_angiogenesis_sets <- c(
    "GOBP_BLOOD_VESSEL_Morphognesis", "GOBP_PLACENTA_BLOOD_VESSEL_DEVELOPMENT",
    "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_DIFFERENTIATION", "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION",
    "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
    "GOBP_BLOOD_VESSEL_MATURATION", "GOBP_BLOOD_VESSEL_MORPHOGENESIS",
    "GOBP_BLOOD_VESSEL_REMODELING", "GOBP_BRANCHING_INVOLVED_IN_BLOOD_VESSEL_MORPHOGENESIS",
    "GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
    "GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION",
    "GOBP_RETINAL_BLOOD_VESSEL_MORPHOGENESIS", "GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT",
    "LU_TUMOR_ANGIOGENESIS_UP", "LU_TUMOR_VASCULATURE_UP", "HALLMARK_ANGIOGENESIS"
  )
  genes_c5 <- C5_gobp %>%
    filter(gs_name %in% blood_angiogenesis_sets) %>%
    select(ensembl_gene) %>%
    as.list() %>%
    unlist() %>%
    unname()

  genes_c2 <- C2_cgp %>%
    filter(gs_name %in% blood_angiogenesis_sets) %>%
    select(ensembl_gene) %>%
    as.list() %>%
    unlist() %>%
    unname()

  genes_hallmark <- Hallmark %>%
    filter(gs_name %in% blood_angiogenesis_sets) %>%
    select(ensembl_gene) %>%
    as.list() %>%
    unlist() %>%
    unname()



  all_angiogenesis_set <- c(list(genes_c2, genes_hallmark, genes_c5) %>% unique()) %>% unlist()

  GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS <- C5_gobp %>%
    filter(gs_name %in% "GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS") %>%
    select(ensembl_gene) %>%
    as.list()

  GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION <- C5_gobp %>%
    filter(gs_name %in% "GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION") %>%
    select(ensembl_gene) %>%
    as.list()
  GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION <- C5_gobp %>%
    filter(gs_name %in% "GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION") %>%
    select(ensembl_gene) %>%
    as.list()
  LU_TUMOR_ANGIOGENESIS_UP <- C2_cgp %>%
    filter(gs_name %in% "LU_TUMOR_ANGIOGENESIS_UP") %>%
    select(ensembl_gene) %>%
    as.list()
  LU_TUMOR_VASCULATURE_UP <- C2_cgp %>%
    filter(gs_name %in% "LU_TUMOR_VASCULATURE_UP") %>%
    select(ensembl_gene) %>%
    as.list()


  intersect_angiogenesis_set <- Reduce(intersect, list(
    genes_hallmark,
    GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS,
    GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION,
    LU_TUMOR_ANGIOGENESIS_UP,
    LU_TUMOR_VASCULATURE_UP,
    GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION
  )) %>% as.list()

  gs_name <- c(
    rep("C5_GOBP", length(genes_c5)),
    rep("C2_CGP", length(genes_c2)),
    rep("Hallmark", length(genes_hallmark)),
    rep("All combined", length(all_angiogenesis_set)),
    rep("Intersection Angiogenesis", length(intersect_angiogenesis_set))
  )

  ensembl_gene <- c(genes_c5, genes_c2, genes_hallmark, all_angiogenesis_set, intersect_angiogenesis_set)
  combined_gs <- data_frame(
    "gs_name" = gs_name,
    "ensembl_gene" = ensembl_gene
  )
  return(combined_gs)
}
