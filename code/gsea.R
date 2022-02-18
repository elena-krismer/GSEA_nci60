setwd("~/Documents/GitHub/nci60_gsea/")

library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(fgsea)

# hallmarks
hallmark_genset = msigdbr(species = "Homo sapiens", category = "H")
# angiogensesis GOBP_BLOOD_VESSL_MORPHOGENESIS
# https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=GOBP_BLOOD_VESSEL_MORPHOGENESIS
c5_gobp_genset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
c5_gomf_genset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
c5_gocc_genset = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")

gene_set_selected <- hallmark_genset

selected_pathway <- "Tyrosine to dopamine"

geneexpression <- read.csv("data/nci_60_expression.txt", sep = "\t") %>% 
  select(-ME.MDA.N) %>% 
  # filter ensembl ids
  dplyr::filter(grepl('ENSG', Identifier)) %>%
  pivot_longer(!Identifier, names_to = "cellline", values_to = "expression") %>%
  dplyr::filter(!grepl('X.', cellline)) # remove obsolet x rows

pathway_activity <- read.csv("data/clustering_pathwayenrichment.txt", sep = "\t")
celllines <- colnames(pathway_activity[4:62,])[4:62] %>% as.vector()


# get activity from model analysis
inactive <- pathway_activity %>% 
  dplyr::filter(pathway_activity$Description == selected_pathway ) %>%
  select(Description,celllines) %>%
  pivot_longer(!Description, names_to = "cellline", values_to = "activity") %>%
  select(cellline, activity) %>%
  group_by(activity) %>%
  filter(activity == 0)

inactive <- inactive$cellline
active <- setdiff(celllines, inactive)

# add column for two groups
geneexpression <- geneexpression %>%
  mutate(group = case_when(
    cellline %in% inactive ~ "INACTIVE",
    cellline %in% active ~ "ACTIVE",
  )) %>%
  na.omit() %>%
  select(-cellline) %>%
  pivot_wider(names_from = group, values_from = expression)

geneexpression$ACTIVE_mean <-  sapply(strsplit(gsub('[c()]', '', 
                                geneexpression$ACTIVE), ","), 
                                 function(x) mean(as.numeric(x), na.rm = TRUE))
geneexpression$INACTIVE_mean <- sapply(strsplit(gsub('[c()]', '', 
                                    geneexpression$INACTIVE), ","), 
                                  function(x) mean(as.numeric(x), na.rm = TRUE))

geneexpression$ACTIVE <- NULL
geneexpression$INACTIVE <- NULL
  
#------------------------------------------------------------------------------
# calculate Fold changes
#------------------------------------------------------------------------------
# extract expression data

FC <- geneexpression %>%
  # Calculate inactive minus active fold change 
  # beause exprssion values are log2, substratction is the same as divison
  # https://www.biostars.org/p/209790/
  #mutate(delta = INACTIVE-ACTIVE) %>%
  # Calculate mean fold change per gene
  group_by(Identifier) %>%
  summarise(mean.delta = log2(ACTIVE_mean/INACTIVE_mean)) %>%
  # arrange by descending fold changes
  arrange(desc(mean.delta))

#------------------------------------------------------------------------------
# foramt for gsea
#----------------------------------------------------------------------------- 

# gsea database
# get all genes into list for one hallmark
# needs to be in list for GSEA
gene_set <- gene_set %>%
  select(gs_name, ensembl_gene) %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(ensembl_gene)) %>%
  deframe()

# get fold changes into vector 
FC.vec <- FC$mean.delta
names(FC.vec) <- FC$Identifier

# set score type
scoreType <- "std"
# if both negative score type is negative
if(min(FC.vex) < 0 && max(FC.vex) < 0){
  scoreType <- "max"
}
if(min(FC.vex) > 0 && max(FC.vex) > 0){
  scoreType <- "neg"
}

#------------------------------------------------------------------------------
# RUN GSEA
#------------------------------------------------------------------------------

gsea <- fgseaSimple(pathways = gene_set,
                     stats = FC.vec,
                     scoreType = scoreType,
                     nperm = 1000) # permutations depending on the p-value

# plot
gsea  %>%
  filter(padj <= 0.05) %>%
  ggplot(aes(x=reorder(pathway, NES),
             y =NES))+
  geom_col() +
  theme_classic() +
  # force equal max min
  lims(y=c(-3.2,3.2))+
  labs(y = "Normalized enrihcmend score (NES)",
       x = "Gene set",
       title = "Hallmark GSEA (FDR < 0.05")


