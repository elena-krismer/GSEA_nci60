# ----------------------------------------------------------------------------
# Linear Combination test
# "GSEA" for continous phenotype
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-212
# ----------------------------------------------------------------------------


T2.like.SAMGS <- function(DATA, cl){
    # DATA : expression data
    #     -> dataframe with  rows=genes,
    #                        columns=samples,
    # weight: weights to T2 statistics of genes. 
    # cl : response vectors for the samples IN THE SAME ORDER AS IN DATA
    #    -> a dataframe with rows=respons,
    #                      colmns=samples
    
    cl<-as.matrix(cl)
    DATA<-as.matrix(DATA)
    Sigma<-cov(t(DATA),t(cl))
    dim.pq<-dim(Sigma)
    if(dim.pq[1]>dim.pq[2]) max(eigen(t(Sigma)%*%Sigma,symmetric=T,only.values=T)[[1]])
    else max(eigen(Sigma%*%t(Sigma),symmetric=T,only.values=T)[[1]])
}



#' @param expressiondata columns - samples, rows - genes
#' @param geneset genset_matrix, columns - genesets, rows - genes, binary matrix
#' @param phenotype_data ID = column - matching sample names, phenotype-column = 
#' data for ranking of samples
#' @param nbPermutations number of permutations default 1000
linear_combination_test <- function(expressiondata, geneset, phenotype_data, nbPermutations = 1000){
  phenotype <- "phenotype"
  ID <- "ID"
  #expressiondata<- geneexpression
  #colnames(expressoniondataexpressiondata$Identifier
  #geneset <- geneset_matrix

  EXPR <- expressiondata %>% t()
  GS <- geneset
  LongData <- phenotype_data
  nb.Samples<-length(unique(LongData$ID)) #Number of Subjects
  model.coefs<-NULL
  
  #Ordering all data
  EXPR<- EXPR[order(rownames(EXPR)),]
  EXPR<- EXPR[,order(colnames(EXPR))]
  GS<-GS[,order(colnames(GS))]
  
  EXPR2 <- EXPR %>%t()
  
  #Step 2: Analyses of Between Subject Variations
  GS=GS[order(rownames(GS)),]
  EXPR2=EXPR2[,order(colnames(EXPR2))]
  GS2= as.data.frame((GS))
  GS.sizes <- sapply(GS2,sum) # size of each gene set
  GS.data<-NULL	
  GS.data <- lapply(data.frame(GS2), function(z) (EXPR2)[z==1, ]) # creat data of each GS (rows=gene trends, columns=samples/families)
  GS.data <- lapply(GS.data,function(z) t(z)) # transfer genes in each GS (columns=gene trends, rows=samples/families)
  
  # (2) Eigen-decompsition of shrinkage pooled covariance matrix for each GS
  
  Cov.Pooled<-lapply(GS.data, function(z) cov.shrink(z,verbose=FALSE, lambda.var=0))  # pooled covariance of genes in each GS
  
  nb.GeneSets<-dim(GS)[2]
  
  for (i in 1:nb.GeneSets){
    EIGEN.decom<-eigen(Cov.Pooled[[i]])
    # eigen decomposition of pooled covariance for each GS
    D<-EIGEN.decom$values;      
    U<-EIGEN.decom$vectors;
    if(sum(D<0)>0){
      class(Cov.Pooled[[i]])="matrix"
      EIGEN.decom<-eigen(nearPD((Cov.Pooled[[i]]))[[1]]);
      D<-EIGEN.decom$values;      
      U<-EIGEN.decom$vectors;
      
    }
    GS.data[[i]]<-t(GS.data[[i]]%*%U)/sqrt(D)
    # adjust data of each GS (rows=genes, columns=samples)
  }
  
  
  phenotype_data_betas3 <- phenotype_data %>% 
    select(!c("Description")) %>% as.vector()
  rownames(phenotype_data_betas3) <- phenotype_data_betas3$ID
  phenotype_data_betas3$ID <- NULL
  betas3 <- phenotype_data_betas3 %>% as.matrix() %>% as.vector()
    #()
  Cov.Pooled<-cov.shrink(betas3,verbose=FALSE, lambda.var=0)
  EIGEN.decom<-eigen(Cov.Pooled);  # eigen decomposition of pooled covariance for beta
  
  D<-EIGEN.decom$values        
  U<-EIGEN.decom$vectors
  cl<-t(betas3%*%U)/sqrt(D)          # adjust beta (rows=respons, columns=samples)
  

  # (3) T-like stats obtained on 'true' data
  sam.sumsquareT.obs  <- sapply(GS.data, function(z) T2.like.SAMGS(z,cl))
  # the T-like statistics obtained on 'true' data
  
  
  # (4) stats obtained on 'permuted' data
  sam.sumsquareT.permut <- matrix(NA,nbPermutations,nb.GeneSets)
  for(i in 1:nbPermutations) {
    ind <- sample(nb.Samples)
    sam.sumsquareT.permut[i,] <- sapply(GS.data, function(z) T2.like.SAMGS(z[,ind],cl))
    # SAMGS statitic for each gene set  - for current permutation
  }
  
  # (5) p-value and q-value
  GeneSets.pval <- apply(t(sam.sumsquareT.permut) >= sam.sumsquareT.obs,1,sum)/nbPermutations
  
  GeneSets.qval <-0
  GeneSets.qval <-p.adjust(GeneSets.pval,method="fdr")
  
  res <- as.data.frame(cbind("GS size"              = GS.sizes,
                             "GS p-value" 	       = GeneSets.pval,
                             "GS q-value"           = GeneSets.qval))
  res$gene_set <- rownames(res)
  return(res)
}