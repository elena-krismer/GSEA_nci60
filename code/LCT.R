# code from: https://sites.ualberta.ca/~yyasui/homepage.html
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-260

#*********************************************************************************
#                               Basic Function: LCT
#                              (maximum correlation)
#                  Gene-set significant Linear Combination Test  
#                           for non-linear association and 
#                          multiple experimental categories
#*********************************************************************************

# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")

library(corpcor)
library(qvalue)
library(splines)

GS.format.dataframe.to.list <- function(GS){
  if(is.data.frame(GS)){
    genes <- rownames(GS)
    L <- NULL
    for(ags in names(GS)){
      w <- which(GS[,ags]==1)
      if(length(w)>0)  {
        L <- c(L,list(genes[w]))
        names(L)[length(L)] <- ags
      }
    }
    L
  }else{
    GS
  }
}

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

ns.ext<-function(x,df=5){
  # x : matrix of basis
  # df: degree of fredom
  
  p<-ncol(x)
  y<-NULL
  for (i in 1:p){ y<-cbind(y,ns(x[,i],5))}
  y
}

LCT <- function(GS, DATA, cl, nbPermutations=1000, nonlinear=FALSE, degree=5, silent=FALSE){
  
  # GS : gene sets
  #      -> a dataframe with rows=genes,
  #                        columns= gene sets,
  #                        GS[i,j]=1 if gene i in gene set j
  #                        GS[i,j]=0 otherwise
  #         OR
  #         a list with each element corresponding to a gene set = a vector of strings (genes identifiers)
  #
  #
  # DATA : expression data
  #     -> a dataframe with  rows=genes,
  #                       columns=samples
  #
  # cl : response data for the samples IN THE SAME ORDER AS IN DATA
  #    -> a dataframe with rows=respons,
  #                      colmns=samples
  #
  #nonlinear: if TRUE, LCT will do a nonlinear test, using nature splines (df=degree) to do nonlinear transformation on GS. 
  
  # (1) pre-treatment of the gene sets and response vector
  
  genes <- rownames(DATA)       # gene names of the microarray data 
  nb.Samples  <- ncol(DATA)     # nb of samples
  nb.GeneSets <- dim(GS)[2]     # nb of gene sets
  
 # GS <-  GS.format.dataframe.to.list(GS)
  GS <-  GS %>% as.list()
  # change format of GS from dataframe to list
  GS <-  lapply(GS,function(z) as.numeric(which(genes %in% z)))
  GS.sizes <- sapply(GS,length) # size of each gene set
  # numericalized index of each GS
  GS.data <- lapply(GS, function(z) as.matrix(DATA[z, ],ncol=nb.Samples))
  # creat data of each GS (rows=genes, columns=samples)
  GS.data <- lapply(GS.data,function(z) t(z))
  # transfer genes in each GS (columns=gene, rows=samples)
  if(nonlinear) GS.data <- lapply(GS.data,function(z) ns.ext(z,df=degree))
  # transformation genes using nature splines (columns=gene, rows=samples)
  cl=t(cl) #transfer response vectors (column=respons, rowa=samples)
  
  
  # (2) Eigen-decompsition of shrinkage pooled covariance matrix for each GS
  
  Cov.Pooled<-lapply(GS.data, function(z) cov.shrink(z,verbose=FALSE, lambda.var=0));
  # # pooled covariance of genes in each GS
  for (i in 1:nb.GeneSets){
     EIGEN.decom<-eigen(Cov.Pooled[[i]]);
  # #eigen decomposition of pooled covariance for each GS
     D<-EIGEN.decom$values;      # shrinkag by adding a possitive constant s0
     U<-EIGEN.decom$vectors;
  #   
      GS.data[[i]]<-t(GS.data[[i]]%*%U)/sqrt(D)
  #   # adjust data of each GS (rows=genes, columns=samples)
  }
  
  Cov.Pooled<-cov.shrink(cl,verbose=FALSE, lambda.var=0);
  EIGEN.decom<-eigen(Cov.Pooled);
  # eigen decomposition of pooled covariance for cl
  D<-EIGEN.decom$values;         
  U<-EIGEN.decom$vectors;
  cl<-t(cl%*%U)/sqrt(D)          # adjust cl (rows=respons, columns=samples)
  
  # (3) T-like stats obtained on 'true' data
  sam.sumsquareT.obs  <- sapply(GS.data, function(z) T2.like.SAMGS(z,cl))
  # the T-like statistics obtained on 'true' data
  
  # (4) stats obtained on 'permuted' data
  sam.sumsquareT.permut <- matrix(NA,nbPermutations,nb.GeneSets)
  for(i in 1:nbPermutations) {
    ind <- sample(nb.Samples)
    sam.sumsquareT.permut[i,] <- sapply(GS.data, function(z) T2.like.SAMGS(z[,ind],cl))
    # SAMGS statitic for each gene set  - for current permutation
    if(!silent & i%%50 == 0)print(paste(i," permutations done."))
  }
  
  # (5) p-value and q-value
  GeneSets.pval <- apply(t(sam.sumsquareT.permut) >= sam.sumsquareT.obs,1,sum)/nbPermutations
  
  if(nb.GeneSets>=2){
    GeneSets.qval <-0; #qvalue(GeneSets.pval)$qvalues
    
    res <- as.data.frame(cbind("GS size"              = GS.sizes,
                               "GS p-value" 	       = GeneSets.pval,
                               "GS q-value"           = GeneSets.qval ))
    res <- cbind(res,"GS name"= names(GS))[c(4,1:3)]
  }
  
  if(nb.GeneSets==1){
    #if there is only one set, no need to calculate q-value.
    res <- as.data.frame(cbind("GS size"              = GS.sizes, ##GeneSets.sizes,
                               "GS p-value" 	       = GeneSets.pval))
    res <- cbind(res,"GS name"= names(GS))[c(3,1:2)]
  }
  rownames(res)<-NULL
  list("GS stats"=res)
}