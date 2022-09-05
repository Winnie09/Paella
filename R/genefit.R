#' Fitting temporal gene expression dynamics
#'
#' Fitting temporal gene expression dynamics along pseudotime for each spatial sub-trajectory.
#'
#' @param mat A gene expression matrix where each row is a gene and each column is a cell/spot. Row names are gene names and column names are cell/spot names. The column names must be present when fitting \code{paella()}.
#' @param paella Direct output from calling \code{paella()}.
#' @param gene A character vector of gene names for which pseudotime pattern will be fitted. Must be a subset of row names of \code{mat}.
#' @param interpoint A positive integer indicating the number of interpolation points.
#' @import Matrix mgcv
#' @export
#' @return A list containing scaled fitted values, original mean, and original standard deviations.
#' @author Wenpin Hou<wh2526@@cumc.columbia.edu>, Zhicheng Ji
#' @examples 
#' pt <- setNames(1:100, paste0('cell',1:100))
#' cd <- matrix(c(rep(1:10,10),rep(c(1:5,11:15),each=10)),ncol=2,dimnames=list(paste0('cell',1:100),c('x','y')))
#' res <- paella(pt,cd)
#' mat <- matrix(rnorm(1000),nrow=10,dimnames=list(paste0('gene',1:10),paste0('cell',1:100)))
#' fit <- genefit(mat,res)

genefit <- function(mat,paella,gene=rownames(mat),interpoint=200) {
  pt <- paella$pseudotime
  cl <- paella$cluster

  mat <- mat[gene,names(cl)]
  
  predptseq <- seq(min(pt),max(pt),length.out=interpoint)
  pl <- list()
  for (i in unique(cl)) {
    tmpe <- mat[,cl==i]
    tmppt <- pt[colnames(tmpe)]
    pl[[i]] <- t(sapply(unname(gene),function(g) {
      predict(mgcv::gam(tmpe[g,]~s(tmppt,k=3)),data.frame(tmppt=predptseq))
    }))
    pl[[i]][,predptseq > max(tmppt) | predptseq < min(tmppt)] <- NA
    colnames(pl[[i]]) <- paste0(i,':',predptseq)
  }
  
  pl <- do.call(cbind,pl)
  ml <- rowMeans(pl,na.rm=T)
  sl <- apply(pl,1,sd,na.rm=T)
  pl <- (pl-ml)/sl
  list(scaledfit=pl,mean=ml,sd=sl)
}
