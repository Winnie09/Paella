#' Cluster genes
#'
#' Assign clusters for genes.
#'
#' @param genefit Direct output from calling \code{genefit()}.
#' @param clunum A positive integer indicating number of clusters. If \code{NULL}, will be automatically selected by \code{findPC}.
#' @param maxclu A positive integer indicating the maximum number of clusters to be considered by \code{findPC}. Not used if \code{clunum} is not \code{NULL}.
#' @import findPC
#' @export
#' @return A list of gene clusters and hclust.
#' @author Wenpin Hou<wh2526@@cumc.columbia.edu>, Zhicheng Ji
#' @examples 
#' pt <- setNames(1:100, paste0('cell',1:100))
#' cd <- matrix(c(rep(1:10,10),rep(c(1:5,11:15),each=10)),ncol=2,dimnames=list(paste0('cell',1:100),c('x','y')))
#' res <- paella(pt,cd)
#' mat <- matrix(rnorm(1000),nrow=10,dimnames=list(paste0('gene',1:10),paste0('cell',1:100)))
#' fit <- genefit(mat,res)
#' clu <- clustergene(fit,maxclu=10)

clustergene <- function(genefit,clunum=NULL,maxclu=20) {
  genefit <- genefit[[1]]
  genefit <- genefit[,colMeans(is.na(genefit)) == 0,drop=F]
  hl <- hclust(dist(genefit))
  if (is.null(clunum)) {
    per <- sapply(1:maxclu,function(clunum) {
      cr <- cutree(hl,clunum)
      sum(sapply(unique(cr),function(i) {
        tmp <- genefit[cr==i,,drop=F]
        sum((t(tmp)-colMeans(tmp))^2)
      }))
    })/sum((t(genefit)-colMeans(genefit))^2)
    clunum <- findPC(sort(per,decreasing = T),number=maxclu)[1,1]  
  }
  list(cluster=cutree(hl,clunum),hclust=hl)
}
