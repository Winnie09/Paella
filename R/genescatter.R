#' Gene expression scatterplot
#'
#' Plotting the temporal gene expression dynamics across spatial sub-trajectories as scatterplots
#'
#' @param mat A gene expression matrix where each row is a gene and each column is a cell/spot. Row names are gene names and column names are cell/spot names. The column names must be present when fitting \code{paella()}.
#' @param genefit Direct output from calling \code{genefit()}.
#' @param paella Direct output from calling \code{paella()}.
#' @param gene Character vector of genes to be plotted
#' @param ncol A positive integer indicating the number of columns for facets.
#' @import ggplot2 reshape2 RColorBrewer
#' @export
#' @return A ggplot2 plot
#' @author Wenpin Hou<wh2526@@cumc.columbia.edu>, Zhicheng Ji
#' @examples 
#' pt <- setNames(1:100, paste0('cell',1:100))
#' cd <- matrix(c(rep(1:10,10),rep(c(1:5,11:15),each=10)),ncol=2,dimnames=list(paste0('cell',1:100),c('x','y')))
#' res <- paella(pt,cd)
#' mat <- matrix(rnorm(1000),nrow=10,dimnames=list(paste0('gene',1:10),paste0('cell',1:100)))
#' fit <- genefit(mat,res)
#' genescatter(mat,fit,res)

genescatter <- function(mat,genefit,paella,gene=rownames(genefit[[1]])[1:6],ncol=2) {
  genefit  <- (genefit[[1]][gene,]*genefit[[3]][gene])+genefit[[2]][gene]
  pt <- paella$pseudotime
  mat <- mat[gene,names(pt)]
  pd <- melt(mat)
  pd$pt <- pt[pd[,2]]
  pd$path <- paella$cluster[pd[,2]]
  
  ld <- melt(genefit[,colMeans(is.na(genefit))==0])
  ld$pt <- as.numeric(sub('.*:','',ld[,2]))
  ld$path <- sub(':.*','',ld$Var2)
  
  ggplot() + geom_point(data=pd,aes(x=pt,y=value,col=path),size=0.1) + geom_line(data=ld,aes(x=pt,y=value,col=path),size=1.5) + facet_wrap(~Var1,scales = 'free_y',ncol=ncol) + theme_classic() + scale_color_manual(values=brewer.pal(length(unique(paella$cluster)),'Set1')) + xlab('Pseudotime') + ylab('Expression') + theme(legend.title=element_blank())
}

