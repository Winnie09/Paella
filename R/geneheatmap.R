#' Gene expression heatmap
#'
#' Plotting the temporal gene expression dynamics across spatial sub-trajectories as a heatmap
#'
#' @param genefit Direct output from calling \code{genefit()}.
#' @param geneclu Direct output from calling \code{geneclu()}.
#' @param diff Direct output from calling \code{difftest()}.
#' @import ComplexHeatmap circlize
#' @export
#' @return A ComplexHeatmap plot
#' @author Wenpin Hou<wh2526@@cumc.columbia.edu>, Zhicheng Ji
#' @examples 
#' pt <- setNames(1:100, paste0('cell',1:100))
#' cd <- matrix(c(rep(1:10,10),rep(c(1:5,11:15),each=10)),ncol=2,dimnames=list(paste0('cell',1:100),c('x','y')))
#' res <- paella(pt,cd)
#' mat <- matrix(rnorm(1000),nrow=10,dimnames=list(paste0('gene',1:10),paste0('cell',1:100)))
#' test <- difftest(mat,res)
#' fit <- genefit(mat,res)
#' clu <- clustergene(fit,maxclu=10)
#' geneheatmap(fit,clu,test)

geneheatmap <- function(genefit,geneclu,diff) {
  ord <- geneclu$hclust$order
  geneclu  <- geneclu$cluster[ord]
  pl <- genefit$scaledfit[ord,]
  
  split <- sub(':.*','',colnames(pl))
  pt <- as.numeric(sub('.*:','',colnames(pl)))
  
  ptcol_fun = colorRamp2(c(0, max(pt)), c("white", "green"))
  set.seed(1)
  clv <- sample(rainbow(length(unique(geneclu)),alpha=0.5))
  names(clv) <- 1:max(geneclu)
  rann <- rowAnnotation(cluster=as.character(geneclu),mean_diff=grepl('Mean',diff[rownames(pl),'Type']),trend_diff=grepl('Trend',diff[rownames(pl),'Type']),col = list(cluster=clv,mean_diff = c('TRUE'='red','FALSE'='grey'),trend_diff = c('TRUE'='red','FALSE'='grey')),show_legend=c(F,T,T))
  cann <- columnAnnotation(pseudotime=pt,col = list(pseudotime=ptcol_fun),annotation_name_side = "left")
  
  col_fun = colorRamp2(c(quantile(pl,0.05,na.rm=T), 0, quantile(pl,0.95,na.rm=T)), c("blue", "white", "red"))
  
  ht <- Heatmap(pl,name='scaled expression',cluster_rows = F,cluster_columns=F,show_row_names = F, show_column_names = F,left_annotation = rann, top_annotation = cann, row_split = geneclu,column_split = split,na_col='grey95',heatmap_legend_param = list(legend_direction = "horizontal"),use_raster=F,col=col_fun)
  draw(ht, heatmap_legend_side="bottom", annotation_legend_side="right",legend_grouping = "original",padding = unit(c(10, 4, 4, 4), "mm"))
}


