#' Testing differential gene expression dynamics
#'
#' Testing if the temporal gene expression dynamics along pseudotime is differential across spatial sub-trajectories for each gene. The tests involve overall, mean, and trend differential.
#'
#' @param mat A gene expression matrix where each row is a gene and each column is a cell/spot. Row names are gene names and column names are cell/spot names. The column names must be present when fitting \code{paella()}.
#' @param paella Direct output from calling \code{paella()}.
#' @param genecut A numeric value between 0 and 1 indicating the gene filtering cutoff. Genes that have expression larger than median of all gene expression values in at least \code{genecut} cells are preserved.
#' @param ncores A positive integer of the number of cores to be used for parallel computing. If 1, no parallel computing will be done.
#' @import Matrix mgcv parallel
#' @export
#' @return A data frame of FDR and pval of overall, mean, and trend differential tests, as well as the type of differential for each gene.
#' @author Wenpin Hou<wh2526@@cumc.columbia.edu>, Zhicheng Ji
#' @examples 
#' pt <- setNames(1:100, paste0('cell',1:100))
#' cd <- matrix(c(rep(1:10,10),rep(c(1:5,11:15),each=10)),ncol=2,dimnames=list(paste0('cell',1:100),c('x','y')))
#' res <- paella(pt,cd)
#' mat <- matrix(rnorm(1000),nrow=10,dimnames=list(paste0('gene',1:10),paste0('cell',1:100)))
#' test <- difftest(mat,res)

difftest <- function(mat,paella,genecut=0.05,ncores=10) {
  clu <- as.factor(paella$cluster)
  pt <- paella$pseudotime
  
  pt <- pt[names(clu)]
  mat <- mat[,names(clu)]
  mede <- median(mat)
  gene <- unique(unlist(sapply(unique(clu),function(sc) {
    names(which(rowMeans(mat[,clu==sc] > mede) > genecut))
  },simplify = F)))
  
  if (ncores==1) {
    pval <- lapply(gene,function(i) {
      print(i)
      g <- mat[i,]
      m1 <- mgcv::gam(g~s(pt,k=3),method='REML')
      m2 <- mgcv::gam(g~clu+s(pt,k=3),method='REML')
      m3 <- mgcv::gam(g~clu+s(pt,by=clu,k=3),method='REML')
      p12 <- anova(m1,m2,test='Chisq')[['Pr(>Chi)']][2]
      p13 <- anova(m1,m3,test='Chisq')[['Pr(>Chi)']][2]
      p23 <- anova(m2,m3,test='Chisq')[['Pr(>Chi)']][2]
      c(p13,p12,p23)
    })
  } else {
    pval <- mclapply(gene,function(i) {
      print(i)
      g <- mat[i,]
      m1 <- mgcv::gam(g~s(pt,k=3),method='REML') # all subtrajectories share a same fitting curve
      m2 <- mgcv::gam(g~clu+s(pt,k=3),method='REML') # mean difference between subtrajectories
      m3 <- mgcv::gam(g~clu+s(pt,by=clu,k=3),method='REML') # mean + trend difference between subtrajectories
      p12 <- anova(m1,m2,test='Chisq')[['Pr(>Chi)']][2]
      p13 <- anova(m1,m3,test='Chisq')[['Pr(>Chi)']][2]
      p23 <- anova(m2,m3,test='Chisq')[['Pr(>Chi)']][2]
      c(p13,p12,p23)
    },mc.cores=ncores)
  }
  names(pval) <- gene
  pval <- do.call(rbind,pval)
  pval[is.na(pval)] <- 1
  colnames(pval) <- c('Overall_pval','Mean_pval','Trend_pval')
  overallfdr <- p.adjust(pval[,'Overall_pval'],method='fdr')
  sig <- names(which(overallfdr < 0.05))
  meanfdr <- trendfdr <- rep(NA,nrow(pval))
  names(meanfdr) <- names(trendfdr) <- rownames(pval)
  meanfdr[sig] <- p.adjust(pval[sig,'Mean_pval'],method='fdr')
  trendfdr[sig] <- p.adjust(pval[sig,'Trend_pval'],method='fdr')
  
  df <- data.frame(Overall_FDR=overallfdr,Mean_FDR=meanfdr,Trend_FDR=trendfdr,pval)
  df$Mean_pval[is.na(df$Mean_FDR)] <- NA
  df$Trend_pval[is.na(df$Trend_FDR)] <- NA
  df <- df[order(df[,'Overall_FDR'],df[,'Overall_pval']),]
  type <- rep('',nrow(df))
  names(type) <- rownames(df)
  type[df$Overall_FDR < 0.05] <- paste0('Overall',type[df$Overall_FDR < 0.05])
  type[rownames(df)[which(df$Mean_FDR < 0.05)]] <- paste0(type[rownames(df)[which(df$Mean_FDR < 0.05)]],',Mean')
  type[rownames(df)[which(df$Trend_FDR < 0.05)]] <- paste0(type[rownames(df)[which(df$Trend_FDR < 0.05)]],',Trend')
  type[type==''] <- 'None'
  df$Type <- paste0(type,' Differential')
  df
}

