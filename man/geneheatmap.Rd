% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneheatmap.R
\name{geneheatmap}
\alias{geneheatmap}
\title{Gene expression heatmap}
\usage{
geneheatmap(genefit, geneclu, diff)
}
\arguments{
\item{genefit}{Direct output from calling \code{genefit()}.}

\item{geneclu}{Direct output from calling \code{geneclu()}.}

\item{diff}{Direct output from calling \code{difftest()}.}
}
\value{
A ComplexHeatmap plot
}
\description{
Plotting the temporal gene expression dynamics across spatial sub-trajectories as a heatmap
}
\examples{
pt <- setNames(1:100, paste0('cell',1:100))
cd <- matrix(c(rep(1:10,10),rep(c(1:5,11:15),each=10)),ncol=2,dimnames=list(paste0('cell',1:100),c('x','y')))
res <- paella(pt,cd)
mat <- matrix(rnorm(1000),nrow=10,dimnames=list(paste0('gene',1:10),paste0('cell',1:100)))
test <- difftest(mat,res)
fit <- genefit(mat,res)
clu <- clustergene(fit,maxclu=10)
geneheatmap(fit,clu,test)
}
\author{
Wenpin Hou<wh2526@cumc.columbia.edu>, Zhicheng Ji
}
