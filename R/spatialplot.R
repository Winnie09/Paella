#' Paella spatial plotting
#'
#' Plotting the spatial sub-trajectories from Paella
#'
#' @param paella Direct output from calling \code{paella()}.
#' @param pointsize Size of the points shown in the plot.
#' @import ggplot2 igraph sp RColorBrewer
#' @export
#' @return A ggplot2 plot
#' @author Wenpin Hou<wh2526@@cumc.columbia.edu>, Zhicheng Ji
#' @examples 
#' pt <- setNames(1:100, paste0('cell',1:100))
#' cd <- matrix(c(rep(1:10,10),rep(c(1:5,11:15),each=10)),ncol=2,dimnames=list(paste0('cell',1:100),c('x','y')))
#' res <- paella(pt,cd)
#' spatialplot(res)

spatialplot <- function(paella,pointsize=1.5) {
  x <- paella$coord[,1]
  y <- paella$coord[,2]
  pt <- paella$pseudotime
  edge <- paella$edge
  edge <- edge[apply(edge,1,order)[1,]==1,]
  clu <- paella$cluster[names(x)]
  
  bound <- do.call(rbind,sapply(unique(clu),function(sc) {
    cell <- names(clu)[clu==sc]
    se <- edge[edge[,1]%in%cell & edge[,2]%in%cell,]
    nei <- sapply(cell,function(i) {
      setdiff(as.vector(se[se[,1]==i|se[,2]==i,]),i)
    })
    int <- apply(se,1,function(i) {
      cm <- intersect(nei[[i[1]]],nei[[i[2]]])
      if (length(cm) <= 1) {
        1
      } else {
        a <- (y[i[2]]-y[i[1]])/(x[i[2]]-x[i[1]])
        b <- y[i[1]]-a*x[i[1]]
        length(unique(sign(y[cm]-a*x[cm]-b)))==1
      }
    })
    se <- se[int==1,]
    
    secenter <- cbind(rowMeans(cbind(x[se[,1]],x[se[,2]])),rowMeans(cbind(y[se[,1]],y[se[,2]])))
    
    # remove enclosed edges using shortest path
    rml <- sapply(1:nrow(se),function(i) {
      if (se[i,1] %in% se[-i,] & se[i,2] %in% se[-i,]) {
        g <- graph_from_data_frame(se[-i,],directed=F)
        p <- names(get.shortest.paths(g,se[i,1],se[i,2])[[1]][[1]])
        if (length(p) > 0) {
          which(point.in.polygon(secenter[,1],secenter[,2],x[p],y[p])==1)
        }  
      }
    },simplify = F)
    rml <- unique(unlist(rml))
    se <- se[setdiff(1:nrow(se),rml),]
    secenter <- cbind(rowMeans(cbind(x[se[,1]],x[se[,2]])),rowMeans(cbind(y[se[,1]],y[se[,2]])))
    
    # remove enclosed edges using all possible paths
    rml <- sapply(1:nrow(se),function(i) {
      if (se[i,1] %in% se[-i,] & se[i,2] %in% se[-i,]) {
        g <- graph_from_data_frame(se[-i,],directed=F)
        p <- all_simple_paths(g,se[i,1],se[i,2])
        if (length(p) > 0) {
          unlist(sapply(p,function(k) {
            k <- names(k)
            which(point.in.polygon(secenter[,1],secenter[,2],x[k],y[k])==1)
          }))
        }
      }
    },simplify = F)
    rml <- unique(unlist(rml))
    se <- se[setdiff(1:nrow(se),rml),]
    
    # fix multi-way
    tab <- table(unlist(se))
    tar <- names(which(tab > 2))
    if (length(tar) > 0) {
      for (star in tar) {
        l <- setdiff(unlist(se[se[,1]==tar|se[,2]==tar,]),star)
        cbn <- combn(length(l),2)
        e <- cbind(l[cbn[1,]],l[cbn[2,]])
        
        center <- cbind(rowMeans(cbind(x[e[,1]],x[e[,2]])),rowMeans(cbind(y[e[,1]],y[e[,2]])))
        hull <- chull(x[l],y[l])
        inout <- point.in.polygon(center[,1],center[,2],x[l][hull],y[l][hull])  
        se <- se[!(se[,1] %in% c(l,star)&se[,2] %in% c(l,star)),]
        se <- rbind(se,e[inout==2,])
      }
      secenter <- cbind(rowMeans(cbind(x[se[,1]],x[se[,2]])),rowMeans(cbind(y[se[,1]],y[se[,2]])))
      
      rml <- sapply(1:nrow(se),function(i) {
        if (se[i,1] %in% se[-i,] & se[i,2] %in% se[-i,]) {
          g <- graph_from_data_frame(se[-i,],directed=F)
          p <- all_simple_paths(g,se[i,1],se[i,2])
          if (length(p) > 0) {
            unlist(sapply(p,function(k) {
              k <- names(k)
              which(point.in.polygon(secenter[,1],secenter[,2],x[k],y[k])==1)
            }))
          }
        }
      },simplify = F)
      rml <- unique(unlist(rml))
      se <- se[setdiff(1:nrow(se),rml),]
    }
    data.frame(se,sc)
  },simplify = F))
    
  cv <- c('#FEFB53','#FEFB53','#F3AF36','#EC7D24','#DE5937','#E74040','#DA3E88','#A643E8','#745FE8','#6D8EF7','#1F66DB','#003291','#050673','#050673')

  cv2 <- brewer.pal(length(unique(bound[,3])),'Set1')
  names(cv2) <- unique(bound[,3])
  ggplot() + geom_point(data=data.frame(x,y,pseudotime=pt),aes(x=x,y=y,fill=pseudotime),shape=21,col='white',size=pointsize,alpha=0.6)+ geom_segment(data=data.frame(x=x[match(bound[,1],names(x))],xend=x[match(bound[,2],names(x))],y=y[match(bound[,1],names(x))],yend=y[match(bound[,2],names(x))],cluster=bound[,3]),aes(x=x,y=y,xend=xend,yend=yend,col=cluster),size=1) + scale_fill_gradientn(colors=cv) + theme_void() + theme(legend.position='bottom',legend.box="vertical") + scale_color_discrete(name = NULL,type=cv2) + guides(colour = guide_legend(order = 1,byrow=T))+ theme(legend.text=element_text(size=11))
}

