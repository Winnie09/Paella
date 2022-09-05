#' Paella: spatial decomposition of cell trajectories
#'
#' Paella takes as input a vector of pseudotime and spatial coordinates of cells or spots. It then identifies multiple spatially distinct while temporally coextending spatial sub-trajectories.
#'
#' @param pt A named numeric vector of pseudotime. The names indicate different cells or spots.
#' @param cd A matrix of spatial coordinates. Each row is a cell or spot. It must have row names indicating the names of cells or spots, and the row names need to agree with the names in \code{pt}.
#' @param ncores A positive integer of the number of cores to be used for parallel computing. If 1, no parallel computing will be done.
#' @param mode A character vector. The elements must be 'both', 'small', or 'large'. If the length of the vector is larger than 1, Paella will automatically select the optimal mode.
#' @param seed A positive integer of random seed
#' @import geometry igraph parallel mgcv
#' @export
#' @return A list of spatial sub-trajectories assignment for each cell/spot, pseudotime, coordinates, and edges.
#' @author Wenpin Hou<wh2526@@cumc.columbia.edu>, Zhicheng Ji
#' @examples 
#' pt <- setNames(1:100, paste0('cell',1:100))
#' cd <- matrix(c(rep(1:10,10),rep(c(1:5,11:15),each=10)),ncol=2,dimnames=list(paste0('cell',1:100),c('x','y')))
#' res <- paella(pt,cd)

paella <- function(pt,cd,ncores=1,mode=c('both','small','large'),seed=10) {
  set.seed(seed)
  
  # set up delaunayn network
  edge <- delaunayn(cd)
  edge <- rbind(edge[ ,c(1,2)],edge[ ,c(1,3)],edge[ ,c(2,3)])
  edge <- unique(edge)
  edge <- rbind(edge,edge[,c(2,1)])
  edge <- cbind(rownames(cd)[edge[,1]],rownames(cd)[edge[,2]])
  
  # filter out isolated spots
  dm <- log(as.matrix(dist(cd)))
  
  dmedge <- dm[cbind(edge[,1],edge[,2])]
  edgvaluelist <- sort(unique(dmedge),decreasing = T)
  cc <- do.call(rbind,sapply(edgvaluelist,function(st) {
    comp <- components(graph_from_data_frame(data.frame(edge[dmedge <= st,]),directed=F))[[1]]
    data.frame(st,comp=sapply(1:max(comp),function(ci) {
      paste0(sort(names(which(comp==ci))),collapse = ':-:')
    }))
  },simplify = F))
  
  cc <- tapply(cc[,1],list(cc[,2]),min)
  cc <- sort(cc,decreasing = T) # important
  
  maxd <- matrix(0,nrow=length(pt),ncol=length(pt),dimnames = list(names(pt),names(pt)))
  for (i in 1:length(cc)) {
    cell <- strsplit(names(cc)[i],':-:')[[1]]
    maxd[cell,cell] <- unname(cc[i])
  }
  diag(maxd) <- 0
  ds <- rowMeans(maxd)  
  
  den <- density(ds)
  difden <- diff(den$y)
  selid <- den$x[which(difden[-length(difden)] < 0 & difden[-1] > 0)]
  selid <- c(min(ds)-1,selid,max(ds)+1)
  group <- cut(ds,selid)
  tab <- table(group)
  keepid <- tab/length(ds)
  if (min(keepid) < 0.01) {
    keepid <- names(keepid)[min(which(keepid < 0.01))-1]
    cut <- max(ds[group==keepid])
    cell <- names(ds)[ds <= cut]
    pt <- pt[cell]
    cd <- cd[cell,]
    edge <- edge[edge[,1] %in% cell & edge[,2] %in% cell,]
  }
  
  # filter out edges too long  
  distedge <- log(sqrt(rowSums((cd[edge[,1],]-cd[edge[,2],])^2)))
  
  den <- density(distedge)
  difden <- diff(den$y)
  selid <- den$x[which(difden[-length(difden)] < 0 & difden[-1] > 0)]
  selid <- c(min(distedge)-1,selid,max(distedge)+1)
  group <- cut(distedge,selid)
  tab <- table(group)
  keepid <- tab/length(distedge)
  if (min(keepid) < 0.01) {
    keepid <- names(keepid)[min(which(keepid < 0.01))-1]
    cut <- max(distedge[group==keepid])
    edge <- edge[distedge <= cut,]
  }
  
  fullpt <- pt
  fullcd <- cd
  fulledge <- edge
  
  # separately fit within each components
  comp <- components(graph_from_edgelist(edge,directed=F))[[1]]
  final <- sapply(unique(comp),function(sccomp) {
    cell <- names(which(comp==sccomp))
    pt <- fullpt[cell]
    cd <- fullcd[cell,]
    edge <- fulledge[edge[,1] %in% cell & edge[,2] %in% cell,]
    x <- cd[,1]
    y <- cd[,2]
    
    # smooth the pseudotime
    oript <- pt
    ptnei <- sapply(cell,function(i) {
      mean(pt[setdiff(as.vector(edge[edge[,1]==i|edge[,2]==i,]),i)])
    })
    
    ptdif <- ptnei-pt
    outcell <- names(which(abs(ptdif) > 2*sd(ptdif)))
    pt[outcell] <- rowMeans(cbind(ptnei[outcell],pt[outcell]))
    
    pt <- predict(gam(pt~te(x,y)))
    names(pt) <- cell
    
    # get tolerance value
    if (ncores==1) {
      ptaddval <- median(sapply(1:100,function(it) {
        permupt <- predict(gam(sample(oript)~te(x,y)))
        names(permupt) <- cell
        2*sd(permupt[edge[,1]]-permupt[edge[,2]])
      })) 
    } else {
      ptaddval <- median(unlist(mclapply(1:100,function(it) {
        permupt <- predict(gam(sample(oript)~te(x,y)))
        names(permupt) <- cell
        2*sd(permupt[edge[,1]]-permupt[edge[,2]])
      },mc.cores=ncores)))  
    }
    
    alldirfit <- sapply(mode,function(edgedir) {
      if (edgedir == 'both') {
        # construct the two graphs
        sg1 <- edge[pt[edge[,1]] < pt[edge[,2]] + ptaddval,]  
        sg2 <- edge[pt[edge[,1]] > pt[edge[,2]] - ptaddval,]  
        sg1 <- graph_from_edgelist(sg1,directed = T)
        sg2 <- graph_from_edgelist(sg2,directed = T)
        dism1 <- distances(sg1,mode='out')
        dism2 <- distances(sg2,mode='out')
        cell <- rownames(dism1)
        dism1 <- dism1[cell,cell]
        dism2 <- dism2[cell,cell]
        
        p <- mean(dism1 < Inf) * mean(dism2 < Inf)
        numcut <- nrow(dism1) * p + 2*sqrt(nrow(dism1) * p*(1-p))
        
        # sequentially identify parts
        k1 <- which(dism1 < Inf,arr.ind=T)
        k1 <- tapply(k1[,2],list(k1[,1]),function(i) i)
        k1 <- k1[sapply(k1,length) >= numcut]
        k2 <- which(dism2 < Inf,arr.ind=T)
        k2 <- tapply(k2[,2],list(k2[,1]),function(i) i)
        k2 <- k2[sapply(k2,length) >= numcut]
        
        exid <- expand.grid(1:length(k1),1:length(k2))
        exid <- cbind(exid[,1],exid[,2])
        if (ncores==1) {
          n <- sapply(1:nrow(exid),function(i) {
            intersect(k1[[exid[i,1]]],k2[[exid[i,2]]])
          },simplify = F)
        } else {
          n <- mclapply(1:nrow(exid),function(i) {
            intersect(k1[[exid[i,1]]],k2[[exid[i,2]]])
          },mc.cores=detectCores())  
        }
      } else {
        # construct the two graphs
        if (edgedir=='small') {
          sg1 <- edge[pt[edge[,1]] < pt[edge[,2]] + ptaddval,]    
        } else if (edgedir=='large') {
          sg1 <- edge[pt[edge[,1]] > pt[edge[,2]] - ptaddval,]   
        }
        sg1 <- graph_from_edgelist(sg1,directed = T)
        dism1 <- distances(sg1,mode='out')
        cell <- rownames(dism1)
        
        p <- mean(dism1 < Inf)
        numcut <- nrow(dism1) * p + 2*sqrt(nrow(dism1) * p*(1-p))
        
        # sequentially identify parts
        k1 <- which(dism1 < Inf,arr.ind=T)
        k1 <- tapply(k1[,2],list(k1[,1]),function(i) i)
        k1 <- k1[sapply(k1,length) >= numcut]
        n <- k1
      }
      
      if (length(n)==0) {
        NULL
      } else {
        len <- sapply(n,length)
        
        part <- list(n[[which.max(len)]])
        n <- n[len >= numcut]
        
        iterdiff <- length(part[[1]])
        while (iterdiff >= numcut) {
          cur <- unlist(part)
          if (ncores==1) {
            un <- sapply(n,function(i) {
              length(union(cur,i))
            })  
          } else {
            un <- unlist(mclapply(n,function(i) {
              length(union(cur,i))
            },mc.cores=ncores))
          }
          part[[length(part)+1]] <- n[[which.max(un)]]
          iterdiff <- length(unique(unlist(part))) - length(unique(unlist(part[-length(part)])))
        }
        part <- part[-length(part)]
        
        # convert to initial clustering
        part <- sapply(part,function(i) cell[i],simplify = F)
        singleid <- table(unlist(part))
        singleid <- names(singleid)[singleid==1]
        clu <- rep(0,length(cell))
        names(clu) <- cell
        for (i in 1:length(part)) clu[intersect(singleid,part[[i]])] <- i
        
        # resolve isolated cells
        for (i in part) {
          ce <- intersect(i,names(clu)[clu > 0])
          te <- edge[edge[,1] %in% ce & edge[,2] %in% ce,]
          cp <- components(graph_from_edgelist(te))[[1]]
          if (length(unique(cp)) > 1) {
            tab <- table(cp)
            smalliso <- setdiff(names(tab),names(tab)[which.max(tab)])
            selcell <- names(cp)[which(cp%in%as.numeric(smalliso))]
            clu[selcell] <- 0
          }
        }
        
        # resolve other cells
        dm <- distances(graph_from_data_frame(data.frame(edge),directed=F))
        outcell <- names(clu)[clu==0]
        singleid <- names(clu)[clu!=0]
        for (i in outcell) {
          clu[i] <- clu[names(which.min(dm[i,singleid]))]
        }
        list(cluster=clu,pseudotime=pt,coord=cd,edge=edge)
      }
    },simplify=F)
    alldirfit <- alldirfit[!sapply(alldirfit,is.null)]
    
    if (length(alldirfit) == 1) {
      alldirfit[[1]]
    } else {
      fromcell <- edge[,1]
      tocell <- edge[,2]
      fromcelllist <- sapply(cell,function(i) which(fromcell==i))
      
      if (ncores==1) {
        paw <- sapply(cell,function(tar) {
          topt <- abs(pt[tocell]-pt[tar])
          inc <- tar
          aed <- NULL
          candid <- which(fromcell == inc & tocell!=inc)
          while(length(inc) < length(cell)) {
            id <- candid[which.min(topt[candid])]
            tmpid <- c(fromcell[id],tocell[id])
            aed <- rbind(aed,tmpid)
            inc <- c(inc,tmpid[2])
            candid <- c(candid,fromcelllist[[tmpid[2]]])
            candid <- candid[!tocell[candid]%in%inc]
          }
          
          gr <- graph_from_data_frame(aed,directed=F)  
          path <- shortest_paths(gr,tar)[[1]]
          dp <- sapply(path,function(i) {
            diff(range(pt[names(i)]))
          })
          na <- sapply(path,function(i) {
            paste0(sort(c(names(i)[length(i)],tar)),collapse = ':-:')
          })
          names(dp) <- na
          dp
        },simplify = F,USE.NAMES = F)  
      } else {
        paw <- mclapply(cell,function(tar) {
          topt <- abs(pt[tocell]-pt[tar])
          inc <- tar
          aed <- NULL
          candid <- which(fromcell == inc & tocell!=inc)
          while(length(inc) < length(cell)) {
            id <- candid[which.min(topt[candid])]
            tmpid <- c(fromcell[id],tocell[id])
            aed <- rbind(aed,tmpid)
            inc <- c(inc,tmpid[2])
            candid <- c(candid,fromcelllist[[tmpid[2]]])
            candid <- candid[!tocell[candid]%in%inc]
          }
          
          gr <- graph_from_data_frame(aed,directed=F)  
          path <- shortest_paths(gr,tar)[[1]]
          dp <- sapply(path,function(i) {
            diff(range(pt[names(i)]))
          })
          na <- sapply(path,function(i) {
            paste0(sort(c(names(i)[length(i)],tar)),collapse = ':-:')
          })
          names(dp) <- na
          dp
        },mc.cores=detectCores())
      }
      
      paw <- do.call('c',paw)
      paw <- tapply(paw,list(names(paw)),min)
      paw <- paw[sub(':-:.*','',names(paw))!=sub('.*:-:','',names(paw))]
      
      ptdif <- abs(pt[sub(':-:.*','',names(paw))]-pt[sub('.*:-:','',names(paw))])
      paw <- paw-ptdif  
      
      oripaw <- paw
      paw <- log2(paw[paw > 0]+1)
      pawc1 <- sub(':-:.*','',names(paw))
      pawc2 <- sub('.*:-:','',names(paw))
      
      scfunc <- function(obj) {
        clu <- obj[[1]]
        edge <- obj[[4]]
        if (length(unique(clu))==1) {
          0
        } else {
          uniclupair <- unique(t(apply(cbind(clu[edge[,1]],clu[edge[,2]]),1,sort)))
          uniclupair <- uniclupair[uniclupair[,1]!=uniclupair[,2],,drop=F]
          sapply(1:nrow(uniclupair),function(i) {
            c1 <- uniclupair[i,1]
            c2 <- uniclupair[i,2]
            within <- paw[(clu[pawc1]==c1&clu[pawc2]==c1)|(clu[pawc1]==c2&clu[pawc2]==c2)]
            cut <- mean(within) + 2*sd(within)
            mean(paw[(clu[pawc1]==c1&clu[pawc2]==c2)|(clu[pawc1]==c2&clu[pawc2]==c1)] > cut)
          })
        }
      }
      
      score <- sapply(alldirfit,scfunc)
      score <- cbind(sapply(score,mean),sapply(score,length),sapply(score,min) > 0.05)
      score <- score[order(-score[,3],-score[,2],-score[,1]),,drop=F]
      alldirfit[[rownames(score)[1]]]
    }
  },simplify = F)
  if (length(final)==1) {
    final <- final[[1]]
    n <- names(final$cluster)
    final$cluster <- paste0('oripart_',final$cluster)
    names(final$cluster) <- n
  } else {
    for (i in 1:length(final)) {
      n <- names(final[[i]]$cluster)
      final[[i]]$cluster <- paste0('region_',i,':oripart_',final[[i]]$cluster)
      names(final[[i]]$cluster) <-  n
    }
    final <- list(cluster=unlist(sapply(final,function(i) i$cluster,simplify = F)),pseudotime=unlist(sapply(final,function(i) i$pseudotime,simplify = F)),coord=do.call(rbind,sapply(final,function(i) i$coord,simplify = F)),edge=do.call(rbind,sapply(final,function(i) i$edge,simplify = F)))
  }
  
  n <- names(sort(final$pseudotime))
  final$cluster <- final$cluster[n]
  final$pseudotime <- final$pseudotime[n]
  final$coord <- final$coord[n,]
  
  if (length(unique(final$cluster))==1) {
    n <- names(final$cluster)
    final$cluster <- rep('part1',length(n))
    names(final$cluster) <- n
  } else {
    fc <- final$coord
    fc[,1] <- fc[,1]-min(fc[,1])
    fc[,1] <- fc[,1]/max(fc[,1])
    fc[,2] <- fc[,2]-min(fc[,2])
    fc[,2] <- fc[,2]/max(fc[,2])
    cc <- aggregate(fc,list(final$cluster),mean)
    rownames(cc) <- cc[,1]
    cc <- as.matrix(cc[,-1])
    dm <- as.matrix(dist(cc))
    ord <- names(which.min(cc[,1]-cc[,2]))
    while(length(ord)<nrow(cc)) {
      tmp <- dm[ord,setdiff(rownames(cc),ord),drop=F]
      ord <- c(ord,colnames(tmp)[which(tmp==min(tmp),arr.ind=T)[1,2]])
    }
    for (i in 1:length(ord)) {
      final$cluster[final$cluster==ord[i]] <- paste0('traj',i)
    }
  }
  
  final
}

