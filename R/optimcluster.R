
# require essential packages
#library(hustlxc)
#library(MASS)
#library(cluster)
#library(randomForest)
#library(clusterSim)
#library(plotrix)
#library(fpc)

#x=iris[,1:4]
#y=iris[,5]

# here I use random forest to calculate proximity matrix
#rf=randomForest(x,proximity=TRUE,importance=TRUE)

# get distance object D, however, the distance matrix can
#+be computed from dist routine.
#D=as.dist(1-rf$proximity)

#load('STAD.Rdata')
#x=tot_num
#D=dist(tot_num)

# Calinski-Harabasz index copied from package fpc
calinhara <- function (x, clustering, cn = max(clustering)) 
{
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    cln <- rep(0, cn)
    W <- matrix(0, p, p)
    for (i in 1:cn) cln[i] <- sum(clustering == i)
    for (i in 1:cn) {
        clx <- x[clustering == i, ]
        cclx <- cov(as.matrix(clx))
        if (cln[i] < 2) 
            cclx <- 0
        W <- W + ((cln[i] - 1) * cclx)
    }
    S <- (n - 1) * cov(x)
    B <- S - W
    out <- (n - cn) * sum(diag(B))/((cn - 1) * sum(diag(W)))
    out
}


optimcluster <- function(x, D=dist(x), plot=FALSE)
{

# pre-defined number of clusters
ncluster=20

# initialize array of metric
sil.idx=rep(0,ncluster)
CI.idx=rep(0,ncluster)
sil.idx[1]=0
CI.idx[1]=0

# calculating clustering metrics, i.e. Silhouete band-wdith and CI index
for (i in 2:ncluster) {
	p=pam(D, k=i)
	sil.idx[i]=mean(silhouette(p)[,3])       # in package cluster
	#CI.idx[i]=index.G1(x, p$clustering, D)   # in package clusterSim
	CI.idx[i]=calinhara(x, p$clustering)     # in package fpc
	message('nclust=',i, ', silhouette index=',sil.idx[i], ', Calinski-Harabasz index=', CI.idx[i])
}

# visualize clustering metrics
if (plot) {
#pdf(fig, height=5, width=6)
twoord.plot(1:20,sil.idx,1:20,CI.idx, ylab='Average Silhouete', rylab='Calinski-Harabasz (CH) Index', xlab='Cluster number', main='Optimal cluster selection')
#dev.off() 
}

# the best partition is derived from the above figure, which is one that
#+maximize CI index and/or silhouette band-width, sometimes arbitrarily.
optim.nclust=which.max(sil.idx)             # select the optimal cluster number via average silhouette
p=pam(D,k=optim.nclust)

# Classical multidimensional scaling (CMDS) of a data matrix, or you can choose
#+to use non-metric multidimensional scaling (NMDS) provided in MASS package
#mds=cmdscale(D, 2, eig=TRUE)                  # CMDS
#mds=isoMDS(D, k=2, maxit=5000, tol=1e-6)     # NMDS

# visualize MDS output
#pdf('clustering_mds_figure.pdf', height=5, width=6)
#s.class2(mds$points, as.factor(p$clustering), grid=FALSE, clabel=0.5, cpoint=0, ellipse.lwd=1)
#dev.off()

list(optim.pam=p, optim.nclust=optim.nclust, silhouette.index=sil.idx, Calinski.Harabasz.index=CI.idx)

}

