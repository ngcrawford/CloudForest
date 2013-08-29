
##################################################################################
#
#    calculate the distance of two sequences on the basis of 
#   number of ancestors between two sequences
#
##################################################################################
nancdist<-function(tree, taxaname)
{
    ntaxa<-length(taxaname)
    nodematrix<-read.tree.nodes(tree,taxaname)$nodes
    if(is.rootedtree(nodematrix)) nodematrix<-unroottree(nodematrix)
    dist<-matrix(0, ntaxa,ntaxa)
    for(i in 1:(ntaxa-1))
        for(j in (i+1):ntaxa)
        {
        anc1<-ancestor(i,nodematrix)
        anc2<-ancestor(j,nodematrix)
        n<-sum(which(t(matrix(rep(anc1,length(anc2)),ncol=length(anc2)))-anc2==0, arr.ind=TRUE)[1,])-3
        if(n==-1) n<-0
        dist[i,j]<-n
        }
    dist<-dist+t(dist)
    z<-list(dist=as.matrix, taxaname=as.vector)
    z$dist<-dist
    z$taxaname<-taxaname
    z
}

NJst<-function(genetrees, taxaname, spname, species.structure)
{
    ntree<-length(genetrees)
    ntaxa<-length(taxaname)
    dist <- matrix(0, nrow = ntree, ncol = ntaxa * ntaxa)
    
    for(i in 1:ntree)
    {
        genetree1 <- read.tree.nodes(genetrees[i])
            thistreetaxa <- genetree1$names
            ntaxaofthistree <- length(thistreetaxa)
            thistreenode <- rep(-1, ntaxaofthistree)
        dist1<-matrix(0,ntaxa,ntaxa)
            for (j in 1:ntaxaofthistree) 
        {
                    thistreenode[j] <- which(taxaname == thistreetaxa[j])
                    if (length(thistreenode[j]) == 0) 
            {
                        print(paste("wrong taxaname", thistreetaxa[j],"in gene", i))
                        return(0)
                    }
            }
        dist1[thistreenode, thistreenode]<-nancdist(genetrees[i],thistreetaxa)$dist
        dist[i,]<-as.numeric(dist1)
    }

    dist[dist == 0] <- NA
        dist2 <- matrix(apply(dist, 2, mean, na.rm = TRUE), ntaxa, ntaxa)
        diag(dist2) <- 0
        if (sum(is.nan(dist2)) > 0) 
    {
            print("missing species!")
            dist2[is.nan(dist2)] <- 10000
        }
        speciesdistance <- pair.dist.mulseq(dist2, species.structure)

    tree<-write.tree(nj(speciesdistance))
    node2name(tree,name=spname)
}

#############################################################################
#
#        generate gene trees from a species tree
#
#############################################################################
# library(phybase)
# sptree<-"(A:0.004,(B:0.003,(C:0.002,(D:0.001,E:0.001):0.001):0.001):0.001);"

# spname<-species.name(sptree)
# nspecies<-length(spname)
# rootnode<-9
# nodematrix<-read.tree.nodes(sptree,spname)$node
# seq<-rep(1,nspecies)
# species.structure<-matrix(0,nspecies,nspecies)
# diag(species.structure)<-1

# ##population size, theta
# nodematrix[,5]<-0.01

# ngene<-5
# genetree<-rep("",ngene)

# ##generate gene trees
# for(i in 1:ngene)
# {
#   genetree[i]<-sim.coaltree.sp(rootnode,nodematrix,nspecies,seq,spname)$gt
# }
    
# ##construct the NJst tree
# NJst(genetree,spname, spname, species.structure)
