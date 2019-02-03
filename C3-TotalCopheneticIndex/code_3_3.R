library(Zseq)
library(gmp)
library(ape)
library(CollessLike)

harm.sum = cumsum(1/1:1000)
EYPhi=function(n){
  return(n*(n+1-2*harm.sum[n]))
}

sapply(3:20,EYPhi)

yule.prob = function(tree){
  if (class(tree)=="phylo") 
    tree=graph.edgelist(tree$edge, directed=TRUE)  
  sp = shortest.paths(tree,mode = "out")
  deg = degree(tree,mode="out")
  leaves = which(deg==0)
  n = length(leaves) 
  k.node = function(node){
    subtree=which(sp[node,]<Inf)
    return(length(intersect(leaves,subtree)))
  } 
  kappas = sapply(which(deg>0), k.node) 
  value = (2^(n-1)/as.numeric(big.factorial(n)))*prod(1/(kappas-1))
  return(value)
}

exp.yule = c()
for(n in 3:8){ 
  trees=read.tree(file=paste("./bintrees_n",n,".txt",sep=""))
  indices = sapply(trees, cophen.index)
  probs=sapply(trees, yule.prob)
  exp.vule[n]=sum(indices*probs)
} 