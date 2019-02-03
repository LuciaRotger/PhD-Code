library(Zseq)
library(gmp)
library(ape)
library(CollessLike)

big.double.factorial = function(n){
  if(n<2) return(1)
  m = (n+2+n%%2)/2
  return(Factorial.Double(m,odd=(n%%2==1))[m])
}

EUPhi = function(n){
  return(as.numeric((n*(n-1)/4)*(big.double.factorial(2*n-2)/
                                   big.double.factorial(2*n-3)-2)))
}

sapply(3:20,EUPhi)

exp.uni = c()
for(n in 3:8){ 
  trees=read.tree(file=paste("./bintrees_n",n,".txt",sep=""))
  indices = sapply(trees, cophen.index) 
  exp.uni[n]=mean(indices)
} 