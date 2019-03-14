library(Zseq)
library(gmp)
library(ape)
library(CollessLike)

# Option 1
tb.ape = read.tree(file = "·/tb_newicks.txt")
# Option 2
load("./treeBASE_database.RData")

bin.tb.ape=tb.ape[sapply(tb.ape,is.rooted)] 
bin.tb.ape=bin.tb.ape[sapply(bin.tb.ape,is.binary)] 
bin.tb.n = sapply(bin.tb.ape,Ntip)
leaves=as.numeric(names(which(table(bin.tb.n)>20)))
bin.tb.mean = c()
indices.tb = list()
for(k in leaves){
  trees = bin.tb.ape[bin.tb.n==k]
  indices.tb[[k]] = sapply(trees, cophen.index)
  value = mean( indices.tb[[k]] )
  bin.tb.mean = rbind(bin.tb.mean,c(k,value))
}

harm.sum = cumsum(1/1:1000)
EYPhi=function(n){
  return(n*(n+1-2*harm.sum[n]))
}
EUPhi = function(n){
  return(as.numeric((n*(n-1)/4)*(big.double.factorial(2*n-2)/
                                   big.double.factorial(2*n-3)-2)))
}
range.exp = 3:130
eyphi.values = sapply(range.exp, EYPhi)
euphi.values = sapply(range.exp, EUPhi)

plot(log(bin.tb.mean[,1]),log(bin.tb.mean[,2]),xlab="Log of number of leaves",
      ylab="Log of means",xlim=c(1.3,4.8),ylim=c(1,11))
lines(log(range.exp),log(eyphi.values),col="red")
lines(log(range.exp),log(euphi.values),col="blue")
legend("topleft", legend=c(expression(E[U]*(Phi[n])),
        expression(E[Y]*(Phi[n]))),col=c( "blue","red"),lty=1,cex=0.8) 
		