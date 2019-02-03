library(Zseq)
library(gmp)
library(ape)
library(CollessLike)

balance.indices2 = function(tree){
  colless.coefficient = (log(0 + exp(1))+log(2 + exp(1)))/2 
  values = balance.indices(tree)
  values[1] = values[1]/colless.coefficient
  return(values)
}

are.tie = function(xx,yy) return(xx==yy) 

exact.ties = function(){
  trees = list()
  all.indices = list()
  num.ties = c(0,0,0) 
  prob.ties = list()
  for(n in 3:7){
    trees[[n]] = read.tree(file = paste("bintrees_n",n,".txt",sep=""))
    total.trees = length(trees[[n]])
    total.pairs = total.trees*(total.trees-1)/2
    all.indices[[n]] = matrix(sapply(trees[[n]],balance.indices2),ncol=3,byrow=T) 
    num.ties[1]=sum(outer(all.indices[[n]][,1],all.indices[[n]][,1],are.tie))
    num.ties[2]=sum(outer(all.indices[[n]][,2],all.indices[[n]][,2],are.tie))
    num.ties[3]=sum(outer(all.indices[[n]][,3],all.indices[[n]][,3],are.tie))
    num.ties = (num.ties-total.trees)/2
    prob.ties[[n]] = num.ties/total.pairs
    print(paste("Ties for n =",n," : ",
                paste(c("p_C=","p_S=","p_Phi"),round(prob.ties[[n]],4),collapse=", "),sep=""))
  }
  return(prob.ties)
}

sim.ties.n = function(n,num.pairs.sim=3000){ 
  num.ties = c(0,0,0)
  for(i in 1:num.pairs.sim){
    t1 = rtree(n,rooted=TRUE)
    continue = TRUE
    while(continue){
      t2 = rtree(n,rooted=TRUE)
      continue = all.equal(t1,t2,use.length=FALSE,use.tip.label=FALSE)
    }
    t1.indices = balance.indices2(t1)
    t2.indices = balance.indices2(t2)
    num.ties = num.ties + (t1.indices==t2.indices)
  } 
  print(paste("n =",n))
  print(paste("Ties :",num.ties))
  prob.ties = num.ties/num.pairs.sim
  print(paste("Prob :",round(prob.ties,4)))
  return(prob.ties)
}

ties.1 = exact.ties()
ties.2 = lapply(8:50, sim.ties.n,num.pairs.sim=3000)
ties = matrix(c(unlist(ties.1),unlist(ties.2)),ncol=3,byrow=TRUE)
colnames(ties)=c("Colless","Sackin","Cophenetic")
rownames(ties)=3:50
ties[1:18,]

plot(log(3:50),log(ties[,3]),type="l",
     xlab="log of the number of leaves",
     ylab="log of the probability of tie", col="blue")
lines(log(3:50),log(ties[,1]),col="red")
lines(log(3:50),log(ties[,2]),col="green")
legend("bottomleft", legend=c("Sackin","Colless","Cophenetic"),
       col=c("green","red", "blue"),lty=1, cex=0.8)
