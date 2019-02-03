library(Zseq)
library(gmp)
library(ape)
library(CollessLike)

big.factorial = function(n){
  if(n<2) return(1)
  return(Factorial(n+1)[n+1])
}

big.double.factorial = function(n){
  if(n<2) return(1)
  m = (n+2+n%%2)/2
  return(Factorial.Double(m,odd=(n%%2==1))[m])
}

big.binomial = function(n,k){
  return(big.factorial(n)/(big.factorial(k)*big.factorial(n-k)))
}

Cknk = function(k,n){
  return(big.binomial(n,k)*((big.double.factorial(2*k-3)*
          big.double.factorial(2*(n-k)-3))/(2*big.double.factorial(2*n-3))))
}

EUS = function(n){
    return(as.numeric(n*(big.double.factorial(2*n-2)/big.double.factorial(2*n-3)-1)))
}

term.S = function(n){
  return((5*n*2^(n-2)*big.factorial(n))/(big.double.factorial(2*n-3))-n*(5*n-2)) 
}

compute.EUS2 = function(n.max=500){
  terms = lapply(2:n.max,term.S)
  terms = c(0,terms)
  exp.values = list(0)
  for(n in 2:n.max){
    sums = 0
    if(n>2){ 
      for(k in 2:(n-1)){ 
        sums = sums + Cknk(k,n)*exp.values[[k]]
      }
      sums = 2*sums
    }
    sums = sums + terms[[n]]
    exp.values[[n]] = sums
    print(n)
  }
  exp.values = sapply(exp.values, as.numeric)
  write.table(exp.values,file = paste("E_U(S2)_",n.max,".txt",sep = ""),
              row.names = F,col.names = F)
  return(exp.values)
}

compute.varUS = function(exp.values, n.max){
  var.form = function(i)return(exp.values[i]-EUS(i)^2)
  var.values = sapply(1:n.max,var.form)
  write.table(var.values,file = paste("var_U(S)_",n.max,".txt",sep = ""),
              row.names = F,col.names = F)
  return(var.values)
}

exp.values.S = compute.EUS2(1000)
var.values.S = compute.varUS(exp.values.S,1000)

exp.values.S[3:20]
var.values.S[3:20]

summary(lm(log(var.values.S[900:1000])~log(900:1000)))

plot(log(1:1000),log(var.values.S),xlab="Log of the number of leaves",
     ylab="Log of the variance")
reg.S=lm(log(var.values.S[500:1000])~log(500:1000))
abline(reg.S,col="red",lwd=2) 
sackin.approx = ((10-3*pi)/3)*(1:1000)^3 
lines(log(1:1000),log(sackin.approx),col="cyan",lty=2,lwd=2)

var.n = function (vec) return(var(vec)*(length(vec)-1)/length(vec))
trees = list()
all.sackin.index = list()
real.var.sackin = c()
for(n in 3:8){
  trees[[n]] = read.tree(file = paste("trees_n",n,".txt",sep = ""))
  all.sackin.index[[n]] = sapply(trees[[n]],sackin.index)
  real.var.sackin[n] = var.n(all.sackin.index[[n]])
  print(paste("var(S_",n,") = ",real.var.sackin[n],sep = ""))
}