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

EUPhi = function(n){
    return(as.numeric((n*(n-1)/4)*(big.double.factorial(2*n-2)/
                                     big.double.factorial(2*n-3)-2)))
}

term.Phi = function(n){
  return(mul.bigq(as.bigq(n*(n-1)/2),(mul.bigq(as.bigq((49*n^3-57*n^2-22*n+24)/48),
                  big.double.factorial(2*n-4)/big.double.factorial(2*n-3))-
                    as.bigq((63*n^2-95*n+28)/30)))) 
}

compute.EUPhi2 = function(n.max=500){
  terms = lapply(2:n.max,term.Phi)
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
  write.table(exp.values,file = paste("E_U(Phi2)_",n.max,".txt",sep = ""),
              row.names = F,col.names = F)
  return(exp.values)
}

compute.varUPhi = function(exp.values,n.max=500){
  var.form = function(i)return(exp.values[i]-EUPhi(i)^2)
  var.values = sapply(1:n.max, var.form)
  write.table(var.values,file = paste("var_U(Phi)_",n.max,".txt",sep = ""),
              row.names = F,col.names = F)
  return(var.values)
}

exp.values.Phi = compute.EUPhi2(1000)
var.values.Phi = compute.varUPhi(exp.values.Phi,1000)

exp.values.Phi[3:20]
var.values.Phi[3:20]

summary(lm(log(var.values.Phi[900:1000])~log(900:1000)))

plot(log(1:1000),log(var.values.Phi),xlab="Log of the number of leaves",
     ylab="Log of the variance")
reg.phi=lm(log(var.values.Phi[500:1000])~log(500:1000))
abline(reg.phi,col="blue",lwd=2)

var.n = function (vec) return(var(vec)*(length(vec)-1)/length(vec))
trees = list()
all.cophen.index = list()
real.var.Phi = c()
for(n in 3:8){
  trees[[n]] = read.tree(file = paste("bintrees_n",n,".txt",sep = ""))
  all.cophen.index[[n]] = sapply(trees[[n]],cophen.index)
  real.var.Phi[n] = var.n(all.cophen.index[[n]])
  print(paste("var(Phi_",n,") = ",real.var.Phi[n],sep = ""))
}