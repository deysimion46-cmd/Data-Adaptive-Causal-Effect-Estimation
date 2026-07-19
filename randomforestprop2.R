library(grf)
library(ranger)
## for other weighting methods
library(WeightIt)

library(fields)

## needed for energy_balance only if solver="cccp"
#library(optiSolve)

## for checking balance
library(cobalt)
library(MASS)

simrep_prop_rf=function(n=250,p=30){
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = F)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob  <- .25*(1+dbeta(x[,1],2,4))
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  y <-2*(x[,1]-1)+.5*(2*a-1)*(1+1/(1+exp(-20*(x[,1]-1/3))))*(1+1/(1+exp(-20*(x[,2]-1/3))))+rnorm(n)
  
  dfxa  <- data.frame(a = a, x = x)
  formula_a <- as.formula(paste0("a ~ ", paste(colnames(dfxa)[-1], collapse = " + ")))
  rf_fit_a <- ranger(formula_a, data = dfxa,
                     num.trees = 1000)
  prop=rf_fit_a$predictions
  prop[a==0]=1-prop[a==0]
  invprop=1/prop
  return(weighted.mean(y[a==1],w=invprop[a==1])-weighted.mean(y[a==0],w=invprop[a==0]))
}
set.seed(12345)
respropinv21=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv21[r]=simrep_prop_rf(500,10)
}
set.seed(12345)
respropinv22=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv22[r]=simrep_prop_rf(500,50)
}
set.seed(12345)
respropinv23=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv23[r]=simrep_prop_rf(500,100)
}
set.seed(12345)
respropinv24=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv24[r]=simrep_prop_rf(1000,10)
}
set.seed(12345)
respropinv25=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv25[r]=simrep_prop_rf(1000,50)
}
set.seed(12345)
respropinv26=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv26[r]=simrep_prop_rf(1000,100)
}
set.seed(12345)
respropinv27=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv27[r]=simrep_prop_rf(2000,10)
}
set.seed(12345)
respropinv28=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv28[r]=simrep_prop_rf(2000,50)
}
set.seed(12345)
respropinv29=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv29[r]=simrep_prop_rf(2000,100)
}

write.table(respropinv21,"respropinv21",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv22,"respropinv22",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv23,"respropinv23",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv24,"respropinv24",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv25,"respropinv25",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv26,"respropinv26",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv27,"respropinv27",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv28,"respropinv28",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv29,"respropinv29",row.names = F,quote = F,col.names = TRUE)
