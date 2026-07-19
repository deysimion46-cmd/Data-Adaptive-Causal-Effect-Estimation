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
  aprob <-1/(1 + exp(-.25*(apply(x[,c(c(1:10),c(21:30))],1,sum))))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  y <- apply(x[,1:20],1,sum) + rnorm(n, sd = sqrt(2))
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
respropinv141=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv141[r]=simrep_prop_rf(500,30)
}
set.seed(12345)
respropinv142=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv142[r]=simrep_prop_rf(500,60)
}
set.seed(12345)
respropinv143=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv143[r]=simrep_prop_rf(500,100)
}
set.seed(12345)
respropinv144=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv144[r]=simrep_prop_rf(1000,30)
}
set.seed(12345)
respropinv145=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv145[r]=simrep_prop_rf(1000,60)
}
set.seed(12345)
respropinv146=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv146[r]=simrep_prop_rf(1000,100)
}
set.seed(12345)
respropinv147=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv147[r]=simrep_prop_rf(2000,30)
}
set.seed(12345)
respropinv148=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv148[r]=simrep_prop_rf(2000,60)
}
set.seed(12345)
respropinv149=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv149[r]=simrep_prop_rf(2000,100)
}

write.table(respropinv141,"respropinv141",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv142,"respropinv142",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv143,"respropinv143",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv144,"respropinv144",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv145,"respropinv145",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv146,"respropinv146",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv147,"respropinv147",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv148,"respropinv148",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv149,"respropinv149",row.names = F,quote = F,col.names = TRUE)
