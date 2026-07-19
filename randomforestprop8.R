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
  aprob <-1/(1 + exp(-((x[,1] > 0) + (x[,2] < -0.5) - 0.5*(x[,3] > 0 & x[,4] < 0) +
                         2 * (x[,4] > 0.5 & x[,5] < -0.5) -
                         2 * (x[,1] > 0.5 & x[,2] < 0.5) + 0.5 * x[,4] - 0.5*x[,5]^2 - 
                         0.5 * x[,6] * (x[,7]>0) + 0.5*x[,8]  + 0.25 * x[,9] ^ 2 - 0.25 * x[,10] ^ 2)))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  y <-5 * (x[,1]>0) + 5 * (x[,2] < -0.5) - 5*(x[,3] > 0 & x[,4] < 0) + 0.5 * x[,4] - 1*x[,5]^2 + 
    5 * (x[,4] > 0.5 & x[,5] < -0.5) - 5 * (x[,1] > 0.5 & x[,2] < 0.5) - 
    5 * x[,6] * (x[,7]>0 + 0.5 * x[,7]) + 0.5*x[,8] + 0.5 * x[,9] ^ 2 - 0.5 * x[,10] ^ 2 + rnorm(n, sd = sqrt(2))
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
respropinv81=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv81[r]=simrep_prop_rf(500,10)
}
set.seed(12345)
respropinv82=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv82[r]=simrep_prop_rf(500,50)
}
set.seed(12345)
respropinv83=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv83[r]=simrep_prop_rf(500,100)
}
set.seed(12345)
respropinv84=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv84[r]=simrep_prop_rf(1000,10)
}
set.seed(12345)
respropinv85=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv85[r]=simrep_prop_rf(1000,50)
}
set.seed(12345)
respropinv86=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv86[r]=simrep_prop_rf(1000,100)
}
set.seed(12345)
respropinv87=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv87[r]=simrep_prop_rf(2000,10)
}
set.seed(12345)
respropinv88=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv88[r]=simrep_prop_rf(2000,50)
}
set.seed(12345)
respropinv89=numeric(40)
nreps=40
for(r in 1:nreps){
  respropinv89[r]=simrep_prop_rf(2000,100)
}

write.table(respropinv81,"respropinv81",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv82,"respropinv82",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv83,"respropinv83",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv84,"respropinv84",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv85,"respropinv85",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv86,"respropinv86",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv87,"respropinv87",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv88,"respropinv88",row.names = F,quote = F,col.names = TRUE)
write.table(respropinv89,"respropinv89",row.names = F,quote = F,col.names = TRUE)
