library(forestBalance)
set.seed(12345)
result_1000_10=numeric(40)
n=1000
p=10
for(i in 1:40){
rho <- -0.25
Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))

x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)

## scramble columns of x so that which pairs of columns
## are correlated is at random
x <- x[,sample.int(p, p, replace = FALSE)]

## true propensity scores..
## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
aprob <- .25*(1+dbeta(x[,1],2,4))

## generate treatment assignments
a <- rbinom(n, 1, aprob)

## generate response with confounding.
## here the true ATE (treatment effect) is 0

y <- 2*(x[,1]-1)+.5*(2*a-1)*(1+1/(1+exp(-20*(x[,1]-1/3))))*(1+1/(1+exp(-20*(x[,2]-1/3))))+rnorm(n)


fit=forest_balance(x,a,y,cross.fitting = TRUE)
result_500_10[i]=fit$ate

}
set.seed(12345)
result_1000_100=numeric(40)
n=1000
p=100
for(i in 1:40){
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- .25*(1+dbeta(x[,1],2,4))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  
  y <- 2*(x[,1]-1)+.5*(2*a-1)*(1+1/(1+exp(-20*(x[,1]-1/3))))*(1+1/(1+exp(-20*(x[,2]-1/3))))+rnorm(n)
  
  
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_500_100[i]=fit$ate
  }
set.seed(12345)
result_2000_10=numeric(40)
n=2000
p=10
for(i in 1:40){
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- .25*(1+dbeta(x[,1],2,4))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  
  y <- 2*(x[,1]-1)+.5*(2*a-1)*(1+1/(1+exp(-20*(x[,1]-1/3))))*(1+1/(1+exp(-20*(x[,2]-1/3))))+rnorm(n)
  
  
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_2000_10[i]=fit$ate
}
set.seed(12345)
result_2000_100=numeric(40)
n=2000
p=100
for(i in 1:40){
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- .25*(1+dbeta(x[,1],2,4))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  
  y <- 2*(x[,1]-1)+.5*(2*a-1)*(1+1/(1+exp(-20*(x[,1]-1/3))))*(1+1/(1+exp(-20*(x[,2]-1/3))))+rnorm(n)
  
  
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_2000_100[i]=fit$ate
}


result_4000_10=numeric(40)
n=4000
p=10
for(i in 1:40){
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- .25*(1+dbeta(x[,1],2,4))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  
  y <- 2*(x[,1]-1)+.5*(2*a-1)*(1+1/(1+exp(-20*(x[,1]-1/3))))*(1+1/(1+exp(-20*(x[,2]-1/3))))+rnorm(n)
  
  
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_4000_10[i]=fit$ate
}
result_4000_100=numeric(40)
n=4000
p=100
for(i in 1:40){
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- .25*(1+dbeta(x[,1],2,4))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  
  y <- 2*(x[,1]-1)+.5*(2*a-1)*(1+1/(1+exp(-20*(x[,1]-1/3))))*(1+1/(1+exp(-20*(x[,2]-1/3))))+rnorm(n)
  
  
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_4000_100[i]=fit$ate
}




res_mat=read.table("resmat21",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv21",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_500_10)
res_mat=res_mat-1.82
g1=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + theme_bw(base_size = 18)+labs(x=NULL)+
  theme(axis.text.x = element_blank())+ggtitle("n=500,p=10")
res_mat=read.table("resmat23",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv23",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_500_100)
res_mat=res_mat-1.82
#2nd plot
g3=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + labs(x=NULL)+ theme_bw(base_size = 18)+
  theme(axis.text.x = element_blank())+labs(y=NULL)+ggtitle("n=500,p=100")

res_mat=read.table("resmat27",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv27",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_2000_10)
res_mat=res_mat-1.82
g7=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle("n=2000,p=10")

res_mat=read.table("resmat29",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv29",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_2000_100)
res_mat=res_mat-1.82
g9=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1) )+labs(y=NULL)+ggtitle("n=2000,p=100")

rplot2=g1+g3+g7+g9
ggsave("Model23.pdf",plot = rplot2,device = "pdf",height = 9,width = 11)

set.seed(12345)
result_2_500_10=numeric(40)
n=500
p=10
for(i in 1:40){
rho <- -0.25
Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))

x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)

## scramble columns of x so that which pairs of columns
## are correlated is at random
x <- x[,sample.int(p, p, replace = FALSE)]

## true propensity scores..
## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
aprob <- 1/(1 + exp(-((x[,1] > 0) + (x[,2] < -0.5) - 0.5*(x[,3] > 0 & x[,4] < 0) +
                        2 * (x[,4] > 0.5 & x[,5] < -0.5) -
                        2 * (x[,1] > 0.5 & x[,2] < 0.5) + 0.5 * x[,4] - 0.5*x[,5]^2 - 
                        0.5 * x[,6] * (x[,7]>0) + 0.5*x[,8]  + 0.25 * x[,9] ^ 2 - 0.25 * x[,10] ^ 2)))

## generate treatment assignments
a <- rbinom(n, 1, aprob)

## generate response with confounding.
## here the true ATE (treatment effect) is 0
y <- 5 * (x[,1]>0) + 5 * (x[,2] < -0.5) - 5*(x[,3] > 0 & x[,4] < 0) + 0.5 * x[,4] - 1*x[,5]^2 + 
  5 * (x[,4] > 0.5 & x[,5] < -0.5) - 5 * (x[,1] > 0.5 & x[,2] < 0.5) - 
  5 * x[,6] * (x[,7]>0 + 0.5 * x[,7]) + 0.5*x[,8] + 0.5 * x[,9] ^ 2 - 0.5 * x[,10] ^ 2 + rnorm(n, sd = sqrt(2))

fit=forest_balance(x,a,y,cross.fitting = TRUE)
result_2_500_10[i]=fit$ate


}


set.seed(12345)
result_2_500_100=numeric(40)
n=500
p=100
for(i in 1:40){
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- 1/(1 + exp(-((x[,1] > 0) + (x[,2] < -0.5) - 0.5*(x[,3] > 0 & x[,4] < 0) +
                          2 * (x[,4] > 0.5 & x[,5] < -0.5) -
                          2 * (x[,1] > 0.5 & x[,2] < 0.5) + 0.5 * x[,4] - 0.5*x[,5]^2 - 
                          0.5 * x[,6] * (x[,7]>0) + 0.5*x[,8]  + 0.25 * x[,9] ^ 2 - 0.25 * x[,10] ^ 2)))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  y <- 5 * (x[,1]>0) + 5 * (x[,2] < -0.5) - 5*(x[,3] > 0 & x[,4] < 0) + 0.5 * x[,4] - 1*x[,5]^2 + 
    5 * (x[,4] > 0.5 & x[,5] < -0.5) - 5 * (x[,1] > 0.5 & x[,2] < 0.5) - 
    5 * x[,6] * (x[,7]>0 + 0.5 * x[,7]) + 0.5*x[,8] + 0.5 * x[,9] ^ 2 - 0.5 * x[,10] ^ 2 + rnorm(n, sd = sqrt(2))
  
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_2_500_100[i]=fit$ate
  
  
}


set.seed(12345)
result_2_2000_10=numeric(40)
n=2000
p=10
for(i in 1:40){
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- 1/(1 + exp(-((x[,1] > 0) + (x[,2] < -0.5) - 0.5*(x[,3] > 0 & x[,4] < 0) +
                          2 * (x[,4] > 0.5 & x[,5] < -0.5) -
                          2 * (x[,1] > 0.5 & x[,2] < 0.5) + 0.5 * x[,4] - 0.5*x[,5]^2 - 
                          0.5 * x[,6] * (x[,7]>0) + 0.5*x[,8]  + 0.25 * x[,9] ^ 2 - 0.25 * x[,10] ^ 2)))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  y <- 5 * (x[,1]>0) + 5 * (x[,2] < -0.5) - 5*(x[,3] > 0 & x[,4] < 0) + 0.5 * x[,4] - 1*x[,5]^2 + 
    5 * (x[,4] > 0.5 & x[,5] < -0.5) - 5 * (x[,1] > 0.5 & x[,2] < 0.5) - 
    5 * x[,6] * (x[,7]>0 + 0.5 * x[,7]) + 0.5*x[,8] + 0.5 * x[,9] ^ 2 - 0.5 * x[,10] ^ 2 + rnorm(n, sd = sqrt(2))
  
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_2_2000_10[i]=fit$ate
  
  
}



set.seed(12345)
result_2_2000_100=numeric(40)
n=2000
p=100
for(i in 1:40){
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- 1/(1 + exp(-((x[,1] > 0) + (x[,2] < -0.5) - 0.5*(x[,3] > 0 & x[,4] < 0) +
                          2 * (x[,4] > 0.5 & x[,5] < -0.5) -
                          2 * (x[,1] > 0.5 & x[,2] < 0.5) + 0.5 * x[,4] - 0.5*x[,5]^2 - 
                          0.5 * x[,6] * (x[,7]>0) + 0.5*x[,8]  + 0.25 * x[,9] ^ 2 - 0.25 * x[,10] ^ 2)))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  y <- 5 * (x[,1]>0) + 5 * (x[,2] < -0.5) - 5*(x[,3] > 0 & x[,4] < 0) + 0.5 * x[,4] - 1*x[,5]^2 + 
    5 * (x[,4] > 0.5 & x[,5] < -0.5) - 5 * (x[,1] > 0.5 & x[,2] < 0.5) - 
    5 * x[,6] * (x[,7]>0 + 0.5 * x[,7]) + 0.5*x[,8] + 0.5 * x[,9] ^ 2 - 0.5 * x[,10] ^ 2 + rnorm(n, sd = sqrt(2))
  
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_2_2000_100[i]=fit$ate
  
  
}





library(ggplot2)
library(tidyverse)
library(patchwork)
res_mat=read.table("resmat81",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv81",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_2_500_10)
g1=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + theme_bw(base_size = 18)+labs(x=NULL)+
  theme(axis.text.x = element_blank())+ggtitle("n=500,p=10")
res_mat=read.table("resmat83",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv83",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_2_500_100)
#2nd plot
g3=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + labs(x=NULL)+ theme_bw(base_size = 18)+
  theme(axis.text.x = element_blank())+labs(y=NULL)+ggtitle("n=500,p=100")

res_mat=read.table("resmat87",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv87",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_2_2000_10)
g7=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle("n=2000,p=10")

res_mat=read.table("resmat89",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv89",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_2_2000_100)
g9=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1) )+labs(y=NULL)+ggtitle("n=2000,p=100")

rplot8=g1+g3+g7+g9
set.seed(12345)
result_3_500_30=numeric(40)
n=500
p=30

for(i in 1:40){

rho <- -0.25
Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))

x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)

## scramble columns of x so that which pairs of columns
## are correlated is at random
x <- x[,sample.int(p, p, replace = FALSE)]

## true propensity scores..
## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
aprob <- 1/(1 + exp(-.25*(apply(x[,c(c(1:10),c(21:30))],1,sum))))

## generate treatment assignments
a <- rbinom(n, 1, aprob)

## generate response with confounding.
## here the true ATE (treatment effect) is 0
y <- apply(x[,1:20],1,sum) + rnorm(n, sd = sqrt(2))
fit=forest_balance(x,a,y,cross.fitting = TRUE)
result_3_500_30[i]=fit$ate
}

set.seed(12345)
result_3_500_100=numeric(40)
n=500
p=100

for(i in 1:40){
  
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- 1/(1 + exp(-.25*(apply(x[,c(c(1:10),c(21:30))],1,sum))))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  y <- apply(x[,1:20],1,sum) + rnorm(n, sd = sqrt(2))
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_3_500_100[i]=fit$ate
}

set.seed(12345)
result_3_2000_30=numeric(40)
n=2000
p=30

for(i in 1:40){
  
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- 1/(1 + exp(-.25*(apply(x[,c(c(1:10),c(21:30))],1,sum))))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  y <- apply(x[,1:20],1,sum) + rnorm(n, sd = sqrt(2))
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_3_2000_30[i]=fit$ate
}

set.seed(12345)
result_3_2000_100=numeric(40)
n=2000
p=100

for(i in 1:40){
  
  rho <- -0.25
  Sig <- (rho) ^ abs(outer(1:p, 1:p, FUN = "-"))
  
  x <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  
  ## scramble columns of x so that which pairs of columns
  ## are correlated is at random
  x <- x[,sample.int(p, p, replace = FALSE)]
  
  ## true propensity scores..
  ## have to make sure the probabilities are not too close to 0 or 1 or else positivity will not hold!
  aprob <- 1/(1 + exp(-.25*(apply(x[,c(c(1:10),c(21:30))],1,sum))))
  
  ## generate treatment assignments
  a <- rbinom(n, 1, aprob)
  
  ## generate response with confounding.
  ## here the true ATE (treatment effect) is 0
  y <- apply(x[,1:20],1,sum) + rnorm(n, sd = sqrt(2))
  fit=forest_balance(x,a,y,cross.fitting = TRUE)
  result_3_2000_100[i]=fit$ate
}

res_mat=read.table("resmat141",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv141",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_3_500_30)
g1=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + theme_bw(base_size = 18)+labs(x=NULL)+
  theme(axis.text.x = element_blank())+ggtitle("n=500,p=30")
res_mat=read.table("resmat143",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv143",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_3_500_100)
#2nd plot
g3=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + labs(x=NULL)+ theme_bw(base_size = 18)+
  theme(axis.text.x = element_blank())+labs(y=NULL)+ggtitle("n=500,p=100")

res_mat=read.table("resmat147",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv147",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_3_2000_30)
g7=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle("n=2000,p=30")

res_mat=read.table("resmat149",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv149",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop,RF_MMD_CF=result_3_2000_100)
g9=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1) )+labs(y=NULL)+ggtitle("n=2000,p=100")

rplot14=g1+g3+g7+g9


# 






