library(osqp)
library(MASS)
library(grf)
library(ranger)
## for other weighting methods
library(WeightIt)

library(fields)

## needed for energy_balance only if solver="cccp"
#library(optiSolve)

## for checking balance
library(cobalt)

library(glmnet)

## this function inputs a fitted random forest
## and outputs predicted leaf node membership 
## for each observation for each tree in the forest:
##
## get predicted leaf node membership for every observation for every tree
## in the random forest. each row is an observation, each column is a tree
forest_2_leaf_node_matrix <- function(forest_object, newdata)
{
  ## check if object is from the grf package or ranger
  if (inherits(forest_object, "grf")) ## grf case
  {
    
    ## this is *really* slow, but I don't know
    ## a better way yet
    ntree <- forest_object[["_num_trees"]]
    nodes_by_trees <- matrix(NA, nrow = NROW(newdata), ncol = ntree)
    
    for (tr in 1:ntree)
    {
      nodes_by_trees[,tr] <- get_leaf_node(get_tree(forest_object, tr), newdata = newdata)
    }
    
  } else if (inherits(forest_object, "ranger")) ## ranger case
  {
    ## get predicted leaf node membership for every observation for every tree
    ## in the random forest.
    ## each row is an observation, each column is a tree
    nodes_by_trees <- predict(forest_object, type = "terminalNodes", data = newdata)$predictions
  }
  
  nodes_by_trees
}

## this function inputs the leaf node membership matrix
## and outputs the random forest kernel matrix
leaf_nodes_2_kernel <- function(nodes_by_trees)
{
  
  n <- nrow(nodes_by_trees)
  
  ## create nxn kernel matrix where the (i,j)th entry
  ## is the kernel similarity between observation i and j.
  ## similarity is the proportion of trees in which
  ## i and j are in the same leaf node
  K_mat <- matrix(0, nrow = n, ncol = n)
  diag(K_mat) <- 1
  
  for (i in 1:(n-1))
  {
    for (j in (i+1):n)
    {
      K_mat[i,j] <- K_mat[j,i] <- mean(nodes_by_trees[i,] == nodes_by_trees[j,])
    }
  }
  
  K_mat
}

## inputs random forest object
## and outputs kernel matrix
forest_2_kernel_matrix <- function(forest_object, newdata)
{
  ## create leaf node matrix
  leaf_node_mat <- forest_2_leaf_node_matrix(forest_object, newdata)
  
  ## return kernel matrix
  leaf_nodes_2_kernel(leaf_node_mat)
}






#CREATE THE DATASET 

n=500
p=10




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


dfxa  <- data.frame(a = a, x = x)
dfxay <- data.frame(a = a, x = x, y = y)

grf_fit_multi <- multi_regression_forest(x, Y = scale(cbind(a, y)),  num.trees = 1000,
                                         min.node.size = 10)


kern_mat_grf_multi <- forest_2_kernel_matrix(grf_fit_multi, newdata = x)

trt=a
s1=which(a==1)
s0=which(a==0)
n1=sum(a)
n0=n-n1
kern=kern_mat_grf_multi


K= trt * t( trt * t(kern)) / sum(trt == 1) ^ 2 +(1 - trt) * t( (1 - trt) * t(kern)) / sum(trt != 1) ^ 2
#chol_decomp <- chol(K_sparse)
#I=Diagonal(nrow(K))
#K_inv=backsolve(chol_decomp, backsolve(chol_decomp, I))
A=matrix(0,nrow=n,ncol=2)
A[,1]=a
A[,2]=1-a


b=numeric(n)
b <-  as.vector(rowSums(trt * kern)) / (sum(trt == 1) * n) +
  as.vector(rowSums((1-trt) * kern)) / (sum(trt != 1) * n)
#for(i in s1){
#b[i]=(1/(n1*n))*sum(kern[i,])
#}
#for(i in s0){
# b[i]=(1/(n0*n))*sum(kern[i,])
#}
#library(Matrix)
current_time=Sys.time()
K_sparse <- Matrix(K, sparse = TRUE)
f1=solve(K_sparse,a)
f2=solve(K_sparse,(1-a))
M1=cbind(f1,f2)
X=t(A)%*%M1
f=c(n1,n0)
y=t(A)%*%solve(K_sparse,b)-f
z=b-A%*%solve(X)%*%y
w=solve(K_sparse,z)
#Weights
print(w)
estimate=weighted.mean(y[a==1],w[a==1])-weighted.mean(y[a==0],w[a==0])


#estimate of ATE
estimate



