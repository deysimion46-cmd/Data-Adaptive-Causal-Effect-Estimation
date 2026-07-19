library(osqp)
library(MASS)
energy_balance <- function(x,         ## matrix of covariates for all RCT and non-RCT patients
                           trt,       ## vector of treatment indicators (either 0 for control or 1 for treated)
                           standardize_x = FALSE, ## standardize the covariates?
                           verbose = FALSE,       ## print out text about quadratic programming algorithm convergence?
                           sigma = 0.5,           ## sigma for gaussian kernel. only used if type = "gaussian_kernel"
                           K_custom = NULL,       ## an n x n kernel matrix (n=nrow(x)). only used if type = "custom"
                           solver = c("osqp", "cccp"), ## which quadratic programming software to use? osqp is fastest
                           type = c("energy",
                                    "gaussian_kernel",
                                    "custom")) ## what kernel/distance to use?
{
  if (standardize_x)
  {
    x <- scale(x)
  }
  
  ## dimension of covariates
  p <- NCOL(x)
  
  
  type    <- match.arg(type)
  solver  <- match.arg(solver)
  
  if (type == "custom")
  {
    if (is.null(K_custom))
    {
      stop("if type set to custom, must specify custom kernel matrix!")
    }
    dist_x <- median_dist <- NULL
  } else
  {
    ## create euclidean distance matrix
    dist_x  <- as.matrix(dist(x))
    
    median_dist <- median(dist_x)
  }
  
  
  
  
  ## we will optimize weighted energy distance
  ## by minimizing f(w) = w'Qw + a'w subject to constraints on w.
  ## constraints are that w >= 0 and Aw = b, where b is a vector
  
  ## the following code constructs the Q matrix and 'a' vector
  ## corresponding to the weighted energy distance.
  
  ## if the treatment is not randomized, we will
  ## balance the treatment arm and the control arm separately
  ## to the target population (instead of the combined trt+control sample).
  ## this makes Q and 'a' slightly different.
  
  
  N <- nrow(x)
  
  if (grepl("gaussian_kernel", type))
  {
    Q_x   <- -exp(-sigma * dist_x ^ 2 / median_dist ^ 2)
  } else if (type == "energy")
  {
    Q_x   <- dist_x
  } else if (type == "custom")
  {
    Q_x <- -K_custom
  }
  
  
  Q_mat <- -trt * t( trt * t(Q_x)) / sum(trt == 1) ^ 2 -
    (1 - trt) * t( (1 - trt) * t(Q_x)) / sum(trt != 1) ^ 2
  
  aa <- 2 * as.vector(rowSums(trt * Q_x)) / (sum(trt == 1) * N) +
    2 * as.vector(rowSums((1-trt) * Q_x)) / (sum(trt != 1) * N)
  
  
  #evals <- eigen(Q_mat)$values
  
  ## only add nugget if we plan to use the solve.QP() function
  # if (min(evals) < 0 & (grepl("gaussian_kernel", type) | grepl("gaussian_anova_kernel", type)))
  # {
  #   nugget <- -min(evals) + 1e-10
  # } else
  # {
  #   nugget <- 0
  # }
  
  nugget <- 0
  
  ## add the smallest value to diagonal such that
  ## the Q matrix is positive definite
  Q_mat <- Q_mat + (nugget) * diag(nrow(Q_mat))
  
  
  ## set up constraints on weights
  AA1           <- matrix(0, ncol = N, nrow = 1)
  AA0           <- matrix(0, ncol = N, nrow = 1)
  
  AA1[,trt == 1] <- 1
  AA0[,trt != 1] <- 1
  
  
  A           <- rbind(AA1, AA0)
  rownames(A) <- paste0("eq", 1:nrow(A))
  sum.constr   <- c(sum(trt == 1), sum(trt != 1))
  
  
  ## A*w = sum.constr forces the weights
  ## to sum to the sample size of each treatment group.
  ## ie w_i : trt_i = 1 will sum to the number of trt == 1
  
  
  if (solver == "osqp")
  {
    ## minimizing f(w) = 0.5*w'Qw - a'w
    
    ## we will set up the quadratic program to be solved by
    ## the solve.QP function accessed from quadprog. this solver
    ## is very fast but cannot handle negative definite Q_mat
    
    
    ## set up the linear (in)equality constraints. The first 1 or 2
    ## constraints are the equality constraints. the remaining are the
    ## inequality constraints (ie making sure all weights are positive)
    
    #Amat <- t( rbind(A, diag(nrow(Q_mat)), -diag(nrow(Q_mat)) ) )
    
    Amat <- rbind(diag(N), AA1, AA0)
    
    lvec <- c(rep(0, N), sum.constr)
    uvec <- c(rep(Inf, N), sum.constr)
    
    ## Amat * w (>)= bvec
    #bvec <- c(sum.constr, rep(1e-14, nrow(Q_mat)), rep(-10*(nrow(Q_mat) ^ (1/3)), nrow(Q_mat)))
    
    #res <- solve.QP(Dmat = Q_mat, dvec = -aa/2, Amat = Amat,
    #                bvec = bvec, meq = length(sum.constr))
    
    
    opt.out <- osqp::solve_osqp(Q_mat,
                                q = aa/2,
                                A = Amat, l = lvec, u = uvec,
                                pars = osqp::osqpSettings(max_iter = 2e5,
                                                          eps_abs = 1e-8,
                                                          eps_rel = 1e-8,
                                                          verbose = FALSE))
    
    # if (!identical(opt.out$info$status, "maximum iterations reached") & !(any(opt.out$x > 1e5)))
    # {
    #   break
    # }
    
    energy_wts <- unname(opt.out$x)
    energy_wts[energy_wts < 0] <- 0 #due to numerical imprecision
    
  } else
  {
    ## minimizing f(w) = w'Qw + a'w
    
    ## we will set up the quadratic program to be solved by
    ## the cccp function accessed from optiSolve. this solver
    ## can handle when Q_mat is negative definite very well!
    
    rownames(Q_mat) <- paste(1:NROW(Q_mat))
    
    qf <- quadfun(Q = Q_mat, a = aa, id = rownames(Q_mat)) #quadratic obj.
    
    lb <- lbcon(0.0, id = rownames(Q_mat)) #lower bound
    
    ub <- ubcon(10*(nrow(Q_mat) ^ (1/3)), id = rownames(Q_mat)) #upper bound. not really needed
    
    ## sup up linear constraints
    lc <- lincon(A   = A,
                 dir = rep("==", length(sum.constr)),
                 val = sum.constr,
                 id  = rownames(Q_mat)) #linear constraint
    
    ## set up the QP
    lcqp <- cop( f = qf, lb = lb, ub = ub, lc = lc )
    
    opt.out <- solvecop(lcqp, solver = "cccp", quiet = !verbose)
    
    energy_wts <- unname(opt.out$x)
    
    lagrangian <- NULL
    
  }
  
  energy_wts[energy_wts < 0] <- 0
  
  list(weights = energy_wts, opt = opt.out)
}


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

## random forests with ranger

## random forests with grf (generalized random forests)
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


### define some functions for generating data, ATE estimates, and the wAMD,
expit = function(x){ 
  pr = ( exp(x) / (1+exp(x)) ) 
  return(pr)
}
ATE_est = function(fY,fw,fA){
  t_ATE = fY*fw
  tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
  return(tt_ATE) 
}
create_weights = function(fp, fA, fw){
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  return(fw)
}

## weighted absolute mean difference, where variables are either
## upweighted or downweighted depending on their outcome model
## coefficient beta
wAMD_function = function(DataM, varlist, trt.var, wgt, beta)
{
  trt = untrt = diff_vec = rep(NA,length(beta)) 
  names(trt) = names(untrt) = names(diff_vec) = varlist
  
  for(jj in 1:length(varlist))
  { 
    this.var         <- paste("w",varlist[jj],sep="") 
    DataM[,this.var] <- DataM[,varlist[jj]] * wgt
    trt[jj]      <- sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(wgt[DataM[,trt.var]==1]) 
    untrt[jj]    <- sum( DataM[DataM[,trt.var]==0, this.var]) / sum(wgt[DataM[,trt.var]==0]) 
    diff_vec[jj] <- abs( trt[jj] - untrt[jj] ) 
  } 
  wdiff_vec <- diff_vec * abs(beta) 
  wAMD      <- c( sum(wdiff_vec))
  ret       <- list( diff_vec  = diff_vec, 
                     wdiff_vec = wdiff_vec, 
                     wAMD      = wAMD )
  return(ret) 
}



outcome_adaptive_lasso <- function(x, trt, y,
                                   nlambda = 100,
                                   gamma_convergence_factor = 2,
                                   log_lam_vec = c(-5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49, 1, 2))
{
  vnames <- colnames(x)
  
  
  p <- NCOL(x)
  n <- NROW(x)
  if (is.null(vnames)) vnames <- paste0("X", 1:p)
  
  lambda_vec <- n ^ log_lam_vec
  names(lambda_vec) = as.character(lambda_vec)
  
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals <- 2*( gamma_convergence_factor - log_lam_vec + 1 )
  names(gamma_vals) = names(lambda_vec)
  
  data <- data.frame(Y = y, A = trt, x)
  
  # Normlize coviarates to have mean 0 and standard deviation 1
  temp.mean = colMeans(data[,vnames])
  Temp.mean = matrix(temp.mean,ncol=length(vnames),nrow=nrow(data),byrow=TRUE)
  data[,vnames] = data[,vnames] - Temp.mean
  temp.sd = apply(data[vnames],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(vnames),nrow=nrow(data),byrow=TRUE)
  data[vnames] = data[,vnames] / Temp.sd
  
  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(vnames,collapse="+")))
  lm.Y = lm(y.form, data=data)
  betaXY = coef(lm.Y)[-c(1,2)]
  
  
  ## Want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(vnames)+1,ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = c("(Intercept)", vnames)
  
  
  
  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda value 
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(vnames,collapse="+")))
  
  
  nlam <- length(lambda_vec)
  
  weight_mat <- matrix(NA, nrow = nrow(data), ncol = nlam)
  
  for( lil in 1:nlam )
  {
    il = unname(lambda_vec[lil])
    ig = gamma_vals[lil]
    
    x_dm <- model.matrix(w.full.form, data = data)[,-1]
    penalty_weights <- abs(betaXY)^(-ig)
    
    
    ### run outcome-adaptive lasso model with appropriate penalty.
    ## fit a sequence of lambdas so that glmnet actually converges for the lambda of interest
    logit_oal <- tryCatch(
      #try to do this
      {
        fit_oal <- glmnet( x = x_dm, y = data$A, penalty.factor = penalty_weights, family = "binomial",
                           nlambda = nlambda,
                           lambda.min.ratio = 1e-8, alpha = 0.999) ## for numerical stability
        fit_oal
      },
      #if an error occurs, tell me the error
      error=function(e) {
        NULL
      },
      warning=function(e) {
        NULL
      }
    )
    
    if (is.null(logit_oal))
    {
      # save propensity score coefficients
      coeff_XA[c("(Intercept)", vnames), lil] = NA
      # create inverse probability of treatment weights
      weight_mat[,lil] = NA
      
      wamd_val <- Inf
    } else
    {
      # generate propensity score
      ps_cur = unname(drop(predict(logit_oal, newx = x_dm, s = il, type = "response")))
      
      # save propensity score coefficients
      coeff_XA[c("(Intercept)", vnames), lil] = drop(as.matrix(predict(logit_oal, type = "coef", s = il)))
      # create inverse probability of treatment weights
      weight_mat[,lil] = create_weights(fp = ps_cur, fA = data$A)
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_cur <- wAMD_function(DataM = data,
                                varlist = vnames,
                                trt.var="A",
                                wgt = weight_mat[,lil], 
                                beta = betaXY)
      
      wamd_val <- wAMD_cur$wAMD
    }
    
    wAMD_vec[lil] = wamd_val
    # save ATE estimate for this lambda value
    ATE[lil] = ATE_est(fY = data$Y,fw = weight_mat[,lil], fA = data$A)
  } # close loop through lambda values
  
  
  # find the lambda value that creates the smallest wAMD
  tt = which.min(wAMD_vec)
  # print out ATE corresponding to smallest wAMD value
  ATE[tt]
  # print out the coefficients for the propensity score that corresponds with smallest wAMD value 
  coeff_XA[,tt]
  
  best_weights <- weight_mat[,tt]
  
  list(ATE = unname(ATE[tt]), weights = best_weights, ATE_vec = ATE,
       coefs = coeff_XA, wAMD = wAMD_vec, gamma_vals = gamma_vals)
}
simRep <- function(n = 250, p = 10)
{
  
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
  Data_X=data.frame(x)
  bA = 0
  
  # set vector of possible lambda's to try
  lambda_vec = c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_vec) = as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) = names(lambda_vec)
  oal_fit=outcome_adaptive_lasso(Data_X,a,y)
  
  
  ## fit a logistic regression model for propensity score
  glmf <- glm(a ~ x, family = binomial())
  
  ## fitted/estimated propensity scores
  fitted_prob <- unname(drop(fitted(glmf)))
  
  OutLasso=oal_fit$ATE
  
  dfxa  <- data.frame(a = a, x = x)
  dfxay <- data.frame(a = a, x = x, y = y)
  
  ## model for treatment
  formula_a <- as.formula(paste0("a ~ ", paste(colnames(dfxa)[-1], collapse = " + ")))
  
  ## model for response
  formula_y <- as.formula(paste0("y ~ ", paste(colnames(dfxa)[-1], collapse = " + ")))
  
  
  ## random forest predicting treatment with the ranger package
  rf_fit_a <- ranger(formula_a, data = dfxa,
                     num.trees = 1000)
  
  ## random forest predicting response with the ranger package
  rf_fit_y <- ranger(formula_y, data = dfxay,
                     num.trees = 1000)
  
  
  ## random forest predicting treatment with the grf package
  grf_fit_a     <- regression_forest(x, Y = a,  num.trees = 1000)
  
  ## random forest predicting response with the grf package
  grf_fit_y     <- regression_forest(x, Y = y,  num.trees = 1000)
  
  ## multi-task random forest predicting treatment and response jointly.
  ## the hope is that this may work better by taking into account *both* the relationship of covariates
  ## with the treatment and their relationship with the response (i.e. hopefully focusin in more on confounders)
  grf_fit_multi <- multi_regression_forest(x, Y = scale(cbind(a, y)),  num.trees = 1000,
                                           min.node.size = 10)
  
  
  
  ## standard energy balancing weights designed to balance distributions of covariates
  ebw_fit <- energy_balance(x, a, solver = "osqp", type="gaussian_kernel",standardize_x = FALSE)
  
  
  ###########################################################################################
  ###########################################################################################
  ##
  ##  construct random forest kernel matrices based on random forest similarity measure
  ##
  ###########################################################################################
  ###########################################################################################
  
  ## kernel matrix from random forest fit on treatment with ranger
  kern_mat_rf_a <- forest_2_kernel_matrix(rf_fit_a, newdata = dfxa)
  
  ## kernel matrix from random forest fit on response with ranger
  kern_mat_rf_y <- forest_2_kernel_matrix(rf_fit_y, newdata = dfxa)
  
  ########################################################################################################
  ## NOTE: making the kernel matrix for GRF forests is extremely slow. will need to make faster somehow
  
  # ## kernel matrix from random forest fit on treatment with grf
  kern_mat_grf_a <- forest_2_kernel_matrix(grf_fit_a, newdata = x)
  # 
  # ## kernel matrix from random forest fit on response with grf
  # kern_mat_grf_y <- forest_2_kernel_matrix(grf_fit_y, newdata = x)
  # 
  # ## kernel matrix from random forest fit on joint model for treatment and response with grf.
  kern_mat_grf_multi <- forest_2_kernel_matrix(grf_fit_multi, newdata = x)
  
  
  ########################################################################################################
  ##
  ##     estimate weights that optimize
  ##     a kernelized version of the weighted energy 
  ##     distance with the random forest kernel
  ##
  ########################################################################################################
  
  ebw_rf_a <- energy_balance(x, a, K_custom = kern_mat_rf_a, solver = "osqp", type = "custom")
  
  ebw_rf_y <- energy_balance(x, a, K_custom = kern_mat_rf_y, solver = "osqp", type = "custom")
  
  ebw_rf_y_sum <- energy_balance(x, a, K_custom = kern_mat_rf_a + kern_mat_rf_y, solver = "osqp", type = "custom")
  
  ebw_grf_a <- energy_balance(x, a, K_custom = kern_mat_grf_a, solver = "osqp", type = "custom")
  # 
  # ebw_grf_y <- energy_balance(x, a, K_custom = kern_mat_grf_y, solver = "osqp", type = "custom")
  # 
  ebw_grf_multi <- energy_balance(x, a, K_custom = kern_mat_grf_multi, solver = "osqp", type = "custom")
  
  
  ## entropy balancing weights designed to exactly balance covariate averages
  ebal_weightit <- weightit(a ~ ., data = dfxa, method = "ebal")
  
  ## estimate propensity score weights by fitting a gradient-boosted tree model
  ## for the propensity score
  gbm_weightit <- weightit(a ~ ., data = dfxa, method = "gbm")
  
  
  ## model for response given x and a (with their interaction as well).
  ## this is for an outcome model approach to estimating the treatment effect
  formul_yxa_int <- as.formula(paste0("y ~ a*(", paste(colnames(dfxa)[-1], collapse = " + "), ")"))
  lm_yxa <- lm(formul_yxa_int, data = dfxay)
  
  dfxa1 <- data.frame(a = 1, x = x)
  dfxa0 <- data.frame(a = 0, x = x)
  
  yxa1 <- predict(lm_yxa, type = "response", newdata = dfxa1)
  yxa0 <- predict(lm_yxa, type = "response", newdata = dfxa0)
  
  ## this is an outcome model approach to estimating the average treatment effect
  outcome_mod_est <- mean(yxa1) - mean(yxa0)
  
  ## this is a doubly robust approach to estimating the average treatment effect.
  ## based on the linear model for the outcome and random forest kernel balancing weights
  dr_est_rf_a <- outcome_mod_est + weighted.mean((y - yxa1)[a == 1], w = ebw_rf_a$weights[a == 1]) - 
    weighted.mean((y - yxa0)[a != 1], w = ebw_rf_a$weights[a != 1])
  
  dr_est_ebw <- outcome_mod_est + weighted.mean((y - yxa1)[a == 1], w = ebw_fit$weights[a == 1]) - 
    weighted.mean((y - yxa0)[a != 1], w = ebw_fit$weights[a != 1])
  
  if (FALSE)
  {
    ## this is a good way to check balance of covariate averages
    bal.tab(a ~ . + x.6:x.7, data = dfxa, weights = list(ebw_rf = ebw_rf_a$weights,
                                                         ebw_grf = ebw_grf_a$weights,
                                                         ebw = ebw_fit$weights,
                                                         ebal = get.w(ebal_weightit),
                                                         gbm = get.w(gbm_weightit)),
            poly = 2, 
            
            un = TRUE,stats=c("mean.diffs","ks.statistics"))
  }
  
  # weighted.mean(y[a == 1], w = ebw_grf$weights[a == 1]) - weighted.mean(y[a != 1], w = ebw_grf$weights[a != 1])
  
  
  true_ps_hajek <- weighted.mean(y[a == 1], w = 1/aprob[a == 1]) - weighted.mean(y[a != 1], w = 1/(1-aprob[a != 1]))
  true_ps_horv_thomp <- mean(y * a / aprob) - mean(y * (1-a) / (1-aprob))
  
  
  est_vec <- c(unweighted = mean(y[a == 1]) - mean(y[a != 1]), 
               ## energy balancing weights estimate
               ebw = weighted.mean(y[a == 1], w = ebw_fit$weights[a == 1]) - weighted.mean(y[a != 1], w = ebw_fit$weights[a != 1]),
               ## RF kernel balancing weights estimates
               rf_bal_a = weighted.mean(y[a == 1], w = ebw_rf_a$weights[a == 1]) - weighted.mean(y[a != 1], w = ebw_rf_a$weights[a != 1]),
               rf_bal_y = weighted.mean(y[a == 1], w = ebw_rf_y$weights[a == 1]) - weighted.mean(y[a != 1], w = ebw_rf_y$weights[a != 1]),
               rf_bal_y_sum = weighted.mean(y[a == 1], w = ebw_rf_y_sum$weights[a == 1]) - weighted.mean(y[a != 1], w = ebw_rf_y_sum$weights[a != 1]),
               grf_bal_a = weighted.mean(y[a == 1], w = ebw_grf_a$weights[a == 1]) - weighted.mean(y[a != 1], w = ebw_grf_a$weights[a != 1]),
               #grf_bal_y = weighted.mean(y[a == 1], w = ebw_grf_y$weights[a == 1]) - weighted.mean(y[a != 1], w = ebw_grf_y$weights[a != 1]),
               grf_bal_multi = weighted.mean(y[a == 1], w = ebw_grf_multi$weights[a == 1]) - weighted.mean(y[a != 1], w = ebw_grf_multi$weights[a != 1]),
               ## entropy balancing weights estimate
               entropy = weighted.mean(y[a == 1], w = get.w(ebal_weightit)[a == 1]) - weighted.mean(y[a != 1], w = get.w(ebal_weightit)[a != 1]),
               ## gbm propensity score weights estimate
               gbm = weighted.mean(y[a == 1], w = get.w(gbm_weightit)[a == 1]) - weighted.mean(y[a != 1], w = get.w(gbm_weightit)[a != 1]),
               ## logistic regression propensity score estimate
               ps = weighted.mean(y[a == 1], w = 1/fitted_prob[a == 1]) - weighted.mean(y[a != 1], w = 1/(1-fitted_prob[a != 1])),
               ## estimates using the TRUE propensity scores
               ## a hajek normalized estimator
               true_ps_hajek = true_ps_hajek,
               ## the classic Horvitz-Thompson estimator (unbiased)
               true_ps_horv_thomp = true_ps_horv_thomp,
               ## outcome model estimate
               outcome_model = outcome_mod_est,
               ## doubly robust estimates
               dr_est_rf_a = dr_est_rf_a,
               dr_est_ebw = dr_est_ebw,
               OutLasso=OutLasso)
  est_vec
}




## number of replications. keeping this
## small for now just so it runs quickly.
## ideally, we would set this to 100 or 1000
nreps <- 40

## sample size:
n <- 500

## dimension of covariates (only 10 impact propensity score and outcome model)

p=30

set.seed(12345)
for (r in 1:nreps)
{
  ## run simulation once with 
  ## sample size and dimension specified
  res_tmp <- simRep(n = n, p = p)
  
  if (r == 1)
  {
    res_mat141 <- res_tmp
  } else
  {
    res_mat141 <- rbind(res_mat141, res_tmp)
  }
}
write.table(res_mat141,"resmat141",row.names = F,quote = F,col.names = TRUE)