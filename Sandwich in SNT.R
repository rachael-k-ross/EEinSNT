
##################################################-
#
# Illustration of sandwich variance estimation 
# for sequential nested trials
#
# Rachael Ross, July 2026
#
##################################################-

library(dplyr)
library(tidyverse)


########################-
#
# Data ----
#
########################-

# set.seed(7)
# n <- 100000
# 
# # Simulate baseline emulated trial (k=1)
# 
# dat_k1 <- tibble(
#   i = seq(1,n),
#   S = 1,
#   W = rbinom(n, 1, 0.2),
#   A = rbinom(n, 1, 0.1 + 0.1*W),
#   Y0 = rbinom(n, 1, 0.3 + 0.1*W - 0.1*0),
#   Y1 = rbinom(n, 1, 0.3 + 0.1*W - 0.1*1),
#   Y = A*Y1 + (1-A)*Y0)
# 
# # Simulate 2nd emulated trial (k=2)
# 
# n_2 <- n - sum(dat_k1$A)
# 
# dat_k2 <- dat_k1 |>
#   filter(A == 0) |> # restrict to those assigned A=1 in first trial
#   select(i,W) |>
#   rename(priorW=W) |>
#   mutate(S = 2, 
#          W = rbinom(n_2, 1, 0.1 + 0.4*priorW),
#          A = rbinom(n_2, 1, 0.1 + 0.1*W),
#          Y0 = rbinom(n_2, 1, 0.3 + 0.1*W - 0.1*0),
#          Y1 = rbinom(n_2, 1, 0.3 + 0.1*W - 0.1*1),
#          Y = A*Y1 + (1-A)*Y0) |>
#   select(-priorW)


# Trial 1
dat_k1 <- tribble(
  ~i,  ~W, ~A, ~Y,
  1,   1,  1,  1,
  2,   0,  0,  0,
  3,   1,  0,  1,
  4,   0,  0,  0,
  5,   0,  1,  1,
  6,   1,  1,  1,
  7,   0,  0,  0,
  8,   1,  0,  1,
  9,   1,  1,  1,
  10,   0,  0,  0,
  11,   0,  0,  0,
  12,   1,  1,  1,
  13,   0,  0,  0,
  14,   1,  0,  1,
  15,   0,  0,  0,
  16,   1,  1,  1,
  17,   0,  0,  0,
  18,   0,  0,  1,
  19,   1,  1,  1,
  20,   1,  0,  0
) |> mutate(S = 1)

# Trial 2: nested among i's with A=0 in trial 1 
dat_k2 <- tribble(
  ~i,  ~W, ~A, ~Y,
  2,   0,  0,  1,
  3,   1,  1,  1,
  4,   0,  0,  0,
  7,   1,  0,  1,
  8,   0,  1,  1,
  10,   0,  0,  0,
  11,   1,  0,  0,
  13,   0,  0,  0,
  14,   1,  1,  1,
  15,   0,  0,  1,
  17,   1,  1,  1,
  18,   0,  0,  0,
  20,   1,  0,  0
) |> mutate(S = 2)
  

# Create long dataset
dat_long <- rbind(dat_k1, dat_k2)



################################-
#
# Sandwich with single trial ----
#
################################-


# Save dataset to be used as df
df <- dat_k1
#df <- dat_k2


#---------------------------------------#
### Step 1 - Point estimates ----
#---------------------------------------#


## Propensity score model

# Save model matrix 
X <- with(df, model.matrix(~ W))

# Estimate logistic model
psmodel <- glm(A ~ X - 1, data = df, binomial(link="logit"))


## Calculate weights
df_wts <- df %>%
  mutate(ps = predict(psmodel, ., type="response"),
         ipw = A/ps + (1-A)/(1-ps))


## Risk difference
rd <- with(df_wts,
           mean(ipw*Y*A - ipw*Y*(1 - A)))

## Store point estimates
ptests <- c(rd,
            psmodel$coefficients)



#---------------------------------------#
### Step 2 - Data setup for sandwich ----
#---------------------------------------#


## Create list of data elements 
dfmat <- with(df,
              list(A = A, # trt vector,
                   Y = Y, # outcome vector,
                   X = X, # model matrix for ps model 
                   nsize = length(A)) # number of unique individuals 
              ) 


#---------------------------------------#
### Step 3 - Stack of ee's ----
#---------------------------------------#


## Create function 
eefx_single <- function(theta, dfmat){ # Inputs: vector of parameters, list of data elements
  
  # Extract alpha parameters from theta
  alpha <- theta[2:length(theta)]
  
  # For ps model 
  pi <- plogis(drop(dfmat$X %*% alpha)) # propensity score
  ee_ps <- as.vector(dfmat$A - pi) * dfmat$X # ee's

  # For risk difference
  ee_rd <- dfmat$A*dfmat$Y/pi - (1-dfmat$A)*dfmat$Y/(1-pi) - theta[1] # ee
  
  # Output: n by p matrix
  return(cbind(ee_rd,
               ee_ps))
}



#---------------------------------------#
### Step 4 - Sandwich estimator ----
#---------------------------------------#


# Save the name of the ee fx created in previous step
eefx <- eefx_single


# Function for summing columns of ee fx 
sumcolee <- function(theta, eefx, dfmat){ 
  
  ee <- eefx(theta, dfmat = dfmat)
  
  return(colSums(ee)) # Output: p-length vector
}


# Sandwich
residuals <- eefx(ptests, dfmat) # residuals at the point estimates
meat <- crossprod(residuals) / dfmat$nsize # meat matrix
bread <- -numDeriv::jacobian(sumcolee, 
                             ptests, 
                             eefx = eefx, 
                             dfmat = dfmat) / dfmat$nsize # bread matrix
bread_inv <- solve(bread) # inverse of bread matrix
sandwich <- (bread_inv %*% meat %*% t(bread_inv)) / dfmat$nsize #putting it all together
ses <- sqrt(diag(sandwich))



#---------------------------------------#
### Results ----
#---------------------------------------#

results <- tibble(
  rd = ptests,
  se = ses,
  lcl = rd - 1.96*se,
  ucl = rd + 1.96*se
)

results


# # To confirm psmodel ses
# ses[2:length(ses)]
# summary(psmodel)$coeff[,2]



################################-
#
# Sandwich with long dataset ----
# results still stratified by trial
#
################################-

# Save dataset to be used as df
df <- dat_long


#---------------------------------------#
### Step 1 - Point estimates ----
#---------------------------------------#


## Propensity score model

# Save model matrix 
X <- with(df, model.matrix(~ W *factor(S)))

# Estimate logistic model
psmodel <- glm(A ~ X - 1, data = df, binomial(link="logit"))


## Calculate weights
df_wts <- df %>%
  mutate(ps = predict(psmodel, ., type="response"),
         ipw = A/ps + (1-A)/(1-ps))


## Risk difference
rd <- df_wts |>
  group_by(S) |>
  summarise(
    r1 = mean(ipw * Y * A),
    r0 = mean(ipw * Y * (1 - A)),
    rd = r1 - r0,
    .groups = "drop"
  ) |>
  arrange(S)



## Store point estimates
ptests <- c(rd$rd,
            psmodel$coefficients)



#---------------------------------------#
### Step 2 - Data setup for sandwich ----
#---------------------------------------#


## Create list of data elements 
dfmat <- with(df,
              list(# same as above
                A = A, # trt vector,
                Y = Y, # outcome vector,
                X = X, # model matrix for ps model
                nsize = sum(S==1), # number of unique individuals
                
                # new elements
                S = S, # trial indicator
                i = i, # unique id for subjects
                K = length(unique(S)) # number of trials
                )) 



#---------------------------------------#
### Step 3 - Stack of ee's ----
#---------------------------------------#


## Create function
eefx_snt <- function(theta, dfmat){ # Inputs: vector of parameters, list of data elements
  
  
  # Extract alpha parameters from theta
  alpha <- theta[(dfmat$K + 1):length(theta)]
  
  
  ## For ps model 
  pi <- plogis(drop(dfmat$X %*% alpha)) # propensity score
  long_ee_ps <- as.vector(dfmat$A - pi) * dfmat$X 
  
  # *Key step*: long_ee_ps is too long, must sum within person -> n by length(alpha)
  ee_ps <- rowsum(long_ee_ps, group = dfmat$i)  

  
  ## For risk differences
  long_ee_rd <- dfmat$A*dfmat$Y/pi - (1-dfmat$A)*dfmat$Y/(1-pi) - theta[dfmat$S] 
  
  # *Key step*: long_ee_rd is too long, must make wide -> n by K
  ee_rd <- matrix(0, nrow = dfmat$nsize, ncol = dfmat$K)
  ee_rd[cbind(dfmat$i, dfmat$S)] <- long_ee_rd
  
      # ## Alternative for ee_rd (less efficient but possibly more intuitive)
      # ee_rd <- matrix(0, nrow = dfmat$nsize, ncol = dfmat$K)
      # 
      # for (k in 1:dfmat$K){
      #   foree <- (dfmat$S == k) * (dfmat$A*dfmat$Y/pi - (1-dfmat$A)*dfmat$Y/(1-pi)  - theta[k]) 
      #   foree_ <- rowsum(foree, group = dfmat$i) 
      #   ee_rd[,k] <- foree_
      # }

    
  ## Output: n by p matrix
  return(cbind(ee_rd,
               ee_ps))
}


#---------------------------------------#
### Step 4 - Sandwich estimator ----
#---------------------------------------#


# Save the name of the ee fx created in previous step
eefx <- eefx_snt


# Function for summing columns of ee fx  - already defined above
 #see sumcolee 


# Sandwich
residuals <- eefx(ptests, dfmat) # residuals at the point estimates
meat <- crossprod(residuals) / dfmat$nsize # meat matrix
bread <- -numDeriv::jacobian(sumcolee, 
                             ptests, 
                             eefx = eefx, 
                             dfmat = dfmat) / dfmat$nsize # bread matrix
bread_inv <- solve(bread) # inverse of bread matrix
sandwich <- (bread_inv %*% meat %*% t(bread_inv)) / dfmat$nsize #putting it all together
ses <- sqrt(diag(sandwich))



#---------------------------------------#
### Results ----
#---------------------------------------#

results <- tibble(
  rd = ptests,
  se = ses,
  lcl = rd - 1.96*se,
  ucl = rd + 1.96*se
)

results

# # To confirm psmodel ses
# ses[3:length(ses)]
# summary(psmodel)$coeff[,2]



################################-
#
# Sandwich with pooling ----
#
################################-

# Save dataset to be used as df
df <- dat_long


#---------------------------------------#
### Step 1 - Point estimates ----
#---------------------------------------#


## Propensity score model

# Save model matrix (dropped interaction term)
X <- with(df, model.matrix(~ W + factor(S)))

# Estimate logistic model
psmodel <- glm(A ~ X - 1, data = df, binomial(link="logit"))


## Calculate weights
df_wts <- df %>%
  mutate(ps = predict(psmodel, ., type="response"),
         ipw = A/ps + (1-A)/(1-ps))


## Risk difference
rd <- with(df_wts,
           mean(ipw*Y*A - ipw*Y*(1 - A)))

## Store point estimates
ptests <- c(rd,
            psmodel$coefficients)



#---------------------------------------#
### Step 2 - Data setup for sandwich ----
#---------------------------------------#


## Create list of data elements 
dfmat <- with(df,
              list(# same as above
                A = A, # trt vector,
                Y = Y, # outcome vector,
                X = X, # model matrix for ps model
                nsize = sum(S==1), # number of unique individuals
                i = i # unique id for subjects
              )) 



#---------------------------------------#
### Step 3 - Stack of ee's ----
#---------------------------------------#


## Create function
eefx_sntpool <- function(theta, dfmat){ # Inputs: vector of parameters, list of data elements
  
  # Extract alpha parameters from theta
  alpha <- theta[2:length(theta)]
  
  
  ## For ps model - same as previous
  pi <- plogis(drop(dfmat$X %*% alpha)) # propensity score
  long_ee_ps <- as.vector(dfmat$A - pi) * dfmat$X 
  
  # *Key step*: long_ee_ps is too long, must sum within person -> n by length(alpha)
  ee_ps <- rowsum(long_ee_ps, group = dfmat$i)  
  
  
  ## For risk difference
  long_ee_rd <- dfmat$A*dfmat$Y/pi - (1-dfmat$A)*dfmat$Y/(1-pi) - theta[1] 
  
  # *Key step*: long_ee_rd is too long, must sum within person -> n by 1
  ee_rd <- rowsum(long_ee_rd, group = dfmat$i)  


  ## Output: n by p matrix
  return(cbind(ee_rd,
               ee_ps))
}


#---------------------------------------#
### Step 4 - Sandwich estimator ----
#---------------------------------------#


# Save the name of the ee fx created in previous step
eefx <- eefx_sntpool


# Function for summing columns of ee fx  - already defined above
#see sumcolee 


# Sandwich
residuals <- eefx(ptests, dfmat) # residuals at the point estimates
meat <- crossprod(residuals) / dfmat$nsize # meat matrix
bread <- -numDeriv::jacobian(sumcolee, 
                             ptests, 
                             eefx = eefx, 
                             dfmat = dfmat) / dfmat$nsize # bread matrix
bread_inv <- solve(bread) # inverse of bread matrix
sandwich <- (bread_inv %*% meat %*% t(bread_inv)) / dfmat$nsize #putting it all together
ses <- sqrt(diag(sandwich))



#---------------------------------------#
### Results ----
#---------------------------------------#

results <- tibble(
  rd = ptests,
  se = ses,
  lcl = rd - 1.96*se,
  ucl = rd + 1.96*se
)

results

