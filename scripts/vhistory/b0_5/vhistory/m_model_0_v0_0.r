################
## analysis.r ##
################

## goal: perform Bayesian data analysis of stream lake dataset 
##       to find relationship between temperature and DBO on network structure

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 20-03-2022 - created v0_0
## 09-05-2022 - renamed model_0_v0_0

###############
## FUNCTIONS ##
###############

#
###

##############
## INITIATE ##
##############

## library
library(Rcpp)
sourceCpp("cpp/DEMCpp_v0.1.cpp")

## load dataset
load("data/dataset_lake_stream.rda")
data = dataset_lake_stream

#
###

##########
## MAIN ##
##########

## visualise data
# plot(data)

#
###

#############
## MODEL 1 ##
#############

## response
Y = data$w_trph_lvl_avg

## explanatory
X = cbind(1,data$temp,data$temp^2)

## predictive model
Yhat = function(X,beta) X%*%beta

## parameters
beta_0 = rnorm(ncol(X),0,0.01)
sd_lik_0 = rnorm(1,0,0.01)

## hyperparameters
sd_pri = c(1.0,1.0)

## log posterior
logLik = function(X,Y,beta,sd_lik)        sum(dnorm(Y-Yhat(X,beta),0,sd_lik,log=T))
logPri = function(beta,sd_lik,sd_pri)     dnorm(log(sd_lik),0,sd_pri[1],log=T) +  sum(dnorm(beta,0,sd_pri[2]),log=T)
logPos = function(X,Y,beta,sd_lik,sd_pri) logLik(X,Y,beta,sd_lik) + logPri(beta,sd_lik,sd_pri)
logPosWrap = function(omega) logPos(X,Y,omega[-1],omega[1],sd_pri)
map = function(omega)   c(log(omega[1]),omega[-1])
unmap = function(omega) c(exp(omega[1]),omega[-1])
dTarget = function(omega) logPosWrap(unmap(omega))

## test 
beta = beta_0
sd_lik = sd_lik_0
omega  = c(sd_lik,beta)
# print(logLik(X,Y,beta,sd_lik))
# print(logPri(beta,sd_lik,sd_pri))
# print(logPos(X,Y,beta,sd_lik,sd_pri))
print(dTarget(omega))

## DEMCpp
chain = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega, "epsilon" = 0.001, "nIt" = 10000))$chainList
chain[,-1] = t(apply(chain[,-1],1,unmap))

## visualise
chain_ = chain[seq(5000,10000,100),]
plot(data.frame(chain_))
# hist(chain_[,4])
# hist(chain_[,5])

## MaP
beta_map = chain_[which.max(chain_[,1]),-1][-1]
plot(X[,2],Y)
s = order(X[,2])
lines(X[s,2],Yhat(X,beta_map)[s,],col="red")

#
###
