#####
## ##
#####

## goal: perform Bayesian data analysis of stream lake dataset 
##       to find relationship between temperature and DBO on network structure

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 07-07-2022 - created v0_0

##############
## INITIATE ##
##############

## library
library(mvtnorm)
library(Rcpp)

## modules
sourceCpp("cpp/DEMCpp_v0.1.cpp")
source("f_HBM_v0_6.r")

## load dataset
load("data/dataset_lake_stream_v0_1.rda")
data = dataset_lake_stream

## path to out
pto = "out_maxTL_model_1"
system(paste("mkdir ",pto,sep=""))

## response
response = "Max. troph. level"

#
###

###############
## FUNCTIONS ##
###############

## std function
std = function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)

## convert string to integer
stringToInteger = function(string)
{
    tmp = 1:length(unique(string))
    names(tmp) = unique(string)
    tmp = tmp[string]
    return(tmp)
}

## sigmoid function
sigmoid = function(x) 1/(1+exp(-x))
logit   = function(x) log(x/(1-x))

#
###

##########
## DATA ##
##########

## remove missing locations
s    = which(is.na(data$CdBH))
data = data[-s,]

## response variable
Y      = std(data$max_troph_lvl)
# Y      = std(data$connectance)
n_data = length(Y)

## explanatory variables
type    = (data$type == "lake")*1
temp    = std(data$temp)
dbo     = std(data$dbo)
rich    = std(data$richness)
year    = std(data$year)
alt     = std(data$altitude)
long    = std(data$long)
latt    = std(data$lat)
locat   = stringToInteger(data$CdBH)

## marix of explanatory variables
X_obs   =                 cbind(1, type,temp,temp^2,type*temp,dbo,type*dbo,temp*dbo,rich,year,alt)
X_mis   =                 cbind(1, type,temp,temp^2,type*temp,  1,type*1  ,temp*1  ,rich,year,alt)
X_pred  = function(x,y,i) cbind(1,(i-1),   x,   x^2,  x*(i-1),  y, y*(i-1),     x*y,   0,   0,  0)
n_beta  = ncol(X_obs)
idx_par_mis = c(6:8)
colnames(X_obs) = c("1","type","temp","temp^2","type*temp","dbo","type*dbo","temp*dbo","rich","year","alt")

## missing values
idx_mis = which(is.na(dbo))
n_mis   = length(idx_mis)

## splitting data into observed and missing
Y_obs   = Y[-idx_mis]
X_obs   = X_obs[-idx_mis,]
Y_mis   = Y[idx_mis]
X_mis   = X_mis[idx_mis,]

## compute distance matrix
x_ = long
y_ = latt
DM = matrix(rep(0,length(x_)^2),ncol=length(x_),nrow=length(x_))
for(i in 1:length(x_))
{
    for(j in 1:length(y_))
    {
        DM[i,j] = sqrt((x_[i] - x_[j])^2 + (y_[i] - y_[j])^2)
    }
}
DM = rbind(DM[-idx_mis,],DM[idx_mis,])
DM = cbind(DM[,-idx_mis],DM[,idx_mis])

## uneven variances (random effects)
n_sd_lik   = length(unique(locat))
idx_sd_lik = c(locat[-idx_mis],locat[idx_mis])

#
###

###########
## MODEL ##
###########

## global variable
## indexing parameters in omega vector
idx_omega_sd_lik  = 1:n_sd_lik
idx_omega_sd_mis  =   n_sd_lik + 1
idx_omega_rho     =   n_sd_lik + 1 + 1:2
idx_omega_mu_mis  =   n_sd_lik + 1 +   2 + 1
idx_omega_beta    =   n_sd_lik + 1 +   2 + 1 + 1:n_beta
idx_omega_xmis    =   n_sd_lik + 1 +   2 + 1 +   n_beta + 1:n_mis
n_param           =   n_sd_lik + 1 +   2 + 1 +   n_beta +   n_mis

## autocorrelation function
ac = function(x,rho) rho[1]*exp(-x/rho[2])

## covariance matrix
Sigma = function(sd_lik,rho,DM)
{
    Sigma = (diag(1,n_data) + (1-diag(1,n_data)) * ac(DM,rho))*(sd_lik%*%t(sd_lik))
	return(Sigma)
}

## predictive model
Yhat = function(X,beta) X%*%beta

## likelihood of response 
logLik = function(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,rho,DM)
{
	## assemble missing data
	X_mis[,idx_par_mis] = X_mis[,idx_par_mis] * cbind(x_mis,x_mis,x_mis)
	
	## combine obs and missing data
	X_ = rbind(X_obs,X_mis)
	Y_ = c(Y_obs,Y_mis)

	## compute likelihood 
	logLik = sum(dmvnorm(x=Y_,mean=Yhat(X_,beta),sigma=Sigma(sd_lik[idx_sd_lik],rho,DM),log=T))
	return(logLik)
}

## probability density of missing variables
logLat_mis = function(x_mis,mu_mis,sd_mis)
{
	return(sum(dnorm(x_mis,mu_mis,sd_mis,log=T)))
}

## log priors
logPri     = function(beta,sd_lik,rho,mu_mis,sd_mis)
{
	logPri = sum( dnorm(  beta,0,10,log=T)) + 
             sum(dlnorm(sd_lik,0, 1,log=T)) + 
                  dbeta(rho[1],2, 5,log=T)  +
                 dlnorm(rho[2],0, 1,log=T)  +
                  dnorm(mu_mis,0,.1,log=T)  + 
                 dlnorm(sd_mis,0,.1,log=T) 
	return(logPri)
}

## random prior sample
rPrior = function()
{
    beta    =  rnorm(  n_beta,0,.1)
    sd_lik  = rlnorm(n_sd_lik,0,1)
    rho     = c(rbeta(1,2,5),rlnorm(1,0,1))
    mu_mis  =  rnorm(1,0,.1)
    sd_mis  = rlnorm(1,0,.1)
    x_mis   = rnorm(n_mis,mu_mis,sd_mis)
    omega   = c(sd_lik,sd_mis,rho,mu_mis,beta,x_mis)
}

## posterior density distribution
logPos     = function(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,rho,mu_mis,sd_mis,DM)
{  
	logPos = logLik(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,rho,DM) + 
             logLat_mis(x_mis,mu_mis,sd_mis) + 
             logPri(beta,sd_lik,rho,mu_mis,sd_mis)
	return(logPos)
}

## wrappers
logPosWrap  = function(omega) 
{
	sd_lik  = omega[idx_omega_sd_lik]
	sd_mis  = omega[idx_omega_sd_mis]
	rho     = omega[idx_omega_rho]
	mu_mis  = omega[idx_omega_mu_mis]
	beta    = omega[idx_omega_beta]
	x_mis   = omega[idx_omega_xmis]
	return(logPos(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,rho,mu_mis,sd_mis,DM))
}
wrap        = function(omega) c(    
                                    log(omega[   idx_omega_sd_lik]),
                                    log(omega[   idx_omega_sd_mis]),
                                sigmoid(omega[   idx_omega_rho[1]]),
                                    log(omega[   idx_omega_rho[2]]),
                                        omega[c(idx_omega_mu_mis,idx_omega_beta,idx_omega_xmis)])
unwrap      = function(omega) c(    
                                    exp(omega[   idx_omega_sd_lik]),
                                    exp(omega[   idx_omega_sd_mis]),
                                  logit(omega[   idx_omega_rho[1]]),
                                    exp(omega[   idx_omega_rho[2]]),
                                        omega[c(idx_omega_mu_mis,idx_omega_beta,idx_omega_xmis)])
dTarget     = function(omega) logPosWrap(unwrap(omega))

#
###

# ###########
# ## CHECK ##
# ###########
# 
# sd_lik  = rlnorm(n_sd_lik,0,.1)
# sd_mis  = rlnorm(       1,0,.1)
# rho     = c(rbeta(1,2,5),
#             rlnorm(1,0,.1))
# mu_mis  = rnorm(      1,0,.1)
# beta    = rnorm( n_beta,0,.1)
# x_mis   = rnorm(  n_mis,0,.1)
# omega   = c(sd_lik,sd_mis,rho,mu_mis,beta,x_mis)
# image(Sigma(sd_lik[idx_sd_lik],rho,DM)[1:100,1:100])
# print(system.time(print(logLik(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,rho,DM))))
# print(system.time(print(logLat_mis(x_mis,mu_mis,sd_mis))))
# print(system.time(print(logPri(beta,sd_lik,rho,mu_mis,sd_mis))))
# print(system.time(print(logPos(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,rho,mu_mis,sd_mis,DM))))
# print(system.time(print(logPosWrap(omega))))
# print(system.time(print(dTarget(wrap(omega)))))
# print(prod(omega == unwrap(wrap(omega))))
# 
# #
# ###
