################
## analysis.r ##
################

## goal: perform Bayesian data analysis of stream lake dataset 
##       to find relationship between temperature and DBO on network structure

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 20-03-2022 - created v0_0
## 27-03-2022 - created v0_1
##            - created simple linear model
## 12-04-2022 - created v0_2
##            - implemened inference of missing values
## 27-04-2022 - created v0_4
##            - implemented spatial autocorrelations
## 09-05-2022 - renamed model_5_v0_0

###############
## FUNCTIONS ##
###############

#
###

##############
## INITIATE ##
##############

## library
library(mvtnorm)
library(Rcpp)
sourceCpp("cpp/DEMCpp_v0.1.cpp")
sourceCpp("cpp/DEMCO_v0.5.cpp")
source("hbm_functions_v0.5.r")

## load dataset
load("data/dataset_lake_stream.rda")
data = dataset_lake_stream

#
###

#################################
## MODEL 2 - MAX TROPHIC LEVEL ##
#################################

## goal: infer missing values in addition to the linear model

## remove missing locations
idx_mis = which(is.na(data$lon))
data = data[-idx_mis,]

## response variable
Y = std(data$max_troph_lvl)
# Y = data$w_trph_lvl_avg
N = length(Y)

## std function
std = function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)

## explanatory variables
year    = std(data$year)
temp    = std(data$temp)
dbo     = std(data$dbo)
type    = (data$type == "lake")*1

## marix of explanatory variables
X_obs   = cbind(1,year,temp,temp^2,type,type*temp,dbo,dbo^2,type*dbo)
X_mis   = cbind(1,year,temp,temp^2,type,type*temp,1,1,type)
n_beta  = ncol(X_obs)

## missing values
idx_mis = which(is.na(dbo))
n_mis   = length(idx_mis)

## splitting data into observed and missing
Y_obs   = Y[-idx_mis]
X_obs   = X_obs[-idx_mis,]
Y_mis   = Y[idx_mis]
X_mis   = X_mis[idx_mis,]

## distance matrix
DM = matrix(0,nrow=N,ncol=N)
for(i in 1:N)
{
  xy_  = as.matrix(data[,c("lat","long")])
  xy_i = t(matrix(xy_[i,],2,N))
	DM[i,] = sqrt(apply((xy_i-xy_)^2,1,sum))
}
DM = rbind(DM[-idx_mis,],DM[idx_mis,])
DM = cbind(DM[,-idx_mis],DM[,idx_mis])

## autocorrelation function
ac = function(x,rho) rho[1]*exp(-rho[2]*x)

## build covariance matrix
Sigma = function(sd_lik,rho,DM)
{
	Sigma = sd_lik^2*ac(DM,rho)
	diag(Sigma) = sd_lik^2
	return(Sigma)
}

## visualise matrix
image(Sigma(1,c(0.5,1),DM)[1:100,1:100])

## predictive model
Yhat = function(X,beta) X%*%beta

## likelihood of response 
logLik = function(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,rho,DM)
{
	## assemble missing data
	X_mis[,7:9] = X_mis[,7:9] * cbind(x_mis,x_mis^2,x_mis)
	
	## combine obs and missing data
	X_ = rbind(X_obs,X_mis)
	Y_ = c(Y_obs,Y_mis)

	## compute likelihood 
	logLik = sum(dmvnorm(x=Y_,mean=Yhat(X_,beta),sigma=Sigma(sd_lik,rho,DM),log=T))
	return(logLik)
}
#
## check
system.time(logLik(Y_obs,X_obs,Y_mis,X_mis,rnorm(n_mis,0,1),rnorm(n_beta,0,0.1),1,c(0.5,0.5),DM))
system.time(Sigma(1,c(1,1),DM))
## => dmvnorm very slow

## probability density of missing variables
logLat_mis = function(x_mis,mu_mis,sd_mis)
{
	return(sum(dnorm(x_mis,mu_mis,sd_mis,log=T)))
}

## log priors
logPri     = function(beta,sd_lik,rho,mu_mis,sd_mis)
{
	logPri = sum(dnorm(beta,         0,10,log=T)) + 
                 dnorm(log(sd_lik),0,1, log=T)  + 
                 dnorm(log(rho[1]),0,1, log=T)  +
	               dnorm(log(rho[2]),0,1, log=T)  +
                 dnorm(mu_mis,     0,1, log=T)  + 
                 dnorm(log(sd_mis),0,1, log=T) 
	return(logPri)
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
	sd_lik  = omega[1]
	rho     = omega[2:3]
	sd_mis  = omega[4]
	mu_mis  = omega[5]
	beta    = omega[5 + 1:n_beta] 
	x_mis   = omega[5 + n_beta  + 1:n_mis] 
	return(logPos(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,rho,mu_mis,sd_mis,DM))
}
wrap        = function(omega) c(log(omega[1:4]),omega[-c(1:4)])
unwrap      = function(omega) c(exp(omega[1:4]),omega[-c(1:4)])
dTarget     = function(omega) logPosWrap(unwrap(omega))

## DEMCpp
chainList   = list()
for(i in 1:3)
{
  check = FALSE
  while(check == FALSE) 
  {
    omega_0    = rnorm(n_mis+n_beta+5,0,0.1)
    check = (dTarget(omega_0) > -Inf)
  }
	chain      = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_0, "epsilon" = 0.001, "nIt" = 100))$chainList
	chain[,-1] = t(apply(chain[,-1],1,unwrap))
	chainList[[i]] = chain
}

## MaP
omega_map = chainList.argmaxPost(chainList)
beta_map  = omega_map[5 + 1:n_beta]

## burn and thin
# chainList_thinned = chainList.thin(chainList.burn(chainList,1:50000))
chainList_thinned = list(chainList[[1]])

## visualise
# chainList.acPlot(chainList_thinned)
# chainList.postPlot(chainList_thinned,1000)
# chainList.bayesPlot(chainList_thinned)

## 
pdf("model_2_results_maxTL.pdf")

## plot predictions temperature
plot(X_obs[,3],Y_obs,pch=16,col=adjustcolor("black",alpha=0.5),xlab="temperature (SU)",ylab="mean trophic level (SU)")
for(i in 1:2)
{
	## predictions
	x        = seq(min(temp),max(temp),0.1)
	X_pred   = cbind(1,0,x,x^2,(i-1),(i-1)*x,0,0,0)
	pred     = chainList.apply(chainList_thinned,function(x_) Yhat(X_pred,x_[-1][5 + 1:n_beta]))
	polygon(x=c(X_pred[,3],rev(X_pred[,3])),y=c(pred$f_q0.05,rev(pred$f_q0.95)),border=NA,col=adjustcolor(i+1,alpha=0.5))
	lines(X_pred[,3],Yhat(X_pred,beta_map),col=i+1)
	lines(X_pred[,3],pred$f_mean,col=i+1)
}
legend("topright",legend=c("stream","lake"),lty=1,col=1:2+1,bty="n")

## plot predictions temperature
plot(X_obs[,7],Y_obs,pch=16,col=adjustcolor("black",alpha=0.5),xlab="dbo (SU)",ylab="mean trophic level (SU)")
for(i in 1:2)
{
	## predictions
	x        = seq(min(dbo,na.rm=T),max(dbo,na.rm=T),0.1)
	X_pred   = cbind(1,0,0,0,(i-1),0,x,x^2,(i-1)*x)
	pred     = chainList.apply(chainList_thinned,function(x_) Yhat(X_pred,x_[-1][4:12]))
	polygon(x=c(X_pred[,7],rev(X_pred[,7])),y=c(pred$f_q0.05,rev(pred$f_q0.95)),border=NA,col=adjustcolor(i+1,alpha=0.5))
	# lines(X_pred[,7],Yhat(X_pred,beta_map),col=i+1)
	lines(X_pred[,7],pred$f_mean,col=i+1)
}
legend("topright",legend=c("stream","lake"),lty=1,col=1:2+1,bty="n")

## 
dev.off()

#
###

# ##################################
# ## MODEL 2 - MEAN TROPHIC LEVEL ##
# ##################################
# 
# ## goal: infer missing values in addition to the linear model
# 
# ## response variable
# # Y = data$max_troph_lvl
# Y = data$w_trph_lvl_avg
# 
# ## std function
# std = function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)
# 
# ## explanatory variables
# year    = std(data$year)
# temp    = std(data$temp)
# dbo     = std(data$dbo)
# type    = (data$type == "lake")*1
# 
# ## marix of explanatory variables
# X_obs   = cbind(1,year,temp,temp^2,type,type*temp,dbo,dbo^2,type*dbo)
# X_mis   = cbind(1,year,temp,temp^2,type,type*temp,1,1,type)
# n_beta  = ncol(X_obs)
# 
# ## missing values
# idx_mis = which(is.na(dbo))
# n_mis   = length(idx_mis)
# 
# ## splitting data into observed and missing
# Y_obs   = Y[-idx_mis]
# X_obs   = X_obs[-idx_mis,]
# Y_mis   = Y[idx_mis]
# X_mis   = X_mis[idx_mis,]
# 
# ## predictive model
# Yhat = function(X,beta) X%*%beta
# 
# ## likelihood of response given observed variables
# logLik_obs = function(Y_obs,X_obs,beta,sd_lik)
# {
# 	return(sum(dnorm(Y_obs,Yhat(X_obs,beta),sd_lik,log=T)))
# }
# 
# ## likelihood of response given missing variables
# logLik_mis = function(Y_mis,X_mis,x_mis,beta,sd_lik)      
# {
# 	X_mis[,7:9] = X_mis[,7:9] * cbind(x_mis,x_mis^2,x_mis)
# 	logLik_mis  = sum(dnorm(Y_mis,Yhat(X_mis,beta),sd_lik,log=T))
# 	return(logLik_mis)
# }
# 
# ## probability density of missing variables
# logLat_mis = function(x_mis,mu_mis,sd_mis)
# {
# 	return(sum(dnorm(x_mis,mu_mis,sd_mis,log=T)))
# }
# 
# ## log priors
# logPri     = function(beta,sd_lik,mu_mis,sd_mis)              
# {
# 	logPri = sum(dnorm(beta,       0,10,log=T)) + 
#                  dnorm(log(sd_lik),0,1, log=T)  + 
#                  dnorm(mu_mis,     0,1, log=T)  + 
#                  dnorm(log(sd_mis),0,1, log=T) 
# 	return(logPri)
# }
# 
# ## posterior density distribution
# logPos     = function(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,mu_mis,sd_mis)
# {  
# 	logPos = logLik_obs(Y_obs,X_obs,beta,sd_lik) + 
#              logLik_mis(Y_mis,X_mis,x_mis,beta,sd_lik) + 
#              logLat_mis(x_mis,mu_mis,sd_mis) + 
#              logPri(beta,sd_lik,mu_mis,sd_mis)
# 	return(logPos)
# }
# 
# ## wrappers
# logPosWrap  = function(omega) 
# {
# 	sd_lik  = omega[1]
# 	sd_mis  = omega[2]
# 	mu_mis  = omega[3]
# 	beta    = omega[3 + 1:n_beta] 
# 	x_mis   = omega[3 + n_beta  + 1:n_mis] 
# 	return(logPos(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,mu_mis,sd_mis))
# }
# wrap       = function(omega) c(log(omega[1:2]),omega[-c(1,2)])
# unwrap     = function(omega) c(exp(omega[1:2]),omega[-c(1,2)])
# dTarget    = function(omega) logPosWrap(unwrap(omega))
# 
# ## DEMCpp
# chainList = list()
# for(i in 1:3)
# {
# 	omega_0    = rnorm(n_mis+n_beta+3,0,1)
# 	chain      = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_0, "epsilon" = 0.001, "nIt" = 100000))$chainList
# 	chain[,-1] = t(apply(chain[,-1],1,unwrap))
# 	chainList[[i]] = chain
# }
# 
# ## MaP
# omega_map = chainList.argmaxPost(chainList)
# beta_map  = omega_map[4:12]
# 
# ## burn and thin
# chainList_thinned = chainList.thin(chainList.burn(chainList,1:50000))
# 
# ## visualise
# # chainList.acPlot(chainList_thinned)
# # chainList.postPlot(chainList_thinned,1000)
# # chainList.bayesPlot(chainList_thinned)
# 
# ## 
# pdf("model_2_results_meanTL.pdf")
# 
# ## plot predictions temperature
# plot(X_obs[,3],Y_obs,pch=16,col=adjustcolor("black",alpha=0.5),xlab="temperature (SU)",ylab="mean trophic level (SU)")
# for(i in 1:2)
# {
# 	## predictions
# 	x        = seq(min(temp),max(temp),0.1)
# 	X_pred   = cbind(1,0,x,x^2,(i-1),(i-1)*x,0,0,0)
# 	pred     = chainList.apply(chainList_thinned,function(x_) Yhat(X_pred,x_[-1][4:12]))
# 	polygon(x=c(X_pred[,3],rev(X_pred[,3])),y=c(pred$f_q0.05,rev(pred$f_q0.95)),border=NA,col=adjustcolor(i+1,alpha=0.5))
# 	# lines(X_pred[,3],Yhat(X_pred,beta_map),col=i+1)
# 	lines(X_pred[,3],pred$f_mean,col=i+1)
# }
# legend("topright",legend=c("stream","lake"),lty=1,col=1:2+1,bty="n")
# 
# ## plot predictions temperature
# plot(X_obs[,7],Y_obs,pch=16,col=adjustcolor("black",alpha=0.5),xlab="dbo (SU)",ylab="mean trophic level (SU)")
# for(i in 1:2)
# {
# 	## predictions
# 	x        = seq(min(dbo,na.rm=T),max(dbo,na.rm=T),0.1)
# 	X_pred   = cbind(1,0,0,0,(i-1),0,x,x^2,(i-1)*x)
# 	pred     = chainList.apply(chainList_thinned,function(x_) Yhat(X_pred,x_[-1][4:12]))
# 	polygon(x=c(X_pred[,7],rev(X_pred[,7])),y=c(pred$f_q0.05,rev(pred$f_q0.95)),border=NA,col=adjustcolor(i+1,alpha=0.5))
# 	# lines(X_pred[,7],Yhat(X_pred,beta_map),col=i+1)
# 	lines(X_pred[,7],pred$f_mean,col=i+1)
# }
# legend("topright",legend=c("stream","lake"),lty=1,col=1:2+1,bty="n")
# 
# ## 
# dev.off()
# 
# #
# ###
