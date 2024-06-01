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
## 05-05-2022 - created v0_3_1
##            - implemented random effects (variances) for different hydrographic regions
##            - added altitude as fixed effect
## 09-05-2022 - renamed file model 4 v0_1
##            - added storage of chains

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
source("hbm_functions_v0.5.r")

## load dataset
load("data/dataset_lake_stream_v0_1.rda")
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

#################################
## MODEL 4 - MAX TROPHIC LEVEL ##
#################################

## goal: infer missing values in addition to the linear model

## remove missing hydrographic bassin
data = data[-which(is.na(data$CdBH)),]

## response variable
Y = data$max_troph_lvl
# Y = data$w_trph_lvl_avg

## std function
std = function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)

## explanatory variables
year    = std(data$year)
temp    = std(data$temp)
dbo     = std(data$dbo)
alt     = std(data$altitude)
type    = (data$type == "lake")*1

## marix of explanatory variables
X_obs   = cbind(1,year,temp,temp^2,type,type*temp,dbo,dbo^2,type*dbo,alt)
X_mis   = cbind(1,year,temp,temp^2,type,type*temp,1,1,type,alt)
n_beta  = ncol(X_obs)

## missing values
idx_mis = which(is.na(dbo))
n_mis   = length(idx_mis)

## uneven variances (random effects)
stringToInteger = function(string)
{
	tmp = 1:length(unique(string))
	names(tmp) = unique(string)
	tmp = tmp[string]
	return(tmp)
}
n_sd_lik       = length(unique(data$CdBH))
idx_sd_lik_obs = stringToInteger(data$CdBH[-idx_mis])
idx_sd_lik_mis = stringToInteger(data$CdBH[idx_mis])

## splitting data into observed and missing
Y_obs   = Y[-idx_mis]
X_obs   = X_obs[-idx_mis,]
Y_mis   = Y[idx_mis]
X_mis   = X_mis[idx_mis,]

## predictive model
Yhat = function(X,beta) X%*%beta

## likelihood of response given observed variables
logLik_obs = function(Y_obs,X_obs,beta,sd_lik)
{
	return(sum(dnorm(Y_obs,Yhat(X_obs,beta),sd_lik[idx_sd_lik_obs],log=T)))
}

## likelihood of response given missing variables
logLik_mis = function(Y_mis,X_mis,x_mis,beta,sd_lik)      
{
	X_mis[,7:9] = X_mis[,7:9] * cbind(x_mis,x_mis^2,x_mis)
	logLik_mis  = sum(dnorm(Y_mis,Yhat(X_mis,beta),sd_lik[idx_sd_lik_mis],log=T))
	return(logLik_mis)
}

## probability density of missing variables
logLat_mis = function(x_mis,mu_mis,sd_mis)
{
	return(sum(dnorm(x_mis,mu_mis,sd_mis,log=T)))
}

## log priors
logPri     = function(beta,sd_lik,mu_mis,sd_mis) 
{
	logPri = sum(dnorm(beta,       0,10,log=T)) + 
             sum(dnorm(log(sd_lik),0,1, log=T)) + 
                 dnorm(mu_mis,     0,1, log=T)  + 
                 dnorm(log(sd_mis),0,1, log=T) 
	return(logPri)
}

## posterior density distribution
logPos     = function(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,mu_mis,sd_mis)
{  
	logPos = logLik_obs(Y_obs,X_obs,beta,sd_lik) + 
             logLik_mis(Y_mis,X_mis,x_mis,beta,sd_lik) + 
             logLat_mis(x_mis,mu_mis,sd_mis) + 
             logPri(beta,sd_lik,mu_mis,sd_mis)
	return(logPos)
}

## wrappers
logPosWrap  = function(omega) 
{
	sd_lik  = omega[1:n_sd_lik]
	sd_mis  = omega[n_sd_lik + 1:1]
	mu_mis  = omega[n_sd_lik + 1 + 1:1]
	beta    = omega[n_sd_lik + 1 + 1 + 1:n_beta] 
	x_mis   = omega[n_sd_lik + 1 + 1 + n_beta + 1:n_mis] 
	return(logPos(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,mu_mis,sd_mis))
}
wrap       = function(omega) c(log(omega[1:(n_sd_lik + 1)]),omega[-c(1:(n_sd_lik + 1))])
unwrap     = function(omega) c(exp(omega[1:(n_sd_lik + 1)]),omega[-c(1:(n_sd_lik + 1))])
dTarget    = function(omega) logPosWrap(unwrap(omega))

# ## pilot DEMCpp
# chainList = list()
# for(i in 1:3)
# {
# 	omega_0        = rnorm(n_sd_lik + 2 + n_beta + n_mis,0,1)
# 	chain          = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_0, "epsilon" = 0.001, "nIt" = 1000000))$chainList
# 	chain[,-1]     = t(apply(chain[,-1],1,unwrap))
# 	chainList[[i]] = chain
# }

## MaP
load("out/model_4_maxTL_chain.RData")
chainList = chainList_thinned
omega_map = wrap(chainList.argmaxPost(chainList))

## full DEMCpp
chainList = list()
for(i in 1:3)
{
	chain          = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_map, "epsilon" = 0.001, "nIt" = 1000000))$chainList
	chain[,-1]     = t(apply(chain[,-1],1,unwrap))
	chainList[[i]] = chain
}

## burn and thin
chainList_thinned = chainList.thin(chainList.burn(chainList,1:500000))
save(file="out/model_4_maxTL_chain.RData",chainList_thinned)
 
## visualise
# chainList.acPlot(chainList_thinned)
# chainList.postPlot(chainList_thinned,1000)
# chainList.bayesPlot(chainList_thinned)

## 
pdf("out/model_4_maxTL.pdf")

## plot predictions temperature
plot(X_obs[,3],Y_obs,pch=16,col=adjustcolor("black",alpha=0.5),xlab="temperature (SU)",ylab="max trophic level (SU)")
for(i in 1:2)
{
	## predictions
	x        = seq(min(temp),max(temp),0.1)
	X_pred   = cbind(1,0,x,x^2,(i-1),(i-1)*x,0,0,0,0)
	pred     = chainList.apply(chainList_thinned,function(x_) Yhat(X_pred,x_[-1][n_sd_lik + 1 + 1 + 1:n_beta]))
	polygon(x=c(X_pred[,3],rev(X_pred[,3])),y=c(pred$f_q0.05,rev(pred$f_q0.95)),border=NA,col=adjustcolor(i+1,alpha=0.5))
	# lines(X_pred[,3],Yhat(X_pred,beta_map),col=i+1)
	lines(X_pred[,3],pred$f_mean,col=i+1)
}
legend("topright",legend=c("stream","lake"),lty=1,col=1:2+1,bty="n")

## plot predictions temperature
plot(X_obs[,7],Y_obs,pch=16,col=adjustcolor("black",alpha=0.5),xlab="dbo (SU)",ylab="max trophic level (SU)")
for(i in 1:2)
{
	## predictions
	x        = seq(min(dbo,na.rm=T),max(dbo,na.rm=T),0.1)
	X_pred   = cbind(1,0,0,0,(i-1),0,x,x^2,(i-1)*x,0)
	pred     = chainList.apply(chainList_thinned,function(x_) Yhat(X_pred,x_[-1][n_sd_lik + 1 + 1 + 1:n_beta]))
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
# ## MODEL 4 - MEAN TROPHIC LEVEL ##
# ##################################
# 
# ## goal: infer missing values in addition to the linear model
# 
# ## remove missing hydrographic bassin
# data = data[-which(is.na(data$CdBH)),]
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
# alt     = std(data$altitude)
# type    = (data$type == "lake")*1
# 
# ## marix of explanatory variables
# X_obs   = cbind(1,year,temp,temp^2,type,type*temp,dbo,dbo^2,type*dbo,alt)
# X_mis   = cbind(1,year,temp,temp^2,type,type*temp,1,1,type,alt)
# n_beta  = ncol(X_obs)
# 
# ## missing values
# idx_mis = which(is.na(dbo))
# n_mis   = length(idx_mis)
# 
# ## uneven variances (random effects)
# stringToInteger = function(string)
# {
# 	tmp = 1:length(unique(string))
# 	names(tmp) = unique(string)
# 	tmp = tmp[string]
# 	return(tmp)
# }
# n_sd_lik       = length(unique(data$CdBH))
# idx_sd_lik_obs = stringToInteger(data$CdBH[-idx_mis])
# idx_sd_lik_mis = stringToInteger(data$CdBH[idx_mis])
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
# 	return(sum(dnorm(Y_obs,Yhat(X_obs,beta),sd_lik[idx_sd_lik_obs],log=T)))
# }
# 
# ## likelihood of response given missing variables
# logLik_mis = function(Y_mis,X_mis,x_mis,beta,sd_lik)      
# {
# 	X_mis[,7:9] = X_mis[,7:9] * cbind(x_mis,x_mis^2,x_mis)
# 	logLik_mis  = sum(dnorm(Y_mis,Yhat(X_mis,beta),sd_lik[idx_sd_lik_mis],log=T))
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
#              sum(dnorm(log(sd_lik),0,1, log=T)) + 
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
# 	sd_lik  = omega[1:n_sd_lik]
# 	sd_mis  = omega[n_sd_lik + 1:1]
# 	mu_mis  = omega[n_sd_lik + 1 + 1:1]
# 	beta    = omega[n_sd_lik + 1 + 1 + 1:n_beta] 
# 	x_mis   = omega[n_sd_lik + 1 + 1 + n_beta + 1:n_mis] 
# 	return(logPos(Y_obs,X_obs,Y_mis,X_mis,x_mis,beta,sd_lik,mu_mis,sd_mis))
# }
# wrap       = function(omega) c(log(omega[1:(n_sd_lik + 1)]),omega[-c(1:(n_sd_lik + 1))])
# unwrap     = function(omega) c(exp(omega[1:(n_sd_lik + 1)]),omega[-c(1:(n_sd_lik + 1))])
# dTarget    = function(omega) logPosWrap(unwrap(omega))
# 
# # ## pilot DEMCpp
# # chainList = list()
# # for(i in 1:3)
# # {
# # 	omega_0        = rnorm(n_sd_lik + 2 + n_beta + n_mis,0,1)
# # 	chain          = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_0, "epsilon" = 0.001, "nIt" = 1000000))$chainList
# # 	chain[,-1]     = t(apply(chain[,-1],1,unwrap))
# # 	chainList[[i]] = chain
# # }
# 
# ## MaP
# load("out/model_4_avgTL_chain.RData")
# chainList = chainList_thinned
# omega_map = wrap(chainList.argmaxPost(chainList))
# 
# ## full DEMCpp
# chainList = list()
# for(i in 1:3)
# {
# 	chain          = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_map, "epsilon" = 0.001, "nIt" = 1000000))$chainList
# 	chain[,-1]     = t(apply(chain[,-1],1,unwrap))
# 	chainList[[i]] = chain
# }
# 
# ## burn and thin
# chainList_thinned = chainList.thin(chainList.burn(chainList,1:500000))
# save(file="out/model_4_avgTL_chain.RData",chainList_thinned)
# 
# ## visualise
# # chainList.acPlot(chainList_thinned)
# # chainList.postPlot(chainList_thinned,1000)
# # chainList.bayesPlot(chainList_thinned)
# 
# ## 
# pdf("out/model_4_avgTL.pdf")
# 
# ## plot predictions temperature
# plot(X_obs[,3],Y_obs,pch=16,col=adjustcolor("black",alpha=0.5),xlab="temperature (SU)",ylab="mean trophic level (SU)")
# for(i in 1:2)
# {
# 	## predictions
# 	x        = seq(min(temp),max(temp),0.1)
# 	X_pred   = cbind(1,0,x,x^2,(i-1),(i-1)*x,0,0,0,0)
# 	pred     = chainList.apply(chainList_thinned,function(x_) Yhat(X_pred,x_[-1][n_sd_lik + 1 + 1 + 1:n_beta]))
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
# 	X_pred   = cbind(1,0,0,0,(i-1),0,x,x^2,(i-1)*x,0)
# 	pred     = chainList.apply(chainList_thinned,function(x_) Yhat(X_pred,x_[-1][n_sd_lik + 1 + 1 + 1:n_beta]))
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
