########
## .r ##
########

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
## 18-06-2022 - created v0_2
##            - implemented adaptive monte-carlo instead of differential-evolution monte-carlo
## 20-05-2022 - renamed file m_avgTL_v0_0
##            - created v0_1
##            - cleaned code

###############
## FUNCTIONS ##
###############

## std function
std = function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)

#
###

##############
## INITIATE ##
##############

## library
library(Rcpp)
sourceCpp("cpp/DEMCpp_v0.1.cpp")
source("f_AMC_v0_2.r")
source("f_HBM_v0_5.r")

## load dataset
load("data/dataset_lake_stream_v0_1.rda")
data = dataset_lake_stream

## path to out
pto = "out_avgTL"
system(paste("mkdir ",pto,sep=""))

## response
response = "Avg. Troph. level"

#
###

############
## MODEL  ##
############

## goal: infer missing values in addition to the linear model

## remove missing hydrographic bassin
data = data[-which(is.na(data$CdBH)),]

## response variable
# Y = data$max_troph_lvl
Y = data$w_trph_lvl_avg

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
                 dnorm(mu_mis,     0,.1, log=T)  + 
                 dnorm(sd_mis,     1,.01, log=T) 
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

#
###

# ###########
# ## TRAIN ##
# ###########
# 
# ## optimisation
# chain = NULL
# for(i in 1:3)
# {
# 	omega_0        = rnorm(n_sd_lik + 2 + n_beta + n_mis,0,0.1)
#     res            = optim(par=omega_0,fn=function(x)-dTarget(x),method="BFGS",control=list(trace=TRUE,maxit=100)) 
#     chain          = rbind(chain,c(res$value,unwrap(res$par)))
# }
# 
# ## MaP
# omega_map = wrap(chain[which.min(chain[,1]),-1])
# # load("out/model_4_maxTL_chain.RData")
# # omega_map = wrap(chainList.argmaxPost(chainList))
# 
# ## pilot chain
# chainList = list()
# for(i in 1:1)
# {
# 	omega_0        = omega_map # rnorm(n_sd_lik + 2 + n_beta + n_mis,0,.1)
# 	# chain          = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_0, "epsilon" = 0.001, "nIt" = 1000000))$chainList
#     chain          = AMC(dTarget   = dTarget,
#                      Theta_0   = omega_0,
#                      scale_p_0 = 4*2.38/sqrt(length(omega_0)),
#                      Sigma_p_0 = diag(1,length(omega_0)),
#                      itVect    =  c(10000,100000),
#                      adaptShape= T,
#                      msg       = T)[-1,]
# 	chain[,-1]     = t(apply(chain[,-1],1,unwrap))
# 	chainList[[i]] = chain
# }
# 
# ## burn and thin
# chainList_thinned = chainList.thin(chainList.burn(chainList,1:75000))
# save(file=paste(pto,"/chain.RData",sep=""),chainList_thinned)
#  
# #
# ###

#############
## FIGURES ##
#############

## goal:

## load chain
load(paste(pto,"/chain.RData",sep=""))

## 
# pdf(paste(pto,"/results.pdf",sep=""))

## PLOT PREDICTIONS TEMPERATURE ##
png(paste(pto,"/fig_1.png",sep=""))
#
plot(X_obs[,3],Y_obs,pch=16,col=adjustcolor("black",alpha=0.5),xlab="Temperature (SU)",ylab=response,main=paste(response," ~ temperature",sep=""))
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
#
legend("topright",legend=c("Stream","Lake"),lty=1,col=1:2+1,bty="n")
#
dev.off()

## PLOT PREDICTIONS TEMPERATURE ##
png(paste(pto,"/fig_2.png",sep=""))
#
plot(X_obs[,7],Y_obs,pch=16,col=adjustcolor("black",alpha=0.5),xlab="DBO (SU)",ylab=response,main=paste(response," ~ DBO",sep=""))
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
#
legend("topright",legend=c("Stream","Lake"),lty=1,col=1:2+1,bty="n")
#
dev.off()

## VISUALISE MISSING VS OBSERVED DBO ##
png(paste(pto,"/fig_3.png",sep=""))
#
x = density(dbo,na.rm=T)$x
y = density(dbo,na.rm=T)$y; y=y/max(y)
plot(x,y,type="l",col="white",xlab="DBO (SU)",ylab="Density (SU)",main="Observed vs missing DBO")
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col=adjustcolor("blue",0.4),border=NA)
#
dbo_mis = chainList.argmaxPost(chainList_thinned)[n_sd_lik + 1 + 1 + n_beta + 1:n_mis]
x = density(dbo_mis,na.rm=T)$x
y = density(dbo_mis,na.rm=T)$y; y=y/max(y)
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col=adjustcolor("red",0.4),border=NA)
#
legend("topright",legend=c("Observed","Missing"),col=adjustcolor(c("blue","red"),0.4),lty=1,bty="n")
#
dev.off()

## VERIFY MODEL ASSUMPTIONS ##
x_mis_      = apply(chainList_thinned[[1]][,-1][,n_sd_lik + 1 + 1 + n_beta + 1:n_mis],2,mean)
X_mis_      = X_mis
X_mis_[,7:9] = X_mis_[,7:9] * cbind(x_mis_,x_mis_^2,x_mis_)
chainList_  = list(chainList_thinned[[1]][,-1][,n_sd_lik + 1 + 1 + 1:n_beta])
Yhat_obs    = chainList.apply(chainList_,function(x)Yhat(X_obs,x))$f_mean
Yhat_mis    = chainList.apply(chainList_,function(x)Yhat(X_mis_,x))$f_mean
res_obs     = Y_obs - Yhat_obs
res_mis     = Y_mis - Yhat_mis

## HISTOGRAM OF RESIDUALS ##
png(paste(pto,"/fig_4.png",sep=""))
#
x = density(res_obs,na.rm=T)$x
y = density(res_obs,na.rm=T)$y; y=y/max(y)
plot(x,y,type="l",col="white",xlab="Residuals",ylab="Density (SU)",main="Observed vs missing DBO residuals")
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col=adjustcolor("blue",0.4),border=NA)
#
x = density(res_mis,na.rm=T)$x
y = density(res_mis,na.rm=T)$y; y=y/max(y)
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col=adjustcolor("red",0.4),border=NA)
#
legend("topright",legend=c("Observed","Missing"),col=adjustcolor(c("blue","red"),0.4),lty=1,bty="n")
#
dev.off()

## QQ plot
png(paste(pto,"/fig_5.png",sep=""))
#
sdVect = apply(chainList_thinned[[1]][,-1][,1:n_sd_lik],2,mean)
par(mfrow=c(3,3))
for(i in 1:n_sd_lik)
{
    x_mis_      = apply(chainList_thinned[[1]][,-1][,n_sd_lik + 1 + 1 + n_beta + 1:n_mis],2,mean)
    X_mis_      = X_mis
    X_mis_[,7:9] = X_mis_[,7:9] * cbind(x_mis_,x_mis_^2,x_mis_)
    chainList_  = list(chainList_thinned[[1]][,-1][,n_sd_lik + 1 + 1 + 1:n_beta])
    Yhat_obs    = chainList.apply(chainList_,function(x)Yhat(X_obs,x))$f_mean
    Yhat_mis    = chainList.apply(chainList_,function(x)Yhat(X_mis_,x))$f_mean
    res_obs     = Y_obs - Yhat_obs
    res_mis     = Y_mis - Yhat_mis
    #
    res_obs_th  = rnorm(length(res_obs),0,sdVect[i])
    res_mis_th  = rnorm(length(res_mis),0,sdVect[i])
    #
    plot(-1:1,xlim=c(-1,1)*1,ylim=c(-1,1)*1,xlab="Theoretical quantiles",ylab="Residuals",main=paste("qqplot for bassin ",i,sep=""),cex=0)
    lines(sort(res_obs_th),sort(res_obs),col=adjustcolor("blue",.4),type="p")
    lines(sort(res_mis_th),sort(res_mis),col=adjustcolor("red",.4),type="p")
    lines(-1:1,-1:1,lty=2)
    #
    legend("bottomright",legend=c("Observed","Missing"),col=adjustcolor(c("blue","red"),0.4),lty=1,bty="n")
}   
par(mfrow=c(1,1))
#
dev.off()

## VISUALISE PARAMETER POSTERIOR DISTRIBUTIONS ##
chain_           = cbind(chainList_thinned[[1]][,1],chainList_thinned[[1]][,-1][,n_sd_lik + 1 + 1 + 1:n_beta])
colnames(chain_) = c("P","1","year","temp","temp^2","type","type*temp","dbo","dbo^2","type*dbo","alt")
png(paste(pto,"/fig_6.png",sep="")); chainList.postPlot(list(chain_),1000); dev.off()
png(paste(pto,"/fig_7.png",sep="")); chainList.bayesPlot(list(chain_)); dev.off()
pdf(paste(pto,"/fig_8.pdf",sep="")); chainList.tracePlot(list(chain_)); dev.off()

## VISUALISE VARIANCES POSTERIOR DISTRIBUTIONS ##
chain_           = cbind(chainList_thinned[[1]][,1],chainList_thinned[[1]][,-1][,1:n_sd_lik])
colnames(chain_) = c("P",paste("sd_",1:n_sd_lik,sep=""))
png(paste(pto,"/fig_9.png",sep="")); chainList.postPlot(list(chain_),1000); dev.off()
png(paste(pto,"/fig_10.png",sep="")); chainList.bayesPlot(list(chain_)); dev.off()
pdf(paste(pto,"/fig_11.pdf",sep="")); chainList.tracePlot(list(chain_)); dev.off()

## VISUALISE MISSING MEAN VARIANCE POSTERIOR DISTRIBUTIONS ##
chain_           = cbind(chainList_thinned[[1]][,1],chainList_thinned[[1]][,-1][,n_sd_lik + 1:2])
colnames(chain_) = c("P","sd_mis","mu_mis")
png(paste(pto,"/fig_12.png",sep="")); chainList.postPlot(list(chain_),1000); dev.off()
png(paste(pto,"/fig_13.png",sep="")); chainList.bayesPlot(list(chain_)); dev.off()
pdf(paste(pto,"/fig_14.pdf",sep="")); chainList.tracePlot(list(chain_)); dev.off()

## VISUALISE MISSING OBSERVATIONS POSTERIOR DISTRIBUTIONS ##
chain_           = cbind(chainList_thinned[[1]][,1],chainList_thinned[[1]][,-1][,n_sd_lik + 1 + 1 + n_beta + 1:n_mis][,1:10])
colnames(chain_) = c("P",paste("mis_",1:10,sep=""))
png(paste(pto,"/fig_15.png",sep="")); chainList.postPlot(list(chain_),1000); dev.off()
png(paste(pto,"/fig_16.png",sep="")); chainList.bayesPlot(list(chain_)); dev.off()
pdf(paste(pto,"/fig_17.pdf",sep="")); chainList.tracePlot(list(chain_)); dev.off()

## 
# dev.off()

#
###
