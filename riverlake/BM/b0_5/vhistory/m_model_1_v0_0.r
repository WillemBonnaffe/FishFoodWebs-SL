################
## analysis.r ##
################

## goal: perform Bayesian data analysis of stream lake dataset 
##       to find relationship between temperature and DBO on network structure

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 20-03-2022 - created v0_0
## 27-03-2022 - created v0_1
## 09-05-2022 - renamed model_1_v0_0

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
load("data/dataset_lake_stream.rda")
data = dataset_lake_stream

## results file
pdf("out/results.pdf")

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
Y = data$max_troph_lvl

## explanatory
temp   = data$temp
type   = (data$type == "lake")*1
year   = data$year - min(data$year)
X      = cbind(1,temp,temp^2,type,type*temp,year)

## predictive model
Yhat = function(X,beta) X%*%beta

## log posterior
logLik     = function(X,Y,beta,sd_lik)        sum(dnorm(Y-Yhat(X,beta),0,sd_lik,log=T))
logPri     = function(beta,sd_lik,sd_pri)     dnorm(log(sd_lik),0,sd_pri[1],log=T) + sum(dnorm(beta,0,sd_pri[2]),log=T)
logPos     = function(X,Y,beta,sd_lik,sd_pri) logLik(X,Y,beta,sd_lik) + logPri(beta,sd_lik,sd_pri)

## wrappers
logPosWrap = function(omega) logPos(X,Y,omega[-1],omega[1],sd_pri)
map        = function(omega) c(log(omega[1]),omega[-1])
unmap      = function(omega) c(exp(omega[1]),omega[-1])
dTarget    = function(omega) logPosWrap(unmap(omega))

## hyperparameters
sd_pri = c(1.0,10.0)

## DEMCpp
chainList = list()
for(i in 1:3)
{
	omega_0    = rnorm(1+ncol(X),0,1)
	chain      = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_0, "epsilon" = 0.001, "nIt" = 100000))$chainList
	chain[,-1] = t(apply(chain[,-1],1,unmap))
	chainList[[i]] = chain
}

## MaP
beta_map = chainList.argmaxPost(chainList)[-1]

## burn and thin
chainList_thinned = chainList.thin(chainList.burn(chainList,1:50000))

## visualise
# chainList.acPlot(chainList_thinned)
chainList.postPlot(chainList_thinned,1000)
chainList.bayesPlot(chainList_thinned)

## plot predictions
plot(X[,2],Y,pch=16,col=adjustcolor("black",alpha=0.5))
for(i in 1:2)
{
	## predictions
	x        = seq(min(temp),max(temp),0.1)
	X_pred   = cbind(1,x,x^2,(i-1),(i-1)*x,0)
	pred     = chainList.apply(chainList_thinned,function(x) Yhat(X_pred,x[-c(1:2)]))
	polygon(x=c(X_pred[,2],rev(X_pred[,2])),y=c(pred$f_q0.05,rev(pred$f_q0.95)),border=NA,col=adjustcolor(i+1,alpha=0.5))
	lines(X_pred[,2],Yhat(X_pred,beta_map),col=i+1)
}
legend("topright",legend=c("stream","lake"),lty=1,col=1:2+1,bty="n")

#
###

###############
## TERMINATE ##
###############

##
dev.off()

#
###
