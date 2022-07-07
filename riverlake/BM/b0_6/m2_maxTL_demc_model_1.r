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

## load data and model
source("m1_maxTL_load_model_1.r")

#
###

###########
## TRAIN ##
###########

## optimisation
chain = NULL
for(i in 1:2)
{
	omega_0        = wrap(rPrior()) # rnorm(n_param,0,0.1)
    res            = optim(par=omega_0,fn=function(x)-dTarget(x),method="BFGS",control=list(trace=TRUE,maxit=100)) 
    chain          = rbind(chain,c(res$value,unwrap(res$par)))
}

## MaP
omega_map = wrap(chain[which.min(chain[,1]),-1])
# load("out/model_4_maxTL_chain.RData")
# omega_map = wrap(chainList.argmaxPost(chainList))

## pilot chain
chainList = list()
for(i in 1:1)
{
	omega_0        = omega_map # rnorm(n_param,0,.1)
	chain          = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_0, "epsilon" = 0.001, "nIt" = 1000))$chainList
	chain[,-1]     = t(apply(chain[,-1],1,unwrap))
	chainList[[i]] = chain
}

## burn and thin
chainList_thinned = chainList.thin(chainList.burn(chainList,1:100))
save(file=paste(pto,"/chain.RData",sep=""),chainList_thinned)
 
#
###
