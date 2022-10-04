#####
## ##
#####

## goal: perform Bayesian data analysis of stream lake dataset 
##       to find relationship between temperature and DBO on network structure

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 08-07-2022 - created v0_0
## 11-07-2022 - created v0_1

###############
## FUNCTIONS ##
###############

#
###

##############
## INITIATE ##
##############

## libraries
library(Rcpp)

## load modules
sourceCpp("cpp/DEMCO_v0.5.cpp")
sourceCpp("cpp/DEMCpp_v0.1.cpp")

## load data and model
# source("m1_con_load.r")
source("m1_mTL_load.r")

## chain name
chainName = "chain_thinned_2"
message(paste(pto,"/",chainName,".RData",sep=""))

#
###

###########
## TRAIN ##
###########

##
## OPTIMISATION

## initiate from prior
k = 0
check = F
while(check == F)
{
	message(paste("attempt ",k))
	omega_0    = rPrior() # omega
	dTarget_   = dTarget(omega_0)
	if (dTarget_ > -7000)
	{
		check = T
		message("chain starting")
	}
	k = k + 1
}

## chain
nIt    = 1000
omega  = omega_0
for(i in 1:30)
{
    omega = DEMCOpp(list("dTarget" = dTarget, "Theta_0" = omega, "epsilon" = 0.001, "lambda" = 100, "nIt" = nIt))$chainList[nIt,-1]
}

##
## DEMC SAMPLING

## initiate from optimisation
omega_0 = omega

# ## initiate omega from pilot chain 
# load(paste(pto,"/chain_thinned.RData",sep=""))
# omega_0 = omega = wrap(chainList.argmaxPost(chainList_thinned))

# ## initiate from prior
# k = 0
# check = F
# while(check == F)
# {
# 	message(paste("attempt ",k))
# 	omega_0    = rPrior() # omega
# 	dTarget_   = dTarget(omega_0)
# 	if (dTarget_ > -Inf)
# 	{
# 		check = T
# 		message("chain starting")
# 	}
# 	k = k + 1
# }

## chain
nIt        = 5000000
chain      = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_0, "epsilon" = 0.001, "nIt" = nIt))$chainList
chain[,-1] = t(apply(chain[,-1],1,unwrap))
chainList  = list(chain)

## burn and thin
# save(file=paste(pto,"/chain.RData",sep=""),chainList)
chainList_thinned = chainList.thin(chainList.burn(chainList,1:round(0.75*nIt)))
save(file=paste(pto,"/",chainName,".RData",sep=""),chainList_thinned)
 
#
###
