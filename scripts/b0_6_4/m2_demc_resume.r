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
sourceCpp("cpp/DEMCpp_v0.2.cpp")

## load data and model
source("m1_con_load.r")
# source("m1_mTL_load.r")

## chain name
chainName = "chain_thinned_1_resumed"
message(paste(pto,"/",chainName,".RData",sep=""))

## Load chains
load(paste(pto,"/chain_thinned_",1,".RData",sep="")); 
primerChain = chainList_thinned[[1]]
dim(primerChain)

## Wrap chain (it was unwrapped before saving)
primerChain[,-1] = t(apply(primerChain[,-1],1,wrap))

# ## Check
# primerChain_unwrapped = primerChain
# primerChain_wrapped = primerChain
# primerChain_wrapped[,-1] = t(apply(primerChain[,-1],1,wrap))
# primerChain_wrapped[,-1] = t(apply(primerChain_wrapped[,-1],1,unwrap))
# tail(primerChain_unwrapped[1,])
# tail(primerChain_wrapped[1,])
# ## OK wrapping unwrapping works

#
###

###########
## TRAIN ##
###########

##
## RESUME DEMC SAMPLING

## chain
nIt = 100000
chain = DEMCpp(list("dTarget" = dTarget, 
                    "epsilon" = 0.001, 
                    "nIt" = nIt, 
                    "primerChain" = primerChain))$chainList
chain[,-1] = t(apply(chain[,-1],1,unwrap))
chainList  = list(chain)

## burn and thin
# save(file=paste(pto,"/chain.RData",sep=""),chainList)
chainList_thinned = chainList.thin(chainList.burn(chainList,1:round(0.75*nIt)))
save(file=paste(pto,"/",chainName,".RData",sep=""),chainList_thinned)
 
#
###
