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
source("m1_con_load.r")

#
###

###########
## TRAIN ##
###########

# ## optimisation
# omega_0 = rPrior() 
# #
# ## initiate omega with previous results
# load("/Users/willembonnaffe/Documents/GitHub/RESOTRO/riverlake/BM/b0_5_1/out_conne_model_0_v0_4/chain_thinned.RData")
# omega_map = chainList.argmaxPost(chainList_thinned)
# omega_0[idx_omega_sd_lik] = omega_map[idx_omega_sd_lik]
# omega_0[idx_omega_sd_mis] = omega_map[idx_omega_sd_mis]
# omega_0[idx_omega_mu_mis] = omega_map[idx_omega_mu_mis]
# omega_0[idx_omega_beta]   = omega_map[idx_omega_beta]
# #
# omega = omega_0 = wrap(omega_0)
# #
# ## chain
# nIt     = 100
# for(i in 1:10)
# {
#     omega = DEMCOpp(list("dTarget" = dTarget, "Theta_0" = omega, "epsilon" = 0.01, "lambda" = 100, "nIt" = nIt))$chainList[nIt,-1]
# }

## DEMC sampling
#
## initiate omega from pilot chain 
load(paste(pto,"/chain_thinned.RData",sep=""))
omega = wrap(chainList.argmaxPost(chainList_thinned))
#
omega_0    = omega
#
## chain
nIt        = 100000
chain      = DEMCpp(list("dTarget" = dTarget, "Theta_0" = omega_0, "epsilon" = 0.001, "nIt" = nIt))$chainList
chain[,-1] = t(apply(chain[,-1],1,unwrap))
chainList  = list(chain)

## burn and thin
# save(file=paste(pto,"/chain.RData",sep=""),chainList)
chainList_thinned = chainList.thin(chainList.burn(chainList,1:round(0.75*nIt)))
save(file=paste(pto,"/chain_thinned.RData",sep=""),chainList_thinned)
 
#
###
