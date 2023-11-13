#####
## ##
#####

## goal:

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

###############
## FUNCTIONS ##
###############

#
###

##############
## INITIATE ##
##############


#
###

#############
## FIGURES ##
#############

## goal:

## load module
source("m1_mTL_load.r")
# source("m1_con_load.r")

## load chains
chainList_ = list()
load(paste(pto,"/chain_thinned_",1,".RData",sep="")); chainList_[[1]] = chainList_thinned[[1]]
load(paste(pto,"/chain_thinned_",2,".RData",sep="")); chainList_[[2]] = chainList_thinned[[1]]
chainList_thinned = chainList_

## format chain
chainList_  = list()
for(i in 1:length(chainList_thinned))
{
    chain_ = cbind(chainList_thinned[[i]][,1], chainList_thinned[[i]][,-1][,idx_omega_beta])
    colnames(chain_) = c("P",colnames(X_obs))
    chainList_[[i]] = chain_    
}
chainList_mtl = chainList_

## load module
# source("m1_mTL_load.r")
source("m1_con_load.r")

## load chains
chainList_ = list()
load(paste(pto,"/chain_thinned_",1,".RData",sep="")); chainList_[[1]] = chainList_thinned[[1]]
load(paste(pto,"/chain_thinned_",2,".RData",sep="")); chainList_[[2]] = chainList_thinned[[1]]
chainList_thinned = chainList_

## format chain
chainList_  = list()
for(i in 1:length(chainList_thinned))
{
  chain_ = cbind(chainList_thinned[[i]][,1], chainList_thinned[[i]][,-1][,idx_omega_beta])
  colnames(chain_) = c("P",colnames(X_obs))
  chainList_[[i]] = chain_    
}
chainList_con = chainList_

## COMBINED BAYES PLOT ##
system("mkdir out")
pdf(paste("out/fig_bayesPlot_beta.pdf",sep=""), width=12, height=6) 

## plot
par(mfrow=c(1,2), mar=c(7,8,2,2), oma=c(2,4,2,2))
chainList.bayesPlot(chainList_mtl)
# legend("topright", "a.", bty="n")
## add the label above the plot
x <- par("usr")[2] - 0.05  # Adjust the x-coordinate as needed
y <- par("usr")[4] + 0.05  # Adjust the y-coordinate to position it above the plot
mtext(text = "a.", side = 3, line = 1, at = x, cex = 1.25)
#
chainList.bayesPlot(chainList_con)
# legend("topright", "b.", bty="n")
## add the label above the plot
x <- par("usr")[2] - 0.05  # Adjust the x-coordinate as needed
y <- par("usr")[4] + 0.05  # Adjust the y-coordinate to position it above the plot
mtext(text = "b.", side = 3, line = 1, at = x, cex = 1.25)
#
par(mfrow=c(1,1))

## terminate
dev.off()

#
###
