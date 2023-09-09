#####
## ##
#####

## goal:

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

###############
## FUNCTIONS ##
###############


## add_axis_and_grid
## Goal: Add custom axis and grid to plot given the vector of x and y data.
## Arguments:
## * x - vector - vector of x coordinates of the data to plot
## * x_ - vector - vector of standardised x coordinates of the data to plot 
## * y - vector - vector of y coordinates of the data to plot
## * y_ - vector - vector of standardised y coordinates of the data to plot
add_axes_and_grid = function(x, y, alpha)
{
  ## format data x
  alpha_x = alpha[1]
  lb_x = floor(min(x)/alpha_x)*alpha_x
  rb_x = ceiling(max(x)/alpha_x)*alpha_x
  dx = 2.5
  x = seq(lb_x, rb_x, dx * alpha_x)

  ## format data y
  alpha_y = alpha[2]
  lb_y = floor(min(y)/alpha_y)*alpha_y
  rb_y = ceiling(max(y)/alpha_y)*alpha_y
  dy = 2.5
  y = seq(lb_y, rb_y, dy * alpha_y)
    
  ## background
  coords = par("usr")
  coords_x = coords[1:2]
  coords_y = coords[3:4]
  polygon(x=c(coords_x, rev(coords_x)), y=c(c(coords_y[1],coords_y[1]), c(coords_y[2],coords_y[2])), col=adjustcolor("lightgrey",alpha=0.2), border=NA)
  
  ## grid guides
  for (l in 1:length(y)) lines(c(x[1]-10,x[length(x)]+10), c(y[l], y[l]), col="white")
  for (l in 1:length(x)) lines(c(x[l], x[l]), c(y[1]-10,y[length(y)]+10), col="white")
  
  ## x axis
  axis(1, label=x, at=x, lwd=0, lwd.ticks=1)
  axis(2, label=y, at=y, lwd=0, lwd.ticks=1)
}

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
legend("topright", "a.", bty="n")
chainList.bayesPlot(chainList_con)
legend("topright", "b.", bty="n")
par(mfrow=c(1,1))

## terminate
dev.off()

#
###
