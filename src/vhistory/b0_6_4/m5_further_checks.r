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

## load module
source("m1_con_load.r")
# source("m1_mTL_load.r")

#
###

#############
## FIGURES ##
#############

## goal:

## load chains
chainList_ = list()
load(paste(pto,"/chain_thinned_",1,".RData",sep="")); chainList_[[1]] = chainList_thinned[[1]]
load(paste(pto,"/chain_thinned_",2,".RData",sep="")); chainList_[[2]] = chainList_thinned[[1]]
chainList_thinned = chainList_

## SUMMARY TABLE ##
summaryTable_ = chainList.summaryTab(chainList_)[[1]]
summaryTable = cbind(rownames(summaryTable_),summaryTable_)
colnames(summaryTable) = c("name",colnames(summaryTable_))
# write.table(summaryTable,file=paste(pto,"/summary.csv",sep=""),sep=",",row.names=F,quote=F)
nscode = (summaryTable$signif[-1]=="*")*1

## VERIFY MODEL ASSUMPTIONS ##
x_mis_ = apply(chainList.unlist(chainList_thinned)[,-1][,idx_omega_xmis],2,mean)
X_mis_ = X_mis_l * X_mis_r(x_mis_)
chainList_ = list(chainList.unlist(chainList_thinned)[,-1][,idx_omega_beta])
Yhat_obs = chainList.apply(chainList_,function(x)Yhat(X_obs,x))$f_mean
Yhat_mis = chainList.apply(chainList_,function(x)Yhat(X_mis_,x))$f_mean
res_obs = Y_obs - Yhat_obs
res_mis = Y_mis - Yhat_mis
res = c(res_obs, res_mis)

## CHECK MULTIPLE SAMPLES PER SITE ##

## Preparing sites
sites = data$station
sites_obs = sites[-idx_mis]
sites_mis = sites[idx_mis]
par(mfrow=c(2,1))
barplot(table(sites_obs))
barplot(table(sites_mis))
par(mfrow=c(1,1))

## Getting multiple sample sites
idx_msample_obs = which(table(sites_obs) > 2)
idx_msample_mis = which(table(sites_mis) > 2)
sites_obs_msampled = names(table(sites_obs))[idx_msample_obs]
sites_mis_msampled = names(table(sites_mis))[idx_msample_mis]

## Fraction of multiple sampled sites
prop_msampled_obs = length(sites_obs_msampled)/length(table(sites_obs))
print(prop_msampled_obs)
prop_msampled_mis = length(sites_mis_msampled)/length(table(sites_mis))
print(prop_msampled_mis)
quantile(table(sites_obs))
quantile(table(sites_mis))

## Looking at residual time series for multiple sites
plot(1:10, cex=0, ylim=c(-3,3), xlim=c(0,12))
for (i in 1:100)
{
  col = rainbow(100)[round(runif(1,1,100))]
  site = sites_obs_msampled[i]
  res_obs_ = res_obs[which(sites_obs == site)]
  lines(1:length(res_obs_), res_obs_, col=col)
}

## CHECK AUTOCORRELATIONS ##

## Computing temporal ACs for multiple sites
cor_vector_obs = NULL
for (i in 1:length(sites_obs_msampled))
{
  ## Get site and observations
  site = sites_obs_msampled[i]
  res_obs_ = res_obs[which(sites_obs == site)]
  
  ## Compute ACs
  res_obs_l = c(res_obs_[length(res_obs_)], res_obs_[-length(res_obs_)]) # left tacking (lag 1)
  res_obs_ll = c(res_obs_l[length(res_obs_l)], res_obs_l[-length(res_obs_l)]) # left left tacking (lag 2)
  res_obs_lll = c(res_obs_ll[length(res_obs_ll)], res_obs_ll[-length(res_obs_ll)])
  res_obs_r = c(res_obs_[-1], res_obs_[1]) # right tacking (lag 1)
  res_obs_rr = c(res_obs_r[-1], res_obs_r[1]) # right right tacking (lag 2)
  res_obs_rrr = c(res_obs_rr[-1], res_obs_rr[1])
  cor_ = 1/3 * cor(res_obs_l, res_obs_r) + 1/3 * cor(res_obs_ll, res_obs_rr) + 1/3 * cor(res_obs_lll, res_obs_rrr)
  
  ## Collect
  cor_vector_obs = c(cor_vector_obs, cor_)
}
## Visualise
barplot(abs(cor_vector_obs), ylim=c(0,1))
lines(c(0, length(cor_vector_obs)*10), c(0.25,0.25), lty=2)
lines(c(0, length(cor_vector_obs)*10), -c(0.25,0.25), lty=2)
print(quantile(abs(cor_vector_obs), probs = c(0.975))) # ~ p(abs(rho) > 0.44) <= 0.025

#
###
