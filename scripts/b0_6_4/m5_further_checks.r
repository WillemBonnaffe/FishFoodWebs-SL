#####
## ##
#####

## goal:

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 07-07-2022 - created v0_0
## 11-07-2022 - created v0_1
##            - split missing covariate matrix into left and right component
##            - added plot of the interaction
## 10-09-2022 - created v0_2
##            - updated figures
## 25-08-2023 - created v0_3
##            - improved figures

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

## VISUALISE PARAMETER POSTERIOR DISTRIBUTIONS ##
chainList_  = list()
for(i in 1:length(chainList_thinned))
{
    chain_ = cbind(chainList_thinned[[i]][,1], chainList_thinned[[i]][,-1][,idx_omega_beta])
    colnames(chain_) = c("P",colnames(X_obs))
    chainList_[[i]] = chain_    
}
pdf(paste(pto,"/fig_postPlot_beta.pdf",sep="")); chainList.postPlot(chainList_, 1000, use_labels=F); dev.off()
pdf(paste(pto,"/fig_bayesPlot_beta.pdf",sep="")); chainList.bayesPlot(chainList_); dev.off()
pdf(paste(pto,"/fig_tracePlot_beta.pdf",sep="")); chainList.tracePlot(chainList_); dev.off()

## SUMMARY TABLE ##
summaryTable_ = chainList.summaryTab(chainList_)[[1]]
summaryTable = cbind(rownames(summaryTable_),summaryTable_)
colnames(summaryTable) = c("name",colnames(summaryTable_))
write.table(summaryTable,file=paste(pto,"/summary.csv",sep=""),sep=",",row.names=F,quote=F)
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

## Fraction of multiply sampled sites
prop_msampled_obs = length(sites_obs_msampled)/length(table(sites_obs))
print(prop_msampled_obs)
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
print(quantile(abs(cor_vector_obs), probs = c(0.975))) # ~ p(abs(rho) > 0.44) <= 0.05

## HISTOGRAM OF RESIDUALS ##
pdf(paste(pto,"/fig_hist_residuals.pdf",sep=""))
#
## plot density observed
x = density(res_obs,na.rm=T)$x
y = density(res_obs,na.rm=T)$y; y=y/max(y)
plot(x,y,type="l",col="white",xlab="Residuals",ylab="Density (SU)", xaxt="n", yaxt="n", bty="l")
add_axes_and_grid(x, y, alpha=c(1, .1))
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col=adjustcolor("blue",0.4),border=NA)
#
## plot density missing
x = density(res_mis,na.rm=T)$x
y = density(res_mis,na.rm=T)$y; y=y/max(y)
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col=adjustcolor("red",0.4),border=NA)
#
## legend
legend("topright", legend=c("Observed BOD", "Missing BOD"), col=adjustcolor(c("blue","red"),0.4), pch=15, bty="n")
#
dev.off()

## QQ PLOT - v0_2 ##
pdf(paste(pto,"/fig_qqplot_residuals.pdf",sep=""))
par(mfrow=c(1,1), mar=c(3,3,2,2), oma=c(2,2,2,2), xpd=NA)

## compute parameters
sdVect = apply(chainList.unlist(chainList_thinned)[,-1][,idx_omega_sd_lik],2,mean)
rho = apply(chainList.unlist(chainList_thinned)[,-1][,idx_omega_rho],2,mean)

## compute distance matrix
long_obs = long[-idx_mis]
long_mis = long[ idx_mis]
latt_obs = latt[-idx_mis]
latt_mis = latt[ idx_mis]
x_       = c(long_obs,long_mis)
y_       = c(latt_obs,latt_mis)
DM        = matrix(rep(0,length(x_)^2),ncol=length(x_),nrow=length(x_))
for(i in 1:length(x_))
{
  for(j in 1:length(y_))
  {
    DM[i,j] = sqrt((x_[i] - x_[j])^2 + (y_[i] - y_[j])^2)
  }
}

## compute sigma
Sigma_ = Sigma(sdVect[idx_sd_lik], rho, DM)

## compute theoretical quantiles
res_th = rmvnorm(n=1, mean=rep(0, n_data), sigma=Sigma_)

## plot
par(xpd=NA)
plot(-1:1, xlim=c(-1,1)*4*sd(res_th), ylim=c(-1,1)*4*sd(res), xlab="Theoretical quantiles", ylab="Residuals", cex=0, xaxt="n", yaxt="n", bty="l")
par(xpd=F)

## grid
x = res_th
y = res
add_axes_and_grid(x, y, alpha=c(1, 1))

## lines
lines(sort(res_th), sort(res), type="p", pch=16, col=gray(runif(n_data, .25, 1), alpha=0.25))
lines((-1:1)*4*sd(res_th),(-1:1)*4*sd(res_th),lty=2)

## legend
# legend("bottomright",legend=c(paste("Bassin ",i,sep=""),"Observed BOD","Missing BOD"),col=adjustcolor(c("white","blue","red"),0.4), pch=1, bty="n")

par(mfrow=c(1,1))
dev.off()

# ## QQ PLOT - v0_1 ##
# pdf(paste(pto,"/fig_qqplot_residuals.pdf",sep=""))
# par(mfrow=c(3,3), mar=c(3,3,2,2), oma=c(2,2,2,2), xpd=NA)
#
# ## compute parameters
# sdVect = apply(chainList.unlist(chainList_thinned)[,-1][,idx_omega_sd_lik],2,mean)
# rho = apply(chainList.unlist(chainList_thinned)[,-1][,idx_omega_rho],2,mean)
#
# ## compute distance matrix
# long_obs = long[-idx_mis]
# long_mis = long[ idx_mis]
# latt_obs = latt[-idx_mis]
# latt_mis = latt[ idx_mis]
# x_       = c(long_obs,long_mis)
# y_       = c(latt_obs,latt_mis)
# DM        = matrix(rep(0,length(x_)^2),ncol=length(x_),nrow=length(x_))
# for(i in 1:length(x_))
# {
#   for(j in 1:length(y_))
#   {
#     DM[i,j] = sqrt((x_[i] - x_[j])^2 + (y_[i] - y_[j])^2)
#   }
# }
#
# ## compute sigma
# Sigma_ = Sigma(sdVect[idx_sd_lik], rho, DM)
#
# ## compute theoretical quantiles
# res_th = rmvnorm(n=1, mean=rep(0, n_data), sigma=Sigma_)
# res_obs_th = res_th[-idx_mis]
# res_mis_th = res_th[idx_mis]
#
# ## plot
# par(xpd=NA)
# plot(-1:1, xlim=c(-1,1)*4*sd(res_obs_th), ylim=c(-1,1)*4*sd(res_obs), xlab="Theoretical quantiles", ylab="Residuals", cex=0, xaxt="n", yaxt="n", bty="l")
# par(xpd=F)
#
# ## grid
# x = res_obs_th
# y = res_obs
# add_axes_and_grid(x, y, alpha=c(1, 1))
#
# ## lines
# lines(sort(res_th),sort(res),col=adjustcolor("black",.4),type="p")
# lines(sort(res_obs_th),sort(res_obs),col=adjustcolor("blue",.4),type="p")
# lines((-1:1)*4*sd(res_obs_th),(-1:1)*4*sd(res_obs),lty=2)
# lines(sort(res_mis_th),sort(res_mis),col=adjustcolor("red",.4),type="p")
# lines((-1:1)*4*sd(res_mis_th),(-1:1)*4*sd(res_mis),lty=2)
#
# ## legend
# legend("bottomright",legend=c(paste("Bassin ",i,sep=""),"Observed BOD","Missing BOD"),col=adjustcolor(c("white","blue","red"),0.4), pch=1, bty="n")
#
# par(mfrow=c(1,1))
# dev.off()

# ## QQ PLOT - v0_0 ##
# pdf(paste(pto,"/fig_qqplot_residuals.pdf",sep=""))
# par(mfrow=c(3,3), mar=c(3,3,2,2), oma=c(2,2,2,2), xpd=NA)
# sdVect = apply(chainList.unlist(chainList_thinned)[,-1][,idx_omega_sd_lik],2,mean)
# for(i in 1:n_sd_lik)
# {
#
#     ## compute theoretical quantiles
#     res_obs_th = rnorm(length(res_obs),0,sdVect[i])
#     res_mis_th = rnorm(length(res_mis),0,sdVect[i])
#
#     ## plot
#     par(xpd=NA)
#     plot(-1:1, xlim=c(-1,1)*4*sd(res_obs_th), ylim=c(-1,1)*4*sd(res_obs), xlab="Theoretical quantiles", ylab="Residuals", cex=0, xaxt="n", yaxt="n", bty="l")
#     par(xpd=F)
#
#     ## grid
#     x = res_obs_th
#     y = res_obs
#     add_axes_and_grid(x, y, alpha=c(1, 1))
#
#     ## lines
#     lines(sort(res_obs_th),sort(res_obs),col=adjustcolor("blue",.4),type="p")
#     lines(sort(res_mis_th),sort(res_mis),col=adjustcolor("red",.4),type="p")
#     lines((-1:1)*4*sd(res_obs_th),(-1:1)*4*sd(res_obs),lty=2)
#
#     ## legend
#     legend("bottomright",legend=c(paste("Bassin ",i,sep=""),"Observed BOD","Missing BOD"),col=adjustcolor(c("white","blue","red"),0.4), pch=1, bty="n")
# }
# par(mfrow=c(1,1))
# dev.off()

## VISUALISE VARIANCES POSTERIOR DISTRIBUTIONS ##
chainList_  = list()
for(i in 1:length(chainList_thinned))
{
    chain_           = cbind(chainList_thinned[[i]][,1],chainList_thinned[[i]][,-1][,idx_omega_sd_lik])
    colnames(chain_) = c("P",paste("sd_",1:n_sd_lik,sep=""))
    chainList_[[i]]  = chain_
}
pdf(paste(pto,"/fig_postPlot_sd_lik.pdf",sep="")); chainList.postPlot(chainList_,1000); dev.off()
pdf(paste(pto,"/fig_bayesPlot_sd_lik.pdf",sep="")); chainList.bayesPlot(chainList_); dev.off()
pdf(paste(pto,"/fig_tracePlot_sd_lik.pdf",sep="")); chainList.tracePlot(chainList_); dev.off()

## VISUALISE MISSING MEAN VARIANCE POSTERIOR DISTRIBUTIONS ##
chainList_  = list()
for(i in 1:length(chainList_thinned))
{
    chain_ = cbind(chainList_thinned[[i]][,1],chainList_thinned[[i]][,-1][,c(idx_omega_mu_mis,idx_omega_sd_mis)])
    colnames(chain_) = c("P","mu_mis","sd_mis")
    chainList_[[i]] = chain_
}
pdf(paste(pto,"/fig_postPlot_sd_mis.pdf",sep="")); chainList.postPlot(chainList_,1000); dev.off()
pdf(paste(pto,"/fig_bayesPlot_sd_mis.pdf",sep="")); chainList.bayesPlot(chainList_); dev.off()
pdf(paste(pto,"/fig_tracePlot_sd_mis.pdf",sep="")); chainList.tracePlot(chainList_); dev.off()

## VISUALISE MISSING OBSERVATIONS POSTERIOR DISTRIBUTIONS ##
chainList_  = list()
for(i in 1:length(chainList_thinned))
{
    chain_           = cbind(chainList_thinned[[i]][,1],chainList_thinned[[i]][,-1][,idx_omega_xmis][,1:10])
    colnames(chain_) = c("P",paste("mis_",1:10,sep=""))
    chainList_[[i]]  = chain_
}
pdf(paste(pto,"/fig_postPlot_missing_bod.pdf",sep="")); chainList.postPlot(chainList_,1000); dev.off()
pdf(paste(pto,"/fig_bayesPlot_missing_bod.pdf",sep="")); chainList.bayesPlot(chainList_); dev.off()
pdf(paste(pto,"/fig_tracePlot_missing_bod.pdf",sep="")); chainList.tracePlot(chainList_); dev.off()

## VISUALISE CORRELATIONS POSTERIOR DISTRIBUTIONS ##
chainList_  = list()
for(i in 1:length(chainList_thinned))
{
    chain_           = cbind(chainList_thinned[[i]][,1],chainList_thinned[[i]][,-1][,idx_omega_rho])
    colnames(chain_) = c("P",paste("rho_",1:2,sep=""))
    chainList_[[i]]  = chain_
}
pdf(paste(pto,"/fig_postPlot_rho.pdf",sep="")); chainList.postPlot(chainList_,1000); dev.off()
pdf(paste(pto,"/fig_bayesPlot_rho.pdf",sep="")); chainList.bayesPlot(chainList_); dev.off()
pdf(paste(pto,"/fig_tracePlot_rho.pdf",sep="")); chainList.tracePlot(chainList_); dev.off()

## COMPUTE SPATIAL CORRELATIONS IN RESIDUALS ##
long_obs = long[-idx_mis]
long_mis = long[ idx_mis]
latt_obs = latt[-idx_mis]
latt_mis = latt[ idx_mis]
x_       = c(long_obs,long_mis)
y_       = c(latt_obs,latt_mis)
#
## compute distance matrix
D        = matrix(rep(0,length(x_)^2),ncol=length(x_),nrow=length(x_))
for(i in 1:length(x_))
{
    for(j in 1:length(y_))
    {
        D[i,j] = sqrt((x_[i] - x_[j])^2 + (y_[i] - y_[j])^2)
    }
}
res_ = c(res_obs,res_mis)
#
## compute correlation between residuals with distance
rho_    = NULL
d_      = NULL
for(i in 1:100)
{
    idx   = order(D[i,])
    res_i = res_[idx]
    x_i   = x_[idx]
    y_i   = y_[idx]
    rho_i = NULL
    d_i   = NULL
    for(j in c(seq(1,10,1),seq(10,100,10),seq(100,2000,100)))
    {
        ## correlation
        res_il  = c(res_i,rep(NA,j))
        res_ir  = c(rep(NA,j),res_i)
        s       = !is.na(res_il*res_ir)
        rho_ij  = cor(res_il[s],res_ir[s])
        #
        ## distance
        x_il   = c(x_i,rep(NA,j))
        x_ir   = c(rep(NA,j),x_i)
        y_il   = c(y_i,rep(NA,j))
        y_ir   = c(rep(NA,j),y_i)
        s      = !is.na(x_il*x_ir)
        d_ij   = mean(sqrt((x_il[s]-x_ir[s])^2 + (y_il[s]-y_ir[s])^2))
        #
        ## concatenate
        rho_i = c(rho_i,rho_ij)
        d_i   = c(  d_i,  d_ij)
    }
    rho_ = rbind(rho_,rho_i)
    d_   = rbind(d_,d_i)
}
rho_mean = apply(rho_,2,mean)
rho_sd   = apply(rho_,2,sd)
d_mean   = apply(  d_,2,mean)
#
## visualise correlation with distance
pdf(paste(pto,"/fig_spatial_autocorrelations.pdf",sep=""));
x = d_mean
y = rho_mean
plot(x, y, xlim=c(min(D),max(D)), ylim=c(0,1), xaxt="n", yaxt="n", bty="l")
add_axes_and_grid(x, y, alpha=c(1, 0.1))
polygon(x=c(d_mean,rev(d_mean)),y=c(rho_mean+2*rho_sd,rev(rho_mean-2*rho_sd)),border=NA,col=grey(0.5,alpha=0.25))
lines(d_mean,rho_mean,col="red")
dev.off()

#
###
