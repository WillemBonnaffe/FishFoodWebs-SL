#####################
## hbm_functions.r ##
#####################

## update: 
## 13-03-2019 - cleaned
## 28-03-2019 - added extra functions (e.g. chainList.write(), ...) + normalized function names
## 08-04-2019 - added extra autocor function
## 20-04-2019 - added/corrected the unlist function
## 21-04-2019 - added a polygon.gamma/truncnorm functions for plotting purposes (might move it to a graphical display thingy)
## 27-03-2022 - created v0_5 
##            - introduced function argmax
## 21-06-2022 - created v0_6
##            - introduced significance in summary tab
## 25-08-2023 - improved chainList.BayesPlot function

## goal: function for bayesian analysis

#####################
## chainList.apply ##
#####################

chainList.apply = function(chainList,f)
{
	chain = chainList.unlist(chainList)
	f_ensemble  = t(apply(chain,1,f))
	f_q0.05 = apply(f_ensemble,2,quantile,p=c(0.05))
	f_q0.95 = apply(f_ensemble,2,quantile,p=c(0.95))
	f_mean  = apply(f_ensemble,2,mean)
	return(list("f_mean"     = f_mean,
				"f_q0.05"    = f_q0.05,
				"f_q0.95"    = f_q0.95,
				"f_ensemble" = f_ensemble))
}

#
###

##########################
## chainList.argmaxPost ##
##########################

## args:
# @ chainList - list - a list of chains

chainList.argmaxPost = function(chainList)
{
	chainTab = chainList.unlist(chainList)
	MaP       = chainTab[which.max(chainTab[,1]),-1]
	return(MaP)
}

#
###

######################
## chainList.unlist ##
######################

## goal: turns a chainList into a tab object

## args:
# @ chainList - list - a list of chains

chainList.unlist <- function(chainList)
{
  
  ##
  chainTab <- NULL
  
  ##
  for(i in 1:length(chainList))
  {
    chainTab <- rbind(chainTab, chainList[[i]])
  }
  
  ##
  colnames(chainTab) <- colnames(chainList[[1]])  
  
  return(chainTab)
  
}

#
###

#####################
## chainList.write ##
#####################

## goal: store a chainList object in a folder

## args:
# @ chainList - list - list of chains
# @ pathToFolder - string - path to folder containing/should contain the chains folder

chainList.write <- function(chainList, pathToFolder)
{
  

  ## create folder
  system(paste("mkdir ", pathToFolder, "/chains", sep=""))
  
  ## path to chains folder
  pathToChainsFolder <- paste(pathToFolder,"/chains",sep="")
    
  ## check if there exists chains already
  chainIndexVect <- 1:length(chainList) + length(list.files(pathToChainsFolder))
  
  ## for each chains
  for(i in 1:length(chainIndexVect))
  {
    
    ##
    system(paste("mkdir ",pathToChainsFolder,"/chain_", chainIndexVect[i], sep=""))
    write.table(chainList[[i]], paste(pathToChainsFolder,"/chain_", chainIndexVect[i],"/chain_", chainIndexVect[i], ".csv", sep=""), sep=";")
    # write(x = paste(names(mcmcList$arguments),": ", mcmcList$arguments, "\n"), file = paste("out/chains/chain_", chainIndexVect[i], "/chain_", chainIndexVect[i], "_args.txt", sep=""))
    
  }     

}

#
###

####################
## chainList.read ##
####################

## goal: read a list of chains stored in chain folder

## args:
# @ pathToFolder - string - path to the folder containing the chain folder

chainList.read <- function(pathToFolder)
{
  
  ## path to chains folder
  pathToChainsFolder <- paste(pathToFolder,"/chains",sep="")
  
  ## initiate
  chainList <- list()
  
  ## discard outlier chains denoted by #
  if(length(grep("#",list.files(pathToChainsFolder)))>0)
  {

    fileVect <- list.files(pathToChainsFolder)[-grep("#",list.files(pathToChainsFolder))]
    
  }else
  {
    fileVect <- list.files(pathToChainsFolder)
  } 
  
  ## read chains
  for(i in 1:length(fileVect))
  {
    
    ##
    chainList[[i]] <- read.table(paste(pathToChainsFolder,"/",fileVect[i],"/",fileVect[i],".csv",sep=""), sep=";", header=T)
    
  }
  
  return(chainList)
  
}

#
###

####################
## chainList.burn ##
####################

## goal: takes a chainlist and returns a burned chain list

## args:
# @ chainList - a list of chains
# @ burnin - a vector of iterations to remove

chainList.burn <- function(chainList, burnin)
{
  
  ## burn
  for(l in 1:length(chainList)){chainList[[l]] <- chainList[[l]][-burnin,]}
  
  ## return chainList
  return(chainList)
  
}

#
###

####################
## chainList.thin ##
####################

## goal: takes a list of chains and return a thinned list of chains

## args:
# @ chainList - a list of chains

chainList.thin <- function(chainList)
{
  
  ## thin
  thinin <- seq(1,nrow(chainList[[1]]),nrow(chainList[[1]])/1000); 
  for(l in 1:length(chainList)){chainList[[l]] <- chainList[[l]][thinin,]}
  
  ## return chainList
  return(chainList)
  
}

#
###

#######################
## chainList.summary ##
#######################

## goal: takes a list of chains and return summary statistics

## args:
# @chainList - list - list of chain objects (tables)

chainList.summary <- function(chainList)
{
  
  ## cut chains in half
  chainList_final <- list() ; j <- 1
  for(i in 1:length(chainList))
  {
    chainList_final[[j]] <- chainList[[i]][1:(nrow(chainList[[i]])/2),]; j <- j+1
    chainList_final[[j]] <- chainList[[i]][(nrow(chainList[[i]])/2):nrow(chainList[[i]]),]; j <- j+1
  }
  chainList <- chainList_final
  
  ## unlist chains
  chains <- NULL ; for(i in 1:length(chainList)){chains <- rbind(chains, chainList[[i]])}
  
  K <- nrow(chainList[[1]])
  w <- apply(X=matrix(data=unlist(lapply(X=chainList, FUN=function(x){apply(X=x, MARGIN=2, FUN=var)})), byrow=T, ncol=ncol(chainList[[1]])), MARGIN=2, FUN=mean)
  b <- K * apply(X=matrix(data=unlist(lapply(X=chainList, FUN=function(x){apply(X=x, MARGIN=2, FUN=mean)})), byrow=T, ncol=ncol(chainList[[1]])), MARGIN=2, FUN=var)
  r_hat <- sqrt(1 + 1/K * (b/w - 1))
  # p <- apply(chains,2,FUN=function(X){length(X[X>0])/length(X[X<0])})
	
	## estimates table
	estimatesTab = cbind(
    MaP=as.vector(t(chains[which.max(chains[,1]),])),
    mean=apply(X=chains, FUN=mean, MARGIN=2),
    sd=apply(X=chains, FUN=sd, MARGIN=2),
    t(apply(X=chains, FUN=quantile, MARGIN=2, probs=c(0.05,0.5,0.95))),
    r_hat)
	# signif =  gsub("1",x=gsub(pattern="0",x=t(apply(estimatesTab[,c(4,6)]>0,1,diff)),replacement="*"),"ns")
	# estimatesTab = data.frame(cbind(round(estimatesTab,4),signif))

    return(list(estimates=estimatesTab,
    covmat=cov(chains[,-1])
  ))
  
}

#
###

##########################
## chainList.summaryTab ##
##########################

## goal: takes a list of chains and return summary statistics (including significance level)

## args:
# @chainList - list - list of chain objects (tables)

chainList.summaryTab <- function(chainList)
{
  
  ## cut chains in half
  chainList_final <- list() ; j <- 1
  for(i in 1:length(chainList))
  {
    chainList_final[[j]] <- chainList[[i]][1:(nrow(chainList[[i]])/2),]; j <- j+1
    chainList_final[[j]] <- chainList[[i]][(nrow(chainList[[i]])/2):nrow(chainList[[i]]),]; j <- j+1
  }
  chainList <- chainList_final
  
  ## unlist chains
  chains <- NULL ; for(i in 1:length(chainList)){chains <- rbind(chains, chainList[[i]])}
  
  K <- nrow(chainList[[1]])
  w <- apply(X=matrix(data=unlist(lapply(X=chainList, FUN=function(x){apply(X=x, MARGIN=2, FUN=var)})), byrow=T, ncol=ncol(chainList[[1]])), MARGIN=2, FUN=mean)
  b <- K * apply(X=matrix(data=unlist(lapply(X=chainList, FUN=function(x){apply(X=x, MARGIN=2, FUN=mean)})), byrow=T, ncol=ncol(chainList[[1]])), MARGIN=2, FUN=var)
  r_hat <- sqrt(1 + 1/K * (b/w - 1))
  # p <- apply(chains,2,FUN=function(X){length(X[X>0])/length(X[X<0])})
	
	## estimates table
	estimatesTab = cbind(
    MaP=as.vector(t(chains[which.max(chains[,1]),])),
    mean=apply(X=chains, FUN=mean, MARGIN=2),
    sd=apply(X=chains, FUN=sd, MARGIN=2),
    t(apply(X=chains, FUN=quantile, MARGIN=2, probs=c(0.05,0.5,0.95))),
    r_hat)
	signif =  gsub("1",x=gsub(pattern="0",x=t(apply(estimatesTab[,c(4,6)]>0,1,diff)),replacement="*"),"ns")
	estimatesTab = data.frame(cbind(round(estimatesTab,4),signif))

    return(list(estimates=estimatesTab,
    covmat=cov(chains[,-1])
  ))
  
}

#
###



#########################
## chainList.bayesPlot ##
#########################

## goal: produce a table summarizing the estimates table obtained from mcmc samples

## args: 
# @chainList - list - of chains
# @logTransform - true/false - whether parameters should be log transformed
chainList.bayesPlot <- function(chainList, logTransform=F, labels=NULL, main="")
{
  ## estimates
  estimatesTab <- chainList.summary(chainList)[["estimates"]]
  estimatesTab <- estimatesTab[-1,]
  s = length(estimatesTab[,"mean"]):1
  # s = rev(order(estimatesTab[,"mean"])) # normal version
  estimatesTab <- estimatesTab[s,]
  
  ## significance
  pvalues = chainList.summaryTab(chainList)[["estimates"]][-1,"signif"][s]
  
  ## chain
  chain = chainList.unlist(chainList)[,-1]
  chain = chain[,s]

  ## plotting parameters
  # par(oma=c(1,4,1,1))

  ## labels
	if(is.null(labels)) labels = rownames(estimatesTab)
  
  ## plot
  lim_x = c(min(c(0,estimatesTab[,"5%"])),max(c(0,estimatesTab[,"95%"])))
  dx = diff(lim_x)
  lim_x = lim_x + 0.25 * c(-1,1) * dx
  lim_y = 0:(nrow(estimatesTab)+1)
  lim_y = c(min(lim_y), max(lim_y))
  plot(0:nrow(estimatesTab), cex=0, bty="l", xlim=lim_x, ylim=lim_y, xlab="Estimated Value", ylab="", xaxt="n", yaxt="n", main=main)
  
  ## background
  coords = par("usr")
  coords_x = coords[1:2]
  coords_y = coords[3:4]
  polygon(x=c(coords_x, rev(coords_x)), y=c(c(coords_y[1],coords_y[1]), c(coords_y[2],coords_y[2])), col=adjustcolor("lightgrey",alpha=0.2), border=NA)
  
  ## x axis
  x = seq(floor(lim_x[1]), ceiling(lim_x[2]), 0.25)
  axis(side = 1, at = x, labels = x, las=1, lwd=0, lwd.ticks=1)
  
  ## y axis
  y = 1:nrow(estimatesTab)
  axis(side = 2, at = 1:nrow(estimatesTab), labels = NA, las=1, lwd=0, lwd.ticks=1)
  text(x=par("usr")[1]-0.1, y=1:nrow(estimatesTab), labels=labels, srt=45, xpd=NA, adj=1.0)
  
  ## grid guides
  for (l in 1:length(y)) lines(c(x[1]-10, x[length(x)]+10), c(y[l], y[l]), col="white")
  for (l in 1:length(x)) lines(c(x[l], x[l]), c(y[1]-10, y[length(y)]+10), col="white")
  lines(c(0,0), c(-1,nrow(estimatesTab)+2), lty=3, lwd=2, col="darkgrey")
  
  ## estimates
  for(i in 1:nrow(estimatesTab))
  {
    
    ## check significance
    if (pvalues[i] == "*") text(estimatesTab[i,"mean"] + 0.05 * dx + 4 * estimatesTab[i,"sd"], i, labels = "*", cex=2, col="red")
    
    ## density
    density_ = density(x=chain[,i])
    density_x = density_$x
    density_y = density_$y / max(density_$y) * 0.33 # / nrow(estimatesTab) * 3
    polygon(x=c(density_x, rev(density_x)), y= c(i, i) + c(density_y, -rev(density_y)), border=NA, col="grey")
    
    ## mean and CI
    lines(estimatesTab[i,c("5%","95%")],c(i,i),col="red", lwd=2);
    points(estimatesTab[i,"mean"],i, pch=16)  
    lines(c(estimatesTab[i,"mean"]-1*estimatesTab[i,"sd"],estimatesTab[i,"mean"]+1*estimatesTab[i,"sd"]),c(i,i),lwd=2, col="black"); 
    
  }
  
  ## plotting parameters
  # par(mar=c(5.1,5.1,4.1,4.1))

}

#
###

#######################
## chainList.postPlot ##
#######################

## goal: visualise posterior distrubtion approximated by chain samples

## args:
# @CMChainList - list - list of CMCchains (i.e. tables)
# @nPoints - int - number of points per graph window

chainList.postPlot <- function(chainList, nPoints, use_labels=T)
{
  
  chains <- NULL ; for(i in 1:length(chainList)){chains <- rbind(chains, chainList[[i]])}
  par(mfrow=c(ncol(chains),ncol(chains)), mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,2,2))
  chains <- chains[order(chains[,1]),]
  
  ## labels
  if(is.null(colnames(chains)) | use_labels == F){colnames(chains) <- paste("X",1:ncol(chains))}
  
  for(i in 1:ncol(chains))
  {
    for(j in 1:ncol(chains))
    {
      
      if(i==ncol(chains)){xaxt <- "s"}else{xaxt <- "n"} ; if(j==1){yaxt <- "s"}else{yaxt <- "n"}
      
      if(i==j)
      {
        plot(1:2, cex=0, xlab="", ylab="", xlim=c(min(chains[,j]),max(chains[,j])), ylim=c(min(chains[,i]),max(chains[,i])), xaxt=xaxt, yaxt=yaxt)
        text(x=(min(chains[,j])+max(chains[,j]))/2, y=(min(chains[,i])+max(chains[,i]))/2 ,labels = paste(colnames(chains)[j]), cex=1)
      }
      
      if(j==i+1)
      {
        hist(chains[,j], freq=F, main="", xaxt="n", yaxt="n") # barplot(height=density(chains[,j])$y, main="", xaxt="n", yaxt="n", horiz = F)
        lines(density(chains[,j]), col="red")
      }
      
      if(i>j)
      {
        plot(1:2, cex=0, xlab="", ylab="", xlim=c(min(chains[,j]),max(chains[,j])), ylim=c(min(chains[,i]),max(chains[,i])), xaxt=xaxt, yaxt=yaxt)            
        subset <- seq(1,nrow(chains),nrow(chains)/nPoints)
        # lines(chains[subset,j], chains[subset,i], col=adjustcolor("black",alpha=0.2), pch=16)
        points(chains[subset,j], chains[subset,i], col=grey(level=0.0+0.9*((chains[subset,1]-min(chains[subset,1]))/(max(chains[subset,1])-min(chains[subset,1]))), alpha=0.75), pch=16)
        # points(x=chains[seq(1,nrow(chains),nrow(chains)/length(chainList)),j], y=chains[seq(1,nrow(chains),nrow(chains)/length(chainList)),i], pch=4, col="blue")
        # points(x=chains[seq(nrow(chains)/length(chainList)-1,nrow(chains),nrow(chains)/length(chainList)),j], y=chains[seq(nrow(chains)/length(chainList)-1,nrow(chains),nrow(chains)/length(chainList)),i], pch=8, col="red")
        points(x=mean(chains[,j]), y=mean(chains[,i]), pch=1, col="red")
        points(x=chains[which.max(chains[,1]),j], y=chains[which.max(chains[,1]),i], pch=8, col="red")
      }
      
      if(j>i+1)
      {
        plot(1:2, cex=0, xlab="", ylab="", xlim=c(min(chains[,j]),max(chains[,j])), ylim=c(min(chains[,i]),max(chains[,i])), xaxt="n", yaxt="n", bty="n")            
      }
    }
  }
  par(mfrow=c(1,1), mar=c(4,4,1,1))
}

#
###

#########################
## chainList.tracePlot ##
#########################

## goal: plot the traces of a chain list

## arg: 
# @chainList - list of matrix of likelihood and parameter values

chainList.tracePlot <- function(chainList)
{
  
  ##
  par_mar_old <- par()$mar
  par_oma_old <- par()$oma
  par(mfrow=c(4,2), mar=c(4,4,1,1), oma=c(0,0,0,0))
  
  ## assemble chains together
  chain <- NULL
  for(j in 1:length(chainList))
  {
    chain <- rbind(chain, cbind(nChain=j, chainList[[j]]),deparse.level=0)
  }
  
  ## visualize chains
  for(i in 2:dim(chain)[2])
  {
    
    ## plots
    plot(chain[,i],ylab=paste(colnames(chain)[i]), xlim=c(0,nrow(chain)/length(chainList)), type="l", col="white")
    # for(j in 1:length(chainList)){lines(1:length(chain[chain[,"nChain"]==j ,i]), chain[chain[,"nChain"]==j,i], col=grey(1-(0.2 + 0.6*j/length(chainList))))}
    for(j in 1:length(chainList)){lines(1:length(chain[chain[,"nChain"]==j ,i]), chain[chain[,"nChain"]==j,i], col=rainbow(length(chainList),start=0.75,end=0.5,alpha=0.5)[j])}
    ## hist
    hist(chain[,i],freq=F,main="", xlab="")
    # lines(density(chain[,i]), col="white")
    # for(j in 1:length(chainList)){lines(density(chain[chain[,"nChain"]==j,i]), col=grey(1-(0.2 + 0.6*j/length(chainList))))}
    for(j in 1:length(chainList)){lines(density(chain[chain[,"nChain"]==j,i]), col=rainbow(length(chainList),start=0.75,end=0.5,alpha=0.5)[j])}
    
  }
  par(mfrow=c(1,1),mar=par_mar_old,oma=par_oma_old)
}

#
###

#################
## acplot.will ##
#################

## goal: plot change in autocorrelation with increasing order

## args:
# @arg1 - X - the data stream
# @arg2 - order - the order up to which ac are calculated

acplot.will <- function(X,order)
{
  
  orderVect <- corVect <- NULL
  for(i in 1:order)
  {
    corVect_local <- NULL; X_local <- X
    for(j in 1:10)
    {
      corVect_local <- c(corVect_local,cor(X_local[-c((length(X_local)-i+1):length(X_local))],X_local[-c(1:i)]))
      X_local <- X_local[-c(1:(length(X)/10))]
    }
    corVect <- c(corVect, mean(corVect_local))
  }
  
  ##
  barplot(height = corVect, ylim=c(-1,1))
}

#
###

######################
## chainAcPlot.will ##
######################

## goal: plots acplots for each elements of the chainList

## args:
# @chainList - list of chains

chainList.acPlot <- function(chainList)
{
  for(l in 1:length(chainList))
  {
    for(i in 2:dim(chainList[[l]])[2])
    {
      acplot.will(X=chainList[[l]][,i],order=50)
    }  
  }
}

#
###


###################
## chainList.ESS ##
###################

## goal: compute the effective sample size

chainList.ESS <- function(chainList)
{
  ESSList <- list()
  l <- 1
  for(l in 1:length(chainList))
  {
    ESSVect <- NULL
    i <- 2
    for(i in 2:dim(chainList[[l]])[2])
    {
      X <- chainList[[l]][,i]
      orderVect <- corVect <- NULL
      o <- 1
      for(o in 1:50)
      {
        corVect_local <- NULL; X_local <- X
        for(j in 1:10)
        {
          corVect_local <- c(corVect_local,cor(X_local[-c((length(X_local)-o+1):length(X_local))],X_local[-c(1:o)]))
          X_local <- X_local[-c(1:(length(X)/10))]
        }
        corVect <- c(corVect, mean(corVect_local))
      }
      ESSVect <- c(ESSVect,dim(chainList[[l]])[1]/(1+2*sum(abs(corVect))))
      
    }
    ESSList[[l]] <- ESSVect
  }
  return(ESSList)
}

#
###

#################################################
## UNIDIMENSIONAL MARGINAL POSTERIOR DENSITIES ##
#################################################

## goal: calculate marginal posterior density over each parameter

## args:
# @postMap - the posterior map (a chain type object)
# @delta - the resolution of the marginalization

marPost <- function(postMap, delta=0.1)
{
  marPostMap <- NULL
  
  ## for each variable
  i <- 2
  for(i in 2:ncol(postMap))
  {
    
    ## initiate
    postMeanVect <- NULL
    meanVect <- NULL
    
    ## define bins for integration
    binVect <- quantile(x=postMap[,i], probs=seq(0, 1, delta))
    
    ## average the value of variable in each bins
    j <- 1
    for(j in 1:(length(binVect)-1))
    {
      
      postMeanVect <- c(postMeanVect, median(postMap[postMap[,i] >= binVect[j] & postMap[,i] <= binVect[j+1] ,1], na.rm=T))
      meanVect <- c(meanVect, median(postMap[postMap[,i] >= binVect[j] & postMap[,i] <= binVect[j+1] ,i], na.rm=T))
      
    }
    
    marPostMap <- cbind(marPostMap, cbind(postMeanVect, meanVect))
    
  }
  colnames(marPostMap) <- as.vector(matrix(c(paste(rep("logP", ncol(postMap)-1),colnames(postMap)[-1], sep="_") , colnames(postMap)[-1]), byrow=T, nrow=2))
  
  return(marPostMap)
}

##############################################
## BIDIMENSIONAL MARGINAL POSTERIOR MAPPING ##
##############################################

## goal: visualize bidimensional marginal posterior densities considering parameters pairwise

## args:
# @postMap - a chain type object
# @d - the resolution of the marginalization

## NOTE: it requires alot of samples => better try it on the final distribution

marPost.bidim <- function(postMap, delta=0.1)
{
  
  marPostMap <- NULL
  
  ## for each variable
  i <- 2
  j <- 3
  k <- l <- 1
  for(i in 2:(ncol(postMap)-1))
  {
    for(j in (i+1):ncol(postMap))
    {
      ## initiate
      yMeanVect <- x_1MeanVect <- x_2MeanVect <-  NULL
      
      ## define bins for integration
      x_1BinVect <- quantile(x=postMap[,i], probs=seq(0, 1, delta))
      x_2BinVect <- quantile(x=postMap[,j], probs=seq(0, 1, delta))
      
      ## average the value of variable in each bin combination
      for(k in 1:(length(x_1BinVect)-1))
      {
        for(l in 1:(length(x_2BinVect)-1))
        {
          subset <- postMap[,i] >= binVect[k] & postMap[,i] <= binVect[k+1] & postMap[,j] >= binVect[l] & postMap[,j] <= binVect[l+1]
          yMeanVect <- c(yMeanVect, median(postMap[subset ,1], na.rm=T))
          x_1MeanVect <- c(x_1meanVect, median(postMap[subset ,i], na.rm=T))
          x_2MeanVect <- c(x_2meanVect, median(postMap[subset ,j], na.rm=T))
        }
      }
      marPostMap <- cbind(marPostMap, cbind(yMeanVect, x_1MeanVect, x_2MeanVect))
    }
  }
  
  colnames(marPostMap) <- as.vector(matrix(c(paste(rep("logP",ncol(postMap)-1), colnames(postMap)[-1], sep="_") , colnames(postMap)[-1], colnames(postMap)[-1]), byrow=T, nrow=3))
  
  return(marPostMap)
}

#
###

