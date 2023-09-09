###########
## AMC.R ##
###########

## update: 

## goal: R implementation of an adaptive MCMCMH

## update log:
## 18-06-2022 - created v0_1
##            - removed naming of argument in target function
##            - reduced the frequency of shape adaptation

##############
## INITIATE ##
##############

##
library(mvtnorm)

#
###

#########
## AMC ##
#########

## goal: MCMCMH with adapting shape

## args:
# @dTarget - function - target function to explore
# @Theta_0 - vector - vector of initial parameters
# @scale_p_0 - float - scaling parameter
# @Sigma_p_0 - matrix - covariance matrix for the target distribution
# @itVect - vector - vector containing the iterations for the adaptsize adapshape and sampling phase
# @msg - bool - whether messages should be printed at each iterations

AMC <- function(dTarget, Theta_0, scale_p_0, Sigma_p_0, itVect, adaptShape=T, msg=T)
{
  
  ####################
  ## INITIATE CHAIN ##
  
  ## initiate param vector
  Theta_local <- Theta_0
  d <- length(Theta_local)
  
  ##
  Sigma_p <- Sigma_p_0
  
  ## evalutate the local value of target at theta_0
  # target_local <- dTarget(Theta=Theta_local, Mu=Mu, Sigma=Sigma)
  target_local <- dTarget(Theta_local)
  
  ## initiate chain
  Theta_chain <- c(dLogPost=target_local, Theta_local); names(Theta_chain) <- c("dLogPost",names(Theta_local))
  
  ################
  ## ADAPT SIZE ##
  
  i = 1 
  k = 1
  count  = 1
  accepted = 0
  idx_accepted = NULL
  for(phase in 1:2)
  {
    for(i in 1:itVect[phase])
    {
      
      Theta_newLocal <- NULL
      
      ## draw new size params
      # scale_p <- scale_p_0 * if (phase == 1) runif(n=1, min=0, max=2) else 1
      scale_p <- scale_p_0 * runif(n=1, min=0, max=2.0)
      
      ## draw mean and dispersion parameters
      Theta_newLocal <- as.vector(rmvnorm(n=1, mean=Theta_local, sigma=scale_p^2*Sigma_p)); names(Theta_newLocal) <- names(Theta_local)
      
      ## calculate new target density
      target_newLocal <- dTarget(Theta_newLocal)
      
      ## calculate metropolis-hasting ratio
      r <- exp(target_newLocal - target_local) 
      if(is.nan(r)){r <- 0}
      
      ## test to accept or reject MH ratio
      if(runif(1, min=0, max=1) < r) 
      {
        Theta_local  = Theta_newLocal
        target_local = target_newLocal
        accepted     = accepted + 1
        idx_accepted = c(idx_accepted,k)
      }
      
      ## update chain
      Theta_chain <- rbind(Theta_chain, c(target_local, Theta_local), deparse.level=0)
      
      ## print current state of chain and acceptance rate
      if(msg==TRUE & count == 1000)
      {
        message(paste(k,
                      " | p: ", round(accepted/k, 2),
                      " | dP: ", round(target_local, 2),
                      " | Theta: ", 
                      round(Theta_local[1], 4),
                      round(Theta_local[2], 4)
        ))

        ## update the covariance matrix
        if(adaptShape==TRUE & phase==2)
        {
          Sigma_p = cov(Theta_chain[(round(k/2)):k,-1])
          # Sigma_p <- cov(Theta_chain[idx_accepted,-1])
          # phase = 3
        }
        count = 0
      }
      count = count + 1
      k     = k + 1
    }
    ## print acceptance rate
    message(paste("p_",phase," = ",round(accepted/k, 2),sep=""))
    
  }
  ## TERMINATE ##
  return(Theta_chain)
}

#
###
