###################################
#Data
###################################

cases = c(0,1,2,3)
households = c(29,9,2,2)
mat = matrix(NA,nrow=4,2)
mat[,1] = cases; mat[,2] = households
mat = as.data.frame(mat)
names(mat) = c('No.of.Cases','No.of.Households')

mcmcSIR.hh <- function(data, iter, sigma) {
  
  # values for the prior hyper-parameters
  lambda.p <- 0.5
  nu.p <- 0.5
  
  lambda.alpha <- 1.0  
  nu.alpha <- 10^(-3)
  
  
  # initial values for the Markov Chain
  alpha.cur <- 0.8;
  p.cur <- 0.5;
  
  # create a matrix to store the values for beta, gamma and I.
  res <- matrix(NA, nrow=iter, ncol=2)
  res[1,] <- c(p.cur, alpha.cur)
  
  ##########################
  # MCMC loop starts here  # 
  ##########################
  
  for (i in 2:iter) {
    
    ########################
    # update p
    ########################
    # propose new value
    p.can = rnorm(1, p.cur, sigma)
    
      if((p.can > 0)&(p.can < 1)){
      # compute the log-densities
      log.pi.cur = log.like(alpha.cur,p.cur) + log(dbeta(p.cur,lambda.p,nu.p))
      log.pi.can = log.like(alpha.cur,p.can) + log(dbeta(p.can,lambda.p,nu.p))
        
      # log-qratio is zero
      log.q.ratio = 0
      
      a.ij = min(log.pi.can - log.pi.cur,log(1))
      u = runif(1)
      if (log(u) < a.ij) {
        p.cur = p.can
      }
    
    }
    ########################
    # update alpha
    ########################
    
    # propose new value
    alpha.can = rnorm(1, alpha.cur, sigma)
    
    # alpha is strictly positive, therefore any proposed value which is negative is automatically rejected.
    if (alpha.can > 0.0) {
      
      # otherwise, if it is positive the accept it with some probability:
      
      # compute the log-densities
      log.pi.cur = log.like(alpha.cur,p.cur) + log(dgamma(alpha.cur,lambda.alpha,1.0)/nu.alpha)
      log.pi.can = log.like(alpha.can,p.cur) + log(dgamma(alpha.can,lambda.alpha,1.0)/nu.alpha)
        
      # log-qratio is zero since we are doing a random Walk Metropolis
      log.q.ratio <- 0;
      
      a.ij = min(log.pi.can - log.pi.cur,log(1))
      u = runif(1)
      if (log(u) < a.ij) {
        alpha.cur = alpha.can
      }
    }
    
    # store the current values of beta, gamma and the sum of the infection times.
    res[i,] <- c(p.cur, alpha.cur)
    
  }
  res = as.data.frame(res)
  names(res) = c('p','alpha')
  res
  
}

out = mcmcSIR.hh(data,10000,0.1)
plot(out$p,type="l")
plot(out$alpha,type="l")
par(mfrow = c(2, 2))
plot(out$p,type='l',xlab="Iteration",ylab="p",ylim=c(0,1))
hist(out$p[1000:10000],breaks=24,xlab="p",main="",col='lightblue',xlim=c(0,1))
plot(out$alpha,type='l',xlab="Iteration",ylab="alpha",ylim=c(0,1))
hist(out$alpha[1000:10000],breaks=24,xlab="alpha",main="",col="salmon",xlim=c(0,1))

plot(out$p[1000:10000],out$alpha[1000:10000],xlab="p",ylab="alpha")
