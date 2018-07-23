##################################################################################################################
##################################################################################################################
# Gamma-distributed infectious period MCMC
##################################################################################################################
##################################################################################################################

#####Functions#####

#########################################################
# count the number of infectives just before time t.
#########################################################
count.no.inf = function(infection.times, removal.times, t) {
  # individual -i- is infected just before time t, if this holds:  I_i < t < R_i
  sum(infection.times < t) - sum(removal.times < t)
}


###############################################################
# compute the integral(S_t*I_t)dt using a double sum
###############################################################
compute.total.pressure = function(data, fs) {
  N = nrow(data);  
  total = 0;
  for (i in 1:fs) { #iterate through each infected individual
    for (j in 1:N) { #iterate through each individual (potential susceptibles)
      #removal[i] must be > infection[i]
      #if infection[j] > removal[i], then infection[j] also > infection[i] and exposure period for j by i is removal[i] - infection[i]
      #if removal[i] > infection[j], then total exposure period for j by i is infection[i] - infection[j] IF i infected first; if j infected first, then no exposure period for j by i
      total = total + min(data$removal[i], data$infection[j]) - min(data$infection[i], data$infection[j])
    }
  }
  total
}

######################################################
# compute the log of the product: prod(# infected at each time, number removed at each time)
######################################################
compute.log.prod = function(data, infection.times, removal.times, fs) {
  
  # find the index of the initial infective.
  time.init.inf = min(data$infection[1:fs]) #time of first infection
  index.init.inf = which(data$infection[1:fs]==time.init.inf) #index of first infection
  
  # create a vector to store the elements of the product
  log.prod.vec = rep(NA, fs);
  
  for (i in 1:fs) { #iterate over all infected individuals (who each have i and r time)
    if (i!=index.init.inf) { #if i does not equal the index of the first individual, then
      #count the number of individuals infected at time point infection.times[i]
      #take log because allows addition of small values rather than multiplication of small values
      log.prod.vec[i] = log(count.no.inf(infection.times, removal.times, infection.times[i]))
    }
    else {
      log.prod.vec[i] = 0
    }
  }
  sum(log.prod.vec)
}


#####MCMC routine#####

mcmcSIR.gamma = function(data, iter, alpha) {
  
  #Data frame "data" should contain two columns; first is the label of the individuals and the second contains the removal times.
  
  # Get the final size and the size of population
  fs = sum(data$removal != Inf)
  N = nrow(data)
  
  # Get a vector for the removal times
  removal.times = data$removal[1:fs]
  
  # First of all, we need to simulate the initial values for the infection
  # times to iniate our MCMC sampler.
  
  # an easy way to place all the infection times just before the first removal time.
  infection.times = seq(0, min(removal.times), len=fs)
  
  # We "augment" the data frame "data" with the infection times.
  data$infection = c(infection.times, rep(Inf, N-fs))
  data = data[,c(1,3,2)] # re order
  
  # Assume  for the parameters beta and gamma
  lambda.beta = 1.0
  nu.beta = 10^(-3)
  
  lambda.gamma = 1.0  
  nu.gamma = 10^(-3)
  
  # compute the current values of the double sum and the log(product).
  double.sum.cur = compute.total.pressure(data, fs);
  log.prod.cur = compute.log.prod(data,infection.times, removal.times, fs);
  log.prod.cur2 = sum(log((removal.times - infection.times)^(alpha-1)))
  sum.R.minus.I.cur = sum(removal.times-infection.times)
  
  # initial values for beta and gamma - draw from the full conditionals.
  beta.cur = rgamma(1, fs - 1 + lambda.beta, 1.0)/(nu.beta + double.sum.cur/N)
  gamma.cur = rgamma(1, fs*alpha + lambda.gamma, 1.0)/(nu.gamma + sum.R.minus.I.cur);
  
  # create a matrix to store the values for beta, gamma and I.
  res = matrix(NA, nrow=iter, ncol=3)
  res[1,] = c(beta.cur, gamma.cur, sum(infection.times))
  
  ##########################
  # MCMC loop starts here  # 
  ##########################
  
  for (i in 2:iter) {
    
    ########################
    # update infection times
    ########################
    
    # first choose an individual whose infection time will be updated.
    choose.ind = sample(1:fs, 1)
    
    # record the current value and propose a candidate value
    inf.time.cur = infection.times[choose.ind]
    inf.time.can = removal.times[choose.ind] - rgamma(1, alpha, 1.0)/gamma.cur;
    
    # compute the (new) value of the double sum, log(prod) and sum.R.minus.I
    
    # first update the data frame data which is used for the computation of the quantities of interest
    data$infection[choose.ind] = inf.time.can;
    infection.times[choose.ind] = inf.time.can;
    
    # compute the log product first. The proposed move should be consistent with the observed data,
    # i.e. results into the same final size. Otherwise, the move is rejected straight away.
    
    log.prod.can = compute.log.prod(data, infection.times, removal.times, fs);
    log.prod.can2 = sum(log((removal.times - infection.times)^(alpha-1)))
    
    if ((log.prod.can != -Inf)&(log.prod.can2 != -Inf)) { 
      
      # compute the other two quantities
      double.sum.can = compute.total.pressure(data, fs);
      sum.R.minus.I.can = sum(removal.times-infection.times);
      
      # compute the log(q_ratio);
      q.num = (gamma.cur^alpha)*((removal.times[choose.ind] - inf.time.cur)^(alpha-1))*exp(-(removal.times[choose.ind] - inf.time.cur)*gamma.cur)/gamma(alpha)
      q.denom = (gamma.cur^alpha)*((removal.times[choose.ind] - inf.time.can)^(alpha-1))*exp(-(removal.times[choose.ind] - inf.time.cur)*gamma.cur)/gamma(alpha)
      
      log.q.ratio = log(q.num) - log(q.denom)
      
      log.pi.can = log.prod.can + log.prod.can2 - (beta.cur/N)*double.sum.can - gamma.cur*sum.R.minus.I.can;
      log.pi.cur = log.prod.cur + log.prod.cur2 - (beta.cur/N)*double.sum.cur - gamma.cur*sum.R.minus.I.cur;
      
      
      # accept/reject move
      u = runif(1);
      
      if (log(u) < log.pi.can - log.pi.cur + log.q.ratio) {
        
        log.prod.cur = log.prod.can;
        log.prod.cur2 = log.prod.can2
        double.sum.cur = double.sum.can;
        sum.R.minus.I.cur = sum.R.minus.I.can;
      }
      else {
        data$infection[choose.ind] = inf.time.cur;
        infection.times[choose.ind] = inf.time.cur;
      }      
    }
    else {
      data$infection[choose.ind] = inf.time.cur;
      infection.times[choose.ind] = inf.time.cur;
    }
    
    #################    
    # update beta
    #################
    beta.cur = rgamma(1, fs - 1 + lambda.beta, 1.0)/(nu.beta + double.sum.cur/N)
    
    #################    
    # update gamma
    #################    
    gamma.cur = rgamma(1, fs*alpha + lambda.gamma, 1.0)/(nu.gamma + sum.R.minus.I.cur)
    
    
    # store the current values of beta, gamma and the sum of the infection times.
    res[i,] = c(beta.cur, gamma.cur, sum(infection.times))
    
  }
  res = as.data.frame(res)
  names(res) = c('beta','gamma','sum.inf')
  res
  
}
par(mfrow = c(2, 2))
data = read.table('data.txt',header=TRUE)
out2 = mcmcSIR.gamma(data,10000,2)
mean(out$beta[1000:10000])
var(out$beta[1000:10000])
plot(out$beta,type='l',xlab="Iteration",ylab="Beta",ylim=c(0,6))
hist(out$beta[1000:10000],breaks=24,xlab="Beta",main="",col='lightblue')

mean(out$gamma[1000:10000])
var(out$gamma[1000:10000])

par(mfrow = c(3, 2))
plot(out$beta,type='l',xlab="Iteration",ylab="Beta",ylim=c(0,6))
hist(out$beta[1000:10000],breaks=24,xlab="Beta",main="",col='lightblue')
plot(out$gamma,type='l',xlab="Iteration",ylab="gamma",ylim=c(0,6))
hist(out$gamma[1000:10000],breaks=24,xlab="gamma",main="",col="salmon")
plot(out$sum.inf,type='l',xlab="Iteration",ylab="Sum. Inf.",ylim=c(45,55))
hist(out$sum.inf[1000:10000],breaks=24,xlab="Sum. Inf.",main="",col="lightgreen")

plot(out2$beta,type='l',xlab="Iteration",ylab="Beta",ylim=c(0,10))
hist(out2$beta[1000:10000],breaks=48,xlab="Beta",main="",col='lightblue')
plot(out2$gamma,type='l',xlab="Iteration",ylab="gamma",ylim=c(0,10))
hist(out2$gamma[1000:10000],breaks=48,xlab="gamma",main="",col="salmon")
plot(out2$sum.inf,type='l',xlab="Iteration",ylab="Sum. Inf.",ylim=c(45,55))
hist(out2$sum.inf[1000:10000],breaks=48,xlab="Sum. Inf.",main="",col="lightgreen")

out.sample = sample(nrow(out[1000:10000,]),2000)
plot(out$beta[out.sample],out$gamma[out.sample],type="p",xlab="Beta",ylab="gamma")
plot(out$gamma[out.sample],out$sum.inf[out.sample],type="p",xlab="gamma",ylab="Sum. Inf.")

hist(out$beta*1/out$gamma,breaks=24,main="",xlab="R0",col="lightyellow")
hist(out2$beta*2/out2$gamma,breaks=24,main="",xlab="R0",col="lightyellow")
