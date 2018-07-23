################################################################################
# R-code to simulate from a Markov and Non-Markov Stochastic epidemic models
################################################################################


########Sampling#######
N = 21
itermax = 100
R0.list = matrix(c(0.9,2,4, 1,1,1),nrow=2, byrow=TRUE)

Num.Dist.Fun = function(N,itermax,R0.list,fun,delta,out.type){
  group.out = as.data.frame(matrix(nrow=dim(R0.list)[2],ncol=itermax))
  
  for (v in 1:dim(R0.list)[2]){
    beta = R0.list[1,v]
    gamma = R0.list[2,v]
    
    out = c()
    for (i in 1:itermax){
      if(out.type=="Num.Inf"){
      if(delta == 0){out.val = fun(N,beta,gamma)$Num.Inf} else {out.val = fun(N,beta,gamma,delta)$Num.Inf}
      } else if(out.type=="Time.Epi"){
        if(delta == 0){out.val = fun(N,beta,gamma)$Time.Epi} else {out.val = fun(N,beta,gamma,delta)$Time.Epi}
      }
      out = append(out, out.val)
    }
    
    group.out[v,] = out
    
  }
  return(group.out)
}

#Basic Markovian SIR
simSIR.Markov.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Markov,0,"Num.Inf")

hist(as.numeric(simSIR.Markov.df[1,]))
hist(as.numeric(simSIR.Markov.df[2,]))
hist(as.numeric(simSIR.Markov.df[3,]))


#Alternative Markovian SIR
simSIR.MarkovAlt.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Markov.alternative,0,"Num.Inf")

hist(as.numeric(simSIR.MarkovAlt.df[1,]))
hist(as.numeric(simSIR.MarkovAlt.df[2,]))
hist(as.numeric(simSIR.MarkovAlt.df[3,]))


#Non-Markovian constant infectious period
simSIR.nonMarkovConst.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Non.Markov.constant,0,"Num.Inf")

hist(as.numeric(simSIR.nonMarkovConst.df[1,]))
hist(as.numeric(simSIR.nonMarkovConst.df[2,]))
hist(as.numeric(simSIR.nonMarkovConst.df[3,]))

#Non-Markovian Gamma infectious period
simSIR.nonMarkovGamma.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Non.Markov.gamma,1,"Num.Inf")

hist(as.numeric(simSIR.nonMarkovGamma.df[1,]))
hist(as.numeric(simSIR.nonMarkovGamma.df[2,]))
hist(as.numeric(simSIR.nonMarkovGamma.df[3,]))

#Non-Markovian Weibull infectious period
simSIR.nonMarkovWeibull.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Non.Markov.Weibull,0,"Num.Inf")

hist(as.numeric(simSIR.nonMarkovWeibull.df[1,]))
hist(as.numeric(simSIR.nonMarkovWeibull.df[2,]))
hist(as.numeric(simSIR.nonMarkovWeibull.df[3,]))

#Basic Markovian SIR
simSIR.Markov.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Markov,0,"Time.Epi")

hist(as.numeric(simSIR.Markov.df[1,]))
hist(as.numeric(simSIR.Markov.df[2,]))
hist(as.numeric(simSIR.Markov.df[3,]))


#Alternative Markovian SIR
simSIR.MarkovAlt.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Markov.alternative,0,"Time.Epi")

hist(as.numeric(simSIR.MarkovAlt.df[1,]))
hist(as.numeric(simSIR.MarkovAlt.df[2,]))
hist(as.numeric(simSIR.MarkovAlt.df[3,]))


#Non-Markovian constant infectious period
simSIR.nonMarkovConst.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Non.Markov.constant,0,"Time.Epi")

hist(as.numeric(simSIR.nonMarkovConst.df[1,]))
hist(as.numeric(simSIR.nonMarkovConst.df[2,]))
hist(as.numeric(simSIR.nonMarkovConst.df[3,]))

#Non-Markovian Gamma infectious period
simSIR.nonMarkovGamma.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Non.Markov.gamma,1,"Time.Epi")

hist(as.numeric(simSIR.nonMarkovGamma.df[1,]))
hist(as.numeric(simSIR.nonMarkovGamma.df[2,]))
hist(as.numeric(simSIR.nonMarkovGamma.df[3,]))

#Non-Markovian Weibull infectious period
simSIR.nonMarkovWeibull.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Non.Markov.Weibull,0,"Time.Epi")

hist(as.numeric(simSIR.nonMarkovWeibull.df[1,]))
hist(as.numeric(simSIR.nonMarkovWeibull.df[2,]))
hist(as.numeric(simSIR.nonMarkovWeibull.df[3,]))

#Non-Markovian constant latent
simSIR.nonMarkovConstantE.df = Num.Dist.Fun(N,itermax,R0.list,simSIR.Non.Markov.constantE,0,"Time.Epi")

hist(as.numeric(simSIR.nonMarkovConstantE.df[1,]))
hist(as.numeric(simSIR.nonMarkovConstantE.df[2,]))
hist(as.numeric(simSIR.nonMarkovConstantE.df[3,]))

simSIR.Markov <- function(N, beta, gamma) {
  
  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;
  
  # recording time;
  t <- 0;
  times <- c(t);
  
  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  while (I > 0) {
    
    # time to next event;
    t <- t + rexp(1, (beta/N)*I*S + gamma*I);
    times <- append(times, t);
    
    if (runif(1) < beta*S/(beta*S + N*gamma)) {
      # infection
      I <- I+1;
      S <- S-1;
      type <- append(type, 1);
    }
    else {
      #removal
      I <- I-1
      type <- append(type, 2);
    }
  }
  
  # record the final size , i.e. the number of initially susceptibles who contracted the disease sometime during the epidemic.
  Num.Inf = length(which(type == 1)) #final number infected
  Time.Epi = length(times)
  
  
  # record the times of events (infections/removals) as well as the type
  
  res <- list("t"=times, "type"=type,"Num.Inf"=Num.Inf,"Time.Epi"=Time.Epi);
  res
}

simSIR.Markov.alternative <- function(N, beta, gamma) {
  
  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;
  
  # recording time;
  t <- 0;
  times <- c(t);
  
  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  while (I > 0) {
    
    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0) {
      t.next.infection <- t +  rexp(1, (beta/N)*I*S)
    }
    else {
      t.next.infection <- Inf;
    }
    
    # time to next removal    
    t.next.removal <- t + rexp(1, gamma*I)
    
    
    # check which of the two events happens first
    if (t.next.infection < t.next.removal) {
      # infection occurs
      I <- I+1;
      S <- S-1;
      type <- append(type, 1);
      times <- append(times, t.next.infection);
      t <- t.next.infection
    }
    else {
      #removal occurs
      I <- I-1
      times <- append(times, t.next.removal);
      type <- append(type, 2);
      t <- t.next.removal
    }
  }
  
  # record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.
  # record the final size , i.e. the number of initially susceptibles who contracted the disease sometime during the epidemic.
  Num.Inf = length(which(type == 1)) #final number infected
  Time.Epi = length(times)
  
  
  # record the times of events (infections/removals) as well as the type
  
  res <- list("t"=times, "type"=type,"Num.Inf"=Num.Inf,"Time.Epi"=Time.Epi);
  res

}

simSIR.Non.Markov.constant <- function(N, beta, k) {
  
  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;
  
  # recording time;
  t <- 0;
  times <- c(t);
  
  # create a vector containing the removal times of all the current infectives.
  r <- k 
  
  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  # a counter for labelling the individuals
  lambda <- 1;
  
  # a vector to store the labels
  labels <- c(1);
  
  while (I > 0) {
    
    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0) {
      T  <- rexp(1, (beta/N)*I*S)
    }
    else {
      T <- Inf;
    }
    
    # time to next removal
    R <- min(r, na.rm=TRUE);
    
    # check which of the two events happens first
    if (t + T < R) {
      # infection occurs
      I <- I+1;
      S <- S-1;
      r <- append(r, t + T + k)
      type <- append(type, 1);
      times <- append(times, t + T);
      
      lambda <- lambda + 1;
      labels <- append(labels, lambda)
      t <- t + T
    }
    else {
      #removal occurs
      I <- I-1
      type <- append(type, 2);
      index.min.r <- which(min(r, na.rm=TRUE)==r)
      r[index.min.r] <- NA
      labels <- append(labels, index.min.r)
      times <- append(times, R);
      t <- R
      
      # update the vector of 
    }
  }
  
  # record the final size , i.e. the number of initially susceptibles who contracted the disease sometime during the epidemic.
  Num.Inf = length(which(type == 1)) #final number infected
  Time.Epi = length(times)

  res <- list("t"=times, "type"=type, "labels"=labels,"Num.Inf"=Num.Inf,"Time.Epi"=Time.Epi);
  res
}

simSIR.Non.Markov.gamma <- function(N, beta, gamma, delta) {
  
  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;
  
  # recording time;
  t <- 0;
  times <- c(t);
  
  # create a vector containing the removal times of all the current infectives.
  k <- rgamma(1, gamma, delta)
  r <- k 
  
  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  # a counter for labelling the individuals
  lambda <- 1;
  
  # a vector to store the labels
  labels <- c(1);
  
  while (I > 0) {
    
    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0) {
      T  <- rexp(1, (beta/N)*I*S)
    }
    else {
      T <- Inf;
    }
    
    # time to next removal
    R <- min(r, na.rm=TRUE);
    
    # check which of the two events happens first
    if (t + T < R) {
      # infection occurs
      I <- I+1;
      S <- S-1;
      k <- rgamma(1, gamma, delta)
      r <- append(r, t + T + k)
      
      lambda <- lambda + 1;
      labels <- append(labels, lambda)
      type <- append(type, 1);
      times <- append(times, t + T);
      t <- t + T
    }
    else {
      #removal occurs
      I <- I-1
      type <- append(type, 2);
      index.min.r <- which(min(r, na.rm=TRUE)==r)
      r[index.min.r] <- NA
      labels <- append(labels, index.min.r)
      times <- append(times, R);
      t <- R      
    }
  }
  
  # record the final size , i.e. the number of initially susceptibles who contracted the disease sometime during the epidemic.
  Num.Inf = length(which(type == 1)) #final number infected
  Time.Epi = length(times)
  
  
  # record the times of events (infections/removals) as well as the type
  
  res <- list("t"=times, "type"=type,"Num.Inf"=Num.Inf,"Time.Epi"=Time.Epi);
  res
}

simSIR.Non.Markov.Weibull <- function(N, beta, gamma) {
  
  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;
  
  # recording time;
  t <- 0;
  times <- c(t);
  
  # create a vector containing the removal times of all the current infectives.
  k <- rweibull(1,gamma,1)
  r <- k 
  
  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  # a counter for labelling the individuals
  lambda <- 1;
  
  # a vector to store the labels
  labels <- c(1);
  
  while (I > 0) {
    
    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0) {
      T  <- rexp(1, (beta/N)*I*S)
    }
    else {
      T <- Inf;
    }
    
    # time to next removal
    R <- min(r, na.rm=TRUE);
    
    # check which of the two events happens first
    if (t + T < R) {
      # infection occurs
      I <- I+1;
      S <- S-1;
      k <- rweibull(1,gamma,1)
      r <- append(r, t + T + k)
      
      lambda <- lambda + 1;
      labels <- append(labels, lambda)
      type <- append(type, 1);
      times <- append(times, t + T);
      t <- t + T
    }
    else {
      #removal occurs
      I <- I-1
      type <- append(type, 2);
      index.min.r <- which(min(r, na.rm=TRUE)==r)
      r[index.min.r] <- NA
      labels <- append(labels, index.min.r)
      times <- append(times, R);
      t <- R      
    }
  }
  
  # record the final size , i.e. the number of initially susceptibles who contracted the disease sometime during the epidemic.
  Num.Inf = length(which(type == 1)) #final number infected
  Time.Epi = length(times)
  
  
  # record the times of events (infections/removals) as well as the type
  
  res <- list("t"=times, "type"=type,"Num.Inf"=Num.Inf,"Time.Epi"=Time.Epi);
  res
}

#Exercise 7

test = simSIR.Non.Markov.constantE(21,4,3,1)
simSIR.Non.Markov.constantE <- function(N, beta, kI,kE) {
  
  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;
  E <- 0;
  
  # recording time;
  t <- 0;
  times <- c(t);
  
  # create a vector containing the removal times of all the current infectives.
  rI <- kI 
  
  # create a vector containing the removal times of all the current latents.
  rE <- c(100); #set to originally be large because no latents 
  
  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  # a counter for labelling the individuals
  lambda <- 1;
  
  # a vector to store the labels
  labels <- c(1);
  
  while (I > 0) {
    
    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0) {
      T  <- rexp(1, (beta/N)*I*S)
    }
    else {
      T <- Inf;
    }
    
    # time to next removal
    RE <- min(rE, na.rm=TRUE);
    RI <- min(rI, na.rm=TRUE);
    R = min(RE,RI)
    # check which of the two events happens first
    if (R == RE){
      if (t + T < RE) {
        # infection occurs
        E <- E+1;
        S <- S-1;
        rE <- append(rE, t + T + kE)
        type <- append(type, 1);
        times <- append(times, t + T);
        
        lambda <- lambda + 1;
        labels <- append(labels, lambda)
        t <- t + T
      } else {
        #removal occurs
        E <- E-1
        I <- I+1
        type <- append(type, 2);
        index.min.r <- which(min(rE, na.rm=TRUE)==RE)
        rE[index.min.r] <- NA
        labels <- append(labels, index.min.r)
        times <- append(times, RE);
        t <- RE
        
        # update the vector of 
      }
    } else if (R == RI){
      if (t + T < RI) {
        # infection occurs
        E <- E+1;
        S <- S-1;
        rE <- append(rE, t + T + kI)
        type <- append(type, 1);
        times <- append(times, t + T);
        
        lambda <- lambda + 1;
        labels <- append(labels, lambda)
        t <- t + T
      } else {
        #removal occurs
        I <- I-1
        type <- append(type, 3);
        index.min.r <- which(min(rI, na.rm=TRUE)==RI)
        rI[index.min.r] <- NA
        labels <- append(labels, index.min.r)
        times <- append(times, RI);
        t <- RI
        
        # update the vector of 
      }
    }
    
  }
  
  # record the final size , i.e. the number of initially susceptibles who contracted the disease sometime during the epidemic.
  Num.Inf = length(which(type == 1)) #final number infected
  Time.Epi = length(times)
  
  res <- list("t"=times, "type"=type, "labels"=labels,"Num.Inf"=Num.Inf,"Time.Epi"=Time.Epi);
  res
}

