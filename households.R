# this is a function to compute the final size
# m initial susceptibles
# a initial infectives
# alpha infection rate
# the length of the fixed infectious period

#P(T = k | Y = y)
compute.conditional.prob <- function(m, j, a, alpha, c){
  if (j==0) {
    res <- (exp(-alpha*m*c))^a
  }
  else {
    
    #output is sum of (m-k  j-k)*p(k)/f(alpha(m-j))^(k+a)
    #f(alpha(m-j)) = exp(-alpha*(m-j)*c)^
    
    part.one <- exp(-alpha*(m-j)*c)^(j+a) 
    part.two <- choose(m, j) #(m j)=factorial(m)/(factorial(j)*factorial(m-j))
    
    sum <- 0;
    for (k in 0:(j-1)) {
      #will calculate p(k) by solving p(x) for x = 0:k and storing the output each time (sum) of p values
      sum <- sum + choose(m-k, j-k)*compute.conditional.prob(m,k,a,alpha,c)/(exp(-alpha*(m-j)*c)^(k+a))
    }
    res <- part.one*(part.two - sum) #part one is p(j=0), so P(T = k) then equals sum of 
    #choose term * p term / f term over all values of possible susceptibles to be infected (0 to m)
  }
  return(res)
}


# We wish to compute P(T = k), for k = 0, ..., n
# "n" is the size of the household
# "alpha" is the person to person infection rate
# "c" is the length of the (fixed) infectious period

#P(T = k)
compute.marginal.prob <- function(k, n, c, alpha, p) {  
  prob <- 0;  
  for (y in 0:n) {
    if (k >= y) {
      #n - y = initital susceptibles
      #j = number individiuals infected in SIR
      #y = initial infectives
      #P(T = k | Y = y)
      cond.prob <- compute.conditional.prob(n - y, k - y, y, alpha, c)
    }
    else {
      cond.prob <- 0
    }
    #calculate over total population size, sum conditional probability (P(T=k|Y=y)) times binomial density of y infectives out of n total with probability of becoming infected 1-p
    prob <- prob + cond.prob*dbinom(y, n, 1 - p)
  }
  prob
}

#Compute log likelihood of data
log.like= function(alpha,p){
  29*log(compute.marginal.prob(0,3,1,alpha,p))+
    9*log(compute.marginal.prob(1,3,1,alpha,p))+
    2*log(compute.marginal.prob(2,3,1,alpha,p))+
    2*log(compute.marginal.prob(3,3,1,alpha,p))
}

