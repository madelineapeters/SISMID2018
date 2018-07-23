####################################
###Compartmental models in pomp
##Boarding-school flu outbreak
####################################

bsflu = read.table('https://kingaa.github.io/sbied/stochsim/bsflu_data.txt')
#B is boys confined to bed, cto boys in convalescence

#Start by treating B as incidence data
#First: estimate parameters of SIR
#Second: decide if SIR model is adequate description of data

#force of infection lambda = muSI = B*I/N
#gamme = muIR
#Assume N is fixed

#Model delta NSI as Binomial(S,1-exp(-lambda*delta_t))
#Model delta NIR  as Binomaial(I,1-exp(-gamma*delta_t))

#Csnippet is piece of C code to specify a model in pomp
#Example: encode SIR model simulator

sir_step = Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
")

#At day zero, assume I = 1 and R = 0, but N may be a parameter to be estimated so S0 = N - 1
#Initialiser C snippet
sir_init = Csnippet("
  S = N - 1;
  I = 1;
  R = 0;
")

#Fold C snippets, with data, into a pomp object:
sir = pomp(bsflu,time="day",t0=0,rprocess=euler.sim(sir_step,delta.t=1/6),
     initializer=sir_init,paramnames=c("N","Beta","gamma"),
     statenames=c("S","I","R"))

#Assume that case reports B result from process by which new infections result in confinement with probability p (prob. that infection is severe enough to be noticed by school)

#Confined cases have much lower transmission rate, treat B as being a count of the number of boys who have moved from I to R over the course of the past day

#Need variable to track daily counts; modify C snippet adding variable H to tally incidence; then replace the rprocess with new one

sir_step = Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init = Csnippet("
  S = N-1;
  I = 1;
  R = 0;
  H = 0;
")

sir = pomp(sir,rprocess=euler.sim(sir_step,delta.t=1/6),initializer=sir_init,
            paramnames=c("Beta","gamma","N"),statenames=c("S","I","R","H"))
#euler.sim is in this case a simulation of state process model

#Model the data B as a binomial process
#Bt ~ Binomial(H(t) - H(t-1),p)
#At time t, variable H we've defined will contain H, not H(t)-H(t-1)
#Overcome by telling pomp that we want H to be set to zero immediately following each observation
#We do this by setting the zeronames argument to pomp

sir = pomp(sir,zeronames="H")

#To include observations in the model, write a dmeasure and an rmeasure component
dmeas = Csnippet('lik = dbinom(B,H,rho,give_log);') #prob density of binomial process
rmeas = Csnippet("B = rbinom(H,rho);") #rho is prob reported

#Then put in pomp object
sir = pomp(sir,rmeasure = rmeas, dmeasure = dmeas,statenames = "H",paramnames = "rho")

#testing the model
#1540 total infections, so N must be larger than that
#final size equation: R0 = -log(1-f)/f where f = R(infinity)/N is the final size of the epidemic
#If the ifnectious period is roughly 1 day, then 1/gamma ~ 1 day and Beta = gamma*R0 ~ 1.5 da
#Can simulate:

sims = simulate(sir,params=c(Beta=2.35,gamma=0.9,rho=0.6,N=2000),nsim=20,as.data.frame=TRUE,include.data=TRUE)
ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=='data'))+geom_line()+guides(color=FALSE)


##################SEIR exercise
sir_step = Csnippet("
                    double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
                    double dN_EI = rbinom(E,1-exp(-alpha*dt));
                    double dN_IR = rbinom(I,1-exp(-gamma*dt));
                    S -= dN_SE;
                    E += dN_SE - dN_EI;
                    I += dN_EI - dN_IR;
                    R += dN_IR;
                    ")

#At day zero, assume I = 1 and R = 0, but N may be a parameter to be estimated so S0 = N - 1
#Initialiser C snippet
sir_init = Csnippet("
                    S = N - 1;
                    E = 0;
                    I = 1;
                    R = 0;
                    ")

#Fold C snippets, with data, into a pomp object:
sir = pomp(bsflu,time="day",t0=0,rprocess=euler.sim(sir_step,delta.t=1/6),
           initializer=sir_init,paramnames=c("N","Beta","alpha","gamma"),
           statenames=c("S","E","I","R"))

#Assume that case reports B result from process by which new infections result in confinement with probability p (prob. that infection is severe enough to be noticed by school)

#Confined cases have much lower transmission rate, treat B as being a count of the number of boys who have moved from I to R over the course of the past day

#Need variable to track daily counts; modify C snippet adding variable H to tally incidence; then replace the rprocess with new one

sir_step = Csnippet("
                    double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
                    double dN_EI = rbinom(E,1-exp(-alpha*dt));
                    double dN_IR = rbinom(I,1-exp(-gamma*dt));
                    S -= dN_SE;
                    E += dN_SE - dN_EI;
                    I += dN_EI - dN_IR;
                    R += dN_IR;
                    H += dN_IR;
                    ")

sir_init = Csnippet("
                    S = N-1;
                    E = 0;
                    I = 1;
                    R = 0;
                    H = 0;
                    ")

sir = pomp(sir,rprocess=euler.sim(sir_step,delta.t=1/6),initializer=sir_init,
           paramnames=c("Beta","alpha","gamma","N"),statenames=c("S","E","I","R","H"))
#euler.sim is in this case a simulation of state process model

#Model the data B as a binomial process
#Bt ~ Binomial(H(t) - H(t-1),p)
#At time t, variable H we've defined will contain H, not H(t)-H(t-1)
#Overcome by telling pomp that we want H to be set to zero immediately following each observation
#We do this by setting the zeronames argument to pomp

sir = pomp(sir,zeronames="H")

#To include observations in the model, write a dmeasure and an rmeasure component
dmeas = Csnippet('lik = dbinom(B,H,rho,give_log);') #prob density of binomial process
rmeas = Csnippet("B = rbinom(H,rho);") #rho is prob reported

#Then put in pomp object
sir = pomp(sir,rmeasure = rmeas, dmeasure = dmeas,statenames = "H",paramnames = "rho")

#testing the model
#1540 total infections, so N must be larger than that
#final size equation: R0 = -log(1-f)/f where f = R(infinity)/N is the final size of the epidemic
#If the ifnectious period is roughly 1 day, then 1/gamma ~ 1 day and Beta = gamma*R0 ~ 1.5 da
#Can simulate:

sims = simulate(sir,params=c(Beta=4.5,alpha=2.5,gamma=0.4,rho=0.7,N=2000),nsim=20,as.data.frame=TRUE,include.data=TRUE)
ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=='data'))+geom_line()+guides(color=FALSE)

##################SEIRRR exercise

sir_step = Csnippet("
                    double N = 763;
                    double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
                    double dN_EI = rbinom(E,1-exp(-mu_E*dt));
                    double dN_IR = rbinom(I,1-exp(-mu_I*dt));
                    double dN_R1R2 = rbinom(R1,1-exp(-mu_R1*dt));
                    double dN_R2R3 = rbinom(R2,1-exp(-mu_R2*dt));
                    S -= dN_SE;
                    E += dN_SE - dN_EI;
                    I += dN_EI - dN_IR;
                    R1 += dN_IR - dN_R1R2;
                    R2 += dN_R1R2 - dN_R2R3;
                    ")

sir_init = Csnippet("
                    S = 762;
                    E = 0;
                    I = 1;
                    R1 = 0;
                    R2 = 0;
                    ")

sir = pomp(sir,rprocess=euler.sim(sir_step,delta.t=1/6),initializer=sir_init,
           paramnames=c("Beta","mu_E","mu_I","mu_R1","mu_R2"),statenames=c("S","E","I","R1","R2"))
#euler.sim is in this case a simulation of state process model

#Model the data B as a binomial process
#Bt ~ Binomial(H(t) - H(t-1),p)
#At time t, variable H we've defined will contain H, not H(t)-H(t-1)
#Overcome by telling pomp that we want H to be set to zero immediately following each observation
#We do this by setting the zeronames argument to pomp

dmeas <- Csnippet("
  lik = dpois(B,rho*R1+1e-6,give_log);
                  ")

rmeas <- Csnippet("
                  B = rpois(rho*R1+1e-6);
                  ")


#Then put in pomp object
flu = pomp(bsflu,times="day",t0=-6,
       rprocess=euler.sim(sir_step,delta.t=1/5),
       initializer=sir_init,rmeasure=rmeas,dmeasure=dmeas,
       statenames=c("S","E","I","R1","R2"),
       paramnames=c("Beta","mu_E","mu_I","mu_R1","mu_R2","rho")
  )

#testing the model
#1540 total infections, so N must be larger than that
#final size equation: R0 = -log(1-f)/f where f = R(infinity)/N is the final size of the epidemic
#If the ifnectious period is roughly 1 day, then 1/gamma ~ 1 day and Beta = gamma*R0 ~ 1.5 da
#Can simulate:

sims = simulate(flu,params=c(Beta=3,mu_E=1,mu_I=0.5,mu_R1=0.25,mu_R2=1/1.8,rho=0.9,N=1275),nsim=20,as.data.frame=TRUE,include.data=TRUE)
ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=='data'))+geom_line()+guides(color=FALSE)


pf = pfilter(flu,Np=5000,params=c(N=1275,Beta=3,mu_I=1/2,mu_E=5,muR1=1/4,mu_R1=1/4,mu_R2=1/1.8,rho=0.9))
logLik(pf)

#graphing the likelihood function: likelihood surface

p = sliceDesign(
  center=c(Beta=2,mu_I=1,mu_E=1,mu_R1=0.24,mu_R2=1/1.8,rho=0.9),
  Beta=rep(seq(0.5,4,length=40),each=3),
  mu_I=rep(seq(0.5,2,length=40),each=3))

library(foreach)
library(doParallel)
registerDoParallel()

set.seed(108028909,kind="L'Ecuyer")

p = foreach(theta=iter(p,'row'),.combine=rbind,
        .inorder=FALSE,
        .options.multicore=list(set.seed=TRUE)
        ) %dopar% {
library(pomp)
pf = pfilter(flu,params=unlist(theta),Np=5000)
theta$loglik = loglike(pf)
theta}

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,
         .options.multicore=list(set.seed=TRUE)
) %dopar% {
  library(pomp)
  pfilter(flu,params=unlist(theta),Np=5000) -> pf
  theta$loglik <- logLik(pf)
  theta
} -> p

bake(file="flu-grid1.rds",seed=421776444,kind="L'Ecuyer",{
  
  expand.grid(Beta=seq(from=1.5,to=5,length=50),
              mu_I=seq(from=0.7,to=4,length=50),mu_E=1,
              mu_R1=1/4,mu_R2=1/1.8,
              rho=0.9) -> p
  
  library(foreach)
  library(doParallel)
  registerDoParallel()
  
  ## Now we do the computation
  foreach (theta=iter(p,"row"),.combine=rbind,.inorder=FALSE,
           .options.multicore=list(set.seed=TRUE)
  ) %dopar% 
  {
    library(pomp)
    pfilter(flu,params=unlist(theta),Np=5000) -> pf
    theta$loglik <- logLik(pf)
    theta
  }
  
})-> p


