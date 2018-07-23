#Load dataframe and rename columns for easier reference
RFD<-read.csv(paste("~/Desktop/Mideo.lab2/Parasite.maturation/Real_Fucking_Data.csv",sep="/"))
names(RFD)[5:7]<-c("ddr","dds","ddt") #donor cell-donor parasites rings and late stages
names(RFD)[9:11]<-c("drr","drs","drt") #donor cell-recipient parasites rings and late stages
names(RFD)[13:15]<-c("rdr","rds","rdt") #recipient cell-donor parasites rings and late stages
names(RFD)[17:19]<-c("rrr","rrs","rrt") #recipient cell-recipient parasites rings and late stages
RFD<-RFD[,2:21] #remove unnecessary columns

#Create new columns in dataframe
RFD<-RFD %>% 
  mutate(., totD = DcDp+DcRp+DcNp) %>% #total donor cell population
  mutate(., totR = RcDp+RcRp+RcNp) %>% #total host cell population
  mutate(.,Dpara=100*DcDp/totD) %>% #parasitaemia of donor cell population w/ donor parasite
  mutate(.,Rpara=100*RcDp/totR) %>% #parasitaemia of host cell population w/ donor parasite
  mutate(., DparaR = 100*DcRp/totD) %>% #parasitaemia of donor cell population w/ host parasite
  mutate(., RparaR = 100*RcRp/totR) #parasitaemia of host cell population w/ host parasite

#Create individual dataframes for host types 
ADcDp<-RFD[51:100,] %>% filter(.,(Mouse=="I1")&(Time<=24)) %>% mutate(., R = ddr*Dpara/100) %>% select(.,Mouse,Time,R)

#############################

rproc = Csnippet("
    double PD = 45;
    double mu_U = 0.025;
    double dU = rbinom(U,1-exp(-mu_U*dt));
    double dL1 = rbinom(L1,1-exp(-mu_L*dt));
    double dL2 = rbinom(L2,1-exp(-mu_L*dt));
    double dL3 = rbinom(L3,1-exp(-mu_L*dt));
    double dL4 = rbinom(L4,1-exp(-mu_L*dt));
    double dL5 = rbinom(L5,1-exp(-mu_L*dt));
    double dL6 = rbinom(L6,1-exp(-mu_L*dt));
    double dL7 = rbinom(L7,1-exp(-mu_L*dt));
    double dL8 = rbinom(L8,1-exp(-mu_L*dt));
    double dL9 = rbinom(L9,1-exp(-mu_L*dt));
    double dL10 = rbinom(L10,1-exp(-mu_L*dt));
    double dL11 = rbinom(L11,1-exp(-mu_L*dt));
    double dL12 = rbinom(L12,1-exp(-mu_L*dt));
    double dR1 = R1;
    double dR2 = R2;
    double dR3 = R3;
    double dR4 = R4;
    double dR5 = R5;
    double dR6 = R6;
    double dR7 = R7;
    double dR8 = R8;
    double dR9 = R9;
    double dR10 = R10;
    double dR11 = R11;
    double dR12 = R12;

    U -= dU - B*dL12;
    R1 += B*dL12 - dR1;
    R2 += dR1 - dR2;
    R3 += dR2 - dR3;
    R4 += dR3 - dR4;
    R5 += dR4 - dR5;
    R6 += dR5 - dR6;
    R7 += dR6 - dR7;
    R8 += dR7 - dR8;
    R9 += dR8 - dR9;
    R10 += dR9 - dR10;
    R11 += dR10 - dR11;
    R12 += dR11 - dR12;
    L1 += dR12 - dL1;
    L2 += dL1 - dL2;
    L3 += dL2 - dL3;
    L4 += dL3 - dL4;
    L5 += dL4 - dL5;
    L6 += dL5 - dL6;
    L7 += dL6 - dL7;
    L8 += dL7 - dL8;
    L9 += dL8 - dL9;
    L10 += dL9 - dL10;
    L11 += dL10 - dL11;
    L12 += dL11 - dL12;
    Tot += R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+U - Tot;

")

#Initialiser C snippet
init = Csnippet("
    double PD = 45;
    U = 1000;
    R1 = PD/24;
    R2 = PD/24;
    R3 = PD/24;
    R4 = PD/24;
    R5 = PD/24;
    R6 = PD/24;
    R7 = PD/24;
    R8 = PD/24;
    R9 = PD/24;
    R10 = PD/24;
    R11 = PD/24;
    R12 = PD/24;
    L1 = PD/24;
    L2 = PD/24;
    L3 = PD/24;
    L4 = PD/24;
    L5 = PD/24;
    L6 = PD/24;
    L7 = PD/24;
    L8 = PD/24;
    L9 = PD/24;
    L10 = PD/24;
    L11 = PD/24;
    L12 = PD/24;
    Tot = R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+U;
")

dmeas <- Csnippet("
  lik = dpois(R,rho*((R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12)/Tot)+1e-6,give_log);
                  ")

rmeas <- Csnippet("
                  R = rpois(rho*((R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12)/Tot)+1e-6);
                  ")

#Fold C snippets, with data, into a pomp object:
pomp.obj = pomp(select(ADcDp,-Mouse),time="Time",
                t0=1,
                rprocess=euler.sim(rproc,delta.t=3),
                rmeasure=rmeas,
                dmeasure=dmeas,
                initializer=init,paramnames=c("B","mu_L","rho"),
                statenames=c("U","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11","R12","L1","L2","L3","L4","L5","L6","L7","L8","L9","L10","L11","L12","Tot")
)


#testing the model
#1540 total infections, so N must be larger than that
#final size equation: R0 = -log(1-f)/f where f = R(infinity)/N is the final size of the epidemic
#If the ifnectious period is roughly 1 day, then 1/gamma ~ 1 day and Beta = gamma*R0 ~ 1.5 da
#Can simulate:

sims = simulate(pomp.obj,params=c(B=0.32,rho=0.6,mu_L=0),nsim=20,as.data.frame=TRUE,include.data=FALSE)
ggplot(sims,mapping=aes(x=time,y=R,group=sim,color=sim=='data'))+geom_line()+guides(color=FALSE)


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


