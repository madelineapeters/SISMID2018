## This script illustrates the ergodic theorem using the Ehrenfest model of diffusion
## Author: Vladimir N. Minin
## last update: 07/11/2018

## this function randomly draws a new state of the Ehrenfest model
## One of your tasks is to finish writing this function

#state space # particles on left
next_state = function(cur.state, num.mol){
  return.value=cur.state+sample(c(-1,1),1,prob=c(cur.state/num.mol,1-(cur.state/num.mol)))
  
  return(return.value)
}


## set the number of molecules and the number of iterations
my_num_mol = 100
sim_size = 100000

## initialize the chain by drawing the initial state uniformly at random from all possible states. R function `sample()' will be handy.
initital.state = sample(1:my_num_mol,1)
my_draws = numeric(sim_size)
my_draws[1] = initital.state

## run the Markov chain

for (i in 2:sim_size){
  ## use next.state function to evolve the Markov chain one step at a time
  
  my_draws[i] = next_state(my_draws[i-1],my_num_mol)
}

## plot the chain

plot(1:sim_size, my_draws, type="l", ylab="Ehhenfest State", xlab="Time Step")

## get the "time averages"

mean(my_draws)
var(my_draws)

## get the "space averages"
my_num_mol/2
my_num_mol/4
