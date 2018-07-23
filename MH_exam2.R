## This script illustrates the Metropolis-Hastings algorithm for
## approximating the posterior distribution of the time of infection
## in a simple SIS model
## Author: Vladimir N. Minin
## last update: 07/11/2018

sis_log_like = function(inf_time, inf_rate, clear_rate, total_time){
  return(log(inf_rate) - inf_rate*inf_time - clear_rate*(total_time-inf_time))
}

# finish this function
sis_proposal = function(cur_state, total_time, delta){
  y = U1.fun(cur_state,delta)
  y = if((y>0)&&(y<total_time)){y} else if(y>total_time){2*total_time-y} else if(y<0){-y}
  
  return(y)
}


inf_time_mcmc = function(start_inf_time, inf_rate, clear_rate, total_time, delta, chain_len){
  
  result_mat = matrix(0, chain_len, 3)
  colnames(result_mat) = c("inf_time", "log_like", "acc_ind")
  
  cur_inf_time = start_inf_time
  result_mat[1,1] = start_inf_time
  result_mat[1,2] = sis_log_like(start_inf_time, inf_rate, clear_rate, total_time)
  
  for (n in 2:chain_len){
    tp = sis_proposal(result_mat[n-1,1],total_time,delta)
    MH.ratio = sis_log_like(tp,inf_rate,clear_rate,total_time)-sis_log_like(result_mat[n-1,1],inf_rate,clear_rate,total_time)
    a.ij = min(MH.ratio,1)
    
    if(log(runif(1))<=a.ij){
      result_mat[n,1] = tp
      result_mat[n,3] = 1
      result_mat[n,2] = sis_log_like(result_mat[n,1],inf_rate,clear_rate,total_time)
    } else {
      result_mat[n,1] = result_mat[n-1,1]
      result_mat[n,3] = 0
      result_mat[n,2] = sis_log_like(result_mat[n,1],inf_rate,clear_rate,total_time)
    }
  }
  
  return(result_mat)
}


## run the above functions
test_sample = inf_time_mcmc(start_inf_time=0.1, inf_rate=0.1, clear_rate=2, total_time=1.0, delta=0.4, chain_len=10000)
mean(test_sample[,3])
plot(test_sample[,1],type='l')

summary(test_sample[1000:10000,])

hist(test_sample[1000:10000,1])

