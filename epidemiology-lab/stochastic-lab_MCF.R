## following along with: Stochastic models in R.Rmd

## Gillepsie Algorithm
gillesp <- function(start,ratefun,trans,pars,times=0:20) {
  t0 <- times[1]                 ## set time to starting time
  ntimes <- length(times) 
  X <- start                     ## set state to starting state
  res <- matrix(nrow=length(times),ncol=length(start),dimnames=list(times,names(start))) 
  ## matrix for results
  for (ctr in 1:(ntimes-1)) {     ## loop over reporting times
    res[ctr,] <- X                ## record current state
    while (t0<times[ctr+1]) {
      rates <- ratefun(X,pars,t0) ## calculate current rates
      if (all(rates==0)) break    ## extinction
      totrate <- sum(rates)       
      elapsed <- rexp(1,totrate)  ## sample elapsed time
      which.trans <- sample(1:nrow(trans),size=1,prob=rates) ## pick transition
      t0 <- t0+elapsed            ## update time
      X <- X+trans[which.trans,]  ## add transition values to current state
    }  }
  cbind(times,res)
}

## setting random seed
set.seed(2789)

#######################################################
####   simulating stochastic SIR model
#######################################################

## initial conditions
start=c(S=100,I=1,R=0)

## definition rates of transition between states
ratefun.SIR = function(X,pars,time)  {
  vals = c(as.list(pars),as.list(X))   ## attach state and pars as lists
  rates = with(vals,                ##allows reference to states and parameters by name 
               c(infection=beta*S*I, death=alpha*I, recovery=gamma*I))
}

## matrix defining states
statenames.SIR = c("S","I","R")       
transnames.SIR = c("infection","death","recovery")

trans.SIR = matrix(c(-1,1,0,0,-1,0, 0,-1,1),
                   byrow=TRUE,        ## default is by column
                   ncol=3,            ## number of columns = number of state variables
                   dimnames=list(transnames.SIR,statenames.SIR))

## vector naming & setting parameter values
pars.SIR= c(beta=0.1,alpha=1,gamma=1)

## vector defining time 
times= seq(0, 5, by = 0.05)

## simulate in silico trials
G.SIR.mult= replicate(100, gillesp(start = c(S = 100, I = 3, R = 0), times = seq(0, 5, by = 0.05),
                                   ratefun = ratefun.SIR, trans = trans.SIR, pars = pars.SIR)[, "I"])
                                # We add [,"I"] to save just the number of infectious individuals. 

## plotting stochastic trajectories
matplot(times,G.SIR.mult,type="l",col="gray",lty=1, xlab="Time", ylab="Number infectious")
lines(times,rowMeans(G.SIR.mult),lwd=2)  ## S1


#######################################################
####   exercises
#######################################################




