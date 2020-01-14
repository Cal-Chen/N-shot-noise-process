### moments of N shot noise processes
library(deSolve)

autoname <- function(prefix,order,N,method='tri'){
  object <- 0
  if(order==0 || order > 2){
    return('Error: order must be 1 or 2')
  }else if(order == 1){
    for(i in 1:N){
      object[i] <- paste(prefix,i,j,sep='')
    }
  }
  else{
    if(method == 'matrix'){
      tmp <- 1
      for(i in 1:N){ 
        for(j in 1:N){
          object[tmp] <- paste(prefix,i,j,sep='')
          tmp <- tmp + 1
        }
      }
    }else{
      tmp <- 1
      for(i in 1:N){ 
        for(j in i:N){
          object[tmp] <- paste(prefix,i,j,sep='')
          tmp <- tmp + 1
        }
      }
    }
  }
  return(object)
}

Moment_1_Nshotnoise <- function(N,c_0,lam,T,delta,rho,alpha){
  
  state <- c(EM = rep(0,N),Elam = lam)
  parameters <- c(c_0=c_0,delta = delta,rho=rho,alpha=alpha)
  
  First_moments_Nshotnoise <- function(t,state,parameters){
    
    EM <- names(state)[1:N]
    Elam <- names(state)[(N+1):(2*N)]
    c_0 <- names(parameters)[1:N]
    
    with(as.list(c(state,parameters)),{
      
      for(i in 1:N){
        assign(paste('d',EM[i],sep=''),eval(parse(text=Elam[i])))
      }
      
      for(i in 1:(N-1)){
        assign(paste('d',Elam[i],sep=''), delta*(eval(parse(text=c_0[i])) -  
                                                   eval(parse(text=Elam[i]))) + 
                 (1/alpha)*  eval(parse(text=Elam[i+1])) ) 
      }
      assign(paste('d',Elam[N],sep=''), delta*(eval(parse(text=c_0[N])) -  
                                                 eval(parse(text=Elam[N]))) + 
               (1/alpha)* rho) 
      
      list(c(sapply(1:N,function(i) eval(parse(text=paste('d',EM[i],sep=''))))
             ,sapply(1:N,function(i) eval(parse(text=paste('d',Elam[i],sep=''))))))
    })
  }
  times <- seq(0,T,0.01)
  out <- ode(y = state,times=times,func = First_moments_Nshotnoise,parms = parameters)
  return(out[length(out[,1]),])
}

Moments_Nshotnoise <- function(N,c_0,lam,T,delta,rho,alpha){
  
  lamij <-0
  tmp <- 1
  for(i in 1:N){
    for(j in i:N){
      lamij[tmp] <- lam[i]*lam[j]
      tmp <- tmp + 1
    }
  }
  
  name_Elam <- autoname('Elam',2,N)
  name_EMlam <- autoname('EMlam',2,N,method='matrix')
  name_EM <- autoname('EM',2,N)
  
  state <- c(a = lamij,b = rep(0,N*N),c=rep(0,N*(N+1)*0.5),EM=rep(0,N),Elam=lam)
  names(state)[1:(N*(N+1)*0.5)] <- name_Elam
  names(state)[(N*(N+1)*0.5 + 1):(N*(3*N+1)*0.5)] <- name_EMlam
  names(state)[(N*(3*N+1)*0.5+1):(N*(2*N+1))] <- name_EM
  parameters = c(c_0=c_0,delta = delta,rho=rho,alpha=alpha)
  
  Second_moments_Nshotnoise <- function(t,state,parameters){
    
    c_0 <- names(parameters)[1:N]
    
    with(as.list(c(state,parameters)),{
      
      for(i in 1:N){
        assign(paste('dEM',i,sep=''),eval(parse(text=paste('Elam',i,sep=''))))
      }
      
      for(i in 1:(N-1)){
        assign(paste('dElam',i,sep=''), delta*(eval(parse(text=c_0[i])) -  
                                        eval(parse(text=paste('Elam',i,sep='')))) + 
                 (1/alpha)*eval(parse(text=paste('Elam',i+1,sep=''))) ) 
      }
      assign(paste('dElam',N,sep=''), delta*(eval(parse(text=c_0[N])) -  
                                        eval(parse(text=paste('Elam',N,sep='')))) + (1/alpha)* rho) 
      
      ### second moments M_ij
      for(i in 1:N){
        for(j in i:N)
          if(j==i){
            assign(paste('dEM',i,i,sep=''), 2*eval(parse(text=paste('EMlam',i,i,sep=''))) + 
                     eval(parse(text=paste('Elam',i,sep=''))))
          }else{
            assign(paste('dEM',i,j,sep=''), eval(parse(text=paste('EMlam',j,i,sep=''))) + 
                     eval(parse(text=paste('EMlam',i,j,sep=''))))
          }
      }
      
      ### moments of M*Lam_ij
      for(j in 1:(N-1)){
        for(i in 1:N){
          if(i == j){
            assign(paste('dEMlam',i,i,sep=''), -delta*eval(parse(text=paste('EMlam',i,i,sep=''))) + 
                     delta*(eval(parse(text=c_0[i])))*eval(parse(text=paste('EM',i,sep='')))  +
                     eval(parse(text=paste('Elam',i,i,sep=''))) + 
                     (1/alpha)*eval(parse(text=paste('EMlam',i,i+1,sep=''))))
            
          }else if(i == N){
            assign(paste('dEMlam',j,N,sep=''),-delta*eval(parse(text=paste('EMlam',j,N,sep=''))) + 
                     delta*(eval(parse(text=c_0[i])))*eval(parse(text=paste('EM',i,sep=''))) +
                     eval(parse(text=paste('Elam',j,N,sep=''))) + 
                     rho*(1/alpha)*eval(parse(text=paste('EM',j,sep=''))) )
            
          }else{
            if(i < j){
              assign(paste('dEMlam',j,i,sep=''),-delta*eval(parse(text=paste('EMlam',j,i,sep=''))) +
                       delta*(eval(parse(text=c_0[i])))*eval(parse(text=paste('EM',i,sep=''))) +
                       eval(parse(text=paste('Elam',i,j,sep=''))) + 
                       (1/alpha)*eval(parse(text=paste('EMlam',j,i+1,sep='')))  )
            }else{
              assign(paste('dEMlam',j,i,sep=''),-delta*eval(parse(text=paste('EMlam',j,i,sep=''))) + 
                       delta*(eval(parse(text=c_0[i])))*eval(parse(text=paste('EM',i,sep=''))) +
                       eval(parse(text=paste('Elam',j,i,sep=''))) + 
                       (1/alpha)*eval(parse(text=paste('EMlam',j,i+1,sep='')))  )
            }
          }
        }
        assign(paste('dEMlam',N,j,sep=''),-delta*eval(parse(text=paste('EMlam',N,j,sep=''))) + 
                 delta*(eval(parse(text=c_0[i])))*eval(parse(text=paste('EM',N,sep=''))) +
                 eval(parse(text=paste('Elam',j,N,sep=''))) + 
                 (1/alpha)*eval(parse(text=paste('EMlam',N,j+1,sep='')))  )
      }
      assign(paste('dEMlam',N,N,sep=''),-delta*eval(parse(text=paste('EMlam',N,N,sep=''))) + 
               delta*(eval(parse(text=c_0[i])))*eval(parse(text=paste('EM',N,sep=''))) +
               eval(parse(text=paste('Elam',N,N,sep=''))) + 
               (rho/alpha)*eval(parse(text=paste('EM',N,sep='')))  )
      
      
      ### moments of lam_i*lam_j
      for(i in 1:(N-1)){
        for(j in i:N){
          if(j == i){
            assign(paste('dElam',j,j,sep=''), -2*delta*eval(parse(text=paste('Elam',j,j,sep='')))+
                     2*(eval(parse(text=c_0[i])))*delta*eval(parse(text=paste('Elam',j,sep='')))+
                     (2/alpha^2)*eval(parse(text=paste('Elam',j+1,sep='')))+ 
                     (2/alpha)*eval(parse(text=paste('Elam',j,j+1,sep=''))) )
            
          }else if(j==N){
            assign(paste('dElam',i,j,sep=''),-2*delta*eval(parse(text=paste('Elam',i,j,sep='')))+
                     eval(parse(text=c_0[i]))*delta*eval(parse(text=paste('Elam',j,sep=''))) + 
                     eval(parse(text=c_0[j]))*delta*eval(parse(text=paste('Elam',i,sep=''))) + 
                     1/alpha * (eval(parse(text=paste('Elam',i+1,j,sep=''))) 
                                + rho*eval(parse(text=paste('Elam',i,sep=''))) )   )
          }else{
            assign(paste('dElam',i,j,sep=''),-2*delta*eval(parse(text=paste('Elam',i,j,sep='')))+
                     eval(parse(text=c_0[i]))*delta*eval(parse(text=paste('Elam',j,sep=''))) + 
                     eval(parse(text=c_0[j]))*delta*eval(parse(text=paste('Elam',i,sep=''))) +
                     1/alpha * (eval(parse(text=paste('Elam',i+1,j,sep=''))) 
                                + eval(parse(text=paste('Elam',i,j+1,sep=''))) )   )
          }
        }
      }
      assign(paste('dElam',N,N,sep=''), -2*delta*eval(parse(text=paste('Elam',N,N,sep='')))+
               2*eval(parse(text=c_0[N]))*delta*eval(parse(text=paste('Elam',N,sep=''))) + 
               (2/alpha^2)*rho + (2/alpha)*rho*eval(parse(text=paste('Elam',N,sep=''))) )
      
      list(c(sapply(1:length(name_Elam),function(i) eval(parse(text=paste('d',name_Elam[i],sep='')))),
             sapply(1:length(name_EMlam),function(i) eval(parse(text=paste('d',name_EMlam[i],sep='')))),
             sapply(1:length(name_EM),function(i) eval(parse(text=paste('d',name_EM[i],sep='')))),
             sapply(1:N,function(i) eval(parse(text=paste('dEM',i,sep='')))),
             sapply(1:N,function(i) eval(parse(text=paste('dElam',i,sep=''))))   ))
    })
  }
  
  times <- seq(0,T,0.01)
  out <- ode(y = state,times=times,func = Second_moments_Nshotnoise,parms = parameters)
  result_tmp <- out[length(out[,1]),]
  result1 <- result_tmp[ sapply(1:N,function(i) paste('EM',i,sep=''))]
  result2 <- result_tmp[ sapply(1:N,function(i) paste('EM',i,i,sep=''))]

  return(list(Moments_1=result1,Moments_2=result2,other=result_tmp))
}


N <- 10
c_0 <- rep(1,N)
lam <- rep(2,N)
T <- 1
delta <- 0.3
alpha <- 2
rho <- 2

name_Elam1 <- names(c(Elam=rep(0,N)))
result <- Moment_1_Nshotnoise(N,c_0,lam,T,delta,rho,alpha)
result[c('Elam1','EM1')]
## theorectical result for delta != 1/alpha
Elambda <- function(t){
  (lam[1]-c_0[1]/(1-(1/(alpha*delta))))*exp(((1/alpha)-delta)*t) + c_0[1]/(1-(1/(alpha*delta)))
}

Elambda(1)
integrate(Elambda,0,1)

# theorectial result for delta = 1/alpha
power_series <- function(x,n){
  result <- sapply(0:n, function(i) (x^i)/gamma(i+1))
  return(sum(result))
}

Reminder <- function(x,n){
  result <- sapply(0:n,function(i) exp(x) - power_series(x,i))
  return(sum(result))
}

Elambda2 <- function(t){
  if(1/alpha != delta){
    return('Error: mu_y and delta must be the same')
  }else{
    result <- lam[1]*power_series(delta*t,N)*exp(-delta*t) + c_0[1]*exp(-delta*t)*Reminder(delta*t,N)
    return(result)
  }
}
Elambda2(T)
0.5*T*sum(weight * sapply(T*0.5*nodes + T*0.5,Elambda2))


### theoretical lower bound for sum lambda_i
sum_lambda_L <- function(t){
  mu <- 1/alpha
  result <- mu*rho / (delta-mu) *(1-exp((mu-delta)*t)) - rho *exp((mu-delta)*t)*(mu^(N+1))/gamma(N+1)* 
  0.5*t*sum(weight * sapply(t*0.5*nodes + t*0.5,function(x) exp(delta-mu*x)*x^N)) - lam[1]*mu*t*exp((mu-delta)*t) + sum(lam)*exp((mu-delta)*t)
  return(result)
}

### theoretical upper bound for sum lambda_i
sum_lambda_U <- function(t){
  mu <- 1/alpha
  result <- mu*rho / (delta-mu) *(1-exp((mu-delta)*t)) - rho *exp((mu-delta)*t)*(mu^(N+1))/gamma(N+1)* 
    0.5*t*sum(weight * sapply(t*0.5*nodes + t*0.5,function(x) exp(-mu*x)*x^N)) - lam[1]*mu*t*exp((mu-delta)*t) + sum(lam)*exp((mu-delta)*t)
  return(result)
}


N <- 10
c_0 <- rep(0,N)
lam <- rep(10/N,N)
T <- 1
delta <- 0.3
alpha <- 2
rho <- 2

name_Elam1 <- names(c(Elam=rep(0,N)))
result <- Moment_1_Nshotnoise(N,c_0,lam,T,delta,rho,alpha)

sum_lambda_U(T)
sum_lambda_L(T)
sum(result[name_Elam1])


##convergence test of moments
N_seq <- 1:20
T <- 1
delta <- 200
alpha <- 0.2
rho <- 5
Moments_seq2 <- rep(0,length(N_seq)-1)


for(i in 2:length(N_seq)){
  N <- N_seq[i]
  c_0 <- rep(2,N)
  lam <- rep(5,N)
  Moments_seq2[i-1] <- Moment_1_Nshotnoise(N,c_0,lam,T,delta,rho,alpha)['EM1']
}


Moments_seq2
par(mfrow=c(1,2))
plot(diff(abs(Moments_seq2)))
plot(Moments_seq2)





  

  
  
  
  