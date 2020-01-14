########################################
####solving the n shot noise processes##
########################################

library(deSolve)

#### fitting a lagrange polynomial(Newton interpolation)
# return the coefficients of (1, x, x^2, ..., x^n)
New_poly <- function(x,y){
  if(length(x) != length(y)){
    return('Error: the length of x and y must be equal')
  }else{
    x_sam <- x
    y_sam <- y
    n <- length(x)
    
    ##calculate the a_i
    a <- rep(0,n)
    a[1] <- y_sam[1]
    a[2] <- (y_sam[2] - a[1])/ (x_sam[2] - x_sam[1])
    
    for(i in 3:n){
      x_sub <- x_sam[1:(i-2)]
      x_sub2 <- x_sam[1:(i-1)]
      nom <- y_sam[i] - a[1] - sum(a[2:(i-1)]*cumprod(x_sam[i]-x_sub))
      deno <- prod(x_sam[i] - x_sub2)
      a[i] <- nom/deno
    }
    #### calculate the coefficients of each x^k
    coef_x <- matrix(0,n,n)
    coef_x[1,1] <- a[1]
    for(i in 2:n){
      x_sub <- x_sam[1:(i-1)]
      vec <- 1:(i-1)
      for(j in 1:i){
        if(j == i){
          coef_x[i,j] <- a[i]
        }else{
          choice <- combn(vec,(i-j))
          choice_x <- matrix(x_sub[choice],nrow(choice),ncol(choice))
          coef_x[i,j] <- a[i]*((-1)^(i-j))*sum(apply(choice_x,2,prod))
        }
      }
    }
    result <- colSums(coef_x)
    return(result)
  }
}

## Richardson Extrapolation for derivatives
Rextraploation <- function(f,x,n,dx=0.1){
  d <- matrix(0,n,n)
  d[1,1] <- (f(x+dx) - f(x-dx))/(2*dx)
  for(i in 2:n){
    d[i,1] <- ( f(x + dx*(0.5^(i-1))) - f(x-dx*(0.5^(i-1))) )/(2*dx*(0.5^(i-1)))
    for(j in 2:i){
      d[i,j] <- (4^(j-1))/(4^(j-1) - 1) * d[i,j-1] -  d[i-1,j-1]/(4^(j-1) - 1)
    }
  }
  return(d[n,n])
}

# laplace transform of the jump distribution
laplace_jump <- function(f,alpha){
  if(f+alpha < 0 ){
    print('Error: lapace transfrom does not exist!')
  }else{
    return( alpha/(alpha + f) )
  }
}

pro_gen2 <- function(position,value,N,c_0,lambda,T,delta,rho,alpha){
  
  theta <- rep(1,N)
  theta[position] <- value
  state = c(a = rep(0,N),R=0)
  parameters = c(theta =theta,c_0=c_0,delta = delta,rho=rho,alpha=alpha)
 
  
  Nshotnoise <- function(t,state,parameters){
    a <- names(state)[1:N]
    theta <- names(parameters)[1:N]
    c_0 <- names(parameters)[(N+1):(2*N)]
    
    with(as.list(c(state,parameters)),{
      
      da1 <- 1 - theta1 - delta*a1
      for(i in 2:N){
        assign(paste('d',a[i],sep=''), 1 - eval(parse(text=theta[i]))*
                 laplace_jump(eval(parse(text=a[i-1])),alpha) - delta*eval(parse(text=a[i]))    )
      }
      dR <- rho*(1- laplace_jump(eval(parse(text=a[N])),alpha)) +
        delta * sum(sapply(1:N,function(i) eval(parse(text=a[i]))) *
                      sapply(1:N,function(i) eval(parse(text=c_0[i])))  )
      list(c(da1,sapply(2:N,function(i) eval(parse(text=paste('d',a[i],sep='')))),dR))
    })
  }
  
  times <- seq(0,T,0.01)
  out <- ode(y = state,times=times,func = Nshotnoise,parms = parameters)
  A <- out[,2:(N+1)]
  R <- out[,length(out[1,])]
  result <- exp(-sum(A[length(A[,1]),]*lam) - R[length(R)])
  return(list(result = result,A=A,R=R))
  
}

prob2 <- function(position,No_points=20){
  result <- 0
  f_tmp <- function(x) pro_gen2(position,x,N,c_0,lam,T,delta,rho,alpha)$result
  x <- seq(-1,1,length.out = No_points) ## |theta| <= 1 in complex plane
  y <- sapply(x, f_tmp)
  result <- New_poly(x,y)
  return(result)
}

N <- 3
c_0 <- rep(1,N)
lam <- rep(2,N)
T <- 1
delta <- 3
alpha <- 0.5
rho <- 2
order <- 15


probability2 <- matrix(0,order,2^N - 1)
Moments <- matrix(0,2,2^N - 1)


tmp <- 1
for(i in 1:N){
  choice <- combn(1:N,i)
  for(j in 1:length(choice[1,])){
    probability2[,tmp] <- prob2(choice[,j],20)[1:order]
    # Moments[1,tmp] <- Rextraploation(function(x) pro_gen2(choice[,j],x,N,c_0,lam,T,delta,rho,alpha)$result,1,3)
    # Moments[2,tmp] <- Moments[1,tmp] + Rextraploation(function(x) Rextraploation(function(x) pro_gen2(choice[,j],x,N,c_0,lam,T,delta,rho,alpha)$result,x,3),1,3)
    # # median2[tmp] <- which(cumsum(probability2[,tmp])>0.5)[1]
    tmp <- tmp +1
  }
}
probability2
Moments




N_seq <- 1:50
T <- 1
delta <- 0.3
alpha <- 2
rho <- 5
Moments_seq1 <- rep(0,length(N_seq)-1)

for(i in 2:length(N_seq)){
  print(i)
  N <- N_seq[i]
  c_0 <- rep(2,N)
  lam <- rep(3,N)
  Moments_seq1[i-1] <- Rextraploation(function(x) pro_gen2(1,x,N,c_0,lam,T,delta,rho,alpha)$result,1,3)
}



