gFlora<-function(M,y,fitf="kopt",k=NULL,mu=40,maxIter=500,maxFitness=0.7,run=20,popSize=200,parallel=TRUE,monitor=FALSE){

  if(missing(k)){k <- ncol(M)}
  
  M = as.matrix(M)
  m <- nrow(M)
  n <- ncol(M)

  if (fitf == "kopt"){
    alpha <- sqrt(.Machine$double.xmax)
    
    penf <- function(x){(rep(1, n) %*% x) - k} # maximal number of taxa
    
    M0 <- apply(M, 2, function(x){x - mean(x)})
    y0 <- y - mean(y)
    P <- t(M0) %*% M0
    Q <- t(M0) %*% y0 %*% t(y0) %*% M0
    Q2 <- t(M0) %*% y0
    
    fitness <- function(x){
      # (t(x) %*% Q2)/sqrt((t(x) %*% P) %*% x)-(alpha*max(penf(x),0))
      # ((t(x) %*% Q) %*% x)/((t(x) %*% P) %*% x)-(alpha*max(penf(x),0))
      cor(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - (alpha*max(penf(x),0))
    }
  }
  if (fitf == "l1"){
    fitness <- function(x){
      cor(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - sum(x)/mu
    }
  }
  
  fitnesses <- data.frame(1:popSize)
  correlations <- data.frame(1:popSize)
  
  if (monitor == "save"){
    monitor <- function(obj) 
    { 
      fitnesses <<- cbind(fitnesses, data.frame(obj@fitness))
      correlations_iter <- c()
      for (i in 1:nrow(obj@population)) {
        individual <- obj@population[i, ]
        abundance <- rowSums(cbind(rep(0,m),rep(0,m), M[,individual==1]))
        r <- cor(abundance, y)
        correlations_iter <- c(correlations_iter, r)
      }
      correlations <<- cbind(correlations, data.frame(correlations_iter))
    }
  }
  
  GA <- GA::ga("binary",
             fitness=fitness,
             nBits=n,
             maxiter=maxIter,
             maxFitness=maxFitness,
             run=run,
             popSize=popSize,
             names=colnames(M),
             monitor=monitor,
             parallel=parallel)
  

  x <- as.numeric(GA@solution)
  fitness <- max(GA@fitness[!is.na(GA@fitness)])
  members <- colnames(M)[GA@solution==1]
  abundance <- rowSums(cbind(rep(0,m),rep(0,m), M[,GA@solution==1]))
  r <- cor(abundance, y)
  return(list(fitness=fitness, x=x, members=members, abundance=abundance, performance=r, fitnesses=fitnesses, correlations=correlations))
}




