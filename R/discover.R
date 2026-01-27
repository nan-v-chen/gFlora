#' graph convolution-based functional co-response group discovery
#' https://doi.org/10.1109/TCBBIO.2025.3560853
#'
#' This function uses a genetic algorithm to discover functional co-response group
#' (columns from matrix \code{M}) that maximizes correlation or association with a functional variable \code{y}.
#'
#' @param M Numeric matrix (samples × taxa). values are topological abundance.
#' @param y Functional variable. Can be continuous, binary or categrial depending on \code{ytype}.
#' @param pk Optional vector of prior knowledge (same length as ncol(M)).
#'        Elements equal to 1 indicate mandatory inclusion, 0 otherwise.
#' @param fitf Fitness function type. Options: \code{"nmax"} or \code{"regul"}.
#' @param ytype Type of functional variable. Options: \code{"continuous"} or \code{"binary"} or \code{"categrial"}.
#' @param funcS Direction of association for binary outcome.
#'        \code{"+"} means True is better, \code{"-"} means False is better.
#' @param Nmax Maximum number of selected taxa.
#' @param alpha Regularization parameter when \code{fitf = "regul"}.
#' @param maxIter Maximum number of GA iterations.
#' @param maxFitness Early stopping threshold for fitness value.
#' @param run Number of consecutive generations without improvement before stopping.
#' @param popSize Population size for genetic algorithm.
#' @param parallel Logical. Whether to enable parallel computing.
#' @param monitor Logical or function. If set to "save", fitness evolution will be stored.
#'
#' @return A list containing:
#' \itemize{
#'   \item fitness: Best fitness value obtained.
#'   \item x: Binary selection vector.
#'   \item members: Names of selected taxa.
#'   \item abundance: Aggregated topological abundance vector.
#'   \item performance: Correlation (continuous y) or eta coefficient (discrete y).
#'   \item fitnesses: Fitness values across generations (if monitored).
#'   \item correlations: Correlations across generations (if monitored).
#' }
#'
#' @details
#' The optimization problem is solved using the \code{GA::ga()} function with binary encoding.
#' For continuous response variables, Pearson correlation is used as fitness.
#' For discrete outcomes, the eta coefficient is computed as fitness.
#'
#' @seealso \code{\link[GA]{ga}}
#'
#' @examples
#' \dontrun{
#' data <- gFlora::load_example_data()
#' M <- data$M
#' y <- data$y
#' out <- gFlora::discover(M, y, Nmax=5)
#' }
#'
#' @importFrom stats cor
#' @export
discover <- function(M, y=NULL, pk=NULL, fitf="nmax", ytype="continuous", funcS="+", Nmax=NULL, alpha=40, maxIter=500, maxFitness=1, run=100, popSize=200, parallel=TRUE, monitor=FALSE){

  if(missing(pk)){pk <- rep(0, ncol(M))}
  if(missing(Nmax)){Nmax <- ncol(M)}

  M = as.matrix(M)
  m <- nrow(M)
  n <- ncol(M)

  eta_f <- function(cont, y){
    grand_mean <- mean(cont)
    ss_total <- sum((cont - grand_mean)^2)
    if (ss_total == 0) return(0)
    ss_between <- sum(tapply(cont, y, function(g) length(g) * (mean(g)-grand_mean)^2))
    eta <- sqrt(ss_between / ss_total)
    if (ytype == "binary"){
      # means of cont related to 0 and 1 in y (in order 0, 1)
      means <- tapply(cont, y, mean)
      if (funcS == "+"){
        # 1 is better than 0
        sign_val <- unname(sign(means[2] - means[1]))
      } else {
        # 0 is better than 1
        sign_val <- unname(sign(means[1] - means[2]))
      }
      eta = sign_val * eta
    }
    return(eta)
  }

  if (fitf == "nmax"){
    pen <- sqrt(.Machine$double.xmax)

    c0 <- function(x){ifelse(min(x-pk) < 0, 1, 0)} # partially known group
    c1 <- function(x){(rep(1, n) %*% x) - Nmax} # maximal number of taxa

    fitness <- function(x){
      if (ytype == "continuous"){
        if (funcS == "+"){
          cor(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - (pen*max(c1(x),0)) - (pen*c0(x))  # gFlora
        } else {
          -cor(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - (pen*max(c1(x),0)) - (pen*c0(x))  #gFlora
        }
      } else {
        eta_f(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - (pen*max(c1(x),0)) - (pen*c0(x))  # discrete y
      }
    }
  }
  if (fitf == "regul"){
    fitness <- function(x){
      cor(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - sum(x)/alpha
    }
  }

  fitnesses <- data.frame(1:popSize)
  correlations <- data.frame(1:popSize)

  if (monitor == "save"){
    monitor <- function(obj)
    {
      fitnesses <<- cbind(fitnesses, data.frame(obj@fitness))
      # print(obj@population)
      correlations_iter <- c()
      for (i in 1:nrow(obj@population)) {
        # print(obj@population[i, ])
        individual <- obj@population[i, ]
        abundance <- rowSums(cbind(rep(0,m),rep(0,m), M[,individual==1]))
        r <- cor(abundance, y)
        correlations_iter <- c(correlations_iter, r)
      }
      correlations <<- cbind(correlations, data.frame(correlations_iter))
      # print(data.frame(correlations_iter))
      # Sys.sleep(20)
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
  solution1 = GA@solution[1,]
  x <- as.numeric(solution1)
  fitness <- max(GA@fitness[!is.na(GA@fitness)])
  members <- colnames(M)[solution1==1]
  abundance <- rowSums(cbind(rep(0,m),rep(0,m), M[,solution1==1]))
  if (ytype == "continuous"){
    r <- cor(abundance, y)
    return(list(fitness=fitness, x=x, members=members, abundance=abundance, performance=r, fitnesses=fitnesses, correlations=correlations))
  } else {
    # discrete y
    eta <- eta_f(abundance, y)
    return(list(fitness=fitness, x=x, members=members, abundance=abundance, performance=eta, fitnesses=fitnesses, correlations=correlations))
  }
}
