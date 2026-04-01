#' graph convolution-based functional co-response group discovery
#' https://doi.org/10.1109/TCBBIO.2025.3560853
#'
#' This function uses a genetic algorithm to discover functional co-response group
#' (columns from matrix \code{M}) that maximizes correlation or association with a functional variable \code{y}.
#'
#' @param M Numeric matrix (samples × taxa). values are topological abundance.
#' @param y Functional variable. Can be continuous, binary or categorical depending on \code{y_type}.
#' @param pk Optional vector of prior knowledge (same length as ncol(M)).
#'        Elements equal to 1 indicate mandatory inclusion, 0 otherwise.
#' @param fit_f Fitness function type. Options: \code{"nmax"} or \code{"regul"}.
#' @param y_type Type of functional variable. Options: \code{"continuous"} or \code{"binary"} or \code{"categorical"}.
#' @param fun_s Direction of association for binary outcome.
#'        \code{"+"} means True is better, \code{"-"} means False is better.
#' @param n_max Maximum number of selected taxa.
#' @param alpha Regularization parameter when \code{fit_f = "regul"}.
#' @param max_iter Maximum number of GA iterations.
#' @param max_fitness Early stopping threshold for fitness value.
#' @param run Number of consecutive generations without improvement before stopping.
#' @param pop_size Population size for genetic algorithm.
#' @param parallel Logical. Whether to enable parallel computing.
#' @param monitor Logical or function. If set to "save", fitness evolution will be stored.
#'
#' @return A list containing:
#' \itemize{
#'   \item fitness: Best fitness value obtained.
#'   \item x: Binary selection vector.
#'   \item members: Names of selected taxa from colomn names of M.
#'   \item abundance: Aggregated topological abundance vector.
#'   \item performance: Correlation (continuous y) or eta coefficient (discrete y).
#'   \item fitnesses: Fitness values across generations (if monitored).
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
#' out <- discover(M, y, n_max=5)
#' }
#'
#' @importFrom stats cor
#' @export
discover <- function (M, y=NULL, pk=NULL, fit_f="nmax", y_type="continuous", fun_s="+", n_max=NULL, alpha=40, max_iter=500, max_fitness=1, run=100, pop_size=200, parallel=TRUE, monitor=FALSE) {
  if (missing(M)) {
    stop("`M` must be provided.")
  }
  if (anyNA(M)) {
    stop("`M` must not contain NA values.")
  }
  M <- as.matrix(M)
  if (!is.numeric(M)) {
    stop("`M` must be a numeric matrix.")
  }

  if (!y_type %in% c("continuous", "binary", "categorical")) {
    stop("`y_type` must be one of 'continuous', 'binary', or 'categorical'.")
  }

  if (is.null(y)) {
    stop("`y` must be provided.")
  }
  if (length(y) != nrow(M)) {
    stop("`y` must have the same length as the number of rows in `M`.")
  }
  if (any(is.na(y))) {
    stop("`y` must not contain NA.")
  }
  if (!is.numeric(y)) {
    stop("`y` must be numeric.")
  }
  if (y_type == "binary") {
    if (!all(y %in% c(0, 1))) {
      stop("When `y_type = 'binary'`, `y` must contain only 0 and 1.")
    }
    if (length(unique(y)) != 2) {
      stop("When `y_type = 'binary'`, `y` must contain both 0 and 1.")
    }
  }
  if (y_type == "categorical") {
    if (any(y != floor(y))) {
      stop("When `y_type = 'categorical'`, numeric `y` must contain integer-coded group labels.")
    }
    if (length(unique(y)) < 2) {
      stop("When `y_type = 'categorical'`, `y` must contain at least two groups.")
    }
  }

  if (!is.null(pk)) {
    if (length(pk) != ncol(M)) {
      stop("`pk` must have length equal to `ncol(M)`.")
    }
    if (!all(pk %in% c(0, 1))) {
      stop("`pk` must contain only 0 and 1.")
    }
  } else {
    pk <- rep(0, ncol(M))
  }

  if (!fit_f %in% c("nmax", "regul")) {
    stop("`fit_f` must be one of 'nmax' or 'regul'.")
  }
  if (!fun_s %in% c("+", "-")) {
    stop("`fun_s` must be '+' or '-'.")
  }
  if (fit_f == "regul") {
    if (y_type != "continuous") {
      stop("When `fit_f = 'regul'`, `y_type` must be 'continuous'.")
    }
    if (fun_s != "+") {
      stop("When `fit_f = 'regul'`, `fun_s` must be '+'.")
    }
  }

  if (is.null(n_max)) {
    n_max <- ncol(M)
  }
  if (!is.numeric(n_max) || length(n_max) != 1 || n_max != floor(n_max)) {
    stop("`n_max` must be a single integer.")
  }
  if (n_max < 1) {
    stop("`n_max` must be at least 1.")
  }
  if (n_max > ncol(M)) {
    stop("`n_max` must not exceed the number of columns of `M`.")
  }

  M = as.matrix(M)
  m <- nrow(M)
  n <- ncol(M)

  eta_f <- function(cont, y) {
    grand_mean <- mean(cont)
    ss_total <- sum((cont - grand_mean)^2)
    if (ss_total == 0) return(0)
    ss_between <- sum(tapply(cont, y, function(g) length(g) * (mean(g)-grand_mean)^2))
    eta <- sqrt(ss_between / ss_total)
    if (y_type == "binary") {
      group_means <- tapply(cont, y, mean)
      mean0 <- unname(group_means["0"])
      mean1 <- unname(group_means["1"])
      if (fun_s == "+") {
        sign_val <- sign(mean1 - mean0)
      } else {
        sign_val <- sign(mean0 - mean1)
      }
      eta = sign_val * eta
    }
    return(eta)
  }

  if (fit_f == "nmax") {
    pen <- sqrt(.Machine$double.xmax)

    c0 <- function(x) {ifelse(min(x-pk) < 0, 1, 0)} # partially known group
    c1 <- function(x) {(rep(1, n) %*% x) - n_max} # maximal number of taxa

    fitness <- function(x) {
      if (y_type == "continuous") {
        if (fun_s == "+") {
          cor(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - (pen*max(c1(x),0)) - (pen*c0(x))  # gFlora
        } else {
          -cor(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - (pen*max(c1(x),0)) - (pen*c0(x))  #gFlora
        }
      } else {
        eta_f(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - (pen*max(c1(x),0)) - (pen*c0(x))  # discrete y
      }
    }
  }
  if (fit_f == "regul") {
    fitness <- function(x) {
      cor(rowSums(cbind(rep(0,m), rep(0,m), M[,x==1])), y) - sum(x)/alpha
    }
  }

  fitnesses <- data.frame(1:pop_size)
  # correlations <- data.frame(1:pop_size)

  if (isTRUE(monitor)) {
    monitor <- function(obj) {
      fitnesses <<- cbind(fitnesses, data.frame(obj@fitness))
      # print(obj@population)
      # correlations_iter <- c()
      # for (i in 1:nrow(obj@population)) {
      #   print(obj@population[i, ])
      #   individual <- obj@population[i, ]
      #   abundance <- rowSums(cbind(rep(0,m),rep(0,m), M[,individual==1]))
      #   r <- cor(abundance, y)
      #   correlations_iter <- c(correlations_iter, r)
      # }
      # correlations <<- cbind(correlations, data.frame(correlations_iter))
      # print(data.frame(correlations_iter))
      # Sys.sleep(20)
    }
  }

  GA <- GA::ga("binary",
               fitness=fitness,
               nBits=n,
               maxiter=max_iter,
               maxFitness=max_fitness,
               run=run,
               popSize=pop_size,
               names=colnames(M),
               monitor=monitor,
               parallel=parallel)
  solution1 = GA@solution[1,]
  x <- as.numeric(solution1)
  fitness <- max(GA@fitness[!is.na(GA@fitness)])
  members <- colnames(M)[solution1==1]
  abundance <- rowSums(cbind(rep(0,m),rep(0,m), M[,solution1==1]))
  if (y_type == "continuous") {
    r <- cor(abundance, y)
    return(list(fitness=fitness, x=x, members=members, abundance=abundance, performance=r, fitnesses=fitnesses))
  } else {
    # discrete y
    eta <- eta_f(abundance, y)
    return(list(fitness=fitness, x=x, members=members, abundance=abundance, performance=eta, fitnesses=fitnesses))
  }
}
