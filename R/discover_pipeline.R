#' Run repeated discovery and build a member co-selection graph
#'
#' Repeatedly calls \code{discover()} for \code{n_iter} iterations, stores the
#' result from each run, aggregates member importance across iterations, and
#' constructs a co-selection graph among selected members.
#'
#' @param n_iter Integer. Number of repeated training iterations.
#'
#' @return A list of length \code{n_iter}. Each element contains:
#' \describe{
#'   \item{iter}{Iteration index.}
#'   \item{taxa}{Selected members of this iteration.}
#'   \item{performance}{Performance value of this iteration.}
#'   \item{result}{Full result returned by \code{discover()}.}
#' }
#'
#' Additional attributes are attached to the returned list:
#' \describe{
#'   \item{nodes}{Node table used for plotting.}
#'   \item{edges}{Edge table used for plotting.}
#'   \item{graph}{A \code{visNetwork} htmlwidget object.}
#' }
#'
#' @details
#' For each iteration, the random seed is set to the iteration index using
#' \code{set.seed(t)} before calling \code{discover()}.
#'
#' Node weights are computed by summing iteration performance values over all
#' iterations in which a taxon is selected. Edge weights are computed by summing
#' the iteration performance over all iterations in which a pair of taxa is
#' co-selected.
#'
#' @importFrom plyr ddply
#' @importFrom reshape2 dcast melt
#' @importFrom visNetwork visNetwork visConfigure visLayout visPhysics visSave
#' @importFrom magrittr %>%
#' @export
discover_pipeline <- function (M, y=NULL, n_iter=10, pk=NULL, fit_f="nmax", y_type="continuous", fun_s="+", n_max=NULL, alpha=40, max_iter=500, max_fitness=1, run=100, pop_size=200, parallel=TRUE) {
  if (!is.numeric(n_iter) || length(n_iter) != 1 || n_iter < 1 || n_iter != as.integer(n_iter)) {
    stop("`n_iter` must be a positive integer.")
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

  res_list <- vector("list", n_iter)
  data_all <- data.frame()
  pairs <- data.frame()
  for (t in seq_len(n_iter)) {
    message(sprintf("Running iteration %d / %d", t, n_iter))
    set.seed(t)

    out <- discover(
      M = M,
      y = y,
      pk = pk,
      fit_f = fit_f,
      y_type = y_type,
      fun_s = fun_s,
      n_max = n_max,
      alpha = alpha,
      max_iter = max_iter,
      max_fitness = max_fitness,
      run = run,
      pop_size = pop_size,
      parallel = parallel
    )

    if (fun_s == "-" && y_type == "continuous"){
      iter_data <- data.table::data.table(taxa=out$members, performance=-out$performance)
    } else {
      iter_data <- data.table::data.table(taxa=out$members, performance=out$performance)
    }
    data_all <- rbind(data_all, iter_data)
    iter_pairs <- expand.grid(iter_data$taxa, iter_data$taxa)
    iter_pairs$value <- iter_data$performance
    pairs <- rbind(pairs, iter_pairs)

    res_list[[t]] <- list(
      iter = t,
      taxa = out$members,
      performance = out$performance,
      result = out
    )
  }

  nodes <- aggregate(performance~taxa, data=data_all, FUN = sum)
  nodes <- nodes[order(nodes$performance, decreasing=TRUE),]
  colnames(nodes) <- c("taxa", "value")
  nodes$value <- nodes$value/max(nodes$value)
  rownames(nodes) <- nodes$taxa
  nodes$id <- rownames(nodes)
  nodes$label <- rownames(nodes)

  weight <- plyr::ddply(pairs, ~Var1+Var2, dplyr::summarise, weight=sum(value))
  weight.filter <- weight
  weight.mat <- reshape2::dcast(weight.filter, Var1~Var2, value.var="weight")
  rownames(weight.mat) <- as.character(weight.mat$Var1)
  weight.mat <- as.matrix(weight.mat[,-1])
  weight.mat[upper.tri(weight.mat)] <- NA
  diag(weight.mat) <- NA
  edges <- reshape2::melt(weight.mat)
  edges <- edges[!is.na(edges$value),]
  colnames(edges) <- c("from","to","width")

  top_ids <- nodes$id[1:n_max]
  nodes$color.background <- ifelse(nodes$id %in% top_ids, "#eca8a6", "#D3D3D3")
  nodes$color.border <- "black"
  edges$color <- "#c5dff4"

  nodes_plot <- nodes[nodes$id %in% top_ids, ]
  edges_plot <- edges[edges$from %in% top_ids & edges$to %in% top_ids, ]

  graph <- visNetwork(nodes_plot, edges_plot) %>%
    visConfigure(enabled=TRUE) %>%
    visLayout(randomSeed=123) %>%
    visPhysics(enabled=FALSE)

  visSave(graph, file = file.path(getwd(), "graph.html"))

  attr(res_list, "nodes") <- nodes_plot
  attr(res_list, "edges") <- edges_plot
  attr(res_list, "graph") <- graph

  return(res_list)
}
