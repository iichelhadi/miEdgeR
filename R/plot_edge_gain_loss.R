#' Plot Edge Gain/Loss Across Consecutive Pseudotime Bins
#'
#' Plots the count of stable, gained, and lost edges between consecutive pseudotime bins
#' and returns edge details (from/to and optional weights).
#'
#' @param pt_results A named list of pseudotime-bin results. Each element can be an igraph
#'   or a list containing $graph (igraph).
#' @param cluster_id Character label used only for plot annotation (default: "C1").
#' @param title Plot title.
#' @return A list with: plot (ggplot or NULL) and edge_details (stable/gained/lost per bin-pair).
#' @export
plot_edge_gain_loss <- function(pt_results,
                                cluster_id = "C1",
                                title = "Edge Gain/Loss/Stable Counts Per Bin") {

  if (!is.list(pt_results) || length(pt_results) < 2) {
    stop("pt_results must be a list with at least two bins.")
  }

  # Helper: extract igraph
  get_graph <- function(x) {
    if (inherits(x, "igraph")) return(x)
    if (is.list(x) && inherits(x$graph, "igraph")) return(x$graph)
    NULL
  }

  # Helper: order bins if names look like PT1, PT2, ...
  bin_names <- names(pt_results)
  if (is.null(bin_names) || any(!nzchar(bin_names))) {
    stop("pt_results must be a *named* list (bin names required).")
  }
  if (all(grepl("^PT\\d+$", bin_names))) {
    ord <- order(as.integer(sub("^PT", "", bin_names)))
    bin_names <- bin_names[ord]
  }

  # Helper: canonical edge keys + optional weights
  edge_df_from_graph <- function(g) {
    if (igraph::ecount(g) == 0) {
      return(data.frame(from = character(0), to = character(0), weight = numeric(0), stringsAsFactors = FALSE))
    }
    g <- g
    if (!("weight" %in% igraph::edge_attr_names(g))) {
      igraph::E(g)$weight <- 1
    }
    df <- igraph::as_data_frame(g, what = "edges")
    if (!("weight" %in% colnames(df))) df$weight <- 1
    df$from <- as.character(df$from)
    df$to   <- as.character(df$to)
    df$w    <- as.numeric(df$weight)
    df$key  <- paste0(pmin(df$from, df$to), "||", pmax(df$from, df$to))
    df
  }

  gl_data <- list()
  edge_details <- list()

  for (i in seq_len(length(bin_names) - 1)) {
    bin1 <- bin_names[i]
    bin2 <- bin_names[i + 1]
    pair_name <- paste(bin1, "to", bin2)

    g1 <- get_graph(pt_results[[bin1]])
    g2 <- get_graph(pt_results[[bin2]])

    if (is.null(g1) || is.null(g2)) {
      message("Skipping ", pair_name, ": bin object is not igraph (or list with $graph).")
      next
    }

    df1 <- edge_df_from_graph(g1)
    df2 <- edge_df_from_graph(g2)

    k1 <- unique(df1$key)
    k2 <- unique(df2$key)

    stable_k <- intersect(k1, k2)
    gained_k <- setdiff(k2, k1)
    lost_k   <- setdiff(k1, k2)

    # Build edge detail tables (include weights if present; weights come from each bin)
    get_edges_by_keys <- function(df, keys) {
      if (length(keys) == 0) return(NULL)
      out <- df[df$key %in% keys, c("from", "to", "w", "key"), drop = FALSE]
      # keep one row per key (in case of duplicates)
      out <- out[!duplicated(out$key), , drop = FALSE]
      colnames(out)[colnames(out) == "w"] <- "weight"
      out
    }

    stable1 <- get_edges_by_keys(df1, stable_k)
    stable2 <- get_edges_by_keys(df2, stable_k)
    stable  <- if (is.null(stable1)) NULL else {
      merge(stable1, stable2, by = "key", suffixes = c("_bin1", "_bin2"), all = TRUE)
    }

    gained <- get_edges_by_keys(df2, gained_k)
    lost   <- get_edges_by_keys(df1, lost_k)

    gl_data[[pair_name]] <- data.frame(
      Bin_Pair = pair_name,
      Cluster  = cluster_id,
      Type     = c("Stable", "Gained", "Lost"),
      Count    = c(length(stable_k), length(gained_k), length(lost_k)),
      stringsAsFactors = FALSE
    )

    edge_details[[pair_name]] <- list(
      stable = stable,   # has weight_bin1 and weight_bin2
      gained = gained,   # has weight (bin2)
      lost   = lost      # has weight (bin1)
    )
  }

  gl_df <- if (length(gl_data) > 0) do.call(rbind, gl_data) else NULL
  if (is.null(gl_df) || nrow(gl_df) == 0 || all(gl_df$Count == 0)) {
    warning("No edge changes detected between bins (or bins were skipped).")
    return(list(plot = NULL, edge_details = edge_details))
  }

  p_gl <- ggplot2::ggplot(gl_df, ggplot2::aes(x = .data$Bin_Pair, y = .data$Count, fill = .data$Type)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(title = title, x = "Bin Transition", y = "Number of Edges") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  list(plot = p_gl, edge_details = edge_details)
}

