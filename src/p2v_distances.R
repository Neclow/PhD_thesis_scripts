library(ape)
library(phangorn)
library(phylo2vec)
library(progress)
library(TreeDist)
library(TreeSearch)
library(TBRDist)

source("extern/rspr/R/rspR.r")

v_move_loop <- function(n_leaves, steps, n = 1) {
  # Start
  v_start <- sample_vector(n_leaves, ordered = FALSE)
  t_start <- read.tree(text = to_newick(v_start))
  v <- v_start
  t <- t_start

  results <- list(
    hamming = integer(steps),
    rf = numeric(steps),
    rrf = numeric(steps),
    nniu = numeric(steps),
    rspr = numeric(steps),
    tbr = numeric(steps)
  )

  pb <- progress_bar$new(
    total = steps,
    format = "  v moves [:bar] :percent eta: :eta"
  )

  for (s in seq_len(steps)) {
    v <- as.integer(v)
    t <- read.tree(text = to_newick(v))

    if (n == 1) {
      # Randomly pick one index to change
      i <- sample(1:(n_leaves - 2), 1)
      idx <- i + 1
      # Pick a direction to change
      delta <- sample(c(-1, 1), 1)

      # Update vector
      new_val <- v[idx] + delta
      if (new_val >= 0 && new_val <= 2 * i) {
        v[idx] <- new_val
      }
    } else {
      subset_size <- sample(1:n, 1)
      indices <- sample(1:(n_leaves - 2), subset_size, replace = FALSE)

      for (i in indices) {
        idx <- i + 1 # convert to R's 1-based indexing
        # Pick direction
        delta <- sample(c(-1, 1), 1)

        # Update vector
        new_val <- v[idx] + delta
        if (new_val >= 0 && new_val <= 2 * i) {
          v[idx] <- new_val
        }
      }
    }

    v <- as.integer(v)
    # Convert to tree
    t <- read.tree(text = to_newick(v))

    results$hamming[s] <- sum(v != v_start)
    results$rf[s] <- RobinsonFoulds(t, t_start)
    results$rrf[s] <- RF.dist(t, t_start, rooted = TRUE)
    results$nniu[s] <- NNIDist(t, t_start)["best_upper"] # / nniNorm
    # results$uspr[s] <- USPRDist(t, t_start)
    results$rspr[s] <- rspr(t, t_start) # / (n_leaves / 2L)
    results$tbr[s] <- TBRDist(t, t_start, exact = TRUE) # / n_leaves

    # Update progress bar
    pb$tick()

    t_start <- t
    v_start <- v
  }

  as.data.frame(results)
}

tree_move_loop <- function(n_leaves, steps) {
  # Start
  v_start <- sample_vector(n_leaves, ordered = FALSE)
  v <- v_start

  # Output
  results <- list(
    nni = numeric(steps),
    spr = numeric(steps),
    tbr = numeric(steps)
  )

  pb <- progress_bar$new(
    total = steps,
    format = "  tree moves [:bar] :percent eta: :eta"
  )

  for (s in seq_len(steps)) {
    # Update starting points
    v_start <- v
    t_start <- read.tree(text = to_newick(v_start))

    # TBR
    t <- RootedTBR(t_start)
    v <- from_newick(remove_parent_labels(write.tree(t)))
    results$tbr[s] <- sum(v != v_start)

    # SPR
    t <- RootedSPR(t_start)
    v <- from_newick(remove_parent_labels(write.tree(t)))
    results$spr[s] <- sum(v != v_start)

    # NNI
    t <- RootedNNI(t_start)
    v <- from_newick(remove_parent_labels(write.tree(t)))
    results$nni[s] <- sum(v != v_start)

    # Update progress bar
    pb$tick()
  }

  as.data.frame(results)
}

n_leaves <- 100
steps <- 1000
# nniNorm <- TreeDist::NNIDiameter(n_leaves)
res_v <- v_move_loop(n_leaves, steps, 1)
write.csv(res_v, "outputs/compare_metrics_from_v_moves_0119.csv")

res_tree <- tree_move_loop(n_leaves, steps)
write.csv(res_tree, "outputs/compare_metrics_from_tree_moves_0119.csv")

n <- 10
res_multi <- v_move_loop(n_leaves, steps, n)
write.csv(
  res_multi,
  paste0("outputs/compare_metrics_from_multi_moves_", n, "_0119.csv")
)
