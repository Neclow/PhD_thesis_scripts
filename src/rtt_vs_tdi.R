# ============================================================================
# Illustration of Root-to-Tip Length and Taxon Discordance Index (TDI)
# ============================================================================

library(ape)
library(phytools)
library(TreeDist)
library(Cairo)

# ----------------------------------------------------------------------------
# 1. Create two gene trees with a discordant taxon
# ----------------------------------------------------------------------------

# Tree for Region 1: ((A,C),((B,D),E))
tree1 <- read.tree(
  text = "((A:0.15,C:0.12):0.2,((B:0.08,D:0.1):0.12,E:0.18):0.15);"
)

# Tree for Region 2: ((A,(B,C)),(D,E))
tree2 <- read.tree(
  text = "((A:0.1,(B:0.45,C:0.12):0.1):0.2,(D:0.14,E:0.16):0.18);"
)

# The discordant taxon
highlight_taxon <- "B"

# ----------------------------------------------------------------------------
# 2. Calculate Root-to-Tip lengths
# ----------------------------------------------------------------------------

get_rtt <- function(tree) {
  root_node <- length(tree$tip.label) + 1
  distances <- dist.nodes(tree)[seq_along(tree$tip.label), root_node]
  names(distances) <- tree$tip.label
  return(distances)
}

rtt1 <- get_rtt(tree1)
rtt2 <- get_rtt(tree2)

# ----------------------------------------------------------------------------
# 3. Calculate MCI and TDI (Taxon Discordance Index)
# ----------------------------------------------------------------------------

mci <- MutualClusteringInfo(tree1, tree2)

# tdi_b = MCI(tree1, tree2) - MCI(tree1 without B, tree2 without B)
tree1_drop <- drop.tip(tree1, highlight_taxon)
tree2_drop <- drop.tip(tree2, highlight_taxon)
mci_reduced <- MutualClusteringInfo(tree1_drop, tree2_drop)
tdi_b <- mci - mci_reduced

# ----------------------------------------------------------------------------
# 4. Visualize as tanglegram with annotations
# ----------------------------------------------------------------------------

# Ladderize trees for cleaner visualization
tree1 <- ladderize(tree1)
tree2 <- ladderize(tree2)

# Create cophylo object
cophy <- cophylo(
  tree1,
  tree2,
  assoc = cbind(tree1$tip.label, tree1$tip.label),
  rotate = TRUE
)

# Create link colors
link_colors <- ifelse(tree1$tip.label == highlight_taxon, "red", "steelblue")

# Function to plot the tanglegram
plot_tanglegram <- function() {
  par(mar = c(4, 2, 3, 2), family = "URWGothic")
  plot(
    cophy,
    link.type = "curved",
    link.lwd = c(2.5, rep(1.5, length(tree1$tip.label) - 1)),
    link.lty = "solid",
    link.col = link_colors,
    fsize = 1.5,
    scale.bar = c(0.1, 0.1),
    pts = FALSE
  )

  # Add branch length labels (positioned above branches)
  edgelabels.cophylo(
    round(cophy$trees[[1]]$edge.length, 2),
    which = "left",
    cex = 1.2,
    frame = "none",
    adj = c(0.5, -0.3)
  )
  edgelabels.cophylo(
    round(cophy$trees[[2]]$edge.length, 2),
    which = "right",
    cex = 1.2,
    frame = "none",
    adj = c(0.5, -0.3)
  )

  mtext("Region 1", side = 3, adj = 0.15, line = 1, font = 2, cex = 1.3)
  mtext("Region 2", side = 3, adj = 0.85, line = 1, font = 2, cex = 1.3)
}

# Save as PDF
CairoPDF("img/rtt_vs_tdi.pdf", family = "URWGothic", width = 10, height = 6)
plot_tanglegram()
dev.off()

# Save as PNG
CairoPNG("img/rtt_vs_tdi.png", width = 1000, height = 600, res = 100)
plot_tanglegram()
dev.off()

# Print metrics to terminal
cat("\n")
cat("MCI(Tree 1, Tree 2) =", round(mci, 3), "\n")
cat(
  "tdi_b = MCI(Tree 1, Tree 2) - MCI(Tree 1 \\ B, Tree 2 \\ B) =",
  round(mci, 3),
  "-",
  round(mci_reduced, 3),
  "=",
  round(tdi_b, 3),
  "\n"
)
cat("RTT_B (Region 1) =", round(rtt1[highlight_taxon], 3), "\n")
cat("RTT_B (Region 2) =", round(rtt2[highlight_taxon], 3), "\n")
cat("\nDone\n")
