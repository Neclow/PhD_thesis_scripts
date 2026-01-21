# Install and load required packages
library(ggtree)
library(ape)
library(ggplot2)

# Read tree from Newick string
tree <- ape::read.tree(file = "animals.phy.treefile")
tree <- chronos(tree)

# Define mammal and reptile clades
mammals <- c(
  "Human",
  "Seal",
  "Cow",
  "Whale",
  "Mouse",
  "Rat",
  "Platypus",
  "Opossum"
)
reptiles <- c("Turtle", "Crocodile", "Bird", "Sphenodon", "Lizard")

# Find MRCA (Most Recent Common Ancestor) nodes
mammal_node <- MRCA(tree, mammals)
reptile_node <- MRCA(tree, reptiles)

# Create the base plot
p <- ggtree(tree, layout = "rectangular", linewidth = 1) +
  geom_tiplab(
    size = 6,
    color = "#2C3E50",
    fontface = "italic",
    offset = 0.01
  ) +
  geom_treescale(x = 0, y = 0, offset = 0.2, color = "#7F8C8D", fontsize = 6) +
  theme_tree2() +
  xlab("Branch length") +
  theme(
    plot.margin = margin(5, 20, 5, 5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20)
  ) +
  xlim(0, max(node.depth.edgelength(tree)) * 1.3)

# Highlight mammal clade
p <- p +
  geom_hilight(node = mammal_node, fill = "#E8F5E9", alpha = 0.5, extend = 0.01)

# Highlight reptile clade
p <- p +
  geom_hilight(
    node = reptile_node,
    fill = "#FFF3E0",
    alpha = 0.5,
    extend = 0.01
  )

# Add support values (bootstrap values)
p <- p +
  geom_nodelab(
    aes(label = label),
    size = 4,
    color = "#E74C3C",
    hjust = -0.2,
    vjust = -0.5
  )

# Add clade labels
p <- p +
  geom_cladelabel(
    node = mammal_node,
    label = "Mammalia",
    offset = 0.15,
    barsize = 2,
    angle = 90,
    offset.text = 0.01,
    hjust = 0.5,
    fontsize = 5
  ) +
  geom_cladelabel(
    node = reptile_node,
    label = "Reptilia",
    offset = 0.15,
    barsize = 2,
    angle = 90,
    offset.text = 0.01,
    hjust = 0.5,
    fontsize = 5
  )

# Save as PNG
ggsave(
  "phylogenetic_tree.png",
  plot = p,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

# Display the plot
print(p)

cat("\nPhylogenetic tree saved as 'phylogenetic_tree.png'\n")
cat("Mammals highlighted in green, Reptiles highlighted in orange\n")
