devtools::install_github(
  "sbhattlab/phylo2vec",
  subdir = "./r-phylo2vec",
  build = FALSE
)

install.packages(
  c("brms", "TreeSearch", "TBRDist"),
  repos = "https://cloud.r-project.org"
)
