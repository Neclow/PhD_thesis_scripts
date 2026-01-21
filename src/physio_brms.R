library(brms)
library(dplyr)
library(janitor)
library(magrittr)
library(readxl)
library(stringr)
library(tidybayes)
library(tidyr)

cond_effects_plot <- function(
  model,
  data,
  y,
  ylabel,
  filename,
  LOD = NULL,
  log_scale = FALSE
) {
  # 1. Extract Data
  ce_list <- conditional_effects(model, effects = "Box")

  plot_data <- ce_list$Box

  # Base plot
  p <- ggplot(plot_data, aes(x = Box, y = estimate__, color = Box))

  # LOD shading
  if (!is.null(LOD)) {
    data$censored <- data[[y]] <= LOD
    p <- p +
      annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = 0,
        ymax = LOD,
        fill = "gray60",
        alpha = 0.3,
        color = NA
      )

    if (log_scale) {
      y_range <- c(data[[y]], plot_data$estimate__)
      y_breaks <- unique(sort(c(
        LOD,
        10^pretty(log10(y_range[y_range > 0]), n = 5)
      )))
      p <- p +
        geom_hline(yintercept = LOD, linetype = "dashed", color = "gray40") +
        scale_y_log10(breaks = y_breaks)
    } else {
      p <- p +
        geom_hline(yintercept = LOD, linetype = "dashed", color = "gray40") +
        scale_y_continuous(
          breaks = c(LOD, pretty(c(data[[y]], plot_data$estimate__)))
        )
    }
  } else {
    data$censored <- FALSE
    if (log_scale) {
      p <- p + scale_y_log10()
    }
  }

  # Raw data
  p <- p +
    geom_jitter(
      data = data,
      aes_string(y = y, x = "Box", shape = "censored"),
      width = 0.15,
      alpha = 0.8,
      size = 2,
      inherit.aes = FALSE,
      color = ifelse(data$censored, "grey60", "black")
    )

  # Model estimates
  p <- p +
    geom_errorbar(
      aes(ymin = lower__, ymax = upper__),
      width = 0.1,
      linewidth = 1
    ) +
    geom_point(size = 4)

  # Formatting
  p <- p +
    labs(
      y = ylabel,
      x = "Box"
    ) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
    )

  # 3. Save as PDF
  # cairo_pdf is often better than standard "pdf" for handling transparency/fonts on Windows
  ggsave(
    filename = filename,
    plot = p,
    device = cairo_pdf,
    width = 4,
    height = 4,
    units = "in"
  )
}

file <- "data/chicken/physiological_data.xlsx"
model_dir <- "data/chicken/brms/models"

boxes <- c("1" = "Control", "4" = "E. coli")

# Sources
# - https://www.eiken.co.jp/uploads/Eiken_RUO%20LZ-SAA%20Leaflet_A4.pdf
# - https://doi.org/10.1186/s13567-022-01040-1
saa_lod <- 1.0 # mg/L

saa <- suppressMessages(
  read_xlsx(file, sheet = "Acute phase protein") |>
    set_colnames(c(
      "SID",
      "OrderDate",
      "TestName",
      "Units",
      "ResultString",
      "Column1",
      "Column2",
      "7",
      "Box"
    )) |>
    select(-c("7", "Column1", "Column2")) |>
    fill(Box) |>
    mutate(
      saa_imputed = ifelse(
        ResultString == "< 0.0",
        saa_lod,
        as.numeric(ResultString)
      ),
      saa_imputed = ifelse(
        saa_imputed < saa_lod,
        saa_lod,
        saa_imputed
      ),
      cens = ifelse(saa_imputed <= saa_lod, "left", "none"),
      Box = boxes[as.character(as.integer(str_extract(Box, "\\d+")))]
    ) |>
    arrange(Box) |>
    drop_na()
)

saa_model <- brm(
  bf(saa_imputed | cens(cens) ~ Box),
  data = saa,
  family = Gamma(link = "log"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(saa_model, file.path(model_dir, "saa.Rds"))

pdf("img/ppcheck_saa.pdf")
pp_check(saa_model)
dev.off()

saa_yrep <- posterior_predict(saa_model, ndraws = 100)
saa_yrep_clipped <- ifelse(saa_yrep < saa_lod, saa_lod, saa_yrep)
pdf("img/ppcheck2_saa.pdf")
bayesplot::ppc_ecdf_overlay(y = saa$saa_imputed, yrep = saa_yrep_clipped) +
  scale_x_log10()
dev.off()

cond_effects_plot(
  saa_model,
  saa,
  "saa_imputed",
  "Serum amyloid A (mg/L)",
  "img/cond2_saa.pdf",
  LOD = saa_lod,
)

#  Family: gamma
#   Links: mu = log
# Formula: saa_imputed | cens(cens) ~ Box
#    Data: saa (Number of observations: 35)
#   Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
#          total post-warmup draws = 8000

# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -6.42      4.65   -19.27    -1.44 1.00     1132      910
# BoxE.coli     7.89      4.68     2.77    20.80 1.00     1143      893

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# shape     0.48      0.20     0.18     0.96 1.00     1591     2329

# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# unique(fitted(saa_model_gamma))
#        Estimate  Est.Error         Q2.5      Q97.5
# [1,] 0.03597115 0.06853652 4.270799e-09  0.2368085
# [2,] 4.82365757 2.61521798 2.133099e+00 11.5079455

les <- suppressMessages(
  read_xlsx(file, sheet = "Lesion scores", skip = 1) |>
    select(-c(2, 3, last_col())) |>
    select(1, "Total score") |>
    slice(-1) |>
    set_colnames(c("Box", "Score")) |>
    filter(!if_all(everything(), is.na)) |>
    mutate(Box = as.integer(Box), Box = recode(Box, !!!boxes)) |>
    arrange(Box) |>
    filter(!is.na(Box))
)

les$MaxScore <- 38

les$Score_prop <- les$Score / les$MaxScore

les_model <- brm(
  bf(Score_prop ~ Box, zi ~ Box),
  data = les,
  family = zero_inflated_beta(),
  seed = 123
)

saveRDS(les_model, file.path(model_dir, "lesion_score.Rds"))

pdf("img/ppcheck_lesion.pdf")
pp_check(les_model)
dev.off()

les_yrep <- posterior_predict(les_model, ndraws = 100)
pdf("img/ppcheck2_lesion.pdf")
bayesplot::ppc_ecdf_overlay(y = les$Score_prop, yrep = les_yrep)
dev.off()

cond_effects_plot(
  les_model,
  les,
  "Score_prop",
  "Lesion Score (normalized)",
  "img/cond2_lesion.pdf"
)

#  Family: zero_inflated_beta
#   Links: mu = logit; zi = logit
# Formula: Score_prop ~ Box
#          zi ~ Box
#    Data: les (Number of observations: 39)
#   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#          total post-warmup draws = 4000

# Regression Coefficients:
#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept       -2.16      0.50    -3.21    -1.24 1.00     2583     2317
# zi_Intercept     1.13      0.54     0.13     2.21 1.00     5073     3205
# BoxE.coli        2.07      0.54     1.04     3.19 1.00     2811     2510
# zi_BoxE.coli    -2.98      0.87    -4.75    -1.41 1.00     3486     2481

# Further Distributional Parameters:
#     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# phi     4.66      1.34     2.46     7.70 1.00     3181     2254

# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

hypothesis(les_model, "BoxE.coli > 0")
# Hypothesis Tests for class b:
#        Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
# 1 (BoxE.coli) 0     3.52      0.58     2.59     4.46        Inf         1
#   Star
# 1    *
# ---
# 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
# '*': For one-sided hypotheses, the posterior probability exceeds 95%;
# for two-sided hypotheses, the value tested against lies outside the 95%-CI.
# Posterior probabilities of point hypotheses assume equal prior probabilities

# Blend lung tissue 1:1 w/v PBS --2mL
# Plated 90 ul total
# Starting from 10^-1 dilution to 10^-7 dilution
# Diluted volume plated after 10^-1 dilution = 0.09 mL
# 1 CFU in 0.09 mL of 10^-1 dilution = 1 CFU in 0.009 mL of original
# cfu/mL = 1 / 0.009 = 111.11
# Since lung was blended 1:1 w/v, multiply by 2 to get cfu/g
# cfu/g = 222.22
cfu_lod <- 222

cfu <- read_xlsx(file, sheet = "cfu per gram lung", skip = 2) |>
  remove_empty("rows") |>
  select(c(1, last_col(offset = 1L))) |>
  set_colnames(c("Box", "cfu")) |>
  filter(Box != 5) |>
  drop_na(Box) |>
  mutate(
    Box = recode(as.integer(Box), !!!boxes),
    cens = ifelse(cfu == 0, "left", "none"),
    cfu_imputed = ifelse(cfu == 0, cfu_lod, cfu)
  ) |>
  arrange(Box)

cfu_model <- brm(
  bf(cfu_imputed | cens(cens) ~ Box),
  data = cfu,
  family = Gamma(link = "log"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(cfu_model, file.path(model_dir, "cfu.Rds"))

#  Family: lognormal
#   Links: mu = identity; sigma = identity
# Formula: cfu_imputed | cens(cens) ~ Box
#    Data: cfu (Number of observations: 39)
#   Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
#          total post-warmup draws = 8000

# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -3.01      2.59    -8.85     1.23 1.00     2622     2952
# BoxE.coli    15.84      2.94    10.84    22.56 1.00     3115     3495

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     5.25      0.91     3.80     7.34 1.00     2892     3691

# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
pdf("img/ppcheck_cfu.pdf")
pp_check(cfu_model)
dev.off()
cfu_yrep <- posterior_predict(cfu_model, ndraws = 100)
cfu_yrep_clipped <- ifelse(cfu_yrep < cfu_lod, cfu_lod, cfu_yrep)
pdf("img/ppcheck2_cfu.pdf")
bayesplot::ppc_ecdf_overlay(y = cfu$cfu_imputed, yrep = cfu_yrep_clipped) +
  scale_x_log10()
dev.off()

cond_effects_plot(
  cfu_model,
  cfu,
  "cfu_imputed",
  "CFU per gram lung (cfu/g)",
  "img/cond2_cfu.pdf",
  log_scale = TRUE,
  LOD = cfu_lod
)
