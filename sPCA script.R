# Spatial Principal Component Analysis (sPCA) Workflow for DArT SNP Data

# This script performs a spatial Principal Component Analysis (sPCA) using individual-level SNP genotype data stored in a DArT-style `dms` object.
# The script is designed to be: 1. reproducible; 2. easy to share with collaborators; 3. clearly annotated for publication-style methods documentation.

# Overview of workflow
# 1. Read setup variables and the filtered `dms` object.
# 2. Extract sample metadata and coordinates.
# 3. Filter samples/loci for sPCA suitability.
# 4. Convert the genotype matrix to a genlight object, then to genind.
# 5. Build a spatial connection network.
# 6. Run sPCA.
# 7. Test for global and local spatial genetic structure.
# 8. Save tables and figures for interpretation.

# This script assumes the `dms` object contains:
#   - dms$gt              : genotype matrix coded as 0 / 1 / 2 / NA
#   - dms$sample_names    : vector of sample names
#   - dms$meta$analyses   : metadata table including latitude and longitude

# The `adegenet::spca()` function requires a `genind` or `genpop` object. Because the filtered DArT genotype data are already stored as a dosage matrix, 
# the script converts:   dms$gt  ->  genlight  ->  genind  ->  spca()

# Positive sPCA eigenvalues indicate global structure (broad-scale spatial autocorrelation), 
# whereas negative eigenvalues indicate local structure (local genetic dissimilarity among neighbours).

# A) Load required packages #####
library(openxlsx)
library(adegenet)
library(ade4)
library(spdep)
library(dartR)
library(ggplot2)
library(dplyr)
library(ozmaps)


# B) User-defined inputs and setup variables #####
#Values are read from `0_setup_variables.xlsx` where possible so that the same script can be reused across species or datasets with minimal editing.
setwd("C:/Users/dimon/rrspecies/Eucalyptus_melliodora/")

setup_variables <- read.xlsx("0_setup_variables.xlsx", colNames = TRUE)
maindir <- setup_variables[1, 2]
species <- setup_variables[2, 2]
dataset <- setup_variables[3, 2]
raw_meta_path <- setup_variables[4, 2]
species_col_name <- setup_variables[5, 2]
site_col_name <- setup_variables[6, 2]
setwd(maindir)


# C) Read the filtered DMS object #####
# Script assumes main QC and filtering have already been completed, and that the saved `dms.RData` object is the version intended for downstream analyses.
dms_path <- paste0(species, "/outputs_", site_col_name, "_", species_col_name, "/r_files/dms.RData")
dms <- readRDS(dms_path)

# D) Create output directory for sPCA results #####
sPCA_dir <- paste0(species, "/outputs_", site_col_name, "_", species_col_name, "/sPCA/")
if (!dir.exists(sPCA_dir)) dir.create(sPCA_dir, recursive = TRUE)


# E) sPCA analysis settings #####
# spca_network_type: Defines how neighbours are connected in space.
# 1 = Delaunay; 2 = Gabriel; 3 = Relative neighbour; 5 = Distance-based graph; 6 = K-nearest neighbours; 7 = Minimum spanning network

# Recommended choice for datasets with duplicate coordinates: Use type = 5 rather than types 1-4. Types 1-4 require unique sample coordinates 
# (you could slightly jitter the duplciate coordinates if youw want)

# spca_d1 / spca_d2:
# Lower and upper distance thresholds for type = 5. If spca_d2 is NULL, it is estimated from the non-zero pairwise distances.

# spca_nfposi / spca_nfnega: Number of positive and negative axes to retain. Positive axes usually capture broad geographic structure.

# spca_nperm: Number of permutations for global/local Monte Carlo tests.

# spca_use_lagged: If TRUE, use lagged scores for some plots. Lagged scores summarise  average spatial signal around a point --> make broad patterns clearer.

spca_network_type <- 5
spca_k <- 10
spca_d1 <- 0
spca_d2 <- NULL
spca_nfposi <- 2
spca_nfnega <- 0
spca_nperm <- 999
spca_use_lagged <- TRUE
spca_point_size <- 2
max_snps_for_spca <- 1000


# F) Extract metadata and retain samples with valid coordinates #####
# sPCA requires geographic coordinates for each retained sample. Samples lacking latitude or longitude are excluded.
spca_meta <- dms$meta$analyses %>% as.data.frame() %>% mutate(sample = dms$sample_names, lat = as.numeric(lat), long = as.numeric(long))
spca_keep <- which(!is.na(spca_meta$lat) & !is.na(spca_meta$long))
if (length(spca_keep) < 3) stop("Not enough samples with valid coordinates for sPCA.")
spca_meta2 <- spca_meta[spca_keep, , drop = FALSE]
spca_gt <- dms$gt[spca_keep, , drop = FALSE]


# G) Remove loci that are completely missing after coordinate filtering #####
# This is a minimum filtering step for the sPCA matrix. Note that loci with partial missingness are retained at this stage; 
# a stricter matrix is created later for the permutation tests.
keep_loci <- colSums(!is.na(spca_gt)) > 0
spca_gt <- spca_gt[, keep_loci, drop = FALSE]
if (ncol(spca_gt) < 2) stop("Not enough loci remaining after filtering for sPCA.")


# H) Optional SNP subsetting to improve runtime #####
# Converting large SNP dosage matrices into `genind` format can produce very large allele tables, which can greatly slow down sPCA.
# For exploratory analyses, subsampling SNPs can reduce runtime substantially.
set.seed(1)
if (ncol(spca_gt) > max_snps_for_spca) {
  keep_snps <- sample(seq_len(ncol(spca_gt)), max_snps_for_spca)
  spca_gt <- spca_gt[, keep_snps, drop = FALSE]
  cat("Subsetting to", ncol(spca_gt), "SNPs for sPCA\n")
} else {
  cat("Using all", ncol(spca_gt), "SNPs for sPCA\n")
}


# I) Convert genotype matrix to genlight, then genind #####
# `adegenet::spca()` operates on `genind` or `genpop` objects rather than `genlight`, so the `gl2gi()` conversion is used as a bridge.
cat("Building genlight object...\n")
spca_gl <- new("genlight", as.matrix(spca_gt))
indNames(spca_gl) <- spca_meta2$sample
locNames(spca_gl) <- colnames(spca_gt)
ploidy(spca_gl) <- 2
spca_gl@other$ind.metrics <- spca_meta2
if (site_col_name %in% colnames(spca_meta2)) pop(spca_gl) <- as.factor(spca_meta2[[site_col_name]])

cat("Converting genlight to genind...\n")
spca_gi <- gl2gi(spca_gl)
spca_xy <- as.matrix(spca_meta2[, c("long", "lat")])
colnames(spca_xy) <- c("x", "y")
spca_gi@other$xy <- spca_xy
if (site_col_name %in% colnames(spca_meta2)) pop(spca_gi) <- as.factor(spca_meta2[[site_col_name]])

cat("genind table dimensions:", dim(spca_gi@tab)[1], "x", dim(spca_gi@tab)[2], "\n")
cat("genind table size:", format(object.size(spca_gi@tab), units = "MB"), "\n")


# J) Build the spatial connection network #####
# The spatial network determines which samples are considered neighbours, 
# as the analysis seeks axes that maximise both genetic variance and spatial autocorrelation.
# For distance-based graphs (type = 5): zero distances caused by duplicate coordinates are ignored when estimating the upper threshold `spca_d2`.
if (spca_network_type == 5) {
  if (is.null(spca_d2)) {
    dmat <- as.matrix(dist(spca_xy))
    dvec <- dmat[upper.tri(dmat)]
    dvec <- dvec[dvec > 0]
    spca_d2 <- quantile(dvec, probs = 0.10, na.rm = TRUE)
  }
  spca_cn <- chooseCN(xy = spca_xy, type = 5, d1 = spca_d1, d2 = spca_d2, plot.nb = FALSE, ask = FALSE)
} else if (spca_network_type == 6) {
  spca_cn <- chooseCN(xy = spca_xy, type = 6, k = spca_k, plot.nb = FALSE, ask = FALSE)
} else {
  spca_cn <- chooseCN(xy = spca_xy, type = spca_network_type, plot.nb = FALSE, ask = FALSE)
}


# K) Run the sPCA #####
# Positive eigenvalues = global structure.
# Negative eigenvalues = local structure.
cat("Running sPCA...\n")
spca_res <- spca(obj = spca_gi, cn = spca_cn, scannf = FALSE, nfposi = spca_nfposi, nfnega = spca_nfnega)
print(spca_res)


# L) Save sPCA eigenvalues table #####
spca_eigs <- data.frame(axis = seq_along(spca_res$eig), eigenvalue = spca_res$eig, structure = ifelse(spca_res$eig > 0, "global", "local"))
write.xlsx(spca_eigs, paste0(sPCA_dir, "sPCA_eigenvalues.xlsx"), overwrite = TRUE)


# M) Global and local Monte Carlo tests #####
# The `global.rtest()` and `local.rtest()` functions do not allow missing values or invariant columns in the input matrix.
#   1. remove columns containing any NA values;
#   2. remove columns with zero variance.

spca_tab_complete <- spca_gi@tab[, colSums(is.na(spca_gi@tab)) == 0, drop = FALSE]
spca_var <- apply(spca_tab_complete, 2, var)
spca_tab_complete <- spca_tab_complete[, spca_var > 0, drop = FALSE]
cat("Final test matrix dimensions:", dim(spca_tab_complete)[1], "x", dim(spca_tab_complete)[2], "\n")
spca_global_test <- global.rtest(spca_tab_complete, spca_res$lw, nperm = spca_nperm)
spca_local_test <- local.rtest(spca_tab_complete, spca_res$lw, nperm = spca_nperm)


# N) Save diagnostic plots #####
pdf(paste0(sPCA_dir, "sPCA_eigenvalues.pdf"), width = 9, height = 5)
barplot(spca_res$eig, col = ifelse(spca_res$eig > 0, "grey35", "grey75"), main = "sPCA eigenvalues", ylab = "Eigenvalue")
abline(h = 0, lty = 2)
dev.off()

pdf(paste0(sPCA_dir, "sPCA_screeplot.pdf"), width = 8, height = 6)
screeplot(spca_res)
dev.off()

pdf(paste0(sPCA_dir, "sPCA_global_test_histogram.pdf"), width = 7, height = 5)
plot(spca_global_test)
dev.off()

pdf(paste0(sPCA_dir, "sPCA_local_test_histogram.pdf"), width = 7, height = 5)
plot(spca_local_test)
dev.off()


# O) Extract sample scores and lagged scores  #####
# `spca_res$li` = raw sPCA scores.
# `spca_res$ls` = lagged sPCA scores.

# Raw scores describe each individual’s position on the sPCA axis. 
# Lagged scores summarise the neighbouring spatial signal and can make broad clinal or regional structure easier to interpret visually.

spca_scores <- as.data.frame(spca_res$li)
colnames(spca_scores) <- paste0("sPC", seq_len(ncol(spca_scores)))
spca_lagged <- as.data.frame(spca_res$ls)
colnames(spca_lagged) <- paste0("lag_sPC", seq_len(ncol(spca_lagged)))

spca_scores_df <- data.frame(sample = rownames(spca_scores), spca_scores, spca_lagged, stringsAsFactors = FALSE) %>% left_join(spca_meta2, by = "sample")
write.xlsx(spca_scores_df, paste0(sPCA_dir, "sPCA_scores.xlsx"), overwrite = TRUE)


# P) Save allele loadings for axis 1  #####
# Squared loadings indicate which alleles contribute most strongly to a given sPCA axis.

if (ncol(spca_res$c1) >= 1) {
  spca_load1 <- spca_res$c1[, 1]^2
  names(spca_load1) <- rownames(spca_res$c1)
  pdf(paste0(sPCA_dir, "sPCA_axis1_loadings.pdf"), width = 12, height = 5)
  loadingplot(spca_load1, threshold = quantile(spca_load1, 0.95, na.rm = TRUE), xlab = "Alleles", ylab = "Squared loading", 
              main = "Contribution of alleles to sPCA axis 1")
  dev.off()
  write.xlsx(data.frame(allele = names(spca_load1), loading_sq = as.numeric(spca_load1)) %>% arrange(desc(loading_sq)), 
             paste0(sPCA_dir, "sPCA_axis1_loadings.xlsx"), overwrite = TRUE)
}


# Q) Decide whether to use raw or lagged scores in downstream figures  #####
score_prefix <- if (spca_use_lagged) "lag_sPC" else "sPC"


# R) Helper function for grouped scatterplots  #####
plot_group_spca <- function(df, group_col_name, group_colours, group_shapes, x_axis = "sPC1", y_axis = "sPC2", x_lab = "sPC1", y_lab = "sPC2", 
                            ellipses = FALSE, hulls = FALSE, point_size = 2, show_legend = TRUE) {
  df <- df[!is.na(df[[group_col_name]]) & !is.na(df[[x_axis]]) & !is.na(df[[y_axis]]), , drop = FALSE]
  group_sizes <- table(df[[group_col_name]])
  ordered <- names(sort(group_sizes, decreasing = TRUE))
  df[[group_col_name]] <- factor(df[[group_col_name]], levels = ordered)
  group_colours <- group_colours[ordered]
  group_shapes <- group_shapes[ordered]
  
  hull_df <- function(dat, x, y, grp) {
    dat %>% group_by(.data[[grp]]) %>% filter(n() >= 3) %>% slice(chull(.data[[x]], .data[[y]])) %>% ungroup()
  }
  
  p <- ggplot(df, aes_string(x = x_axis, y = y_axis, colour = group_col_name, fill = group_col_name, shape = group_col_name)) +
    {if (hulls) geom_polygon(data = hull_df(df, x_axis, y_axis, group_col_name), aes_string(x = x_axis, y = y_axis, fill = group_col_name), 
                             alpha = 0.15, colour = NA)} +
    geom_point(size = point_size, stroke = 0.7) +
    {if (ellipses) stat_ellipse(aes_string(colour = group_col_name), level = 0.95, linewidth = 0.5)} +
    xlab(x_lab) + ylab(y_lab) + scale_colour_manual(values = group_colours) + scale_fill_manual(values = group_colours) + 
    scale_shape_manual(values = group_shapes) +
    theme_bw() + theme(panel.grid = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 10, face = "italic"), 
                       legend.position = if (show_legend) "right" else "none")
  return(p)
}


# S) Fallback colour and shape palettes  #####
# If these objects were not created earlier in a larger workflow, generate defaults so the script can run as a standalone file.

if (!exists("sp_colours") || !exists("sp_shapes")) {
  cat("sp_colours/sp_shapes not found — generating defaults\n")
  uq_species <- unique(dms$meta$analyses[, species_col_name])
  uq_species <- uq_species[order(uq_species)]
  sp_colours <- setNames(grDevices::rainbow(length(uq_species)), uq_species)
  sp_shapes <- setNames(rep(c(21:25, 3, 4, 8, 12), length.out = length(uq_species)), uq_species)
}

if (!exists("site_colours") || !exists("site_shapes")) {
  cat("site_colours/site_shapes not found — generating defaults\n")
  uq_sites <- unique(dms$meta$analyses[, site_col_name])
  uq_sites <- uq_sites[order(uq_sites)]
  site_colours <- setNames(grDevices::rainbow(length(uq_sites)), uq_sites)
  site_shapes <- setNames(rep(c(21:25, 3, 4, 8, 12), length.out = length(uq_sites)), uq_sites)
}

cat("Species levels:", length(sp_colours), "\n")
cat("Site levels:", length(site_colours), "\n")


# T) Plot sPCA scores by species #####
spca_species_plot <- plot_group_spca(df = spca_scores_df, group_col_name = species_col_name, group_colours = sp_colours, group_shapes = sp_shapes, 
                                     x_axis = paste0(score_prefix, "1"), y_axis = paste0(score_prefix, "2"), 
                                     x_lab = paste0(score_prefix, "1"), y_lab = paste0(score_prefix, "2"), ellipses = FALSE, hulls = FALSE, 
                                     point_size = spca_point_size, show_legend = TRUE)

ggsave(paste0(sPCA_dir, "sPCA_species_", score_prefix, "12.png"), spca_species_plot, width = 16, height = 12, dpi = 600, units = "cm", bg = "white")


# U) Plot sPCA scores by site #####
spca_site_plot <- plot_group_spca(df = spca_scores_df, group_col_name = site_col_name, group_colours = site_colours, group_shapes = site_shapes, 
                                  x_axis = paste0(score_prefix, "1"), y_axis = paste0(score_prefix, "2"), x_lab = paste0(score_prefix, "1"), 
                                  y_lab = paste0(score_prefix, "2"), ellipses = TRUE, hulls = TRUE, point_size = spca_point_size, show_legend = TRUE)

ggsave(paste0(sPCA_dir, "sPCA_sites_", score_prefix, "12.png"), spca_site_plot, width = 20, height = 14, dpi = 600, units = "cm", bg = "white")


# V) Plot latitude-coloured sPCA scatterplot #####
# Latitude is coloured with a blue–white–red gradient:
#   - more negative latitude values (further south) = blue
#   - higher latitude values (further north)        = red
spca_lat_plot <- ggplot(spca_scores_df, aes_string(x = paste0(score_prefix, "1"), y = paste0(score_prefix, "2"), colour = "lat")) +
  geom_point(size = 2) +
  scale_colour_gradient(low = "blue", high = "red", name = "Latitude", na.value = "grey30") +
  xlab(paste0(score_prefix, "1")) + ylab(paste0(score_prefix, "2")) + theme_bw() + theme(panel.grid = element_blank())

ggsave(paste0(sPCA_dir, "sPCA_latitude_", score_prefix, "12.png"), spca_lat_plot, width = 14, height = 12, dpi = 600, units = "cm", bg = "white")


# W) Helper function for mapping sPCA axes in geographic space #####
plot_spca_map <- function(df, axis_col, axis_lab) {
  x_pad <- diff(range(df$long, na.rm = TRUE)) * 0.08
  y_pad <- diff(range(df$lat, na.rm = TRUE)) * 0.08
  if (is.na(x_pad) || x_pad == 0) x_pad <- 0.2
  if (is.na(y_pad) || y_pad == 0) y_pad <- 0.2
  
  xlims <- c(min(df$long, na.rm = TRUE) - x_pad, max(df$long, na.rm = TRUE) + x_pad)
  ylims <- c(min(df$lat, na.rm = TRUE) - y_pad, max(df$lat, na.rm = TRUE) + y_pad)
  
  ggplot() +
    geom_sf(data = ozmaps::abs_ste, fill = NA, colour = "grey35", linewidth = 0.4) +
    geom_point(data = df, aes(x = long, y = lat, fill = .data[[axis_col]]), shape = 21, colour = "black", size = 3, stroke = 0.3) +
    scale_fill_gradient(low = "red", high = "blue", name = axis_lab, na.value = "grey80") +
    coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) +
    theme_bw() + labs(x = "Longitude", y = "Latitude") + theme(panel.grid = element_blank())
}


# X) Map the first two sPCA axes #####
spca_map1 <- plot_spca_map(spca_scores_df, paste0(score_prefix, "1"), paste0(score_prefix, "1"))
ggsave(paste0(sPCA_dir, "sPCA_map_", score_prefix, "1.png"), spca_map1, width = 16, height = 12, dpi = 600, units = "cm", bg = "white")

if (paste0(score_prefix, "2") %in% colnames(spca_scores_df)) {
  spca_map2 <- plot_spca_map(spca_scores_df, paste0(score_prefix, "2"), paste0(score_prefix, "2"))
  ggsave(paste0(sPCA_dir, "sPCA_map_", score_prefix, "2.png"), spca_map2, width = 16, height = 12, dpi = 600, units = "cm", bg = "white")
}


# Main outputs written to `sPCA_dir`:
#   - sPCA_eigenvalues.xlsx
#   - sPCA_scores.xlsx
#   - sPCA_axis1_loadings.xlsx
#   - sPCA diagnostic PDFs
#   - grouped scatterplots and geographic maps

# After running:
#   1. Examine the eigenvalue distribution.
#   2. Confirm whether global and/or local tests are significant.
#   3. Inspect whether lagged and raw scores tell a consistent story.
#   4. Re-run with alternative SNP subsets if using exploratory subsampling.
