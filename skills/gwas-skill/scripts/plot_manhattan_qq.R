#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
})

options(scipen = 9999)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && (is.na(x) || identical(x, ""))) return(y)
  x
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    trait = "GWAS trait",
    gwas_table = "",
    chrom_sizes = "",
    out_prefix = "gwas_plot",
    point_size = "0.8"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Bad argument: ", key)
    name <- gsub("-", "_", substring(key, 3))
    if (!(name %in% names(defaults))) stop("Unknown argument: ", key)
    defaults[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  defaults$point_size <- as.numeric(defaults$point_size)
  defaults
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

standardize_gwas <- function(df) {
  name_map <- c(
    snp_id = "snp_id", SNP = "snp_id",
    chrom = "chrom", CHR = "chrom",
    bp = "bp", BP = "bp",
    p_value = "p_value", P = "p_value"
  )
  for (src in names(name_map)) {
    if (src %in% colnames(df)) colnames(df)[colnames(df) == src] <- name_map[[src]]
  }
  required <- c("snp_id", "chrom", "bp", "p_value")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("GWAS table is missing required columns: ", paste(missing, collapse = ", "))
  df$chrom <- as.character(df$chrom)
  df$bp <- as.numeric(df$bp)
  df$p_value <- as.numeric(df$p_value)
  min_nonzero_p <- min(df$p_value[df$p_value > 0], na.rm = TRUE)
  if (!is.finite(min_nonzero_p)) min_nonzero_p <- 1e-300
  df$p_value[df$p_value == 0] <- min_nonzero_p
  df <- df[!is.na(df$p_value) & !is.na(df$bp) & !is.na(df$chrom), , drop = FALSE]
  df$mlog10_p <- -log10(pmax(df$p_value, .Machine$double.xmin))
  df
}

infer_chrom_sizes <- function(df) {
  out <- aggregate(bp ~ chrom, data = df, FUN = max)
  colnames(out) <- c("chrom", "length_bp")
  out
}

load_chrom_sizes <- function(path, df) {
  if (is.null(path) || identical(path, "") || !file.exists(path)) {
    chrom_sizes <- infer_chrom_sizes(df)
  } else {
    chrom_sizes <- read_tsv(path)
    required <- c("chrom", "length_bp")
    missing <- setdiff(required, colnames(chrom_sizes))
    if (length(missing) > 0) stop("chrom_sizes table is missing required columns: ", paste(missing, collapse = ", "))
  }
  chrom_sizes$chrom <- as.character(chrom_sizes$chrom)
  chrom_sizes$length_bp <- as.numeric(chrom_sizes$length_bp)
  chrom_sizes <- chrom_sizes[order(chrom_sizes$chrom), , drop = FALSE]
  chrom_sizes
}

build_manhattan <- function(df, chrom_sizes) {
  chrom_sizes$offset <- c(0, cumsum(chrom_sizes$length_bp)[-nrow(chrom_sizes)])
  chrom_sizes$center <- chrom_sizes$offset + chrom_sizes$length_bp / 2
  df <- merge(df, chrom_sizes[, c("chrom", "offset")], by = "chrom", all.x = TRUE, sort = FALSE)
  df$x <- df$bp + df$offset
  df$chr_group <- as.integer(as.factor(df$chrom)) %% 2
  df
}

build_qq <- function(df) {
  p <- sort(df$p_value)
  exp <- -log10((seq_along(p) - 0.5) / length(p))
  obs <- -log10(p)
  data.frame(expected = exp, observed = obs)
}

main <- function() {
  args <- parse_args()
  gwas <- standardize_gwas(read_tsv(args$gwas_table))
  chrom_sizes <- load_chrom_sizes(args$chrom_sizes, gwas)
  manhattan <- build_manhattan(gwas, chrom_sizes)
  qq <- build_qq(gwas)

  greyLine <- -log10(1 / nrow(gwas))
  redLine  <- -log10(0.05 / nrow(gwas))

  write.table(gwas[, c("snp_id", "chrom", "bp", "p_value", "mlog10_p")], paste0(args$out_prefix, ".manhattan_points.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(qq, paste0(args$out_prefix, ".qq_points.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  manhattan_plot <- ggplot(manhattan, aes(x = x, y = mlog10_p, color = factor(chr_group))) +
    geom_point(size = args$point_size, alpha = 0.85) +
    scale_color_manual(values = c("#4C78A8", "#F58518")) +
    geom_hline(yintercept = greyLine, color = "grey50", linetype = "dashed") +
    geom_hline(yintercept = redLine, color = "red3", linetype = "dashed") +
    scale_x_continuous(
      breaks = chrom_sizes$center,
      labels = chrom_sizes$chrom,
      expand = c(0.01, 0.01)
    ) +
    labs(title = args$trait, x = "Chromosome", y = expression(-log[10](P))) +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

  qq_plot <- ggplot(qq, aes(x = expected, y = observed)) +
    geom_point(size = 0.9, color = "#4C78A8", alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red3") +
    labs(title = paste0(args$trait, " QQ plot"), x = "Expected -log10(P)", y = "Observed -log10(P)") +
    theme_cowplot(font_size = 12)

  combined <- plot_grid(manhattan_plot, qq_plot, ncol = 1, rel_heights = c(1.6, 1))

  ggsave(paste0(args$out_prefix, ".png"), combined, width = 12, height = 8, dpi = 180)
  ggsave(paste0(args$out_prefix, ".pdf"), combined, width = 12, height = 8)
}

main()
