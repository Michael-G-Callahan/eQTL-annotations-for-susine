#!/usr/bin/env Rscript

required_packages <- c("arrow", "dplyr", "readr", "susieR", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
  library(tibble)
})

parse_cli_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    key <- sub("^--", "", key)
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[parts[[1]]]] <- paste(parts[-1], collapse = "=")
      i <- i + 1L
    } else {
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
        out[[key]] <- TRUE
        i <- i + 1L
      } else {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      }
    }
  }
  out
}

scalar_or_na <- function(x, default = NA_real_) {
  if (length(x) == 0 || is.null(x) || all(is.na(x))) default else x[[1]]
}

normalize_prob_vec <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  s <- sum(x)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  x / s
}

prob_entropy <- function(prob) {
  p <- normalize_prob_vec(prob)
  if (!any(p > 0)) return(NA_real_)
  -sum(p[p > 0] * log(p[p > 0]))
}

prob_k_eff <- function(prob) {
  h <- prob_entropy(prob)
  if (!is.finite(h)) return(NA_real_)
  exp(h)
}

get_credible_set <- function(alpha, rho = 0.95) {
  a <- as.numeric(alpha)
  a[!is.finite(a)] <- 0
  a[a < 0] <- 0
  if (sum(a) <= 0) return(integer(0))
  o <- order(-a, seq_along(a))
  k <- which(cumsum(a[o]) >= rho)[1]
  if (is.na(k)) integer(0) else o[seq_len(k)]
}

prob_entropy_core <- function(prob, rho = 0.95) {
  p <- normalize_prob_vec(prob)
  idx <- get_credible_set(p, rho = rho)
  if (!length(idx)) return(NA_real_)
  prob_entropy(p[idx])
}

prob_k_eff_core <- function(prob, rho = 0.95) {
  h <- prob_entropy_core(prob, rho = rho)
  if (!is.finite(h)) return(NA_real_)
  exp(h)
}

cs_purity_min_abs <- function(R, idx) {
  idx <- as.integer(idx)
  if (!length(idx)) return(NA_real_)
  if (length(idx) == 1L) return(1)
  sub_r <- abs(R[idx, idx, drop = FALSE])
  vals <- sub_r[upper.tri(sub_r)]
  if (!length(vals)) return(NA_real_)
  min(vals, na.rm = TRUE)
}

safe_quantile <- function(x, prob) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, probs = prob, na.rm = TRUE))
}

z_score_metrics <- function(z, top_k = 10L) {
  z <- as.numeric(z)
  z2 <- z^2
  z4 <- z2^2
  ord <- order(z2, decreasing = TRUE, na.last = NA)
  top_idx <- head(ord, max(1L, as.integer(top_k)))
  rest_idx <- setdiff(ord, top_idx)
  top_sum <- sum(z2[top_idx], na.rm = TRUE)
  rest_sum <- sum(z2[rest_idx], na.rm = TRUE)
  tibble(
    z_max_abs = max(abs(z), na.rm = TRUE),
    z_topk_ratio = if (rest_sum > 0) top_sum / rest_sum else NA_real_,
    z_count_abs_gt_3 = sum(abs(z) > 3, na.rm = TRUE),
    z_eff_signals = if (sum(z4, na.rm = TRUE) > 0) {
      (sum(z2, na.rm = TRUE)^2) / sum(z4, na.rm = TRUE)
    } else {
      NA_real_
    }
  )
}

build_dense_ld <- function(ld_long_df, p, locus_id) {
  required <- c("snp_index_1", "snp_index_2", "r")
  missing <- setdiff(required, names(ld_long_df))
  if (length(missing)) {
    stop("LD long table missing columns: ", paste(missing, collapse = ", "))
  }
  idx1 <- as.integer(ld_long_df$snp_index_1)
  idx2 <- as.integer(ld_long_df$snp_index_2)
  r <- as.numeric(ld_long_df$r)
  if (any(!is.finite(r))) stop("LD long table contains non-finite r values.")
  if (any(idx1 < 0 | idx1 >= p | idx2 < 0 | idx2 >= p, na.rm = TRUE)) {
    stop("LD indices out of range for ", locus_id)
  }
  expected_pairs <- p * (p - 1) / 2
  if (nrow(ld_long_df) != expected_pairs) {
    stop(
      "LD long table has ", nrow(ld_long_df), " rows; expected ",
      expected_pairs, " upper-triangle pairs."
    )
  }
  pair_key <- paste(pmin(idx1, idx2), pmax(idx1, idx2), sep = ":")
  if (any(duplicated(pair_key))) stop("Duplicate LD pairs detected.")

  R <- matrix(0, nrow = p, ncol = p)
  diag(R) <- 1
  R[cbind(idx1 + 1L, idx2 + 1L)] <- r
  R[cbind(idx2 + 1L, idx1 + 1L)] <- r
  if (max(abs(R - t(R))) > 1e-8) stop("Reconstructed LD matrix is not symmetric.")
  if (max(abs(diag(R) - 1)) > 1e-8) stop("Reconstructed LD diagonal is not all ones.")
  R
}

load_locus_inputs <- function(row, max_variants) {
  for (field in c("master_path", "order_path", "ld_long_path")) {
    path <- as.character(row[[field]][[1]])
    if (!file.exists(path)) stop("Missing required input: ", path)
  }

  master_df <- readr::read_csv(row$master_path[[1]], show_col_types = FALSE)
  if ("ld_included" %in% names(master_df)) {
    master_df <- master_df %>% dplyr::filter(.data$ld_included)
  }
  master_required <- c("variant_id", "z_score", "sample_size")
  missing_master <- setdiff(master_required, names(master_df))
  if (length(missing_master)) {
    stop("Master variant file missing columns: ", paste(missing_master, collapse = ", "))
  }
  optional_master_cols <- intersect(c("chrom", "pos", "ref", "alt"), names(master_df))
  master_df <- master_df %>%
    dplyr::transmute(
      variant_id = as.character(.data$variant_id),
      z_score = as.numeric(.data$z_score),
      sample_size = as.numeric(.data$sample_size),
      dplyr::across(dplyr::all_of(optional_master_cols))
    )

  order_df <- readr::read_tsv(row$order_path[[1]], show_col_types = FALSE)
  order_required <- c("id", "index")
  missing_order <- setdiff(order_required, names(order_df))
  if (length(missing_order)) {
    stop("LD order file missing columns: ", paste(missing_order, collapse = ", "))
  }
  order_df <- order_df %>%
    dplyr::mutate(index = as.integer(.data$index), id = as.character(.data$id)) %>%
    dplyr::arrange(.data$index)

  aligned_df <- order_df %>%
    dplyr::left_join(master_df, by = c("id" = "variant_id")) %>%
    dplyr::arrange(.data$index)
  p <- nrow(aligned_df)
  if (p < 1L) stop("No aligned variants.")
  if (is.finite(max_variants) && p > max_variants) {
    stop("too_large_for_baseline_susie: p=", p, " exceeds cap ", max_variants)
  }
  if (!all(aligned_df$index == seq_len(p) - 1L)) {
    stop("LD order index is not contiguous from 0 to p-1.")
  }
  if (any(is.na(aligned_df$z_score))) stop("Aligned z_score vector contains NA values.")
  if (any(!is.finite(aligned_df$z_score))) stop("Aligned z_score vector contains non-finite values.")
  if (any(is.na(aligned_df$sample_size)) || any(!is.finite(aligned_df$sample_size))) {
    stop("Aligned sample_size vector contains missing or non-finite values.")
  }
  if (any(aligned_df$sample_size <= 2, na.rm = TRUE)) {
    stop("Aligned sample_size vector contains values <= 2.")
  }

  ld_long_df <- arrow::read_parquet(row$ld_long_path[[1]]) %>% tibble::as_tibble()
  R <- build_dense_ld(ld_long_df, p = p, locus_id = row$locus_id[[1]])

  list(
    aligned_df = aligned_df,
    R = R,
    z = aligned_df$z_score,
    n_rss = round(mean(aligned_df$sample_size, na.rm = TRUE))
  )
}

empty_locus_metrics <- function() {
  tibble(
    locus_id = character(),
    gene_name = character(),
    n_variants = integer(),
    n_rss = integer(),
    status = character(),
    error_message = character(),
    elapsed_seconds = double(),
    elbo_last = double(),
    n_iter = integer(),
    max_pip = double(),
    pip_sum = double(),
    pip_entropy = double(),
    pip_k_eff = double(),
    pip_top1_share = double(),
    pip_top2_gap = double(),
    n_pip_gt_0_01 = integer(),
    n_pip_gt_0_05 = integer(),
    n_pip_gt_0_1 = integer(),
    n_pip_gt_0_5 = integer(),
    n_pip_gt_0_9 = integer(),
    z_max_abs = double(),
    z_topk_ratio = double(),
    z_count_abs_gt_3 = integer(),
    z_eff_signals = double(),
    n_active_effects = integer(),
    min_active_alpha_max = double(),
    median_active_alpha_max = double(),
    max_active_alpha_max = double(),
    min_active_alpha_k_eff_core95 = double(),
    median_active_alpha_k_eff_core95 = double(),
    max_active_alpha_k_eff_core95 = double()
  )
}

empty_effect_metrics <- function() {
  tibble(
    locus_id = character(),
    gene_name = character(),
    effect_l = integer(),
    alpha_max = double(),
    alpha_entropy = double(),
    alpha_entropy_core95 = double(),
    alpha_k_eff = double(),
    alpha_k_eff_core95 = double(),
    cs_size_raw = integer(),
    cs_purity = double()
  )
}

empty_top_variants <- function() {
  tibble(
    locus_id = character(),
    gene_name = character(),
    rank_pip = integer(),
    variant_id = character(),
    ld_matrix_index = integer(),
    chrom = character(),
    pos = integer(),
    ref = character(),
    alt = character(),
    z_score = double(),
    pip = double(),
    posterior_mean = double()
  )
}

write_task_outputs <- function(output_dir, locus_metrics, effect_metrics, top_variants) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(locus_metrics, file.path(output_dir, "locus_metrics.csv"))
  readr::write_csv(effect_metrics, file.path(output_dir, "effect_metrics.csv"))
  readr::write_csv(top_variants, file.path(output_dir, "top_variants.csv"))
}

summarise_effects <- function(alpha, R, locus_id, gene_name) {
  if (is.null(alpha) || !length(alpha)) return(empty_effect_metrics())
  alpha <- as.matrix(alpha)
  rows <- lapply(seq_len(nrow(alpha)), function(l) {
    a <- as.numeric(alpha[l, ])
    cs <- get_credible_set(normalize_prob_vec(a), rho = 0.95)
    tibble(
      locus_id = locus_id,
      gene_name = gene_name,
      effect_l = as.integer(l),
      alpha_max = max(a, na.rm = TRUE),
      alpha_entropy = prob_entropy(a),
      alpha_entropy_core95 = prob_entropy_core(a, rho = 0.95),
      alpha_k_eff = prob_k_eff(a),
      alpha_k_eff_core95 = prob_k_eff_core(a, rho = 0.95),
      cs_size_raw = length(cs),
      cs_purity = cs_purity_min_abs(R, cs)
    )
  })
  dplyr::bind_rows(rows)
}

summarise_locus <- function(row, inputs, fit, effect_metrics, elapsed_seconds,
                            active_alpha_max_threshold) {
  pip <- as.numeric(fit$pip)
  pip_sum <- sum(pip, na.rm = TRUE)
  pip_sorted <- sort(pip, decreasing = TRUE, na.last = NA)
  z_metrics <- z_score_metrics(inputs$z)
  active_effects <- effect_metrics %>%
    dplyr::filter(.data$alpha_max >= active_alpha_max_threshold)

  tibble(
    locus_id = row$locus_id[[1]],
    gene_name = row$gene_name[[1]],
    n_variants = length(pip),
    n_rss = inputs$n_rss,
    status = "completed",
    error_message = NA_character_,
    elapsed_seconds = elapsed_seconds,
    elbo_last = tail(fit$elbo, 1),
    n_iter = length(fit$elbo),
    max_pip = max(pip, na.rm = TRUE),
    pip_sum = pip_sum,
    pip_entropy = prob_entropy(pip),
    pip_k_eff = prob_k_eff(pip),
    pip_top1_share = if (pip_sum > 0) max(pip, na.rm = TRUE) / pip_sum else NA_real_,
    pip_top2_gap = if (length(pip_sorted) >= 2L) pip_sorted[[1]] - pip_sorted[[2]] else NA_real_,
    n_pip_gt_0_01 = sum(pip > 0.01, na.rm = TRUE),
    n_pip_gt_0_05 = sum(pip > 0.05, na.rm = TRUE),
    n_pip_gt_0_1 = sum(pip > 0.1, na.rm = TRUE),
    n_pip_gt_0_5 = sum(pip > 0.5, na.rm = TRUE),
    n_pip_gt_0_9 = sum(pip > 0.9, na.rm = TRUE),
    z_max_abs = z_metrics$z_max_abs[[1]],
    z_topk_ratio = z_metrics$z_topk_ratio[[1]],
    z_count_abs_gt_3 = z_metrics$z_count_abs_gt_3[[1]],
    z_eff_signals = z_metrics$z_eff_signals[[1]],
    n_active_effects = nrow(active_effects),
    min_active_alpha_max = if (nrow(active_effects)) min(active_effects$alpha_max, na.rm = TRUE) else NA_real_,
    median_active_alpha_max = if (nrow(active_effects)) median(active_effects$alpha_max, na.rm = TRUE) else NA_real_,
    max_active_alpha_max = if (nrow(active_effects)) max(active_effects$alpha_max, na.rm = TRUE) else NA_real_,
    min_active_alpha_k_eff_core95 = if (nrow(active_effects)) min(active_effects$alpha_k_eff_core95, na.rm = TRUE) else NA_real_,
    median_active_alpha_k_eff_core95 = if (nrow(active_effects)) median(active_effects$alpha_k_eff_core95, na.rm = TRUE) else NA_real_,
    max_active_alpha_k_eff_core95 = if (nrow(active_effects)) max(active_effects$alpha_k_eff_core95, na.rm = TRUE) else NA_real_
  )
}

build_top_variants <- function(row, inputs, fit, top_n = 20L) {
  pip <- as.numeric(fit$pip)
  posterior_mean <- if (!is.null(fit$alpha) && !is.null(fit$mu)) {
    colSums(as.matrix(fit$alpha) * as.matrix(fit$mu))
  } else {
    rep(NA_real_, length(pip))
  }
  df <- inputs$aligned_df %>%
    dplyr::mutate(
      pip = pip,
      posterior_mean = posterior_mean,
      rank_pip = dplyr::row_number(dplyr::desc(.data$pip))
    ) %>%
    dplyr::arrange(dplyr::desc(.data$pip), .data$index) %>%
    dplyr::slice_head(n = top_n)

  tibble(
    locus_id = row$locus_id[[1]],
    gene_name = row$gene_name[[1]],
    rank_pip = seq_len(nrow(df)),
    variant_id = as.character(df$id),
    ld_matrix_index = as.integer(df$index),
    chrom = as.character(if ("chrom" %in% names(df)) df$chrom else NA_character_),
    pos = as.integer(if ("pos" %in% names(df)) df$pos else NA_integer_),
    ref = as.character(if ("ref" %in% names(df)) df$ref else NA_character_),
    alt = as.character(if ("alt" %in% names(df)) df$alt else NA_character_),
    z_score = as.numeric(df$z_score),
    pip = as.numeric(df$pip),
    posterior_mean = as.numeric(df$posterior_mean)
  )
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  for (required in c("task-manifest", "task-id")) {
    if (is.null(args[[required]])) stop("Missing required argument --", required)
  }
  task_manifest_path <- args[["task-manifest"]]
  task_id <- as.integer(args[["task-id"]])
  susie_l <- as.integer(args[["L"]] %||% 10L)
  max_variants <- as.numeric(args[["max-variants"]] %||% Inf)
  active_alpha_max_threshold <- as.numeric(args[["active-alpha-max-threshold"]] %||% 0.01)

  manifest <- readr::read_csv(task_manifest_path, show_col_types = FALSE)
  if (!is.finite(task_id) || task_id < 1L || task_id > nrow(manifest)) {
    stop("task-id must be in 1..", nrow(manifest))
  }
  row <- manifest[task_id, , drop = FALSE]
  output_dir <- row$output_dir[[1]]
  start <- proc.time()[["elapsed"]]

  tryCatch({
    inputs <- load_locus_inputs(row, max_variants = max_variants)
    fit <- susieR::susie_rss(
      z = inputs$z,
      R = inputs$R,
      n = inputs$n_rss,
      L = susie_l,
      estimate_residual_variance = TRUE,
      check_prior = FALSE
    )
    elapsed <- proc.time()[["elapsed"]] - start
    effect_metrics <- summarise_effects(
      alpha = fit$alpha,
      R = inputs$R,
      locus_id = row$locus_id[[1]],
      gene_name = row$gene_name[[1]]
    )
    locus_metrics <- summarise_locus(
      row = row,
      inputs = inputs,
      fit = fit,
      effect_metrics = effect_metrics,
      elapsed_seconds = elapsed,
      active_alpha_max_threshold = active_alpha_max_threshold
    )
    top_variants <- build_top_variants(row, inputs, fit, top_n = 20L)
    write_task_outputs(output_dir, locus_metrics, effect_metrics, top_variants)
  }, error = function(e) {
    elapsed <- proc.time()[["elapsed"]] - start
    status <- if (grepl("^too_large_for_baseline_susie:", conditionMessage(e))) {
      "too_large_for_baseline_susie"
    } else {
      "failed"
    }
    locus_metrics <- empty_locus_metrics() %>%
      tibble::add_row(
        locus_id = row$locus_id[[1]],
        gene_name = row$gene_name[[1]],
        n_variants = as.integer(scalar_or_na(row$n_variants, NA_integer_)),
        n_rss = NA_integer_,
        status = status,
        error_message = conditionMessage(e),
        elapsed_seconds = elapsed
      )
    write_task_outputs(output_dir, locus_metrics, empty_effect_metrics(), empty_top_variants())
  })
}

`%||%` <- function(x, y) if (is.null(x)) y else x

main()
