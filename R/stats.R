# stats.R
# -----------------------------------------------------------------------------
# Statistical testing functions:
#   1. run_global_test()         - global PERMANOVA / ANOSIM
#   2. run_pairwise_test()       - pairwise PERMANOVA / ANOSIM
#   3. run_alpha_lmm()           - LMM for one alpha diversity metric
#   4. run_alpha_lmm_all_metrics() - loops over all metrics
#   5. run_alpha_lm_simple()     - simple LM (no random effects, fallback)
#
# Compatible with:
#   config.R  -> padj_method, subject_id_col, is_baseline_col
#   data.R    -> pcoa_obj$meta, $dist, $meta_enriched, $config
#   app.R     -> called via eventReactive on button press
#   mapping   -> CatID, Timepoint, Treatment, Sequence, IsBaseline
# -----------------------------------------------------------------------------

# -- Internal helper: reject the same column in multiple statistical roles --
# Picking one column as primary variable + covariate + strata risks
# rank-deficiency in adonis2 and is almost always a UI mistake.
.check_var_conflicts <- function(v, covariates = NULL, strata = NULL) {
  if (!is.null(covariates) && v %in% covariates)
    stop(sprintf(
      "'%s' cannot be used as both the primary grouping variable and a covariate. Remove it from the covariate list.",
      v
    ))
  if (!is.null(strata) && identical(v, strata))
    stop(sprintf(
      "'%s' cannot be used as both the primary grouping variable and strata. Choose a different strata column (e.g. a subject ID) or a different grouping variable.",
      v
    ))
  if (!is.null(strata) && !is.null(covariates) && strata %in% covariates)
    stop(sprintf(
      "'%s' cannot be used as both a covariate and strata. Remove it from one.",
      strata
    ))
}

# -- Internal helper: reject a continuous column as a grouping variable -----
# PERMANOVA/ANOSIM group samples by factor level; a numeric column with
# many unique values (e.g. Shannon diversity) produces near-singleton
# "groups" and a meaningless test, not an error, so it needs an explicit
# guard rather than relying on the test itself to fail loudly.
.check_not_continuous <- function(meta_f, v, max_levels = 15) {
  vals <- meta_f[[v]]
  if (is.numeric(vals)) {
    n_unique <- length(unique(vals[!is.na(vals)]))
    if (n_unique > max_levels)
      stop(sprintf(
        paste0(
          "'%s' looks like a continuous numeric variable (%d unique values), ",
          "not a grouping factor. PERMANOVA/ANOSIM require a categorical ",
          "grouping variable (e.g. Treatment, Timepoint) - a continuous ",
          "column such as an alpha-diversity metric belongs in the ",
          "'covariates' list instead, not as the primary variable."
        ),
        v, n_unique
      ))
  }
}

# -- Internal helper: build a restricted-permutation block factor safely ---
# Only restricts permutations when there are >=2 blocks AND every block has
# >=2 observations. permute::how() can error (or silently do the wrong
# thing) on singleton blocks - e.g. a subject present in only one phase -
# so in that case we fall back to unrestricted permutations with a warning
# rather than letting the test crash.
.safe_strata_blocks <- function(meta_f, strata) {
  if (is.null(strata) || !strata %in% colnames(meta_f)) return(NULL)

  blocks <- factor(meta_f[[strata]])
  tbl    <- table(blocks)

  if (length(tbl) < 2) {
    warning(sprintf(
      "Strata column '%s' has fewer than 2 blocks after filtering - permutations left unrestricted.",
      strata
    ))
    return(NULL)
  }

  if (any(tbl < 2)) {
    warning(sprintf(
      "Strata column '%s' has block(s) with a single observation (%s) - permutations left unrestricted to avoid permute::how() errors.",
      strata, paste(names(tbl)[tbl < 2], collapse = ", ")
    ))
    return(NULL)
  }

  blocks
}

# -- Internal helper: filter baseline + missing --------------------------------
# Removes baseline rows and rows with NA in any requested variable
.filter_meta <- function(meta, vars_needed,
                         is_baseline_col = "IsBaseline",
                         exclude_baseline = TRUE) {

  # Exclude baseline samples from treatment comparisons
  if (isTRUE(exclude_baseline) &&
      is_baseline_col %in% colnames(meta)) {
    meta <- meta[meta[[is_baseline_col]] == FALSE, , drop = FALSE]
  }

  # Complete cases across all needed variables
  sub <- meta[, intersect(vars_needed, colnames(meta)), drop = FALSE]
  keep <- stats::complete.cases(sub)

  if (sum(keep) == 0)
    stop(sprintf(
      "No complete cases remaining after filtering. Check columns: %s",
      paste(vars_needed, collapse = ", ")
    ))

  if (sum(!keep) > 0)
    message(sprintf(
      "%d sample(s) dropped due to NA in: %s",
      sum(!keep),
      paste(vars_needed, collapse = ", ")
    ))

  meta[keep, , drop = FALSE]
}

# -- 1. run_global_test() ------------------------------------------------------

#' Run global PERMANOVA or ANOSIM on a PCoA object
#'
#' For crossover designs, supply \code{covariates = c("Timepoint",
#' "Sequence")} and \code{strata = "CatID"} to restrict permutations
#' within individual animals and adjust for period/sequence effects.
#'
#' @param obj         A \code{pcoa_obj} from \code{prepare_pcoa()}.
#' @param v           Primary grouping variable (metadata column).
#' @param method      Either \code{"permanova"} or \code{"anosim"}.
#' @param covariates  Character vector of covariate column names to add
#'   before \code{v} in the PERMANOVA formula. Ignored for ANOSIM.
#' @param strata      Metadata column for restricted permutations
#'   (e.g. \code{"CatID"} for crossover design). \code{NULL} = unrestricted.
#' @param permutations Number of permutations (default 999).
#' @param exclude_baseline Logical; exclude \code{IsBaseline == TRUE} rows.
#' @return A data frame with the global test result.
#' @importFrom vegan adonis2 anosim
#' @importFrom permute how
#' @importFrom stats as.formula complete.cases
#' @export
run_global_test <- function(obj, v,
                            method           = c("permanova", "anosim"),
                            covariates       = NULL,
                            strata           = NULL,
                            permutations     = 999,
                            exclude_baseline = TRUE) {

  method <- match.arg(method)
  cfg    <- obj$config

  # -- Guard against conflicting variable roles ------------------
  .check_var_conflicts(v, covariates, strata)

  # -- Determine columns needed for complete-case filter --------
  vars_needed <- unique(c(v, covariates, strata))

  # -- Filter metadata -------------------------------------------
  meta_f <- .filter_meta(
    meta             = obj$meta,
    vars_needed      = vars_needed,
    is_baseline_col  = cfg$is_baseline_col,
    exclude_baseline = exclude_baseline
  )

  # -- Guard against a continuous column as the grouping variable -
  .check_not_continuous(meta_f, v)

  # -- Subset distance matrix to filtered samples ---------------
  d <- stats::as.dist(
    as.matrix(obj$dist)[rownames(meta_f), rownames(meta_f)]
  )

  # -- Validate primary grouping variable ------------------------
  grp <- factor(meta_f[[v]])
  if (nlevels(grp) < 2)
    stop(sprintf(
      "Variable '%s' has fewer than 2 levels after filtering: %s",
      v, paste(levels(grp), collapse = ", ")
    ))

  # -- PERMANOVA ------------------------------------------------
  if (method == "permanova") {

    # Build formula: covariates entered BEFORE primary variable
    # (sequential SS - order matters in adonis2 with by="terms")
    rhs <- paste(c(covariates, v), collapse = " + ")
    fml <- stats::as.formula(paste("d ~", rhs))

    # Restricted permutations within strata - only if every block
    # has enough observations (.safe_strata_blocks warns + falls back
    # to unrestricted permutations otherwise).
    perm_ctrl <- permutations
    blocks    <- .safe_strata_blocks(meta_f, strata)
    if (!is.null(blocks)) {
      perm_ctrl <- permute::how(
        nperm  = permutations,
        blocks = blocks
      )
      message(sprintf(
        "Permutations restricted within strata: %s (%d blocks)",
        strata, length(unique(blocks))
      ))
    }

    fit   <- vegan::adonis2(fml, data = meta_f,
                            permutations = perm_ctrl,
                            by = "terms")

    # Extract the row for the primary variable v
    row_v <- which(rownames(fit) == v)
    if (length(row_v) == 0)
      stop(sprintf("Variable '%s' not found in adonis2 output.", v))

    data.frame(
      Test         = "PERMANOVA (global)",
      Variable     = v,
      Covariates   = if (is.null(covariates)) NA_character_
      else paste(covariates, collapse = " + "),
      Strata       = if (is.null(strata)) NA_character_ else strata,
      N_samples    = nrow(meta_f),
      Statistic_R2 = round(fit$R2[row_v], 4),
      F_value      = round(fit$F[row_v], 4),
      p_value      = fit$`Pr(>F)`[row_v],
      Permutations = permutations,
      stringsAsFactors = FALSE
    )

    # -- ANOSIM ---------------------------------------------------
  } else {

    if (!is.null(covariates))
      warning("ANOSIM does not support covariates. 'covariates' ignored.")

    strata_vec <- .safe_strata_blocks(meta_f, strata)

    fit <- vegan::anosim(d, grp,
                         permutations = permutations,
                         strata       = strata_vec)

    data.frame(
      Test         = "ANOSIM (global)",
      Variable     = v,
      Covariates   = NA_character_,
      Strata       = if (is.null(strata)) NA_character_ else strata,
      N_samples    = nrow(meta_f),
      Statistic_R  = round(fit$statistic, 4),
      F_value      = NA_real_,
      p_value      = fit$signif,
      Permutations = permutations,
      stringsAsFactors = FALSE
    )
  }
}

# -- 2. run_pairwise_test() ----------------------------------------------------

#' Run pairwise PERMANOVA or ANOSIM with p-value adjustment
#'
#' @inheritParams run_global_test
#' @return A data frame with pairwise results and BH-adjusted p-values.
#' @importFrom vegan adonis2 anosim
#' @importFrom permute how
#' @importFrom stats as.dist p.adjust complete.cases as.formula
#' @importFrom utils combn
#' @export
run_pairwise_test <- function(obj, v,
                              method           = c("permanova", "anosim"),
                              covariates       = NULL,
                              strata           = NULL,
                              permutations     = 999,
                              exclude_baseline = TRUE) {

  method    <- match.arg(method)
  cfg       <- obj$config
  padj_meth <- cfg$padj_method   # from pcoa_config() (config.R)

  # -- Guard against conflicting variable roles ------------------
  .check_var_conflicts(v, covariates, strata)

  # -- Filter metadata -------------------------------------------
  vars_needed <- unique(c(v, covariates, strata))

  meta_f <- .filter_meta(
    meta             = obj$meta,
    vars_needed      = vars_needed,
    is_baseline_col  = cfg$is_baseline_col,
    exclude_baseline = exclude_baseline
  )

  # -- Guard against a continuous column as the grouping variable -
  .check_not_continuous(meta_f, v)

  full_mat <- as.matrix(obj$dist)[rownames(meta_f), rownames(meta_f)]
  grp      <- factor(meta_f[[v]])
  lvls     <- levels(grp)

  if (length(lvls) < 2)
    stop(sprintf("Variable '%s' must have at least 2 levels.", v))

  pairs <- utils::combn(lvls, 2, simplify = FALSE)

  # -- Loop over all pairs ---------------------------------------
  res <- lapply(pairs, function(pr) {

    idx    <- grp %in% pr
    md_sub <- meta_f[idx, , drop = FALSE]
    d_sub  <- stats::as.dist(full_mat[idx, idx])
    g_sub  <- factor(md_sub[[v]])

    if (method == "permanova") {

      rhs <- paste(c(covariates, v), collapse = " + ")
      fml <- stats::as.formula(paste("d_sub ~", rhs))

      # Restricted permutations within strata, but only when every block
      # in THIS pair's subset has enough observations (e.g. a subject
      # present in only one phase would otherwise error inside
      # permute::how() for this particular pairwise comparison).
      perm_ctrl <- permutations
      blocks    <- .safe_strata_blocks(md_sub, strata)
      if (!is.null(blocks)) {
        perm_ctrl <- permute::how(
          nperm  = permutations,
          blocks = blocks
        )
      }

      fit   <- vegan::adonis2(fml, data = md_sub,
                              permutations = perm_ctrl,
                              by = "terms")
      row_v <- which(rownames(fit) == v)

      data.frame(
        Group1    = pr[1],
        Group2    = pr[2],
        N         = nrow(md_sub),
        Statistic = round(fit$R2[row_v], 4),
        F_value   = round(fit$F[row_v],  4),
        p_value   = fit$`Pr(>F)`[row_v],
        stringsAsFactors = FALSE
      )

    } else {

      strata_vec <- .safe_strata_blocks(md_sub, strata)

      fit <- vegan::anosim(d_sub, g_sub,
                           permutations = permutations,
                           strata       = strata_vec)

      data.frame(
        Group1    = pr[1],
        Group2    = pr[2],
        N         = nrow(md_sub),
        Statistic = round(fit$statistic, 4),
        F_value   = NA_real_,
        p_value   = fit$signif,
        stringsAsFactors = FALSE
      )
    }
  })

  out <- do.call(rbind, res)

  # -- p-value adjustment ----------------------------------------
  out$p_adjusted    <- stats::p.adjust(out$p_value, method = padj_meth)
  out$Test          <- ifelse(method == "permanova",
                              "PERMANOVA (pairwise)", "ANOSIM (pairwise)")
  out$Variable      <- v
  out$Covariates    <- if (is.null(covariates)) NA_character_
  else paste(covariates, collapse = " + ")
  out$Adjust_method <- padj_meth

  # Reorder columns for readability
  out[, c("Test", "Variable", "Group1", "Group2", "N",
          "Statistic", "F_value", "p_value",
          "p_adjusted", "Adjust_method", "Covariates")]
}

# -- 3. run_alpha_lmm() --------------------------------------------------------

#' Fit a linear mixed model for one alpha diversity metric
#'
#' Designed for the 2x2 crossover design:
#' \code{metric ~ Treatment + Baseline_metric + Timepoint + Sequence
#'               + (1|CatID)}
#'
#' Any subset of fixed effects can be included via \code{fixed_effects}.
#' The random effect (\code{random_effect}) is always \code{CatID} by
#' default to account for repeated measures within each cat.
#' Baseline samples are automatically excluded.
#'
#' @param meta_enriched Data frame from \code{pcoa_obj$meta_enriched}
#'   (output of \code{prepare_all()} in data.R).
#' @param metric        Alpha diversity column to model (e.g.
#'   \code{"Shannon"}).
#' @param fixed_effects Character vector of fixed effect column names.
#'   Default: \code{c("Treatment", "Timepoint", "Sequence")}.
#'   Add \code{"Baseline_Shannon"} (or relevant metric) here if present.
#' @param random_effect Metadata column for the random intercept
#'   (default \code{"CatID"}).
#' @param treatment_ref Reference level for the treatment variable
#'   (default \code{"CON"}).
#' @param treatment_col Name of the treatment column (default
#'   \code{"Treatment"}).
#' @param is_baseline_col Column flagging baseline rows (default
#'   \code{"IsBaseline"}).
#' @param padj_method   p-value adjustment method (default \code{"BH"}).
#' @return A data frame with one row per model term containing:
#'   Metric, Term, Estimate, SE, df, t_value, p_value, p_adjusted,
#'   and the model formula used.
#' @importFrom lmerTest lmer
#' @importFrom stats as.formula p.adjust relevel complete.cases
#' @export
run_alpha_lmm <- function(meta_enriched,
                          metric          = "Shannon",
                          fixed_effects   = c("Treatment",
                                              "Timepoint",
                                              "Sequence"),
                          random_effect   = "CatID",
                          treatment_ref   = "CON",
                          treatment_col   = "Treatment",
                          is_baseline_col = "IsBaseline",
                          padj_method     = "BH") {

  # -- Validate --------------------------------------------------
  if (!metric %in% colnames(meta_enriched))
    stop(sprintf("Metric '%s' not found in meta_enriched.", metric))

  if (!random_effect %in% colnames(meta_enriched))
    stop(sprintf(
      "Random effect column '%s' not found. Check subject_id_col.",
      random_effect
    ))

  # -- Filter baseline rows ----------------------------------
  md <- meta_enriched
  if (is_baseline_col %in% colnames(md)) {
    md <- md[md[[is_baseline_col]] == FALSE, , drop = FALSE]
    message(sprintf("Baseline samples excluded. N remaining = %d",
                    nrow(md)))
  }

  # -- Only keep fixed effects that exist in the data ------------
  available_fe <- intersect(fixed_effects, colnames(md))
  dropped_fe   <- setdiff(fixed_effects, colnames(md))

  if (length(dropped_fe) > 0)
    warning(sprintf(
      "Fixed effect(s) not found and skipped: %s",
      paste(dropped_fe, collapse = ", ")
    ))

  if (length(available_fe) == 0)
    stop("No fixed effects available in meta_enriched.")

  # -- Set treatment reference level -----------------------------
  if (treatment_col %in% colnames(md)) {
    md[[treatment_col]] <- stats::relevel(
      factor(md[[treatment_col]]),
      ref = treatment_ref
    )
  }

  # -- Complete cases for all model variables --------------------
  model_vars <- c(metric, available_fe, random_effect)
  cc         <- stats::complete.cases(md[, model_vars, drop = FALSE])

  if (sum(!cc) > 0)
    message(sprintf(
      "%d row(s) dropped due to NA in model variables (e.g. MONA).",
      sum(!cc)
    ))

  md <- md[cc, , drop = FALSE]

  if (nrow(md) < 5)
    stop("Too few observations after filtering to fit LMM.")



  # -- Build formula ---------------------------------------------
  # e.g.: Shannon ~ Treatment + Baseline_Shannon + Timepoint +
  #                 Sequence + (1|CatID)
  # Treatment is always LAST in fixed effects so its p-value
  # is estimated after accounting for all covariates

  # Separate treatment from other covariates for clear ordering:
  # order: baseline covariate -> period -> sequence -> treatment
  covariate_terms <- setdiff(available_fe, treatment_col)
  ordered_fe      <- c(covariate_terms, treatment_col)

  fml_str <- sprintf(
    "%s ~ %s + (1|%s)",
    metric,
    paste(ordered_fe, collapse = " + "),
    random_effect
  )
  fml <- stats::as.formula(fml_str)
  message(sprintf("Fitting LMM: %s", fml_str))

  # -- Fit model -------------------------------------------------
  # Fit with lmerTest::lmer() (not lme4::lmer()) - this returns a
  # lmerModLmerTest object, which is required for summary() below to
  # dispatch to lmerTest's S3 method and report Satterthwaite df/p-values.
  # Fitting with plain lme4::lmer() would silently skip Satterthwaite
  # correction (or error, since lmerTest does not export a standalone
  # `summary` function - only an S3 method for its own model class).
  fit <- tryCatch(
    lmerTest::lmer(fml, data = md, REML = FALSE),
    error = function(e) {
      stop(sprintf(
        "LMM failed for metric '%s': %s\nFormula: %s",
        metric, e$message, fml_str
      ))
    }
  )

  # -- Extract coefficients ----------------------------------------
  # summary() dispatches to lmerTest:::summary.lmerModLmerTest here
  # because fit was created with lmerTest::lmer(), giving Satterthwaite
  # df and p-values.
  sum_fit  <- summary(fit)
  coef_tbl <- as.data.frame(sum_fit$coefficients)

  # Rename columns consistently
  colnames(coef_tbl) <- make.names(colnames(coef_tbl))
  coef_tbl$Term      <- rownames(coef_tbl)
  coef_tbl$Metric    <- metric
  coef_tbl$Formula   <- fml_str

  # Standardise column names across R versions
  col_map <- c(
    "Estimate"   = "Estimate",
    "Std..Error" = "SE",
    "df"         = "df",
    "t.value"    = "t_value",
    "Pr...t.."   = "p_value"
  )
  for (old in names(col_map)) {
    if (old %in% colnames(coef_tbl))
      colnames(coef_tbl)[colnames(coef_tbl) == old] <- col_map[old]
  }

  # -- Keep only needed columns (handle missing df gracefully) ---
  keep_cols <- intersect(
    c("Metric", "Term", "Estimate", "SE", "df", "t_value",
      "p_value", "Formula"),
    colnames(coef_tbl)
  )
  coef_tbl <- coef_tbl[, keep_cols, drop = FALSE]
  rownames(coef_tbl) <- NULL

  # -- FDR correction across terms -------------------------------
  if ("p_value" %in% colnames(coef_tbl)) {
    coef_tbl$p_adjusted <- stats::p.adjust(
      coef_tbl$p_value,
      method = padj_method
    )
  }

  coef_tbl
}

# -- 4. run_alpha_lmm_all_metrics() --------------------------------------------

#' Run LMM for all available alpha diversity metrics
#'
#' Loops over Shannon, Simpson, Observed_ASVs, and Chao1 (any present
#' in \code{meta_enriched}), fits a separate LMM for each, and returns
#' a single combined data frame.
#'
#' For the crossover design, the function automatically detects
#' and includes the corresponding baseline covariate column
#' (e.g. \code{Baseline_Shannon} for the Shannon metric) if present
#' in \code{meta_enriched}.
#'
#' @param meta_enriched Data frame from \code{pcoa_obj$meta_enriched}.
#' @param fixed_effects Character vector of fixed effects to include
#'   (excluding the auto-detected baseline covariate).
#'   Default: \code{c("Timepoint", "Sequence", "Treatment")}.
#' @param random_effect Subject ID column for random intercept
#'   (default \code{"CatID"}).
#' @param treatment_ref Reference level for treatment
#'   (default \code{"CON"}).
#' @param treatment_col Treatment column name (default \code{"Treatment"}).
#' @param is_baseline_col Baseline flag column (default \code{"IsBaseline"}).
#' @param padj_method p-value adjustment method (default \code{"BH"}).
#' @return A combined data frame with one row per model term per metric,
#'   with FDR correction applied across Treatment terms only.
#' @export
run_alpha_lmm_all_metrics <- function(
    meta_enriched,
    fixed_effects   = c("Timepoint", "Sequence", "Treatment"),
    random_effect   = "CatID",
    treatment_ref   = "CON",
    treatment_col   = "Treatment",
    is_baseline_col = "IsBaseline",
    padj_method     = "BH") {

  # -- Detect available metrics in the data ----------------------
  # ALPHA_METRICS is the single source of truth (config.R) - do not
  # re-hardcode the metric list here, it will drift.
  avail_metrics <- intersect(ALPHA_METRICS, colnames(meta_enriched))

  if (length(avail_metrics) == 0)
    stop(paste0(
      "No alpha diversity metrics found in meta_enriched.\n",
      "Expected one or more of: ",
      paste(ALPHA_METRICS, collapse = ", ")
    ))

  message(sprintf(
    "Running LMM for %d metric(s): %s",
    length(avail_metrics),
    paste(avail_metrics, collapse = ", ")
  ))

  # -- Loop over each metric -------------------------------------
  results <- lapply(avail_metrics, function(metric) {

    # Auto-detect baseline covariate for this metric
    # e.g. Baseline_Shannon for Shannon metric
    baseline_col <- paste0("Baseline_", metric)

    fe_for_this_metric <- if (baseline_col %in% colnames(meta_enriched))  {
      # Insert baseline covariate FIRST (before period/sequence/treatment)
      # so it is partialled out first in the sequential model
      unique(c(baseline_col, fixed_effects))
    } else {
      warning(sprintf(
        "Baseline covariate '%s' not found - fitting without it.",
        baseline_col
      ))
      fixed_effects
    }

    tryCatch(
      run_alpha_lmm(
        meta_enriched   = meta_enriched,
        metric          = metric,
        fixed_effects   = fe_for_this_metric,
        random_effect   = random_effect,
        treatment_ref   = treatment_ref,
        treatment_col   = treatment_col,
        is_baseline_col = is_baseline_col,
        padj_method     = padj_method
      ),
      error = function(e) {
        warning(sprintf("LMM failed for '%s': %s", metric, e$message))
        # Return empty placeholder so other metrics still run
        data.frame(
          Metric      = metric,
          Term        = "ERROR",
          Estimate    = NA_real_,
          SE          = NA_real_,
          df          = NA_real_,
          t_value     = NA_real_,
          p_value     = NA_real_,
          p_adjusted  = NA_real_,
          Formula     = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    )
  })

  out <- do.call(rbind, results)

  # -- Cross-metric FDR: correct Treatment p-values across metrics -
  # This adjusts for testing the same hypothesis (Treatment effect)
  # across multiple alpha diversity metrics simultaneously
  treat_pattern <- paste0("^", treatment_col)
  treat_rows    <- grepl(treat_pattern, out$Term)

  if (any(treat_rows, na.rm = TRUE)) {
    out$p_adjusted_cross_metric        <- NA_real_
    out$p_adjusted_cross_metric[treat_rows] <- stats::p.adjust(
      out$p_value[treat_rows],
      method = padj_method
    )
    message(sprintf(
      "Cross-metric FDR applied to %d Treatment term(s) using %s",
      sum(treat_rows), padj_method
    ))
  }

  rownames(out) <- NULL
  out
}

# -- 5. run_alpha_lm_simple() -------------------------------------------------

#' Simple linear model for alpha diversity (no random effects)
#'
#' Fallback when only one observation per subject exists, or for
#' exploratory analysis without mixed modelling.
#' For the crossover design, \code{run_alpha_lmm()} is preferred.
#'
#' @inheritParams run_alpha_lmm
#' @return A data frame with one row per model term.
#' @importFrom stats lm p.adjust as.formula relevel complete.cases
#' @export
run_alpha_lm_simple <- function(
    meta_enriched,
    metric          = "Shannon",
    fixed_effects   = c("Treatment", "Timepoint", "Sequence"),
    treatment_ref   = "CON",
    treatment_col   = "Treatment",
    is_baseline_col = "IsBaseline",
    padj_method     = "BH") {

  if (!metric %in% colnames(meta_enriched))
    stop(sprintf("Metric '%s' not found.", metric))

  # -- Filter baseline ---------------------------------------
  md <- meta_enriched
  if (is_baseline_col %in% colnames(md))
    md <- md[md[[is_baseline_col]] == FALSE, , drop = FALSE]

  # -- Set reference level ---------------------------------------
  if (treatment_col %in% colnames(md))
    md[[treatment_col]] <- stats::relevel(
      factor(md[[treatment_col]]), ref = treatment_ref
    )

  # -- Available fixed effects ------------------------------------
  available_fe <- intersect(fixed_effects, colnames(md))

  # -- Complete cases --------------------------------------------
  model_vars <- c(metric, available_fe)
  cc         <- stats::complete.cases(md[, model_vars, drop = FALSE])
  md         <- md[cc, , drop = FALSE]

  # -- Build and fit formula -------------------------------------
  fml_str <- sprintf("%s ~ %s",
                     metric,
                     paste(available_fe, collapse = " + "))
  fml     <- stats::as.formula(fml_str)
  message(sprintf("Fitting LM: %s", fml_str))

  fit     <- stats::lm(fml, data = md)
  coef_tbl <- as.data.frame(summary(fit)$coefficients)

  colnames(coef_tbl) <- c("Estimate", "SE", "t_value", "p_value")
  coef_tbl$Term      <- rownames(coef_tbl)
  coef_tbl$Metric    <- metric
  coef_tbl$Formula   <- fml_str
  coef_tbl$p_adjusted <- stats::p.adjust(coef_tbl$p_value,
                                         method = padj_method)
  rownames(coef_tbl) <- NULL

  coef_tbl[, c("Metric", "Term", "Estimate", "SE",
               "t_value", "p_value", "p_adjusted", "Formula")]
}
