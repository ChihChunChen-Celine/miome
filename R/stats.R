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

# -- Internal helper: parse a free-text fixed-effects formula fragment -----
# Accepts a formula RHS like "Treatment + Timepoint" or
# "Treatment * Timepoint" or "Treatment + Timepoint + Treatment:Sequence".
# Returns the set of underlying variable names referenced (for existence
# checks and complete-case filtering) - unlike the old character-vector
# API, this correctly picks up variables that only appear inside an
# interaction term.
.parse_fixed_formula <- function(fixed_formula, meta) {
  if (is.null(fixed_formula) || !nzchar(trimws(fixed_formula)))
    stop(paste0(
      "Fixed effects formula is empty. Please enter at least one term, ",
      "e.g. 'Treatment + Timepoint' or 'Treatment * Timepoint'."
    ), call. = FALSE)

  fml <- tryCatch(
    stats::as.formula(paste("~", fixed_formula)),
    error = function(e)
      stop("Cannot parse fixed-effects formula '", fixed_formula, "':\n  ",
           e$message,
           "\n(Tip: wrap column names containing spaces or special ",
           "characters in backticks, e.g. `Time Point`.)",
           call. = FALSE)
  )

  vars_used <- unique(all.vars(fml))
  missing_cols <- setdiff(vars_used, colnames(meta))
  if (length(missing_cols) > 0)
    stop(sprintf(
      "Fixed-effects formula term(s) not found in metadata: %s\nAvailable columns: %s",
      paste(missing_cols, collapse = ", "),
      paste(colnames(meta), collapse = ", ")
    ), call. = FALSE)

  list(fml = fml, vars_used = vars_used)
}

# -- Internal helper: rebuild a fixed-effects RHS with every bare variable
# name backtick-quoted, so column names containing spaces or other
# non-syntactic characters never break stats::as.formula() downstream.
# Uses terms()'s term.labels (not the raw user text) so "*" has already
# been expanded into main effects + interaction(s), and interaction terms
# come back as "VarA:VarB" - each side is quoted independently.
.build_quoted_fixed_rhs <- function(fml) {
  term_labels <- attr(stats::terms(fml), "term.labels")
  if (length(term_labels) == 0)
    stop("Fixed-effects formula resolved to zero terms.", call. = FALSE)

  quoted_terms <- vapply(term_labels, function(term) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    parts <- vapply(parts, function(p) {
      # Leave function-wrapped terms (e.g. I(Age^2)) unquoted - backticks
      # only make sense around bare variable names.
      if (grepl("[()]", p)) p else sprintf("`%s`", p)
    }, character(1))
    paste(parts, collapse = ":")
  }, character(1))

  paste(quoted_terms, collapse = " + ")
}

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
                            exclude_baseline = TRUE,
                            formula_rhs      = NULL) {

  method <- match.arg(method)
  cfg    <- obj$config

  # -- Guard against conflicting variable roles ------------------
  .check_var_conflicts(v, covariates, strata)

  # -- Determine columns needed for complete-case filter --------
  # If a free-text formula is supplied, extract its variables so
  # complete-case filtering covers everything the model references
  # (including variables that only appear inside an interaction).
  if (!is.null(formula_rhs) && nzchar(trimws(formula_rhs))) {
    fml_tmp     <- stats::as.formula(paste("~", formula_rhs))
    formula_vars <- all.vars(fml_tmp)               # base::all.vars (NOT stats::)
    missing_cols <- setdiff(formula_vars, colnames(obj$meta))
    if (length(missing_cols) > 0)
      stop(sprintf(
        "Formula term(s) not found in metadata: %s\nAvailable: %s",
        paste(missing_cols, collapse = ", "),
        paste(colnames(obj$meta), collapse = ", ")
      ), call. = FALSE)
    vars_needed <- unique(c(formula_vars, strata))
  } else {
    vars_needed <- unique(c(v, covariates, strata))
  }
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

    # Build formula. If a free-text RHS is supplied (supports "*" and ":"
    # for interactions) use it verbatim; otherwise fall back to the
    # covariates-then-primary-variable additive construction.
    if (!is.null(formula_rhs) && nzchar(trimws(formula_rhs))) {
      rhs <- formula_rhs
    } else {
      rhs <- paste(c(covariates, v), collapse = " + ")
    }
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
    ft <- as.data.frame(fit)
    ft$Term <- rownames(ft)

    if (!is.null(formula_rhs) && nzchar(trimws(formula_rhs))) {
      # Free-text formula: return EVERY model term (main effects +
      # interaction[s] + Residual/Total) so the interaction row is visible.
      keep <- !ft$Term %in% c("Residual", "Total")
      out_tbl <- data.frame(
        Test         = "PERMANOVA (global)",
        Variable     = ft$Term[keep],
        Covariates   = NA_character_,
        Strata       = if (is.null(strata)) NA_character_ else strata,
        N_samples    = nrow(meta_f),
        Statistic_R2 = round(ft$R2[keep], 4),
        F_value      = round(ft$F[keep], 4),
        p_value      = ft$`Pr(>F)`[keep],
        Permutations = permutations,
        stringsAsFactors = FALSE
      )
      rownames(out_tbl) <- NULL
      out_tbl
    } else {
      # Original behaviour: single row for the primary variable v.
      row_v <- which(ft$Term == v)
      if (length(row_v) == 0)
        stop(sprintf("Variable '%s' not found in adonis2 output.", v))
      data.frame(
        Test         = "PERMANOVA (global)",
        Variable     = v,
        Covariates   = if (is.null(covariates)) NA_character_
        else paste(covariates, collapse = " + "),
        Strata       = if (is.null(strata)) NA_character_ else strata,
        N_samples    = nrow(meta_f),
        Statistic_R2 = round(ft$R2[row_v], 4),
        F_value      = round(ft$F[row_v], 4),
        p_value      = ft$`Pr(>F)`[row_v],
        Permutations = permutations,
        stringsAsFactors = FALSE
      )
    }
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
#' Accepts a free-text fixed-effects formula fragment, so interactions are
#' supported directly, e.g. \code{fixed_formula = "Treatment * Timepoint"}
#' or \code{"Treatment + Timepoint + Treatment:Sequence"}.
#' The random effect (\code{random_effect}) is always a single random
#' intercept, e.g. \code{(1|CatID)}. Baseline samples are automatically
#' excluded.
#'
#' @param meta_enriched Data frame from \code{pcoa_obj$meta_enriched}.
#' @param metric        Alpha diversity column to model (e.g. \code{"Shannon"}).
#' @param fixed_formula Free-text fixed-effects formula fragment (no
#'   leading \code{~}, no random effect term). Supports \code{+}, \code{*},
#'   \code{:}. Column names with spaces or special characters must be
#'   backtick-quoted, e.g. \code{`Time Point`}.
#'   Default: \code{"Treatment + Timepoint + Sequence"}.
#' @param random_effect Metadata column for the random intercept
#'   (default \code{"CatID"}).
#' @param treatment_ref Reference level for the treatment variable
#'   (default \code{"CON"}).
#' @param treatment_col Name of the treatment column (default
#'   \code{"Treatment"}).
#' @param is_baseline_col Column flagging baseline rows (default
#'   \code{"IsBaseline"}).
#' @param padj_method   p-value adjustment method (default \code{"BH"}).
#' @return A data frame with one row per model term.
#' @importFrom lmerTest lmer
#' @importFrom stats as.formula p.adjust relevel complete.cases terms
#' @export
run_alpha_lmm <- function(meta_enriched,
                          metric          = "Shannon",
                          fixed_formula   = "Treatment + Timepoint + Sequence",
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

  parsed <- .parse_fixed_formula(fixed_formula, meta_enriched)

  # -- Filter baseline rows ----------------------------------
  md <- meta_enriched
  if (is_baseline_col %in% colnames(md)) {
    md <- md[md[[is_baseline_col]] == FALSE, , drop = FALSE]
    message(sprintf("Baseline samples excluded. N remaining = %d",
                    nrow(md)))
  }

  # -- Set treatment reference level (only if treatment_col is actually
  #    used somewhere in the formula) -----------------------------
  if (treatment_col %in% colnames(md) &&
      treatment_col %in% parsed$vars_used) {
    md[[treatment_col]] <- stats::relevel(
      factor(md[[treatment_col]]),
      ref = treatment_ref
    )
  }

  # -- Complete cases for all model variables --------------------
  model_vars <- unique(c(metric, parsed$vars_used, random_effect))
  cc         <- stats::complete.cases(md[, model_vars, drop = FALSE])

  if (sum(!cc) > 0)
    message(sprintf(
      "%d row(s) dropped due to NA in model variables (e.g. MONA).",
      sum(!cc)
    ))

  md <- md[cc, , drop = FALSE]

  if (nrow(md) < 5)
    stop("Too few observations after filtering to fit LMM.")

  # -- Build formula -----------------------------------------------
  # Every bare variable name (including each side of an interaction)
  # is backtick-quoted, so column names with spaces never break
  # as.formula() - see .build_quoted_fixed_rhs().
  quoted_rhs <- .build_quoted_fixed_rhs(parsed$fml)
  fml_str <- sprintf("`%s` ~ %s + (1|`%s`)", metric, quoted_rhs, random_effect)
  fml <- tryCatch(
    stats::as.formula(fml_str),
    error = function(e)
      stop(sprintf("Could not build model formula: %s\nFormula string: %s",
                   e$message, fml_str), call. = FALSE)
  )
  message(sprintf("Fitting LMM: %s", fml_str))

  # -- Fit model -------------------------------------------------
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
  sum_fit  <- summary(fit)
  coef_tbl <- as.data.frame(sum_fit$coefficients)

  colnames(coef_tbl) <- make.names(colnames(coef_tbl))
  coef_tbl$Term      <- rownames(coef_tbl)
  coef_tbl$Metric    <- metric
  coef_tbl$Formula   <- fml_str

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

  keep_cols <- intersect(
    c("Metric", "Term", "Estimate", "SE", "df", "t_value",
      "p_value", "Formula"),
    colnames(coef_tbl)
  )
  coef_tbl <- coef_tbl[, keep_cols, drop = FALSE]
  rownames(coef_tbl) <- NULL

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
    fixed_formula   = "Timepoint + Sequence + Treatment",
    random_effect   = "CatID",
    treatment_ref   = "CON",
    treatment_col   = "Treatment",
    is_baseline_col = "IsBaseline",
    padj_method     = "BH") {

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

  results <- lapply(avail_metrics, function(metric) {

    baseline_col <- paste0("Baseline_", metric)

    fml_for_this_metric <- if (baseline_col %in% colnames(meta_enriched)) {
      # Prepend baseline covariate so it's partialled out first
      paste0(baseline_col, " + ", fixed_formula)
    } else {
      warning(sprintf(
        "Baseline covariate '%s' not found - fitting without it.",
        baseline_col
      ))
      fixed_formula
    }

    tryCatch(
      run_alpha_lmm(
        meta_enriched   = meta_enriched,
        metric          = metric,
        fixed_formula   = fml_for_this_metric,
        random_effect   = random_effect,
        treatment_ref   = treatment_ref,
        treatment_col   = treatment_col,
        is_baseline_col = is_baseline_col,
        padj_method     = padj_method
      ),
      error = function(e) {
        warning(sprintf("LMM failed for '%s': %s", metric, e$message))
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

  # -- Cross-metric FDR: correct Treatment MAIN-EFFECT p-values across
  # metrics only - excludes interaction terms like "TreatmentA:TimepointB"
  # so they aren't mixed into the same adjustment set as the main effect.
  treat_pattern <- paste0("^", treatment_col)
  treat_rows    <- grepl(treat_pattern, out$Term) & !grepl(":", out$Term)

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
# -- run_alpha_lmm_interaction_all() ------------------------------------------
#' Fit interaction LMMs for all alpha metrics WITHOUT baseline-covariate injection
#'
#' Unlike run_alpha_lmm_all_metrics(), this wrapper does NOT prepend a
#' Baseline_<metric> covariate. It fits the same user-supplied fixed-effects
#' formula (typically an interaction model such as "Treatment * Timepoint")
#' to every available metric, then applies cross-metric FDR correction to a
#' user-chosen term of interest (default: the Treatment:Timepoint interaction).
#'
#' Intended for a simple two-timepoint parallel design where day-1 (initial)
#' IS the baseline and is modelled directly via the Timepoint factor, so no
#' separate ANCOVA baseline covariate is wanted.
#'
#' @param meta_enriched Data frame from pcoa_obj$meta_enriched.
#' @param fixed_formula Fixed-effects fragment (no ~, no random term).
#'   Default: "Treatment * Timepoint".
#' @param random_effect Subject ID column for the random intercept
#'   (default "CatID").
#' @param treatment_ref Reference level for treatment (default "CON").
#' @param treatment_col Treatment column name (default "Treatment").
#' @param is_baseline_col Baseline flag column (default "IsBaseline").
#' @param padj_method p-value adjustment method (default "BH").
#' @param fdr_term_pattern Regex identifying the term(s) to FDR-correct
#'   across metrics. Default "(?=.*Treatment)(?=.*:)" (perl) = any term that
#'   contains "Treatment" AND a colon, i.e. the treatment interaction.
#' @return Combined data frame: one row per model term per metric, plus a
#'   p_adjusted_cross_metric column populated for the targeted term.
#' @export
run_alpha_lmm_interaction_all <- function(
    meta_enriched,
    fixed_formula    = "Treatment * Timepoint",
    random_effect    = "CatID",
    treatment_ref    = "CON",
    treatment_col    = "Treatment",
    is_baseline_col  = "IsBaseline",
    padj_method      = "BH",
    fdr_term_pattern = "(?=.*Treatment)(?=.*:)") {

  # -- Discover available metrics (uses the same ALPHA_METRICS constant) ----
  avail_metrics <- intersect(ALPHA_METRICS, colnames(meta_enriched))
  if (length(avail_metrics) == 0)
    stop(paste0(
      "No alpha diversity metrics found in meta_enriched.\n",
      "Expected one or more of: ",
      paste(ALPHA_METRICS, collapse = ", ")
    ))

  message(sprintf(
    "Fitting interaction LMM (no baseline covariate) for %d metric(s): %s",
    length(avail_metrics), paste(avail_metrics, collapse = ", ")
  ))
  message(sprintf("Fixed-effects formula: %s", fixed_formula))

  # -- Fit one model per metric, reusing the vetted run_alpha_lmm() ---------
  results <- lapply(avail_metrics, function(metric) {
    tryCatch(
      run_alpha_lmm(
        meta_enriched   = meta_enriched,
        metric          = metric,
        fixed_formula   = fixed_formula,   # NOTE: passed through unchanged
        random_effect   = random_effect,
        treatment_ref   = treatment_ref,
        treatment_col   = treatment_col,
        is_baseline_col = is_baseline_col,
        padj_method     = padj_method
      ),
      error = function(e) {
        warning(sprintf("LMM failed for '%s': %s", metric, e$message))
        data.frame(
          Metric = metric, Term = "ERROR",
          Estimate = NA_real_, SE = NA_real_, df = NA_real_,
          t_value = NA_real_, p_value = NA_real_, p_adjusted = NA_real_,
          Formula = NA_character_, stringsAsFactors = FALSE
        )
      }
    )
  })

  out <- do.call(rbind, results)

  # -- Cross-metric FDR on the interaction term(s) only ---------------------
  # perl=TRUE so the look-ahead pattern "(?=.*Treatment)(?=.*:)" works:
  # it matches any term that contains BOTH "Treatment" and a colon,
  # e.g. "TreatmentBG:TimepointFinal".
  target_rows <- grepl(fdr_term_pattern, out$Term, perl = TRUE) &
    out$Term != "ERROR"

  out$p_adjusted_cross_metric <- NA_real_
  if (any(target_rows, na.rm = TRUE)) {
    out$p_adjusted_cross_metric[target_rows] <- stats::p.adjust(
      out$p_value[target_rows],
      method = padj_method
    )
    message(sprintf(
      "Cross-metric FDR applied to %d interaction term(s) using %s.",
      sum(target_rows), padj_method
    ))
  } else {
    warning(paste0(
      "No interaction terms matched pattern '", fdr_term_pattern,
      "'. Cross-metric FDR column left as NA.\n",
      "Check your fixed_formula actually contains an interaction, and that ",
      "term names look like e.g. 'TreatmentBG:TimepointFinal'."
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
    fixed_formula   = "Treatment + Timepoint + Sequence",
    treatment_ref   = "CON",
    treatment_col   = "Treatment",
    is_baseline_col = "IsBaseline",
    padj_method     = "BH") {

  if (!metric %in% colnames(meta_enriched))
    stop(sprintf("Metric '%s' not found.", metric))

  parsed <- .parse_fixed_formula(fixed_formula, meta_enriched)

  md <- meta_enriched
  if (is_baseline_col %in% colnames(md))
    md <- md[md[[is_baseline_col]] == FALSE, , drop = FALSE]

  if (treatment_col %in% colnames(md) &&
      treatment_col %in% parsed$vars_used)
    md[[treatment_col]] <- stats::relevel(
      factor(md[[treatment_col]]), ref = treatment_ref
    )

  model_vars <- unique(c(metric, parsed$vars_used))
  cc         <- stats::complete.cases(md[, model_vars, drop = FALSE])
  md         <- md[cc, , drop = FALSE]

  quoted_rhs <- .build_quoted_fixed_rhs(parsed$fml)
  fml_str <- sprintf("`%s` ~ %s", metric, quoted_rhs)
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
