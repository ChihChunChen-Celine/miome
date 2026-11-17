# stats.R
# -----------------------------------------------------------------------------
# Statistical testing:
#   1. run_global_test()            - global PERMANOVA / ANOSIM
#   2. run_pairwise_test()          - pairwise PERMANOVA / ANOSIM
#   3. run_alpha_lmm()              - LMM for one alpha diversity metric
#   4. run_alpha_lmm_all_metrics()  - LMM for every metric (+ cross-metric FDR)
#   5. run_alpha_lm_simple()        - LM without random effects
# -----------------------------------------------------------------------------

# -- Internal: parse a fixed-effects RHS ("Treatment * Timepoint") -----------
.parse_fixed_formula <- function(fixed_formula, meta) {
  if (is.null(fixed_formula) || !nzchar(trimws(fixed_formula)))
    stop("The formula is empty. Enter at least one term, e.g. ",
         "'Treatment + Timepoint' or 'Treatment * Timepoint'.", call. = FALSE)
  fml <- tryCatch(
    stats::as.formula(paste("~", fixed_formula)),
    error = function(e)
      stop("Cannot read formula '", fixed_formula, "': ", conditionMessage(e),
           "\n(Tip: wrap column names containing spaces in backticks, ",
           "e.g. `Time Point`.)", call. = FALSE))
  vars <- unique(all.vars(fml))
  miss <- setdiff(vars, colnames(meta))
  if (length(miss))
    stop(sprintf("Formula term(s) not found in metadata: %s\nAvailable columns: %s",
                 paste(miss, collapse = ", "), paste(colnames(meta), collapse = ", ")),
         call. = FALSE)
  list(fml = fml, vars_used = vars)
}

# -- Internal: rebuild the RHS with every variable backtick-quoted -----------
.build_quoted_fixed_rhs <- function(fml) {
  labs <- attr(stats::terms(fml), "term.labels")
  if (!length(labs)) stop("Formula resolved to zero terms.", call. = FALSE)
  quote1 <- function(p) if (grepl("[()`]", p)) p else sprintf("`%s`", p)
  paste(vapply(labs, function(t)
    paste(vapply(strsplit(t, ":", fixed = TRUE)[[1]], quote1, ""), collapse = ":"),
    ""), collapse = " + ")
}

.quote_terms <- function(vars) paste(sprintf("`%s`", vars), collapse = " + ")

# -- Internal: reject the same column in several roles -----------------------
.check_var_conflicts <- function(v, covariates = NULL, strata = NULL) {
  if (v %in% covariates)
    stop(sprintf("'%s' cannot be both the grouping variable and a covariate.", v),
         call. = FALSE)
  if (identical(v, strata))
    stop(sprintf("'%s' cannot be both the grouping variable and strata. Choose a subject ID as strata.", v),
         call. = FALSE)
  if (!is.null(strata) && strata %in% covariates)
    stop(sprintf("'%s' cannot be both a covariate and strata.", strata),
         call. = FALSE)
}

# -- Internal: reject a continuous column as grouping variable ---------------
.check_not_continuous <- function(meta_f, v, max_levels = 15) {
  x <- meta_f[[v]]
  n <- length(unique(x[!is.na(x)]))
  if (is.numeric(x) && n > max_levels)
    stop(sprintf(paste0(
      "'%s' looks like a continuous number (%d unique values), not a group. ",
      "PERMANOVA/ANOSIM need a categorical grouping variable such as Treatment; ",
      "put continuous columns in 'Covariates' instead."), v, n), call. = FALSE)
}

# -- Internal: restricted permutations -------------------------------------------
# Samples are only shuffled within their own block (e.g. within each animal).
# A block with a single sample (e.g. an animal that dropped out) is fine: that
# sample simply stays in place. The restriction is only impossible when NO
# block has two or more samples.
# Returns: blocks (factor or NULL), label for the results table, the
# `permutations` argument for vegan, the number of permutations actually
# available, and a note for the user.
.strata_plan <- function(meta_f, strata, permutations) {
  plan <- list(blocks = NULL, label = NA_character_, perm = permutations,
               n_perm = permutations, note = "")
  if (is.null(strata) || !strata %in% colnames(meta_f)) return(plan)

  blocks <- factor(meta_f[[strata]])
  size   <- table(blocks)
  if (!any(size >= 2)) {
    warning(sprintf(paste0("Strata '%s': every block has a single sample, so there is ",
                           "nothing to permute within blocks - permutations left unrestricted."),
                    strata), call. = FALSE)
    plan$label <- sprintf("%s (NOT applied: every block has one sample)", strata)
    return(plan)
  }
  n_single <- sum(size == 1)
  possible <- prod(factorial(as.numeric(size)))        # distinct within-block shuffles
  plan$blocks <- blocks
  plan$perm   <- permute::how(nperm = permutations, blocks = blocks)
  plan$label  <- if (n_single)
    sprintf("%s (%d single-sample block%s kept fixed)", strata, n_single,
            if (n_single > 1) "s" else "") else strata
  # When the complete set of within-block shuffles is small enough, use ALL of
  # them: the p-value is then exact and identical on every run.
  if (is.finite(possible) && possible <= max(10000, permutations + 1)) {
    plan$perm   <- permute::how(nperm = possible, blocks = blocks,
                                maxperm = possible + 1, minperm = possible + 1)
    plan$n_perm <- possible - 1
    plan$note   <- sprintf("exact test: all %d possible permutations within %s used (smallest possible p = %s)",
                           possible, strata, signif(1 / possible, 2))
  }
  plan
}

# TRUE for model terms whose variables never change within a block: such
# between-subject terms cannot be tested by shuffling within blocks.
.constant_within <- function(terms, meta_f, blocks) {
  if (is.null(blocks)) return(rep(FALSE, length(terms)))
  const <- function(v) v %in% colnames(meta_f) &&
    all(tapply(as.character(meta_f[[v]]), blocks, function(x) length(unique(x)) <= 1))
  vapply(strsplit(gsub("`", "", terms), ":", fixed = TRUE),
         function(vs) all(vapply(vs, const, logical(1))), logical(1))
}

# -- Internal: drop baseline rows (NA-safe) and incomplete cases --------------
.filter_meta <- function(meta, vars_needed, is_baseline_col = "IsBaseline",
                         exclude_baseline = TRUE) {
  if (isTRUE(exclude_baseline))
    meta <- meta[!.is_baseline(meta, is_baseline_col), , drop = FALSE]
  keep <- stats::complete.cases(meta[, intersect(vars_needed, colnames(meta)),
                                     drop = FALSE])
  if (!any(keep))
    stop(sprintf("No complete samples left after filtering. Check columns: %s",
                 paste(vars_needed, collapse = ", ")), call. = FALSE)
  if (any(!keep))
    message(sprintf("%d sample(s) dropped due to NA in: %s", sum(!keep),
                    paste(vars_needed, collapse = ", ")))
  meta[keep, , drop = FALSE]
}

.sub_dist <- function(obj, ids) {
  m <- if (!is.null(obj$dist_mat)) obj$dist_mat else as.matrix(obj$dist)
  stats::as.dist(m[ids, ids])
}

# -- 1. run_global_test() ------------------------------------------------------
#' Run a global PERMANOVA or ANOSIM on a PCoA object
#'
#' For crossover designs supply \code{covariates = c("Timepoint", "Sequence")}
#' and \code{strata = "CatID"} to adjust for period/sequence and restrict
#' permutations within subjects.
#'
#' @param obj A \code{pcoa_obj} from \code{prepare_all()} / \code{prepare_pcoa()}.
#' @param v Primary grouping variable (metadata column).
#' @param method \code{"permanova"} or \code{"anosim"}.
#' @param covariates Covariate columns entered before \code{v} in the
#'   PERMANOVA model. Ignored for ANOSIM.
#' @param strata Column for restricted permutations (e.g. a subject ID), or
#'   \code{NULL}.
#' @param permutations Number of permutations.
#' @param exclude_baseline Logical; drop baseline samples first.
#' @param formula_rhs Optional full PERMANOVA model, e.g.
#'   \code{"Treatment * Timepoint"}. When given it replaces
#'   \code{covariates} + \code{v} and every model term is reported.
#' @return A data frame with one row per fitted model term (covariates first,
#'   in model order; terms are tested sequentially). \code{Strata} says how
#'   the restriction was applied (blocks with one sample stay fixed).
#'   \code{Note} flags terms that do not vary within the strata - their
#'   p-value is \code{NA} because shuffling within blocks cannot test them -
#'   and cases where few distinct permutations exist.
#' @export
run_global_test <- function(obj, v,
                            method           = c("permanova", "anosim"),
                            covariates       = NULL,
                            strata           = NULL,
                            permutations     = 999,
                            exclude_baseline = TRUE,
                            formula_rhs      = NULL) {
  method   <- match.arg(method)
  cfg      <- obj$config
  use_fml  <- method == "permanova" && !is.null(formula_rhs) &&
              nzchar(trimws(formula_rhs))

  if (use_fml) {
    parsed      <- .parse_fixed_formula(formula_rhs, obj$meta)
    if (!is.null(strata) && strata %in% parsed$vars_used)
      stop(sprintf("'%s' is used as strata and also appears in the formula.", strata),
           call. = FALSE)
    vars_needed <- unique(c(parsed$vars_used, strata))
  } else {
    .check_var_conflicts(v, covariates, strata)
    vars_needed <- unique(c(v, covariates, strata))
  }

  meta_f <- .filter_meta(obj$meta, vars_needed, cfg$is_baseline_col, exclude_baseline)
  d      <- .sub_dist(obj, rownames(meta_f))

  if (!use_fml) {
    .check_not_continuous(meta_f, v)
    grp <- factor(meta_f[[v]])
    if (nlevels(grp) < 2)
      stop(sprintf("'%s' has fewer than 2 groups after filtering: %s",
                   v, paste(levels(grp), collapse = ", ")), call. = FALSE)
  }

  plan   <- .strata_plan(meta_f, strata, permutations)
  blocks <- plan$blocks
  notes  <- function(terms) {
    fixed <- .constant_within(terms, meta_f, blocks)
    n <- ifelse(fixed, sprintf("does not vary within %s: cannot be tested with restricted permutations",
                               strata), "")
    list(fixed = fixed, text = trimws(paste(n, ifelse(nzchar(n) & nzchar(plan$note), ";", ""),
                                            plan$note)))
  }

  if (method == "anosim") {
    if (!is.null(covariates))
      warning("ANOSIM does not support covariates; they were ignored.", call. = FALSE)
    message(sprintf("ANOSIM: %s | strata: %s | n = %d", v,
                    if (is.null(blocks)) "none" else strata, nrow(meta_f)))
    fit <- suppressMessages(vegan::anosim(d, grp, permutations = plan$perm))
    nt  <- notes(v)
    return(data.frame(Test = "ANOSIM (global)", Variable = v,
                      Covariates = NA_character_, Strata = plan$label,
                      N_samples = nrow(meta_f),
                      Statistic_R = round(fit$statistic, 4), F_value = NA_real_,
                      p_value = if (nt$fixed) NA_real_ else fit$signif,
                      Permutations = plan$n_perm, Note = nt$text))
  }

  rhs <- if (use_fml) .build_quoted_fixed_rhs(parsed$fml)
         else .quote_terms(c(covariates, v))
  message(sprintf("PERMANOVA: distance ~ %s | strata: %s | n = %d | %d permutations",
                  gsub("`", "", rhs), if (is.null(blocks)) "none" else plan$label,
                  nrow(meta_f), plan$n_perm))
  fit <- suppressMessages(
    vegan::adonis2(stats::as.formula(paste("d ~", rhs)), data = meta_f,
                   permutations = plan$perm, by = "terms"))
  ft  <- as.data.frame(fit)
  ft$Term <- gsub("`", "", rownames(ft))
  # Every fitted term (covariates first, in model order) is reported, so the
  # table itself shows which terms PERMANOVA actually used.
  rows <- which(!ft$Term %in% c("Residual", "Total"))
  if (!use_fml && !v %in% ft$Term)
    stop(sprintf("'%s' not found in adonis2 output.", v), call. = FALSE)
  nt <- notes(ft$Term[rows])

  out <- data.frame(
    Test         = "PERMANOVA (global)",
    Variable     = ft$Term[rows],
    Role         = if (use_fml) "model term"
                   else ifelse(ft$Term[rows] == v, "primary variable", "covariate"),
    Covariates   = if (use_fml) paste("formula:", formula_rhs)
                   else if (length(covariates)) paste(covariates, collapse = " + ")
                   else NA_character_,
    Strata       = plan$label,
    N_samples    = nrow(meta_f),
    Statistic_R2 = round(ft$R2[rows], 4),
    F_value      = round(ft$F[rows], 4),
    p_value      = ifelse(nt$fixed, NA_real_, ft$`Pr(>F)`[rows]),
    Permutations = plan$n_perm,
    Note         = nt$text
  )
  rownames(out) <- NULL
  out
}

# -- 2. run_pairwise_test() ----------------------------------------------------
#' Run pairwise PERMANOVA or ANOSIM with p-value adjustment
#'
#' The adjustment method is taken from \code{obj$config$padj_method}.
#'
#' @inheritParams run_global_test
#' @return A data frame with one row per pair of groups.
#' @export
run_pairwise_test <- function(obj, v,
                              method           = c("permanova", "anosim"),
                              covariates       = NULL,
                              strata           = NULL,
                              permutations     = 999,
                              exclude_baseline = TRUE) {
  method <- match.arg(method)
  cfg    <- obj$config
  .check_var_conflicts(v, covariates, strata)
  meta_f <- .filter_meta(obj$meta, unique(c(v, covariates, strata)),
                         cfg$is_baseline_col, exclude_baseline)
  .check_not_continuous(meta_f, v)

  grp  <- factor(meta_f[[v]])
  if (nlevels(grp) < 2)
    stop(sprintf("'%s' must have at least 2 groups.", v), call. = FALSE)
  rhs  <- .quote_terms(c(covariates, v))

  out <- do.call(rbind, lapply(utils::combn(levels(grp), 2, simplify = FALSE), function(pr) {
    md  <- meta_f[grp %in% pr, , drop = FALSE]
    d   <- .sub_dist(obj, rownames(md))
    plan <- suppressWarnings(.strata_plan(md, strata, permutations))
    if (method == "permanova") {
      fit <- suppressMessages(
        vegan::adonis2(stats::as.formula(paste("d ~", rhs)), data = md,
                       permutations = plan$perm, by = "terms"))
      r <- which(gsub("`", "", rownames(fit)) == v)
      stat <- fit$R2[r]; fval <- fit$F[r]; p <- fit$`Pr(>F)`[r]
    } else {
      fit  <- suppressMessages(vegan::anosim(d, factor(md[[v]]), permutations = plan$perm))
      stat <- fit$statistic; fval <- NA_real_; p <- fit$signif
    }
    # a grouping that never changes within a block cannot be tested within blocks
    if (.constant_within(v, md, plan$blocks)) p <- NA_real_
    data.frame(Group1 = pr[1], Group2 = pr[2], N = nrow(md),
               Statistic = round(stat, 4), F_value = round(fval, 4), p_value = p,
               Strata = plan$label, Permutations = plan$n_perm)
  }))

  out$p_adjusted    <- stats::p.adjust(out$p_value, method = cfg$padj_method)
  out$Test          <- if (method == "permanova") "PERMANOVA (pairwise)" else "ANOSIM (pairwise)"
  out$Variable      <- v
  out$Adjust_method <- cfg$padj_method
  out$Covariates    <- if (length(covariates)) paste(covariates, collapse = " + ") else NA_character_
  out[, c("Test", "Variable", "Group1", "Group2", "N", "Statistic", "F_value",
          "p_value", "p_adjusted", "Adjust_method", "Covariates", "Strata",
          "Permutations")]
}

# -- Internal: shared model-data preparation for LMM / LM ---------------------
.prep_model_data <- function(meta_enriched, metric, fixed_formula, extra_vars,
                             treatment_ref, treatment_col, is_baseline_col) {
  if (!metric %in% colnames(meta_enriched))
    stop(sprintf("Metric '%s' not found.", metric), call. = FALSE)
  parsed <- .parse_fixed_formula(fixed_formula, meta_enriched)
  md <- meta_enriched[!.is_baseline(meta_enriched, is_baseline_col), , drop = FALSE]

  if (treatment_col %in% parsed$vars_used) {
    md[[treatment_col]] <- factor(md[[treatment_col]])
    if (!is.null(treatment_ref) && treatment_ref %in% levels(md[[treatment_col]]))
      md[[treatment_col]] <- stats::relevel(md[[treatment_col]], ref = treatment_ref)
  }
  vars <- unique(c(metric, parsed$vars_used, extra_vars))
  cc   <- stats::complete.cases(md[, vars, drop = FALSE])
  if (any(!cc))
    message(sprintf("%d row(s) dropped due to NA in model variables.", sum(!cc)))
  list(md = md[cc, , drop = FALSE], rhs = .build_quoted_fixed_rhs(parsed$fml))
}

.coef_table <- function(fit, metric, fml_str, padj_method) {
  ct <- as.data.frame(summary(fit)$coefficients)
  out <- data.frame(Metric = metric, Term = rownames(ct), Estimate = ct[, 1],
                    SE = ct[, 2], stringsAsFactors = FALSE)
  if (ncol(ct) == 5) out$df <- ct[, 3]          # lmerTest: Satterthwaite df
  out$t_value    <- ct[, ncol(ct) - 1]
  out$p_value    <- ct[, ncol(ct)]
  out$p_adjusted <- stats::p.adjust(out$p_value, method = padj_method)
  out$N          <- stats::nobs(fit)            # samples actually used by this model
  out$Formula    <- fml_str
  rownames(out)  <- NULL
  out
}

# -- 3. run_alpha_lmm() --------------------------------------------------------
#' Fit a linear mixed model for one alpha diversity metric
#'
#' Baseline samples are excluded from the model rows. Interactions are
#' supported in \code{fixed_formula}, e.g. \code{"Treatment * Timepoint"}.
#'
#' @param meta_enriched Data frame from \code{pcoa_obj$meta_enriched}.
#' @param metric Alpha diversity column to model (e.g. \code{"Shannon"}).
#' @param fixed_formula Fixed-effects formula without \code{~} or random
#'   effect. Backtick-quote names with spaces.
#' @param random_effect One or more columns, each added as its own random
#'   intercept, e.g. \code{c("CatID", "Litter")} gives
#'   \code{(1|CatID) + (1|Litter)}. Nested designs work as crossed terms as
#'   long as IDs are unique across groups (e.g. litter "L1" is never reused in
#'   another house).
#' @param treatment_ref Reference level of the treatment column.
#' @param treatment_col Treatment column name.
#' @param is_baseline_col Baseline flag column.
#' @param padj_method p-value adjustment method for \code{p_adjusted}.
#' @return A data frame with one row per model term.
#' @export
run_alpha_lmm <- function(meta_enriched,
                          metric          = "Shannon",
                          fixed_formula   = "Treatment + Timepoint + Sequence",
                          random_effect   = "CatID",
                          treatment_ref   = "CON",
                          treatment_col   = "Treatment",
                          is_baseline_col = "IsBaseline",
                          padj_method     = "BH") {
  random_effect <- unique(random_effect[!is.na(random_effect) & nzchar(random_effect)])
  if (!length(random_effect))
    stop("Choose at least one random effect column.", call. = FALSE)
  miss <- setdiff(random_effect, colnames(meta_enriched))
  if (length(miss))
    stop(sprintf("Random effect column(s) not found: %s", paste(miss, collapse = ", ")),
         call. = FALSE)
  both <- intersect(random_effect, .parse_fixed_formula(fixed_formula, meta_enriched)$vars_used)
  if (length(both))
    stop(sprintf("'%s' is in both the fixed effects and the random effects - use it in only one.",
                 paste(both, collapse = "', '")), call. = FALSE)

  p <- .prep_model_data(meta_enriched, metric, fixed_formula, random_effect,
                        treatment_ref, treatment_col, is_baseline_col)
  if (nrow(p$md) < 5) stop("Too few observations to fit an LMM.", call. = FALSE)
  for (re in random_effect) {
    n_lv <- length(unique(p$md[[re]]))
    if (n_lv < 2 || n_lv >= nrow(p$md))
      stop(sprintf(paste0("Random effect '%s' has %d groups for %d samples - it needs at ",
                          "least 2 groups and fewer groups than samples (repeated samples ",
                          "per group)."), re, n_lv, nrow(p$md)), call. = FALSE)
  }

  rand    <- paste(sprintf("(1|`%s`)", random_effect), collapse = " + ")
  fml_str <- sprintf("`%s` ~ %s + %s", metric, p$rhs, rand)
  message("Fitting LMM: ", fml_str)
  fit <- tryCatch(
    lmerTest::lmer(stats::as.formula(fml_str), data = p$md, REML = FALSE),
    error = function(e)
      stop(sprintf("LMM failed for '%s': %s", metric, conditionMessage(e)), call. = FALSE))
  .coef_table(fit, metric, fml_str, padj_method)
}

# -- 4. run_alpha_lmm_all_metrics() --------------------------------------------
#' Fit the LMM for every alpha diversity metric
#'
#' Loops over every metric in \code{ALPHA_METRICS} present in
#' \code{meta_enriched}. With no \code{random_effect} (\code{NULL} or
#' \code{character(0)}) an ordinary linear model is fitted instead, for
#' studies with one sample per animal. With \code{add_baseline = TRUE} each metric's own
#' baseline value (\code{Baseline_<metric>}) is added as the first covariate
#' when available. A cross-metric adjusted p-value is computed for the
#' treatment terms chosen with \code{fdr_target}.
#'
#' @inheritParams run_alpha_lmm
#' @param add_baseline Logical; prepend \code{Baseline_<metric>} when present.
#' @param fdr_target Which treatment terms get the cross-metric adjustment:
#'   \code{"main"} (treatment main effects, default), \code{"interaction"}
#'   (interaction terms involving the treatment) or \code{"auto"} (the
#'   interaction when the model has one, otherwise the main effect). Terms are matched exactly
#'   on the treatment column's levels, so another column whose name merely
#'   starts the same (e.g. \code{Group2} next to \code{Group}) is not picked up.
#' @param fdr_term_pattern Optional Perl regex that overrides
#'   \code{fdr_target} for custom selections.
#' @return One combined data frame, one row per term per metric. Besides
#'   \code{p_adjusted} (within one metric's model) it carries two
#'   cross-metric adjustments of the treatment terms, identical when the
#'   treatment has two groups:
#'   \describe{
#'     \item{p_adj_across_metrics}{each comparison adjusted over the metrics
#'       only (e.g. the MUO-vs-Control p-values of all metrics).}
#'     \item{p_adj_across_metrics_and_groups}{every treatment term of every
#'       metric adjusted as one family. Usually the larger p-value, though
#'       Benjamini-Hochberg is not monotone in family size, so not always.}
#'   }
#' @export
run_alpha_lmm_all_metrics <- function(meta_enriched,
                                      fixed_formula    = "Timepoint + Sequence + Treatment",
                                      random_effect    = "CatID",
                                      treatment_ref    = "CON",
                                      treatment_col    = "Treatment",
                                      is_baseline_col  = "IsBaseline",
                                      padj_method      = "BH",
                                      add_baseline     = TRUE,
                                      fdr_target       = c("main", "interaction", "auto"),
                                      fdr_term_pattern = NULL) {
  fdr_target <- match.arg(fdr_target)
  metrics <- intersect(ALPHA_METRICS, colnames(meta_enriched))
  empty   <- metrics[vapply(metrics, function(m) all(is.na(meta_enriched[[m]])), logical(1))]
  if (length(empty))
    message("Skipping metric(s) with no values (e.g. Chao1 needs raw counts): ",
            paste(empty, collapse = ", "))
  metrics <- setdiff(metrics, empty)
  if (!length(metrics))
    stop("No alpha diversity metrics found. Expected: ",
         paste(ALPHA_METRICS, collapse = ", "), call. = FALSE)
  .fit_models_all(meta_enriched, metrics, fixed_formula, random_effect, treatment_ref,
                  treatment_col, is_baseline_col, padj_method, add_baseline,
                  fdr_target, fdr_term_pattern)
}

# -- Internal engine: one model per response column + cross-response adjustment --
# Shared by alpha diversity (responses = metrics) and targeted measurements
# (responses = the features the user selected). Returns the alpha-style table
# (column "Metric" = response name).
.fit_models_all <- function(meta_enriched, metrics, fixed_formula, random_effect,
                            treatment_ref, treatment_col, is_baseline_col,
                            padj_method = "BH", add_baseline = TRUE,
                            fdr_target = "main", fdr_term_pattern = NULL) {
  random_effect <- unique(random_effect[!is.na(random_effect) & nzchar(random_effect)])
  mixed   <- length(random_effect) > 0
  if (!mixed)
    message("No random effect chosen: fitting an ordinary linear model ",
            "(appropriate when every animal has one sample).")

  out <- do.call(rbind, lapply(metrics, function(metric) {
    bcol <- paste0("Baseline_", metric)
    fml  <- if (add_baseline && bcol %in% colnames(meta_enriched))
      paste(bcol, "+", fixed_formula) else fixed_formula
    tryCatch(
      if (mixed)
        run_alpha_lmm(meta_enriched, metric, fml, random_effect, treatment_ref,
                      treatment_col, is_baseline_col, padj_method)
      else
        cbind(run_alpha_lm_simple(meta_enriched, metric, fml, treatment_ref,
                                  treatment_col, is_baseline_col, padj_method),
              df = NA_real_)[, c("Metric", "Term", "Estimate", "SE", "df", "t_value",
                                 "p_value", "p_adjusted", "N", "Formula")],
      error = function(e) {
        warning(conditionMessage(e), call. = FALSE)
        data.frame(Metric = metric, Term = "ERROR", Estimate = NA_real_,
                   SE = NA_real_, df = NA_real_, t_value = NA_real_,
                   p_value = NA_real_, p_adjusted = NA_real_, N = NA_integer_,
                   Formula = conditionMessage(e))
      })
  }))

  if (!is.null(fdr_term_pattern)) {
    target <- grepl(fdr_term_pattern, out$Term, perl = TRUE)
  } else {
    lv    <- meta_enriched[[treatment_col]][!.is_baseline(meta_enriched, is_baseline_col)]
    tset  <- paste0(treatment_col, unique(as.character(lv[!is.na(lv)])))
    comp  <- strsplit(gsub("`", "", out$Term), ":", fixed = TRUE)
    has_t <- vapply(comp, function(x) any(x %in% tset), logical(1))
    if (fdr_target == "auto")
      fdr_target <- if (any(has_t & lengths(comp) > 1)) "interaction" else "main"
    target <- has_t & (if (fdr_target == "interaction") lengths(comp) > 1 else lengths(comp) == 1)
  }
  target <- target & out$Term != "ERROR" & !is.na(out$p_value)

  # Two adjustments across metrics, because both are defensible and they differ
  # only when the treatment has 3+ groups (i.e. several comparisons per metric):
  #   p_adj_across_metrics             - each comparison is its own family: the
  #     same term (e.g. GroupMUO) is adjusted over the metrics only.
  #   p_adj_across_metrics_and_groups  - one family for everything: every
  #     treatment term of every metric is adjusted together.
  out$p_adj_across_metrics <- NA_real_
  if (any(target))
    out$p_adj_across_metrics[target] <- stats::ave(
      out$p_value[target], out$Term[target],
      FUN = function(p) stats::p.adjust(p, method = padj_method))
  out$p_adj_across_metrics_and_groups <- NA_real_
  out$p_adj_across_metrics_and_groups[target] <-
    stats::p.adjust(out$p_value[target], method = padj_method)

  attr(out, "fdr_target") <- if (is.null(fdr_term_pattern)) fdr_target else "custom"
  attr(out, "fdr_terms")  <- sort(unique(out$Term[target]))
  if (!any(target))
    warning(sprintf("No treatment %s term found for '%s'; cross-metric columns left empty.",
                    if (fdr_target == "interaction") "interaction" else "main-effect",
                    treatment_col), call. = FALSE)
  rownames(out) <- NULL
  out
}

#' Interaction LMMs for all metrics without a baseline covariate
#'
#' Shortcut for \code{run_alpha_lmm_all_metrics(add_baseline = FALSE,
#' fdr_target = "interaction")}. Kept so existing scripts keep working.
#'
#' @inheritParams run_alpha_lmm_all_metrics
#' @return See \code{run_alpha_lmm_all_metrics()}.
#' @export
run_alpha_lmm_interaction_all <- function(meta_enriched,
                                          fixed_formula    = "Treatment * Timepoint",
                                          random_effect    = "CatID",
                                          treatment_ref    = "CON",
                                          treatment_col    = "Treatment",
                                          is_baseline_col  = "IsBaseline",
                                          padj_method      = "BH",
                                          fdr_term_pattern = NULL) {
  run_alpha_lmm_all_metrics(meta_enriched, fixed_formula, random_effect,
                            treatment_ref, treatment_col, is_baseline_col,
                            padj_method, add_baseline = FALSE,
                            fdr_target = "interaction",
                            fdr_term_pattern = fdr_term_pattern)
}

# -- 5. run_alpha_lm_simple() -------------------------------------------------
#' Linear model for alpha diversity (no random effects)
#'
#' Fallback when each subject has a single observation.
#'
#' @inheritParams run_alpha_lmm
#' @return A data frame with one row per model term.
#' @export
run_alpha_lm_simple <- function(meta_enriched,
                                metric          = "Shannon",
                                fixed_formula   = "Treatment + Timepoint + Sequence",
                                treatment_ref   = "CON",
                                treatment_col   = "Treatment",
                                is_baseline_col = "IsBaseline",
                                padj_method     = "BH") {
  p <- .prep_model_data(meta_enriched, metric, fixed_formula, NULL,
                        treatment_ref, treatment_col, is_baseline_col)
  if (nrow(p$md) < 5) stop("Too few observations to fit a linear model.", call. = FALSE)
  fml_str <- sprintf("`%s` ~ %s", metric, p$rhs)
  message("Fitting LM: ", fml_str)
  .coef_table(stats::lm(stats::as.formula(fml_str), data = p$md),
              metric, fml_str, padj_method)
}
