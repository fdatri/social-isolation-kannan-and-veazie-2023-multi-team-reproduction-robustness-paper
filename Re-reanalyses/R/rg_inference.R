#Random-groups inference####
library(data.table)

#' Partition units into s random-groups
#'
#' The random-groups method (interpenetrating subsamples, Mahalanobis 1946;
#' modern treatment in Wolter 2007, ch. 2): a random, roughly-equal-size,
#' disjoint partition of the unit identifiers into s groups. If
#' \code{balance} is given, units are ordered on the balancing variables
#' before assignment (a serpentine allocation), so each group is closer to a
#' miniature of the design instead of a purely random split.
#'
#' @param units vector of unique unit identifiers (e.g. respondents or PSUs).
#' @param s number of groups.
#' @param balance optional list of vectors, same length and order as \code{units},
#'   to allocate groups so that they are balanced on these variables.
#' @return An integer vector, same length as \code{units}, with the group index
#'   (1 to s) of each unit.
rg_partition <- function(units, s, balance = NULL) {
  n <- length(units)
  r <- runif(n)
  ord <- if (is.null(balance)) order(r) else
    do.call(order, c(as.list(balance), list(r)))
  grp <- integer(n)
  grp[ord] <- rep_len(c(seq_len(s), rev(seq_len(s))), n)
  grp
}

#' Bind a list of estimate vectors into a matrix, keeping only valid draws
#'
#' Internal helper for \code{rg_inference}: drops failed fits (non-numeric or
#' unnamed results) and row-binds the rest, aligning columns by name so that a
#' fit that dropped a term (e.g. a missing factor level in a small group)
#' doesn't break the alignment.
#'
#' @param x a list of named numeric vectors, one per group (possibly with
#'   \code{NULL} entries for failed fits).
#' @return A matrix with one row per valid group and one column per estimand
#'   name seen anywhere in \code{x}, or \code{NULL} if nothing was valid.
rg_bind <- function(x) {
  x <- Filter(function(z) is.numeric(z) && length(z) && !is.null(names(z)), x)
  if (!length(x)) return(NULL)
  nm <- unique(unlist(lapply(x, names)))
  do.call(rbind, lapply(x, function(z) z[nm]))
}

#' Random-groups variance for an arbitrary estimator
#'
#' Design-based variance estimation that works for any estimator : partition
#' the sample into s disjoint random-groups, recompute \code{estimate} inside
#' each one, and read the variance off the dispersion between groups. The
#' point estimate reported is the full-sample one (passed in as
#' \code{theta_full}, or computed here if missing); only the groups are used
#' for the variance, which avoids the small-sample bias of averaging the
#' group-level estimates. Standard errors use t-quantiles with \code{s - 1}
#' degrees of freedom (or \code{n_valid - 1} if some groups failed).
#'
#' Using the full-sample estimate as the point estimate, with the between-group
#' dispersion only as its variance, is Wolter's v1 (2007, sec. 2.2, eq. 2.2.2 as
#' applied to \eqn{\hat\theta} , justified by \eqn{\hat\theta}
#' being the more efficient estimator for nonlinear statistics.
#' Averaging the group estimates instead, as \code{test_on = "groups"} does, is the more
#' bias-prone route for our post-stratified estimands, since a group only
#' has m = n/s subjects to estimate the model's variance components from.
#'
#' @param data data.table with all rows the estimator needs.
#' @param estimate function(d) -> named numeric vector. Called once on the full
#'   data (unless \code{theta_full} is supplied) and once per random-group. Should
#'   return the same set of names every time it succeeds; use
#'   \code{tryCatch(..., error = function(e) NULL)} internally if a group can fail
#'   to fit (e.g. a missing level) - failures are dropped automatically.
#' @param unit column name identifying the partitioning unit (e.g. \code{"subj"}
#'   or \code{"psu"}). This must be the level at which the design clusters, or the
#'   variance recovered will understate the true one.
#' @param s number of random-groups. Ten is the classical default (Wolter);
#'   higher is fine and often cheaper once the model itself is fast, as long as
#'   each group still has enough data for the model to behave (variance
#'   components, rare factor levels).
#' @param theta_full named numeric vector: the full-sample point estimate,
#'   normally from a fit you already have. If \code{NULL}, computed by calling
#'   \code{estimate} on all of \code{data}.
#' @param balance optional balancing variables passed to \code{rg_partition}.
#' @param weights optional name of a weight column to renormalise within each
#'   group before calling \code{estimate} (to \code{sum(w) == nrow(group))});
#'   necessary in GLMMs, where the scale of the weights affects the
#'   likelihood/penalty balance and not just the point estimate.
#' @param prepare optional function(d) -> d, applied after weight
#'   renormalisation, for anything else that needs recomputing inside each
#'   group (e.g. \code{droplevels}, group-specific poststratification shares).
#' @param B number of repeated partitions, averaged to stabilise the variance
#'   estimate (reduces the noise from which particular partition was drawn,
#'   not the true sampling variance). Degrees of freedom stay at \code{s - 1}
#'   (or \code{n_valid - 1}), since repeated partitions share the same data
#'   and are not independent draws.
#' @param null_value value to test against for the t-statistic and p-value.
#' @param conf_level confidence level for the interval.
#' @param min_valid minimum fraction of groups (of \code{s}) that must return a
#'   finite value for a given estimand, or it is dropped.
#' @param test_on \code{"full"} (default) centres the interval and test on
#'   \code{theta_full}; \code{"groups"} centres on the mean of the group estimates
#'   instead.
#' @param seed optional seed, set once before partitioning (and before any
#'   repeated partitions in \code{B > 1}).
#' @param cl optional parallel cluster (as returned by
#'   \code{parallel::makeCluster}) to fit the groups in parallel via
#'   \code{parLapplyLB}.
#' @param n_workers Number of local PSOCK workers when cl is NULL.
#' @param worker_objects Named list of estimator dependencies to export to workers.
#' @param full_models Optional fitted model or named list passed to estimate(models, data, ...).
#' @param partition_files One RDS per partition, with assignment (unit and group columns)
#'   and fits (one entry per group, with a models element). Files are read sequentially.
#' @param estimate_args Additional named arguments for the post-fit estimator.
#' @param keep_groups if \code{TRUE}, attaches the per-group estimates as an
#'   attribute (\code{"groups"}) for diagnostics (e.g. checking whether the
#'   variance estimate is stable, or whether groups are failing).
#' @return A data.table with one row per estimand: \code{term, est_full,
#'   est_groups, se, df, n_valid, n_part, ci_lo, ci_hi, t, p}.
rg_inference <- function(data, estimate, unit, s = 10, theta_full = NULL,
                         balance = NULL, weights = NULL, prepare = NULL,
                         B = 1, null_value = 0, conf_level = 0.95,
                         min_valid = 2 / 3, test_on = c("full", "groups"),
                         seed = NULL, cl = NULL, keep_groups = TRUE,
                         n_workers = 1L, worker_objects = list(),
                         full_models = NULL, partition_files = NULL, estimate_args = list()) {
  test_on <- match.arg(test_on)
  data <- as.data.table(data)
  stopifnot(unit %in% names(data))
  if (is.null(partition_files)) {
    stopifnot(is.null(full_models), !length(estimate_args))
    if (!is.null(seed)) set.seed(seed)
    #u, units: length n, one entry per row / one per unique partitioning unit
    u <- as.character(data[[unit]])
    units <- unique(u)
    bal <- if (is.null(balance)) NULL else
      as.list(data[match(units, u), ..balance])

    run <- function(d) {
      tryCatch({
        d <- data.table::copy(d)
        if (!is.null(weights)) d[, (weights) := get(weights) * .N / sum(get(weights))]
        if (!is.null(prepare)) d <- prepare(d)
        estimate(d)
      }, error = function(e) NULL)
    }
    #Ship only the estimator and its settings in the worker closure.
    environment(run) <- list2env(list(estimate = estimate, prepare = prepare, weights = weights),
                                 parent = environment(rg_inference))

    #theta_full: named vector, one entry per estimand
    if (is.null(theta_full)) theta_full <- run(data)

    #jobs: list of length B*s, each element a data subset for one (partition, group) pair
    jobs <- unlist(lapply(seq_len(B), function(b) {
      g <- rg_partition(units, s, bal)
      split(data, factor(g[match(u, units)], levels = seq_len(s)))
    }), recursive = FALSE)
    stopifnot(length(n_workers) == 1L, is.finite(n_workers), n_workers >= 1,
              n_workers == as.integer(n_workers))
    if (is.null(cl) && n_workers > 1L) {
      cl <- parallel::makePSOCKcluster(min(n_workers, length(jobs)))
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, library(data.table))
    }
    if (!is.null(cl) && length(worker_objects)) {
      parallel::clusterExport(cl, names(worker_objects), envir = list2env(worker_objects))
    }
    res <- if (is.null(cl)) lapply(jobs, run) else parallel::parLapplyLB(cl, jobs, run)
  } else {
    #Use the saved assignments and refit-free estimates for each partition
    stopifnot(length(partition_files) > 0L, !is.null(full_models),
              is.null(weights), is.null(prepare), is.null(balance))
    if (!missing(B)) stopifnot(B == length(partition_files))
    B <- length(partition_files)
    check_s <- !missing(s)
    if (is.null(theta_full))
      theta_full <- do.call(estimate, c(list(full_models, data), estimate_args))
    stopifnot(is.numeric(theta_full), length(theta_full) > 0L,
              !is.null(names(theta_full)), !anyDuplicated(names(theta_full)))
    res <- list()
    for (b in seq_along(partition_files)) {
      partition <- readRDS(partition_files[b])
      assignment <- partition$assignment
      if (b == 1L && !check_s) s <- length(partition$fits)
      stopifnot(length(partition$fits) == s, unit %in% names(assignment),
                !anyDuplicated(assignment[[unit]]),
                setequal(assignment[[unit]], data[[unit]]),
                setequal(assignment$group, seq_len(s)))
      group_id <- assignment$group[match(data[[unit]], assignment[[unit]])]
      for (g in seq_len(s)) {
        d <- droplevels(data[which(group_id == g)])
        z <- tryCatch(do.call(estimate, c(list(partition$fits[[g]]$models, d), estimate_args)),
                      error = function(e) {
                        warning(sprintf("Partition %d, group %d: %s", b, g, conditionMessage(e)))
                        NULL
                      })
        res[length(res) + 1L] <- list(z)
      }
      rm(partition)
    }
  }

  #th: list of length B, each element an (up to s) x (n estimands) matrix from rg_bind
  th <- lapply(split(res, rep(seq_len(B), each = s)), rg_bind)
  th <- Filter(Negate(is.null), th)
  if (!length(th)) return(NULL)
  #nm: length = number of distinct estimand names seen across all groups/partitions
  nm <- unique(unlist(lapply(th, colnames)))
  
  pick <- function(m, k) if (k %in% colnames(m)) m[, k] else NA_real_
  #out: one row per estimand in nm
  out <- rbindlist(lapply(nm, function(k) {
    #v: 3 x B, one column per partition (mean group variance, mean group estimate, mean group count)
    v <- vapply(th, function(m) {
      z <- pick(m, k); z <- z[is.finite(z)]
      if (length(z) < max(2, ceiling(min_valid * s))) c(NA, NA, NA)
      else c(var(z) / length(z), mean(z), length(z))
    }, numeric(3))
    v <- v[, is.finite(v[1, ]), drop = FALSE]
    if (!ncol(v)) return(NULL)
    se <- sqrt(mean(v[1, ]))
    nv <- mean(v[3, ])
    df <- round(nv) - 1
    est <- if (test_on == "groups" || is.null(theta_full[k]) ||
               !k %in% names(theta_full)) mean(v[2, ]) else unname(theta_full[k])
    tq <- qt(1 - (1 - conf_level) / 2, df)
    data.table(term = k, est_full = unname(theta_full[k]), est_groups = mean(v[2, ]),
               se = se, df = df, n_valid = nv, n_part = ncol(v),
               ci_lo = est - tq * se, ci_hi = est + tq * se,
               t = (est - null_value) / se,
               p = 2 * pt(-abs((est - null_value) / se), df))
  }))
  if (keep_groups) {
    attr(out, "group_details") <- lapply(res, attr, which = "details", exact = TRUE)
    attr(out, "groups") <- th
    attr(out, "s") <- s
    attr(out, "B") <- B
  }
  out[]
}
