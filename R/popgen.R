#' Allele frequency after one generation of selection
#'
#' Computes the frequency of allele A after one generation of selection at a
#' single biallelic locus, assuming random mating and a fitness of
#' `1 + h*s` for heterozygotes and `1 + s` for AA homozygotes.
#'
#' @param q Current frequency of allele A.
#' @param h Dominance coefficient of allele A.
#' @param s Selection coefficient of allele A.
#'
#' @return The allele frequency of A in the next generation.
#' @export
#'
#' @examples
#' fitfreq(q = 0.1, h = 0.5, s = 0.01)
fitfreq <- function(q, h, s) {
  p <- 1 - q
  wbar <- (1) * p^2 + (1 + h * s) * 2 * p * q + (1 + s) * q^2
  (q^2 * (1 + s) + p * q * (1 + h * s)) / wbar
}

#' Allele frequency trajectory over time: apply selection, then drift comes from rbinom
#'
#' Simulates a Wright-Fisher trajectory of a biallelic locus under selection
#' and genetic drift. Each generation, the deterministic allele frequency
#' after selection is computed with [fitfreq()], then the number of A alleles
#' in the next generation is drawn from a binomial distribution.
#'
#' @param N Population size (number of diploid individuals).
#' @param q Initial frequency of allele A.
#' @param h Dominance coefficient of allele A.
#' @param s Selection coefficient of allele A. Either a scalar (recycled for
#'   every generation) or a vector of length `G - 1` giving the selection
#'   coefficient applied in each generation.
#' @param G Number of generations to simulate.
#'
#' @return A numeric vector of length `G` with the frequency of allele A in
#'   each generation.
#' @export
#'
#' @examples
#' WF_sel(N = 100, q = 0.1, h = 0.5, s = 0.01, G = 50)
WF_sel <- function(N, q, h, s, G) {
  if (length(s) == 1) {
    s <- rep(s, G - 1)
  } else if (length(s) != G - 1) {
    stop("s must be a scalar or a vector of length G - 1.")
  }

  counts <- numeric(G)
  counts[1] <- N * q
  for (i in 2:G) {
    counts[i] <- stats::rbinom(1, N, fitfreq(counts[i - 1] / N, h, s[i - 1]))
  }
  counts / N
}

#' Standardized variance in allele-frequency change
#'
#' Estimates Fc, the standardized variance in allele-frequency change between
#' two time points (Waples 1989): the sum of squared allele-frequency shifts
#' across loci, normalized by the sum of `p0 * (1 - p0)` weighting terms.
#'
#' @param p0 Numeric vector of initial allele frequencies, between 0 and 1.
#' @param pt Numeric vector of allele frequencies after `t` generations,
#'   same length as `p0`.
#'
#' @return A single numeric value, the standardized variance in allele
#'   frequency change.
#' @export
#'
#' @examples
#' fc(p0 = 0.5, pt = 0.25)
fc <- function(p0, pt) {
  stopifnot(
    length(p0) == length(pt),
    all(p0 >= 0, na.rm = TRUE),
    all(pt >= 0, na.rm = TRUE),
    all(p0 <= 1, na.rm = TRUE),
    all(pt <= 1, na.rm = TRUE)
  )

  fsum_num <- (p0 - pt)^2
  fsum_denom <- p0 * (1 - p0)

  sum(fsum_num, na.rm = TRUE) / sum(fsum_denom, na.rm = TRUE)
}

#' Estimate Ne from the temporal method, corrected for selfing
#'
#' Applies the Waples (1989) correction for biased N/Ne estimates to convert
#' an estimate of N/Ne derived from the selfing rate into an estimate of the
#' effective population size Ne. Takes lower and upper bounds of the N/Ne
#' ratio and returns the corresponding bounds of Ne, using the standardized
#' variance in allele-frequency change computed with [fc()].
#'
#' @param r_low Lower bound of the N/Ne ratio estimated from the selfing rate.
#' @param r_high Upper bound of the N/Ne ratio estimated from the selfing rate.
#' @param p0 Numeric vector of initial allele frequencies, between 0 and 1,
#'   passed to [fc()].
#' @param pt Numeric vector of final allele frequencies, between 0 and 1,
#'   passed to [fc()].
#' @param S0 Sample size (number of diploid individuals) at time 0.
#' @param St Sample size (number of diploid individuals) at time t.
#' @param t Number of generations between the two samples.
#'
#' @return A numeric vector of length 2: Ne estimated with `r_low` and with
#'   `r_high`.
#' @export
#'
#' @examples
#' waples_ne(r_low = 2, r_high = 2, p0 = 0.5, pt = 0.25, S0 = 20, St = 20, t = 1)
waples_ne <- function(r_low, r_high, p0, pt, S0, St, t) {
  r <- c(r_low, r_high)
  Fc <- fc(p0, pt)
  ne <- (r * t - 2) / (2 * r * (Fc - 1 / S0 - 1 / St))
  return(ne)
}

#' Expected linkage disequilibrium between loci a given distance apart
#'
#' Computes the expected r^2 between two loci separated by distance `d`,
#' accounting for finite sample size and population-scaled recombination
#' rate (Hill and Weir 1988).
#'
#' @param d Distance between two loci in bp.
#' @param n Sample size.
#' @param C Population-scaled recombination rate: 4Nc.
#'
#' @return The expected r^2.
#' @export
#'
#' @examples
#' hill_weir_r2(d = 1000, n = 30, C = 0.001)
#' hill_weir_r2(d = c(0, 1000, 10000), n = 30, C = 0.001)
hill_weir_r2 <- function(d, n, C) {
  part1 <- (10 + C * d) / ((2 + C * d) * (11 + C * d))
  part2 <- 1 + ((3 + C * d) * (12 + 12 * C * d + (C * d)^2)) /
    (n * (2 + C * d) * (11 + C * d))
  r2 <- part1 * part2
  return(r2)
}

# ---------------------------------------------------------------------------
# Dataset structures for temporal-replicate allele frequency studies
# ---------------------------------------------------------------------------

#' Allele frequency matrix (L variants x S samples)
#'
#' Constructs a validated allele frequency matrix: rows are variants (loci),
#' columns are samples. Values must lie in \code{[0, 1]} (NA is allowed for
#' missing genotypes). Works with the other dataset objects via
#' [validate_af_dataset()] and [extract_samples()].
#'
#' @param x A numeric matrix with L rows (variants) and S columns (samples).
#'
#' @return A `freq_matrix` object: a numeric matrix of allele frequencies.
#' @export
#'
#' @examples
#' m <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
#' freq_matrix(m)
freq_matrix <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix.")
  }
  if (nrow(x) < 1 || ncol(x) < 1) {
    stop("x must have at least one row and one column.")
  }
  if (any(x < 0, na.rm = TRUE) || any(x > 1, na.rm = TRUE)) {
    stop("All allele frequencies must be between 0 and 1.")
  }
  class(x) <- c("freq_matrix", class(x))
  x
}

#' Variant coordinates (L variants x 2)
#'
#' Constructs a validated variant coordinate table with one row per variant
#' and the columns `CHROM` (character) and `POS` (numeric position). Column
#' names are matched case-insensitively and normalized to uppercase, so the
#' object plugs directly into [rventures::plot_manhattan()].
#'
#' @param x A data.frame with L rows and 2 columns: chromosome and position.
#'
#' @return A `snp_coords` object: a data.frame with columns `CHROM` and `POS`.
#' @export
#'
#' @examples
#' d <- data.frame(chrom = c("1", "1", "2"), pos = c(100, 250, 80))
#' snp_coords(d)
snp_coords <- function(x) {
  if (!is.data.frame(x)) {
    stop("x must be a data.frame.")
  }
  if (ncol(x) != 2) {
    stop("x must have exactly 2 columns: chromosome and position.")
  }
  names(x) <- toupper(names(x))
  if (!all(c("CHROM", "POS") %in% names(x))) {
    stop("x must have columns 'chrom' and 'pos'.")
  }
  x$CHROM <- as.character(x$CHROM)
  x$POS <- as.numeric(x$POS)
  if (any(is.na(x$CHROM)) || any(is.na(x$POS))) {
    stop("CHROM and POS must not contain NA values.")
  }
  if (any(x$POS < 1)) {
    stop("POS must be at least 1.")
  }
  class(x) <- c("snp_coords", class(x))
  x
}

#' Initial allele frequencies (L variants)
#'
#' Constructs a validated vector of initial allele frequencies, one per
#' variant. A single-column matrix is accepted and coerced to a vector.
#'
#' @param x A numeric vector (or single-column matrix) of frequencies in
#'   \code{[0, 1]}.
#'
#' @return A `p0_vec` object: a numeric vector of initial allele frequencies.
#' @export
#'
#' @examples
#' p0_vec(c(0.1, 0.5, 0.9))
p0_vec <- function(x) {
  if (is.matrix(x) && ncol(x) == 1) {
    x <- as.vector(x)
  }
  if (!is.numeric(x) || is.matrix(x)) {
    stop("x must be a numeric vector.")
  }
  if (length(x) < 1) {
    stop("x must contain at least one frequency.")
  }
  if (any(x < 0, na.rm = TRUE) || any(x > 1, na.rm = TRUE)) {
    stop("All allele frequencies must be between 0 and 1.")
  }
  class(x) <- c("p0_vec", class(x))
  x
}

#' Sample metadata (S samples x 3)
#'
#' Constructs a validated sample metadata table with one row per sample and
#' the columns `population` (character), `time_point` (integer), and
#' `replicate` (character).
#'
#' @param x A data.frame with S rows and 3 columns, in the order population,
#'   time point, replicate (or with those column names in any order).
#'
#' @return A `sample_info` object: a data.frame with columns `population`,
#'   `time_point`, and `replicate`.
#' @export
#'
#' @examples
#' d <- data.frame(
#'   population = c("AA", "AA", "BB"),
#'   time_point = c(0, 10, 0),
#'   replicate = c("R1", "R1", "R1")
#' )
#' sample_info(d)
sample_info <- function(x) {
  if (!is.data.frame(x)) {
    stop("x must be a data.frame.")
  }
  if (ncol(x) != 3) {
    stop("x must have exactly 3 columns: population, time_point, replicate.")
  }
  if (!setequal(names(x), c("population", "time_point", "replicate"))) {
    stop("x must have columns 'population', 'time_point', and 'replicate'.")
  }
  x <- x[, c("population", "time_point", "replicate"), drop = FALSE]
  x$population <- as.character(x$population)
  x$time_point <- as.integer(x$time_point)
  x$replicate <- as.character(x$replicate)
  if (any(is.na(x$population)) || any(is.na(x$time_point)) ||
      any(is.na(x$replicate))) {
    stop("sample metadata must not contain NA values.")
  }
  if (any(x$time_point < 0)) {
    stop("time_point must be non-negative.")
  }
  class(x) <- c("sample_info", class(x))
  x
}

#' Print a freq_matrix
#'
#' @param x A `freq_matrix` object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The object, invisibly.
#' @export
#'
#' @examples
#' m <- freq_matrix(matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2))
#' print(m)
print.freq_matrix <- function(x, ...) {
  cat(sprintf("freq_matrix: %d variants x %d samples\n", nrow(x), ncol(x)))
  if (!is.null(colnames(x))) {
    cat("Samples:", paste(colnames(x), collapse = ", "), "\n")
  }
  print(utils::head(unclass(round(x, 4))), ...)
  invisible(x)
}

#' Print a snp_coords object
#'
#' @param x A `snp_coords` object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The object, invisibly.
#' @export
#'
#' @examples
#' d <- snp_coords(data.frame(chrom = c("1", "2"), pos = c(100, 200)))
#' print(d)
print.snp_coords <- function(x, ...) {
  cat(sprintf("snp_coords: %d variants on %d chromosome(s)\n",
              nrow(x), length(unique(x$CHROM))))
  print(utils::head(unclass(x)), ...)
  invisible(x)
}

#' Print a p0_vec object
#'
#' @param x A `p0_vec` object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The object, invisibly.
#' @export
#'
#' @examples
#' p <- p0_vec(c(0.1, 0.5, 0.9))
#' print(p)
print.p0_vec <- function(x, ...) {
  cat(sprintf("p0_vec: %d initial allele frequencies\n", length(x)))
  print(utils::head(unclass(round(x, 4))), ...)
  invisible(x)
}

#' Print a sample_info object
#'
#' @param x A `sample_info` object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The object, invisibly.
#' @export
#'
#' @examples
#' d <- sample_info(data.frame(
#'   population = c("AA", "BB"), time_point = c(0, 10), replicate = c("R1", "R1")
#' ))
#' print(d)
print.sample_info <- function(x, ...) {
  cat(sprintf("sample_info: %d samples\n", nrow(x)))
  print(utils::head(unclass(x)), ...)
  invisible(x)
}

#' Validate a set of temporal allele frequency dataset objects
#'
#' Cross-checks the four dataset objects for compatible dimensions: the same
#' number of variants (rows of the frequency matrix and coordinates, length
#' of the initial frequencies) and the same number of samples (columns of the
#' frequency matrix, rows of the metadata). Optionally checks that frequency
#' matrix column names follow the `pop_time_rep` naming scheme.
#'
#' @param freq_mat A `freq_matrix` object.
#' @param snp_coords A `snp_coords` object.
#' @param p0 A `p0_vec` object.
#' @param sample_meta A `sample_info` object.
#' @param check_colnames If `TRUE` and `freq_mat` has column names, check
#'   they match the `pop_time_rep` scheme derived from `sample_meta`.
#'
#' @return `TRUE`, invisibly, if all checks pass; an error otherwise.
#' @export
#'
#' @examples
#' m <- matrix(runif(10), nrow = 5)
#' fm <- freq_matrix(m)
#' coords <- snp_coords(data.frame(chrom = rep("1", 5), pos = 1:5 * 100))
#' p0 <- p0_vec(m[, 1])
#' meta <- sample_info(data.frame(
#'   population = c("AA", "BB"), time_point = c(0, 10), replicate = c("R1", "R1")
#' ))
#' validate_af_dataset(fm, coords, p0, meta)
validate_af_dataset <- function(freq_mat, snp_coords, p0, sample_meta,
                                check_colnames = TRUE) {
  stopifnot(
    inherits(freq_mat, "freq_matrix"),
    inherits(snp_coords, "snp_coords"),
    inherits(p0, "p0_vec"),
    inherits(sample_meta, "sample_info")
  )

  n_variants <- nrow(freq_mat)
  if (nrow(snp_coords) != n_variants || length(p0) != n_variants) {
    stop(sprintf(paste(
      "Variant count mismatch: freq_matrix has %d rows,",
      "snp_coords has %d rows, p0_vec has %d entries."
    ), n_variants, nrow(snp_coords), length(p0)))
  }

  n_samples <- ncol(freq_mat)
  if (nrow(sample_meta) != n_samples) {
    stop(sprintf(paste(
      "Sample count mismatch: freq_matrix has %d columns,",
      "sample_info has %d rows."
    ), n_samples, nrow(sample_meta)))
  }

  if (check_colnames && !is.null(colnames(freq_mat))) {
    expected <- paste(sample_meta$population, sample_meta$time_point,
                      sample_meta$replicate, sep = "_")
    if (!identical(colnames(freq_mat), expected)) {
      stop("freq_matrix colnames do not match the sample_info metadata ",
           "(expected pop_time_rep names in sample order).")
    }
  }

  invisible(TRUE)
}

#' Extract samples from an allele frequency matrix
#'
#' Subsets the columns of a `freq_matrix` by population, time point, and/or
#' replicate (using `%in%` matching against the sample metadata), rebuilds
#' column names as `pop_time_rep`, and attaches the subset metadata as the
#' attribute `"sample_info"`. The rebuilt names let the result feed directly
#' into [covmat_pop_pair()] and provide replicate/time labels for
#' [replicate_gt()].
#'
#' @param freq_mat A `freq_matrix` object.
#' @param sample_meta The `sample_info` object for `freq_mat`.
#' @param population Optional character vector of populations to keep.
#' @param time Optional integer vector of time points to keep.
#' @param replicate Optional character vector of replicates to keep.
#'
#' @return A classed `freq_matrix` with columns in metadata order, column
#'   names `pop_time_rep`, and attribute `"sample_info"` holding the subset
#'   metadata. Errors if no sample matches.
#' @export
#'
#' @examples
#' m <- matrix(runif(8), nrow = 4, ncol = 4)
#' fm <- freq_matrix(m)
#' meta <- sample_info(data.frame(
#'   population = c("AA", "AA", "BB", "BB"),
#'   time_point = c(0, 10, 0, 10),
#'   replicate = c("R1", "R1", "R1", "R1")
#' ))
#' extract_samples(fm, meta, population = "AA")
extract_samples <- function(freq_mat, sample_meta, population = NULL,
                            time = NULL, replicate = NULL) {
  stopifnot(
    inherits(freq_mat, "freq_matrix"),
    inherits(sample_meta, "sample_info"),
    ncol(freq_mat) == nrow(sample_meta)
  )

  keep <- rep(TRUE, nrow(sample_meta))
  if (!is.null(population)) {
    keep <- keep & sample_meta$population %in% population
  }
  if (!is.null(time)) {
    keep <- keep & sample_meta$time_point %in% time
  }
  if (!is.null(replicate)) {
    keep <- keep & sample_meta$replicate %in% replicate
  }
  if (!any(keep)) {
    stop("No samples match the given population/time/replicate filters.")
  }

  sub_meta <- sample_meta[keep, , drop = FALSE]
  out <- freq_mat[, keep, drop = FALSE]
  colnames(out) <- paste(sub_meta$population, sub_meta$time_point,
                         sub_meta$replicate, sep = "_")
  attr(out, "sample_info") <- sub_meta
  class(out) <- c("freq_matrix", class(out))
  out
}

# ---------------------------------------------------------------------------
# Temporal-replicate popgen: frequency changes, covariance, and G statistics
# ---------------------------------------------------------------------------

#' Sum of heterozygosity
#'
#' Computes `sum(2 * x * (1 - x))` over a vector of allele frequencies, one
#' per locus, skipping NAs.
#'
#' @param x Numeric vector of allele frequencies, one per locus.
#'
#' @return A single number: the summed heterozygosity.
#' @export
#'
#' @examples
#' sum_of_het(c(0.5, 0.25, 1))
sum_of_het <- function(x) {
  sum(2 * x * (1 - x), na.rm = TRUE)
}

#' Pairwise geometric mean
#'
#' Computes the mean of `sqrt(x_i * x_j)` over all unordered pairs of
#' elements in a numeric vector.
#'
#' @param my_vector Numeric vector of length at least 2.
#'
#' @return A single number: the mean of the pairwise geometric means.
#' @export
#'
#' @examples
#' geom_pairwise_mean(c(1, 4, 9))
geom_pairwise_mean <- function(my_vector) {
  stopifnot(length(my_vector) > 1, is.numeric(my_vector))
  n <- length(my_vector)
  var_sum <- 0
  for (i in 1:(n - 1)) {
    var_a <- my_vector[i]
    for (j in (i + 1):n) {
      var_b <- my_vector[j]
      var_sum <- var_sum + sqrt(var_a * var_b)
    }
  }
  var_mean <- var_sum / ((n^2 - n) / 2)
  return(var_mean)
}

#' Allele frequency changes between adjacent generations
#'
#' Computes the differences between adjacent columns of an allele frequency
#' matrix: columns are time points in ascending order, rows are variants.
#'
#' @param pmat Numeric matrix of allele frequencies; time points ascending
#'   from the left column to the right column.
#'
#' @return A numeric matrix with the same number of rows and one fewer
#'   column than `pmat`: the allele frequency change in each interval.
#' @export
#'
#' @examples
#' pmat <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.2, 0.3, 0.4, 0.5), nrow = 2, byrow = TRUE)
#' freq_increments(pmat)
freq_increments <- function(pmat) {
  pmat[, -1] - pmat[, -ncol(pmat)]
}

#' Mark entries after first NA as also NA
#'
#' Given a time series of allele frequencies ordered from time 0, 1, 2, ...,
#' t, sets every entry from the first NA onward to NA, truncating the
#' trajectory at its first missing value.
#'
#' @param x Numeric vector assumed to be ordered in time.
#'
#' @return The input vector with all entries from the first NA onward set
#'   to NA. Unchanged if there is no NA.
#' @export
#'
#' @examples
#' rm_na_after_na(c(0.1, 0.2, NA, 0.4, 0.5))
rm_na_after_na <- function(x) {
  y <- which(is.na(x))
  if (length(y) > 0) {
    first_na <- min(y) # get first na in sequence
    x[first_na:length(x)] <- NA # mark instances after first na as also na
  }
  return(x)
}

#' Arcsine-square root transform of allele frequencies
#'
#' Applies the variance-stabilizing arcsine-square root transformation,
#' `asin(sqrt(x))`, commonly used before temporal analyses of allele
#' frequencies.
#'
#' @param x Numeric vector of values in \code{[0, 1]}; NA entries pass
#'   through as NA.
#'
#' @return The transformed numeric vector, in \code{[0, pi/2]}.
#' @export
#'
#' @examples
#' arcsin_sqrt(c(0, 0.5, 1))
arcsin_sqrt <- function(x) {
  stopifnot(
    is.numeric(x),
    all(x >= 0, na.rm = TRUE),
    all(x <= 1, na.rm = TRUE)
  )
  asin(sqrt(x))
}

#' Sign permute a block of allele frequency changes
#'
#' Randomly multiplies allele frequency changes by +1 or -1 at one of several
#' granularities, for use in randomization tests of linked selection.
#' Because covariance is signed, flipping signs breaks real shared-direction
#' structure while preserving per-element magnitudes.
#'
#' @param pdiff Data.frame or matrix of allele frequency changes (rows =
#'   variants, columns = time intervals).
#' @param procedure One of `"none"` (return input unchanged), `"window"`
#'   (one random sign per genomic window, applied to all variants in the
#'   window), `"cell"` (one random sign per matrix entry), or `"genome"`
#'   (one random sign per column). Invalid values raise an error.
#' @param windows Numeric or character vector of window IDs, one per row of
#'   `pdiff`; required when `procedure = "window"`.
#'
#' @return The sign-permuted allele frequency changes: same type as the
#'   input for `"none"`, a data.frame for `"window"` and `"cell"`, and a
#'   matrix for `"genome"`.
#' @export
#'
#' @examples
#' set.seed(1)
#' pdiff <- data.frame(d1 = c(-0.1, 0.2, -0.3, 0.4), d2 = c(0.1, -0.2, 0.3, -0.4))
#' sign_permute_increments(pdiff, procedure = "genome")
sign_permute_increments <- function(pdiff, procedure = "none", windows = NULL) {
  stopifnot(ncol(pdiff) >= 1, nrow(pdiff) >= 1)
  pdiff <- as.data.frame(pdiff)
  if (procedure == "window") {
    # within a window, randomly multiply allele frequency change by +1 or -1
    stopifnot(!is.null(windows), length(windows) == nrow(pdiff))
    pdiff$window <- windows
    pdiff <- dplyr::group_by(pdiff, window)
    pdiff <- dplyr::mutate(
      pdiff, dplyr::across(dplyr::where(is.numeric), ~ .x * sample(c(1, -1), size = 1))
    )
    pdiff <- dplyr::ungroup(pdiff)
    pdiff <- pdiff[, setdiff(names(pdiff), "window"), drop = FALSE]
  } else if (procedure == "cell") {
    # randomly multiply every cell frequency change by +1 or -1
    pdiff <- dplyr::mutate(
      pdiff, dplyr::across(
        dplyr::where(is.numeric), ~ .x * sample(c(1, -1), dplyr::n(), replace = TRUE)
      )
    )
  } else if (procedure == "genome") {
    # randomly multiply entire columns by +1 or -1
    sign_permutations <- sample(c(-1, 1), size = ncol(pdiff), replace = TRUE)
    # Create the diagonal matrix from vector
    P <- diag(sign_permutations)
    # Perform matrix multiplication to flip all increments at once
    pdiff <- as.matrix(pdiff) %*% P
  } else if (procedure == "none") {
    # no flipping
  } else {
    stop("Procedure not valid. Choose from: none, window, cell, genome.")
  }
  return(pdiff)
}

#' Standardize covariance matrix by heterozygosity
#'
#' Divides each variance and covariance of a temporal covariance matrix by
#' the summed heterozygosity of the corresponding interval, so that entries
#' are comparable across intervals with different amounts of standing
#' variation. Covariances are divided by `half_het_sums[min(i, j)]`, the
#' heterozygosity of the earlier of the two intervals.
#'
#' @param pmat Numeric matrix of allele frequencies excluding the last time
#'   point: one column per interval, holding the frequencies at the start of
#'   that interval.
#' @param covmat Square covariance matrix of allele frequency changes, with
#'   one row/column per interval.
#'
#' @return The standardized covariance matrix.
#' @export
#'
#' @examples
#' pmat <- matrix(c(0.5, 0.2, 0.5, 0.2), nrow = 2)
#' covmat <- cov(t(pmat))
#' standard_cov_by_het(pmat, covmat)
standard_cov_by_het <- function(pmat, covmat) {
  stopifnot(ncol(pmat) == ncol(covmat))
  half_het_sums <- 0.5 * apply(pmat, MARGIN = 2, FUN = sum_of_het)
  stopifnot(all(half_het_sums > 0))
  # standaridize variances
  diag(covmat) <- diag(covmat) / half_het_sums
  # standardize covariances
  for (i in 1:nrow(covmat)) {
    for (j in 1:ncol(covmat)) {
      if (i != j) {
        covmat[i, j] <- covmat[i, j] / half_het_sums[min(c(i, j))]
      }
    }
  }
  return(covmat)
}

#' Calculate covariances from allele frequency matrix
#'
#' Computes the covariance matrix of allele frequency changes between
#' adjacent generations (columns = time intervals), optionally with a sample
#' size correction, heterozygosity standardization, arcsine-transformed
#' input, and sign randomization for null comparisons.
#'
#' @param pmat Numeric matrix of allele frequencies; the first column is the
#'   initial generation and time points ascend left to right.
#' @param n Optional numeric vector of sample sizes (number of chromosomes
#'   sampled per time point); required if `correct_for_n = TRUE`. Must have
#'   one entry per time point (i.e., `ncol(pmat)`).
#' @param correct_for_n Boolean; whether to correct covariances and variances
#'   for finite sample size.
#' @param standard_by_het Boolean; whether to standardize the covariance
#'   matrix by heterozygosity via [standard_cov_by_het()].
#' @param input_asin_trans Boolean; whether `pmat` holds arcsine-square-root
#'   transformed frequencies (changes the sample size correction formulas).
#' @param procedure Sign randomization procedure passed to
#'   [sign_permute_increments()]; use `"none"` (default) for no
#'   randomization.
#' @param windows Optional vector of window IDs, one per row of `pmat`,
#'   passed to [sign_permute_increments()] when `procedure = "window"`.
#'
#' @return The (optionally corrected/standardized) covariance matrix of
#'   allele frequency changes.
#' @export
#'
#' @examples
#' set.seed(1)
#' pmat <- matrix(runif(15), nrow = 5, ncol = 3)
#' covmat_from_pmat(pmat, n = c(50, 50, 50))
covmat_from_pmat <- function(pmat, n = NULL, correct_for_n = TRUE,
                             standard_by_het = FALSE, input_asin_trans = FALSE,
                             procedure = "none", windows = NULL) {
  # allele frequency changes between adjacent generations
  pdiff <- freq_increments(pmat)
  # randomly swap signs of allele frequency changes
  pdiff <- sign_permute_increments(pdiff, procedure = procedure, windows = windows)
  # ensure pdiff is a matrix even with a single increment (2-column pmat)
  if (is.null(dim(pdiff))) {
    pdiff <- matrix(pdiff, ncol = 1)
  }
  # covariance matrices
  covmat <- stats::cov(pdiff, use = "pairwise.complete.obs")
  # correct for sample size, if needed
  if (correct_for_n) {
    if (any(n < 2)) {
      stop("All n must be >= 2")
    }
    if (length(n) != (ncol(covmat) + 1)) {
      stop("Should be a sample size for every time point.")
    }
    ### ARCSIN SQRT TRANSFORMED FREQUENCIES ###
    if (input_asin_trans) {
      # check transformed data are actually input
      if (all(pmat[, 1] < pi)) {
        warning("Are you sure values are arcsin transformed?")
      }
      # correct overlapping covariances
      # seq_len() yields integer(0) when ncol(covmat) == 1, so this is skipped
      for (i in seq_len(ncol(covmat) - 1)) {
        corrected_cov <- covmat[i, i + 1] + (1 / n[i + 1])
        covmat[i, i + 1] <- corrected_cov
        covmat[i + 1, i] <- corrected_cov
      }
      # correct variances
      for (i in seq_len(ncol(covmat))) {
        corrected_var <- covmat[i, i] - (1 / n[i]) - (1 / n[i + 1])
        if (corrected_var < 0) {
          warning("Sample size correction makes variance negative. Setting variance to zero")
          covmat[i, i] <- 0
        } else {
          covmat[i, i] <- corrected_var
        }
      }
      ### RAW FREQUENCIES ###
    } else {
      # correct overlapping covariances
      # seq_len() yields integer(0) when ncol(covmat) == 1, so this is skipped
      for (i in seq_len(ncol(covmat) - 1)) {
        corrected_cov <- mean(pdiff[, i] * pdiff[, i + 1], na.rm = TRUE) +
          mean(pmat[, i + 1] * (1 - pmat[, i + 1]) / (n[i + 1] - 1), na.rm = TRUE)
        covmat[i, i + 1] <- corrected_cov
        covmat[i + 1, i] <- corrected_cov
      }
      # correct variances
      for (i in seq_len(ncol(covmat))) {
        corrected_var <- mean((pdiff[, i])^2, na.rm = TRUE) -
          mean(pmat[, i] * (1 - pmat[, i]) / (n[i] - 1), na.rm = TRUE) -
          mean(pmat[, i + 1] * (1 - pmat[, i + 1]) / (n[i + 1] - 1), na.rm = TRUE)
        if (corrected_var < 0 | is.na(corrected_var)) {
          warning("Sample size correction makes variance negative. Setting variance to zero")
          covmat[i, i] <- 0
        } else {
          covmat[i, i] <- corrected_var
        }
      }
    }
  }

  if (standard_by_het) {
    covmat <- standard_cov_by_het(pmat[, -ncol(pmat), drop = FALSE], covmat)
  }

  return(covmat)
}

#' Rolling sum of matrix, starting from top left
#'
#' Sums the elements of nested top-left squares of a square matrix: element
#' 1 is `mat[1, 1]`, element 2 is `sum(mat[1:2, 1:2])`, and so on up to the
#' full matrix sum.
#'
#' @param mat A square numeric matrix.
#'
#' @return A numeric vector of length `ncol(mat)`: the sums of the nested
#'   top-left squares.
#' @export
#'
#' @examples
#' rolling_matrix_sum(matrix(1:9, nrow = 3))
rolling_matrix_sum <- function(mat) {
  stopifnot(nrow(mat) == ncol(mat))
  rolling_sum <- c()
  for (i in 1:ncol(mat)) {
    rolling_sum <- c(rolling_sum, sum(mat[1:i, 1:i]))
  }
  stopifnot(length(rolling_sum) == ncol(mat))
  return(rolling_sum)
}

#' G statistic, corrected for selection increasing variance in allele
#' frequency change
#'
#' Computes G' (Buffalo and Coop style): 1 minus the expected drift-only
#' variance `times * E[p0(1-p0)] / (2N)`, scaled by the observed rolling sums
#' of the (optionally absolute-valued) temporal covariance matrix. Values
#' near 1 indicate covariance consistent with linked selection; values near
#' or below 0 indicate drift-like, uncorrelated change.
#'
#' @param times Numeric vector of time points, one per interval of `pmat`
#'   (length `ncol(pmat) - 1`), used to order the rolling sums.
#' @param pmat Numeric matrix of allele frequencies, rows = variants,
#'   columns = time points ascending; the first column is the initial
#'   generation.
#' @param N Effective population size.
#' @param n Numeric vector of sample sizes (chromosomes) per time point,
#'   passed to [covmat_from_pmat()] for the sample size correction.
#' @param take_abs Boolean; whether to accumulate the absolute values of the
#'   covariance matrix instead of the raw values.
#'
#' @return A numeric vector, one entry per time point: the corrected G
#'   statistic accumulated up to that time.
#' @export
#'
#' @examples
#' set.seed(1)
#' pmat <- matrix(runif(15), nrow = 5, ncol = 3)
#' g_prime(times = 2:3, pmat = pmat, N = 100, n = c(50, 50, 50))
g_prime <- function(times, pmat, N, n, take_abs = FALSE) {
  stopifnot(all(times >= 0), N > 2)
  covmat <- covmat_from_pmat(pmat, n)
  ep0 <- mean(pmat[, 1] * (1 - pmat[, 1]), na.rm = TRUE)
  if (take_abs) {
    var_roll <- rolling_matrix_sum(abs(covmat))
  } else {
    var_roll <- rolling_matrix_sum(covmat)
  }
  stopifnot(length(times) == length(var_roll))
  return(1 - (times * ep0) / (2 * N * var_roll))
}

#' Calculate linked selection ratios from covariance matrix
#'
#' Accumulates variances (diagonal) and covariances (strict upper triangle)
#' of a temporal covariance matrix over nested top-left squares and derives
#' several G-style ratios of shared to total variance in allele frequency
#' change.
#'
#' @param covmat Square covariance matrix of allele frequency changes, one
#'   row/column per time interval.
#'
#' @return A data.frame with one row per time point (`gen`) and columns for
#'   the G ratios (`gt`, `gt_abs`, `gt_abs_2`, `gt_pos`), the accumulated
#'   sums (`sum_cov`, `sum_abs_cov`, `sum_pos_cov`, `sum_neg_cov`,
#'   `sum_var`), and the positive-covariance fraction (`ratio`).
#' @export
#'
#' @examples
#' covmat <- matrix(c(1, 0.5, 0.2, 0.5, 2, 0.3, 0.2, 0.3, 3), nrow = 3)
#' gt_from_covmat(covmat)
gt_from_covmat <- function(covmat) {
  covmat[lower.tri(covmat)] <- NA

  covariances <- covmat[upper.tri(covmat)]
  variances <- diag(covmat)

  sums_var <- cumsum(variances)

  sums_cov <- c(0)
  sums_abs_cov <- c(0)
  sums_pos_cov <- c(0)
  sums_neg_cov <- c(0)
  for (j in 2:nrow(covmat)) {
    cov_sub <- covariances[1:(((j - 1) * (j)) / 2)]

    sum_abs_cov <- 2 * sum(abs(cov_sub), na.rm = TRUE)
    sum_cov <- 2 * sum(cov_sub, na.rm = TRUE)
    sum_neg_cov <- 2 * sum(abs(cov_sub[(cov_sub < 0)]), na.rm = TRUE)
    sum_pos_cov <- 2 * sum(cov_sub[(cov_sub > 0)], na.rm = TRUE)

    sums_cov <- c(sums_cov, sum_cov)
    sums_abs_cov <- c(sums_abs_cov, sum_abs_cov)
    sums_pos_cov <- c(sums_pos_cov, sum_pos_cov)
    sums_neg_cov <- c(sums_neg_cov, sum_neg_cov)
  }

  result <- data.frame(
    gen = 1:nrow(covmat),
    gt = sums_cov / (sums_var + sums_cov),
    gt_abs = sums_abs_cov / (sums_var + sums_abs_cov),
    gt_abs_2 = sums_abs_cov / (sums_var + sums_cov),
    gt_pos = sums_pos_cov / (sums_var + sums_pos_cov),
    sum_cov = sums_cov,
    sum_abs_cov = sums_abs_cov,
    sum_pos_cov = sums_pos_cov,
    sum_neg_cov = sums_neg_cov,
    sum_var = sums_var,
    ratio = sums_pos_cov / (sums_pos_cov + sums_neg_cov)
  )
  return(result)
}

#' Calculate convergence correlation within environment
#'
#' Estimates the correlation of allele frequency changes among independent
#' replicate populations within the same environment: the mean off-diagonal
#' covariance divided by the geometric mean of the replicate variances.
#'
#' @param pdiff Numeric matrix of allele frequency changes, rows = variants,
#'   columns = replicate time series (all replicates of one site, one time
#'   point), with at least 2 rows and 2 columns.
#'
#' @return A numeric vector of length 2: `c(numerator, denominator)` — the
#'   mean between-replicate covariance and the geometric mean within-replicate
#'   variance.
#' @export
#'
#' @examples
#' set.seed(1)
#' pdiff <- matrix(rnorm(40, sd = 0.05), nrow = 20, ncol = 2)
#' conv_cor_wn_env(pdiff)
conv_cor_wn_env <- function(pdiff) {
  stopifnot(ncol(pdiff) > 1, nrow(pdiff) > 1)
  covmat <- stats::cov(pdiff, use = "pairwise.complete.obs")
  rep_var <- diag(covmat)
  rep_cov <- covmat[upper.tri(covmat)]
  numerator <- sum(2 * rep_cov, na.rm = TRUE) / (2 * length(rep_cov))
  denominator <- geom_pairwise_mean(rep_var)
  stopifnot(denominator > 0, abs(numerator) <= denominator)
  return(c(numerator, denominator))
}

#' Replicate G(t): proportion of variance in allele frequency change due to
#' shared selection pressure
#'
#' Splits the covariance matrix of allele frequency changes into
#' within-replicate variances and between-replicate covariances using the
#' supplied replicate and time labels, and returns means and totals of each
#' (raw, positive-only, and absolute-valued), for estimating the fraction of
#' variance in allele frequency change attributable to shared selection.
#'
#' @param pdiff Numeric matrix of allele frequency changes, rows = variants,
#'   columns = replicate time series.
#' @param rep_labels Vector of replicate labels, one per column of `pdiff`.
#' @param time_labels Vector of time point labels, one per column of
#'   `pdiff`; must include at least 2 distinct values, as must `rep_labels`.
#'
#' @return A one-row data.frame with mean and total (raw, positive, absolute,
#'   and trace) within-replicate variances and between-replicate covariances,
#'   plus the counts of replicates (`n_var`) and replicate pairs
#'   (`n_covar`).
#' @export
#'
#' @examples
#' set.seed(1)
#' pdiff <- matrix(rnorm(60, sd = 0.05), nrow = 15, ncol = 4)
#' replicate_gt(pdiff, rep_labels = c(1, 1, 2, 2), time_labels = c(1, 2, 1, 2))
replicate_gt <- function(pdiff, rep_labels, time_labels) {
  stopifnot(ncol(pdiff) > 1,
            nrow(pdiff) > 1,
            length(rep_labels) == ncol(pdiff),
            length(time_labels) == ncol(pdiff),
            length(unique(rep_labels)) >= 2,
            length(unique(time_labels)) >= 2)

  rep_eq_bool_mat <- outer(rep_labels, rep_labels, "==")

  time_eq_bool_mat <- outer(time_labels, time_labels, "==")

  covmat <- stats::cov(pdiff, use = "pairwise.complete.obs")

  # total variance within replicates
  var_only_mat <- covmat[rep_eq_bool_mat]
  total_var <- sum(var_only_mat)
  total_var_pos <- sum(var_only_mat[(var_only_mat > 0)])
  total_var_abs <- sum(abs(var_only_mat))
  var_trace <- covmat[(rep_eq_bool_mat & time_eq_bool_mat)]
  total_var_trace <- sum(var_trace)
  n_var <- length(unique(rep_labels))

  # total covariance across different replicates AND time points
  cov_only_mat <- covmat[(!rep_eq_bool_mat) & (!time_eq_bool_mat)]
  total_covar <- sum(cov_only_mat)
  total_covar_abs <- sum(abs(cov_only_mat))
  total_covar_pos <- sum(cov_only_mat[(cov_only_mat > 0)])
  n_covar <- (n_var^2 - n_var) / 2

  mean_var <- total_var / n_var
  mean_var_pos <- total_var_pos / n_var
  mean_var_abs <- total_var_abs / n_var
  mean_var_trace <- total_var_trace / n_var
  mean_covar <- total_covar / n_covar
  mean_covar_pos <- total_covar_pos / n_covar
  mean_covar_abs <- total_covar_abs / n_covar

  stopifnot(total_covar <= total_covar_pos,
            mean_covar <= mean_covar_pos,
            total_var_pos <= total_var_abs,
            total_covar < sum(abs(covmat)),
            total_var < sum(abs(covmat)),
            total_covar <= total_covar_abs,
            isTRUE(all.equal(sum(var_only_mat) + sum(covmat[!rep_eq_bool_mat]),
                             sum(covmat))),
            mean_var >= 0,
            total_var >= 0,
            total_var_trace >= 0,
            mean_var_trace >= 0)

  return(data.frame(
    mean_var = mean_var,
    mean_var_pos = mean_var_pos,
    mean_var_abs = mean_var_abs,
    mean_var_trace = mean_var_trace,
    mean_covar = mean_covar,
    mean_covar_pos = mean_covar_pos,
    mean_covar_abs = mean_covar_abs,
    total_var = total_var,
    total_var_pos = total_var_pos,
    total_var_abs = total_var_abs,
    total_var_trace = total_var_trace,
    total_covar = total_covar,
    total_covar_pos = total_covar_pos,
    total_covar_abs = total_covar_abs,
    n_var = n_var,
    n_covar = n_covar
  ))
}

#' Calculate covariances across a pair of environments
#'
#' Computes the mean standardized covariance of allele frequency changes
#' between two populations: the average covariance between the two sets of
#' replicate time series, normalized by the mean pairwise product of
#' within-population standard deviations.
#'
#' @param pmat1 Numeric matrix of allele frequencies in population 1, rows =
#'   variants, columns = replicate time series. Column names must contain
#'   `pop1`.
#' @param pmat2 Numeric matrix of allele frequencies in population 2, same
#'   layout as `pmat1`. Column names must contain `pop2`.
#' @param pop1 Character string identifying population 1; used to grep the
#'   columns of `pmat1`.
#' @param pop2 Character string identifying population 2; used to grep the
#'   columns of `pmat2`.
#' @param gen0 Numeric vector of allele frequencies at generation 0, one per
#'   variant (shared by both populations).
#'
#' @return A single number: the mean standardized between-population
#'   covariance.
#' @export
#'
#' @examples
#' gen0 <- c(0.5, 0.4, 0.3)
#' p1 <- cbind(AA_rep1 = gen0 + 0.01, AA_rep2 = gen0 - 0.01)
#' p2 <- cbind(BB_rep1 = gen0 - 0.01, BB_rep2 = gen0 + 0.01)
#' covmat_pop_pair(p1, p2, pop1 = "AA", pop2 = "BB", gen0 = gen0)
covmat_pop_pair <- function(pmat1, pmat2, pop1, pop2, gen0) {
  stopifnot(
    !is.null(colnames(pmat1)),
    !is.null(colnames(pmat2))
  )
  # allele frequency changes between adjacent generations
  pdiff1 <- pmat1 - gen0
  pdiff2 <- pmat2 - gen0

  covmat <- stats::cov(cbind(pdiff1, pdiff2), use = "pairwise.complete.obs")

  # mean product of std of between population pairs
  rep_var <- diag(covmat)
  rep_var_1 <- rep_var[1:ncol(pmat1)]
  rep_var_2 <- rep_var[(ncol(pmat1) + 1):length(rep_var)]
  mean_pair_std <- mean(sqrt(rep_var_1 %o% rep_var_2))

  bw_pop_cov <- covmat[base::grep(pop1, rownames(covmat)),
                       base::grep(pop2, colnames(covmat))]

  return(mean(bw_pop_cov) / mean_pair_std)
}
