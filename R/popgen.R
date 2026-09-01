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
