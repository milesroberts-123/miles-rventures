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

#' Deprecated: Wright-Fisher simulation under selection
#'
#' [WF.sel()] is deprecated; use [WF_sel()] instead. Behaves identically.
#'
#' @inheritParams WF_sel
#'
#' @return See [WF_sel()].
#' @export
#' @keywords internal
WF.sel <- function(N, q, h, s, G) {
  warning("WF.sel() is deprecated; use WF_sel() instead.", call. = FALSE)
  WF_sel(N = N, q = q, h = h, s = s, G = G)
}
