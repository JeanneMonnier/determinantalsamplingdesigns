
#' periodicdsd
#'
#' @description \code{periodicdsd} constructs a complex matrix, kernel of a determinantal sampling designs such that first order inclusion probabilities are equal and the sampling design presents periodic behaviour. It is also possible to return a matrix \code{V} such that \code{VV*} is the kernel matrix.
#'
#' @param n a positive number, the number of items to be drawn, such that \code{n <= N}.
#' @param r a positive number, such that \code{r} and \code{N} are two relatively prime integers
#' @param N a positive number, the number of items to choose from.
#' @param K a logical indicating which model of matrix should be returned. See 'Details.'
#'
#' @details \code{periodicdsd}, for periodic determinantal sampling design, constructs a toeplitz matrix constructed upon primitive Nth roots of the unity. This matrix then presents periodic behaviour.
#'
#' @details Given that the first order inclusion probabilities are equal, these are equal \code{n/N}.
#'
#' @return For \code{K = TRUE}, the kernel matrix of the determinantal sampling design.
#' @return For \code{K = FALSE}, the matrix \code{V} such that \code{VV*} is the kernel matrix of the determinantal sampling design.
#'
#' @importFrom rlang abort
#' @importFrom numbers coprime
#' @export
#'
#' @references V., Loonis, X., Mary, (2017) "Determinantal Sampling Designs".
#'
#' @examples
#' periodicdsd(4, 2, 7, TRUE)
#' #returns a matrix K, kernel of the sampling design corresponding
#'
#' periodicdsd(15, 1, 365)
#' #returns a matrix V
#'
#' @seealso \link{lavancier}
#'
periodicdsd <- function(n, r, N, K = FALSE) {

  if (N < n) {
    rlang::abort("N should be bigger than n")
  }

  if (N <= r) {
    rlang::abort("N should be bigger than r")
  }

  if (!numbers::coprime(r, N)) {
    rlang::abort("r and N should be coprime integers")
  }

  V <- 1 / sqrt(N) * .periodicdsd_inter(n, r, N, N)

  if (!K){
    return(V)
  }
  else{
    return(V %*% t(Conj(V)))
  }
}

#' periodicdsd_inter
#'
#' This function is called by \code{periodicdsd}.
#'
#' @param n a positive number, the number of items to be drawn, such that \code{n <= N}.
#' @param r a positive number, such that \code{r} and \code{N} are two relatively prime integers
#' @param N a positive number, the number of items to choose from.
#' @param i a positive integer, between 0 and n. It is the increment of the recursion.
#'
#' @export
#'
#' @return The ith step of the construction of the matrix v.
#'
#'
#'
.periodicdsd_inter <- function(n, r, N, i) {

  if (i == 1) {

    z <- complex(modulus = 1, argument = 2*pi*r/N)
    v <- matrix(1, nrow = 1, ncol = n)

    if (n > 2) {
      for (k in 2:n) {
        v[1,k] <- v[1, (k-1)] * z
      }
    }

    return(v)
  }

  else {

    R <- .periodicdsd_inter(n, r, N, i-1)
    return(rbind(R, R[1, ] ^ i))

  }
}
