#' lavancier
#'
#' @description \code{lavancier} samples from fixed size determinantal sampling designs.
#'
#' @param v a matrix. See 'Details.'
#' @param s a positive number, the number of samples to take.
#' @param B logical indicating which model of vector should be returned. See 'Details.'
#' @param C a logical vector of length \code{n} (the number of elements sampled and the number of columns of v), indicating for each draw if the probability weights for the next draw should be plotted. See 'Details.'
#'
#' @details \code{v} has to be a matrix of size \code{Nxn} where \code{N} is the number of elements to choose from and \code{n} the number of elements sampled such that \code{vv*} is the kernel matrix of your sampling design. Two functions of this package provides such matrices. The function loonismary provides real matrices with a vector of first order inclusion probabilities as only argument. The function \code{periodicdsd} provides complex matrices with fixed size, equal first order probabilities, that exhibit some periodic behavior. Both models suit the function. See the examples.
#'
#' If \code{s = 1}, the function returns a vector of length \code{n} of elements sampled.
#'
#' If \code{s > 1}, the functions returns a dataframe of s columns in which the ith column contains the ith sample.
#'
#' For \code{lavancier}, the default for \code{s} is \code{1}.
#'
#'
#' @return For \code{B = TRUE}, a vector of ones and zeros of length \code{N} (the number of elements to choose from), where ones indicate the elements drawn.
#' @return For \code{B = FALSE}, a vector of length \code{n} (the number of elements sampled) with elements drawn.
#'
#' @export
#'
#' @importFrom stats runif
#' @import ggplot2
#' @importFrom rlang abort
#'
#' @examples
#' lavancier(loonismary(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8)))
#' #returns a vector of 4 elements sampled from 1:7.
#' lavancier(loonismary(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8)), 5)
#' #returns a dataframe with 5 columns of 4 elements sampled from 1:7.
#' lavancier(loonismary(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8)), B = TRUE)
#' #returns a vector of 4 elements sampled from 1:7
#' #and plots the evolution of the probability weights at every draw.
#'
#' @references V., Loonis, X., Mary, (2017) "Determinantal Sampling Designs".
#'
#' @seealso \link{loonismary} and \link{periodicdsd}.
#'
lavancier <- function(v, s = 1, B = FALSE, C = rep(FALSE, 5000)){

  if (length(C) < ncol(v)) {
    rlang::abort("Length of the vector C is too small. There are more draws than logicals in C indicating if their probabilities should be plotted or not.")
  }

  if (is.numeric(v)) {
    return(.lavancier_mult(v, s, B, C))
  }
  else{ return(.lavancier_mult_complex(v, s, B, C))}
}














#' .lavancier_mult_complex
#'
#' This function is called by \code{lavancier}.
#'
#' @param v a matrix.
#' @param s a positive number, the number of samples to take.
#' @param B logical indicating which model of vector should be returned.
#' @param C a logical vector of length \code{n} (the number of elements sampled and the number of columns of v), indicating for each draw if the probability weights for the next draw should be plotted.
#'
#' @export
#'
#' @return s samples from DSD(vv*)
#'

.lavancier_mult_complex <- function(v, s, B, C){

  if (s == 1) {
    return(.lavancier_01_B_C_complex(v, B, C))
  }

  else{
    echant <- replicate(s, .lavancier_01_B_C_complex(v, B, C))
    colnames(echant) <- paste("Sample", 1:s)
    return(echant)
  }
}














#' .lavancier_mult
#'
#' @param v a real matrix.
#' @param s a positive number, the number of samples to take.
#' @param B logical indicating which model of vector should be returned.
#' @param C a logical vector of length \code{n} (the number of elements sampled and the number of columns of v), indicating for each draw if the probability weights for the next draw should be plotted.
#'
#' @export
#'
#' @return s samples from DSD(vv*)
#'
.lavancier_mult <- function(v, s, B, C){

  if (s == 1) {
    return(.lavancier_01_B_C(v, B, C))
  }
  else{
    echant <- replicate(s, .lavancier_01_B_C(v, B, C))
    colnames(echant) <- paste("Sample", 1:s)
    return(echant)
  }
}














#' .lavancier_01_B_C_complex
#'
#' This function is called by \code{.lavancier_mult_complex}.
#'
#' @param v a matrix.
#' @param B logical indicating which model of vector should be returned.
#' @param C a logical vector of length \code{n} (the number of elements sampled and the number of columns of v), indicating for each draw if the probability weights for the next draw should be plotted.
#'
#' @export
#'
#' @return A sample from DSD(vv*)

.lavancier_01_B_C_complex <- function(v, B = TRUE, C = rep(FALSE, nrow(v))){

  indices <- .data <- NULL

  N <- nrow(v)
  n <- ncol(v)
  echant <- rep(0, N)

  #Step 1: Sampling the first element
  w <- v

  ref <- stats::runif(1)
  total <- 0
  i <- 0
  pi1 <- Re( diag( v %*% t(Conj(v)) ) )

  while (total < ref) {
    i <- i + 1
    total <- total + ( pi1[i] / n )
  }
  echant[i] <- 1


  #Graphic representation -------------------------


  if (C[1]) {

    pin <- pi1 / n
    d <- data.frame(c(1:N), pin)
    colnames(d)[1] <- "indices"

    p <- ggplot2::ggplot(data = d, aes(x = indices))
    print(p +
            ggplot2::geom_point (aes(y = .data$pin), colour = "#999999") +
            ggplot2::geom_area(aes(y = .data$pin), colour = "#999999", fill = "#999999", alpha = 0.8) +
            ggplot2::geom_vline(xintercept = i) +
            ggplot2::theme_bw() +
            ggplot2::ylab("Pik,n") +
            ggplot2::xlab("k in U") +
            ggplot2::ggtitle(paste("Inclusion probability weights at 1st draw")))
    #----------------------------------------------------

  }
  else {
    d <- c()
  }

  M <- v[i,]
  e1 <- M / c(Re (sqrt (t(M) %*% Conj(M)) ) )


  #Step 2: Sampling the n-1 others elements
  for (j in 1:(n-1)) {

    r <- n-j
    inter <- v %*% Conj(e1)
    pi1 <- pi1 - t(inter * Conj(inter))
    pi2 <- Re( 1 / r*pi1 )

    ref <- stats::runif(1)
    total <- 0
    i <- 0

    while (total < ref) {
      i <- i + 1
      total <- total + pi2[i]
    }
    echant[i] <- 1


    #Graphic representation -----------------


    if (C[j+1]) {

      pk <- t(pi2)

      if(length(d) != 0){
        d <- cbind(d, pk)
        l <- length(d)
      }

      else{
        d <- data.frame(c(1:N), pk)
        colnames(d)[1] <- "indices"
        l <- 2
      }

      colnames(d)[l] <- paste("pi",r)

      p <- ggplot2::ggplot(data = d, aes(x = indices))
      print(p +
              ggplot2::geom_point (aes(y = d[[l]]), colour = "#999999") +
              ggplot2::geom_area(aes(y = d[[l]]), colour = "#999999", fill = "#999999", alpha = 0.8) +
              ggplot2::geom_vline(xintercept = i) +
              ggplot2::theme_bw() +
              ggplot2::ylab(paste("Pik,",r)) +
              ggplot2::xlab("k in U") +
              ggplot2::ggtitle(paste("Inclusion probability weights at",j+1 ,"th draw")))
      #------------------------------------------
    }

    w <- w - t( t(Conj(e1)) %*% t(w) ) %*% t(e1)
    L <- w[i, ]
    e1 <- L / c(Re(sqrt (t(L) %*% Conj(L) )))

  }
  if(B) {
    return(echant)
  }
  else {
    return((1:N)[echant==1])
  }

}














#' .lavancier_01_B_C
#'
#' This function is called by \code{.lavancier_mult}.
#'
#' @param v a real matrix.
#' @param B logical indicating which model of vector should be returned.
#' @param C a logical vector of length \code{n} (the number of elements sampled and the number of columns of v), indicating for each draw if the probability weights for the next draw should be plotted.
#'
#' @export
#'
#' @return A sample from DSD(vv*)
#'
.lavancier_01_B_C <- function(v, B = TRUE, C = rep(FALSE, nrow(v))){

  indices <- .data <- NULL

  N <- nrow(v)
  n <- ncol(v)
  echant <- rep(0, N)

  #First step: Sampling the first element
  w <- v

  ref <- stats::runif(1)
  total <- 0
  i <- 0
  pi1 <- diag(v  %*% t(v))

  while (total < ref) {
    i <- i + 1
    total <- total + ( pi1[i] / n )
  }
  echant[i] <- 1


  #Graphical representation -------------------------


  if (C[1]) {

    pin <- pi1 / n
    d <- data.frame(c(1:N), pin)
    colnames(d)[1] <- "indices"

    p <- ggplot2::ggplot(data = d, aes(x = indices))
    print(p +
            ggplot2::geom_point (aes(y = .data$pin), colour = "#999999") +
            ggplot2::geom_area(aes(y = .data$pin), colour = "#999999", fill = "#999999", alpha = 0.8) +
            ggplot2::geom_vline(xintercept = i) +
            ggplot2::theme_bw() +
            ggplot2::ylab("Pik,n") +
            ggplot2::xlab("k in U") +
            ggplot2::ggtitle(paste("Inclusion probability weights at 1st draw")))
    #----------------------------------------------------

  }
  else{
    d <- c()
  }

  l <- v[i,]
  e1 <- l / as.numeric( sqrt( t(l) %*% l ) )


  #Step 2: Sampling the n-1 others elements
  for (j in 1:(n-1)) {

    r <- n-j
    inter <- (v %*% e1)
    pi1 <- pi1 - t( inter * inter )
    pi2 <- 1 / r * pi1

    ref <- stats::runif(1)
    total <- 0
    i <- 0

    while (total < ref) {
      i <- i + 1
      total <- total + pi2[i]
    }
    echant[i] <- 1


    #Graphic representation -----------------


    if (C[j+1]) {

      pk <- t(pi2)

      if (length(d) != 0) {
        d <- cbind(d, pk)
        l <- length(d)
      }

      else{
        d <- data.frame(c(1:N), pk)
        colnames(d)[1] <- "indices"
        l <- 2
      }

      colnames(d)[l] <- paste("pi",r)

      p <- ggplot2::ggplot(data = d, aes(x = indices))
      print(p +
              ggplot2::geom_point (aes(y = d[[l]]), colour = "#999999") +
              ggplot2::geom_area(aes(y = d[[l]]), colour = "#999999", fill = "#999999", alpha = 0.8) +
              ggplot2::geom_vline(xintercept = i) +
              ggplot2::theme_bw() +
              ggplot2::ylab(paste("Pik,",r)) +
              ggplot2::xlab("k in U") +
              ggplot2::ggtitle(paste("Inclusion probability weights at ",j+1 ,"th draw")))
      #------------------------------------------
    }

    w <- w - (w %*% e1) %*% t(e1)
    L <- w[i,]
    e1 <- L / as.numeric( sqrt( t(L) %*% L ))

  }
  if(B) {
    return(echant)
  }
  else {
    return((1:N)[echant==1])
  }

}


