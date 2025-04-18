#' @title Adaptive Double Sparse Iterative Hard Thresholding Algorithm (ADSIHT)
#'
#' @description An implementation of the sparse group selection in linear regression model via ADSIHT.
#'
#' @param x Input matrix, of dimension \eqn{n \times p}; each row is an observation
#' vector and each column is a predictor.
#' @param y The response variable of \code{n} observations.
#' @param group A vector indicating which group each variable belongs to
#' For variables in the same group, they should be located in adjacent columns of \code{x}
#' and their corresponding index in \code{group} should be the same.
#' Denote the first group as \code{1}, the second \code{2}, etc.
#' @param s0 A vector that controls the degrees with group.
#' Default is \eqn{d^((l-1)/(L-1))} : \eqn{1 \leq l \leq L}, where d is the maximum group size.
#' @param kappa A parameter that controls the rapid of the decrease of threshold. Default is 0.9.
#' @param ic.type The type of criterion for choosing the support size.
#' Available options are \code{"dsic"}, \code{"loss"}.
#' Default is \code{"dsic"}.
#' @param ic.scale A non-negative value used for multiplying the penalty term
#' in information criterion. Default: \code{ic.scale = 3}.
#' @param ic.coef A non-negative value used for multiplying the penalty term
#' for choosing the optimal stopping time. Default: \code{ic.coef = 3}.
#' @param L The length of the sequence of s0. Default: \code{L = 5}.
#' @param weight The weight of the samples, with the default value set to 1 for each sample.
#' @param coef1 A positive value to control the sub-optimal stopping time.
#' @param coef2 A positive value to control the overall stopping time. A small value leads to larger search range.
#' @param eta A parameter controls the step size in the gradient descent step.
#' Default: \code{eta = 0.8}.
#' @param max_iter A parameter that controls the maximum number of line search, ignored if \code{OLS} is employed.
#' @param method Whether \code{ols} (default) or \code{linesearch} method should be employed.
#'
#' @return A \code{list} object comprising:
#' \item{beta}{A \eqn{p}-by-\code{length(s0)} matrix of coefficients, stored in column format.}
#' \item{intercept}{A \code{length(s0)} vector of intercepts}.
#' \item{lambda}{A \code{length(s0)} vector of threshold values}
#' \item{A_out}{The selected variables given threshold value in \code{lambda}.}
#' \item{ic}{The values of the specified criterion for each fitted model given threshold \code{lamdba}.}
#'
#' @author Yanhang Zhang, Zhifan Li, Shixiang Liu, Jianxin Yin.
#'
#' @export
#'
#' @examples
#'
#' n <- 200
#' m <- 100
#' d <- 10
#' s <- 5
#' s0 <- 5
#' data <- gen.data(n, m, d, s, s0)
#' fit <- ADSIHT(data$x, data$y, data$group)
#' fit$A_out[which.min(fit$ic)]

ADSIHT <- function(x, y, group,
                  s0,
                  kappa = 0.9,
                  ic.type = c("dsic", "loss"),
                  ic.scale = 3.0,
                  ic.coef = 3.0,
                  L = 5,
                  weight = rep(1, nrow(x)),
                  coef1 = 1,
                  coef2 = 1,
                  eta = 0.8,
                  max_iter = 20,
                  method = "ols")
{
  if(missing(group)) group <- 1:ncol(x)
  p <- ncol(x)
  n <- nrow(x)
  N <- length(unique(group))
  if (length(group)!= ncol(x)) stop("The length of group should be the same with ncol(x)")
  if(!is.matrix(x)) x <- as.matrix(x)
  vn <- colnames(x)
  orderGi <- order(group)
  x <- x[, orderGi]
  vn <- vn[orderGi]
  group <- group[orderGi]
  gi <- unique(group)
  index <- match(gi, group)-1
  if (missing(s0)) {
    s0 <- max(table(group))^(seq(1, L-1, length.out = L)/(L-1))
  }
  ic.type <- match.arg(ic.type)
  ic_type <- switch(ic.type,
                    "loss" = 0,
                    "dsic" = 1
                    )
  if (method == "ols") {
    method = TRUE
  } else {
    method = FALSE
  }
  res <- DSIHT_Cpp(x, y, weight = weight, sequence = s0, ic_type = ic_type, ic_scale = ic.scale, kappa = kappa, g_index = index, ic_coef = ic.coef, coef1 = coef1, coef2 = coef2, eta = eta, max_iter = max_iter, method = method, nor = TRUE)
  return(res)
}
