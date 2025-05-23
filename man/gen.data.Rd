% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen.data.R
\name{gen.data}
\alias{gen.data}
\title{Generate simulated data}
\usage{
gen.data(
  n,
  m,
  d,
  s,
  s0,
  cor.type = 1,
  beta.type = 1,
  rho = 0.5,
  sigma1 = 1,
  sigma2 = 1,
  seed = 1
)
}
\arguments{
\item{n}{The number of observations.}

\item{m}{The number of groups of interest.}

\item{d}{The group size of each group. Only even group structure is allowed here.}

\item{s}{The number of important groups in the underlying regression model.}

\item{s0}{The number of important variables in each important group.}

\item{cor.type}{The structure of correlation.
\code{cor.type = 1} denotes the independence structure,
where the covariance matrix has \eqn{(i,j)} entry equals \eqn{I(i \neq j)}.
\code{cor.type = 2} denotes the exponential structure,
where the covariance matrix has \eqn{(i,j)} entry equals \eqn{rho^{|i-j|}}.
\code{cor.type = 3} denotes the constant structure,
where the non-diagonal entries of covariance
matrix are \eqn{rho} and diagonal entries are 1.}

\item{beta.type}{The structure of coefficients.
\code{beta.type = 1} denotes the homogenous setup,
where each entry has the same magnitude.
\code{beta.type = 2} denotes the heterogeneous structure,
where the coefficients are drawn from a normal distribution.}

\item{rho}{A parameter used to characterize the pairwise correlation in
predictors. Default is \code{0.5}..}

\item{sigma1}{The value controlling the strength of the gaussian noise. A large value implies strong noise. Default \code{sigma1 = 1}.}

\item{sigma2}{The value controlling the strength of the coefficients. A large value implies large coefficients. Default \code{sigma2 = 1}.}

\item{seed}{random seed. Default: \code{seed = 1}.}
}
\value{
A \code{list} object comprising:
\item{x}{Design matrix of predictors.}
\item{y}{Response variable.}
\item{beta}{The coefficients used in the underlying regression model.}
\item{group}{The group index of each variable.}
\item{true.group}{The important groups in the sparse group linear model.}
\item{true.variable}{The important variables in the sparse group linear model.}
}
\description{
Generate simulated data for sparse group linear model.
}
\examples{

# Generate simulated data
n <- 200
m <- 100
d <- 10
s <- 5
s0 <- 5
data <- gen.data(n, m, d, s, s0)
str(data)
}
\author{
Yanhang Zhang, Zhifan Li, Jianxin Yin.
}
