\name{autoDI_methods}
\alias{AIC.autoDI}
\alias{AICc.autoDI}
\alias{BIC.autoDI}
\alias{BICc.autoDI}
\alias{anova.autoDI}
\alias{coef.autoDI}
\alias{print.autoDI}
\alias{summary.autoDI}
\alias{formula.autoDI}
\alias{hnp.autoDI}
\alias{logLik.autoDI}
\alias{model.matrix.autoDI}
\alias{plot.autoDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for \code{autoDI} Objects
}
\description{
Different methods that can be used with objects of class \code{autoDI}.
}
\usage{
\method{AIC}{autoDI}(object, ...)
\method{AICc}{autoDI}(obj)
\method{BIC}{autoDI}(object, ...)
\method{AICc}{autoDI}(obj)
\method{anova}{autoDI}(object, ...)
\method{coef}{autoDI}(object, ...)
\method{print}{autoDI}(x, ...)
\method{summary}{autoDI}(object, ...)
\method{formula}{autoDI}(x, ...)
\method{hnp}{autoDI}(object, ...)
\method{logLik}{autoDI}(object, ...)
\method{model.matrix}{autoDI}(object, ...)
\method{plot}{autoDI}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an \code{autoDI} model object}
  \item{obj}{an \code{autoDI} model object}
\item{x}{an \code{autoDI} model object}
\item{...}{further arguments passed to \code{anova}, \code{coef}, \code{print}, \code{summary}, \code{formula}, \code{hnp}, \code{model.matrix}, or \code{plot}}
}
\author{
Rafael A. Moral, John Connolly and Caroline Brophy
}

\seealso{
\code{\link{autoDI}}
}
%\examples{
%## examples
%}