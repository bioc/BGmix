\name{readBGX}
\alias{readBGX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Reads output from BGX package, for input to BGmix.}
\description{
Reads posterior mean parameters from BGX, and outputs objects suitable for input to BGmix.
}
\usage{
readBGX(path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{ directory containing BGX output}
}
\value{
  \item{ybar }{ybar object (see \code{\link{BGmix}} help for details)}
  \item{ss }{ss object (see \code{\link{BGmix}} help for details)}
  ...
}
\author{ Ernest Turro}
\keyword{ hplot }
