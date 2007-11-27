\name{plotCompare}
\alias{plotCompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Scatter plot with equal axes.}
\description{
Plots a scatter plot of two variables with equal scales for the axes.
}
\usage{
plotCompare(var1, var2, limi = 0, xlab = substitute(var1), ylab = substitute(var2), log = "", title = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var1}{ data to plot (x co-ordinate)}
  \item{var2}{ data to plot (y co-ordinate)}
  \item{limi}{ limits of axes. If not specified, axes limits are determined from input data.}
  \item{xlab}{ x-axis label}
  \item{ylab}{ y-axis label}
  \item{log}{ specifies if axes are on the log scale (as argument to 'par')}
  \item{title}{ title of plot}
  \item{\dots}{ other parameters input to plot}
}
\value{
Outputs the limits used in the plot (the input 'limi' argument if specified).
}
\author{ Alex Lewin}
\keyword{ hplot }
\examples{
x <- runif(100)
y <- rbeta(100,0.5,0.5)
plotCompare(x,y)
}
