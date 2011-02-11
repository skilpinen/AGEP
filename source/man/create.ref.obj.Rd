\name{create.ref.obj}
\Rdversion{1.1}
\alias{create.ref.obj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create gene distributions from reference data  
}
\description{
Creates gene, or other biological feature, distributions from reference data, to be used with \code{do.match}.
}
\usage{
create.ref.obj(data, annotation, min.data = 5, min.bw = -Inf, max.bw = Inf, resolution = 512, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  numeric data matrix where columns are samples and rows are genes, or other quantifiable features.
}
  \item{annotation}{
  character vector whose length is the same as the number of columns in \code{data}. Indicates for each sample the tissue (or other category) it belongs to.
}
  \item{min.data}{
  the minimum number datapoints a gene (or other feature) must have per class for the distribution to be calculated. Default \code{5}.
}
  \item{min.bw}{
  sets the minimum bandwidth used to calculated the densities. Default \code{-Inf}.
}
  \item{max.bw}{
  sets the maximum bandwidth used to calculated the densities. Default \code{Inf}.
}
  \item{resolution}{
  the number of evaluation points used in the density calculation. Use multiples of 512. Default \code{512}.
}
  \item{verbose}{
  if \code{TRUE}, prints to the screen the progress per tissue (category) and gene (feature).
}
}
\details{
This function creates density estimates for each row of a numeric matrix whose columns are set into categories by \code{annotation}. It generates by 512 (adjustable by \code{resolution}) evaluation points that span the range of values for that entire row, and then calculates the densities for each category for that row within the range. The return object is rather large, and this function is primarily meant to be used only with the \code{do.match} function.
}
\value{
A list with the following components:

\item{data.densities}{
contains a list with a matrix for each tissue (category) whose rows are evaluation points and columns are genes (features), the value of the elements indicating the density value of the gene in that evaluation point in that tissue.
}
\item{froms}{
a vector containing the expression value (or other measure) of the first evaluation point for each gene (feature). Currently always 0.
}
\item{tos}{
a vector containing the expression value (or other measure) of the last evaluation point for each gene (feature). Currently the maximum observed expression plus a certain percentage.
}
\item{scaling.factors}{
contains a list with a vector for each tissue (category), the elements of the vectors correspond to genes (features), the numbers indicate the multiplier for that gene's density values so as to reach an area of 1 for the whole density graph.
}
\item{resolution}{
a number indicating the amount of evaluation points used.
}
\item{bws}{
contains a list with a vector for each tissue (category), the elements of the vectors correspond to genes (features), the numbers indicate the bandwidth used to calculate the density for the gene in the tissue.
}
\item{ref.data}{
copy of the \code{data} argument.
}
\item{annotation}{
copy of the \code{annotation} argument.
}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sami Kilpinen, Kalle Ojala
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
