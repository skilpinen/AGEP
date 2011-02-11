\name{do.match}
\Rdversion{1.1}
\alias{do.match}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare a sample by gene distributions
}
\description{
Compares a sample to a number of tissues (or other categories) for whose gene (or other feature) density estimates have been pre-calculated.
}
\usage{
do.match(query.sample, query.annotation = NULL, query.id = NULL, ref.data.obj, verbose = TRUE, verbose.res = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query.sample}{
  named numeric vector representing the sample to be matched. The names of the elements must match the genes (features) whose densities are found in \code{ref.data.obj}.
}
  \item{query.annotation}{
  optional string that denotes the a priori known tissue (category) the matched sample belongs to.
}
  \item{query.id}{
  optional string containing some identification of the matched sample.
}
  \item{ref.data.obj}{
  the object returned by the \code{create.ref.obj} function.
}
  \item{verbose}{
  determines whether the progress of the process will be printed to the screen. Default \code{FALSE}.
}
  \item{verbose.res}{
  if \code{TRUE} the functions returns all the scores for all genes (features). If \code{FALSE}, the function only returns the final matching scores to all tissues (categories). Default \code{FALSE}.
}
}
\details{
This function performs the AGEP matching of a single sample to a reference data. The reference data must have been processed with the \code{create.ref.obj} function, which calculates the densities of all genes (features) for all tissues (categories) in the data. Has two optional outputs, controlled by the \code{verbose} argument. If set to \code{TRUE}, will return all scores from all steps along with the final similarity scores.
}
\value{
A list with the following components:

\item{em.scores}{
a matrix with the expression match scores for each gene (feature) and each tissue (category). Returned only when \code{verbose.res} is set to \code{TRUE}.
}
\item{ts.scores}{
a matrix with the tissue specificity scores for each gene (feature) and each tissue (category). Returned only when \code{verbose.res} is set to \code{TRUE}.
}
\item{query.id}{
the \code{query.id} argument, if not given then the default \code{NULL}.
}
\item{query.annotation}{
the \code{query.annotation} argument, if not given then the default \code{NULL}.
}
\item{tissue.scores}{
a numeric named vector that contains tissue similarity scores for all tissues (categories) used in matching. Elements names contain the tissues (categories). 
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
