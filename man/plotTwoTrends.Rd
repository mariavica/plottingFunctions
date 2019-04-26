\name{plotTwoTrends}
\alias{plotTwoTrends}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot two trends}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotTwoTrends (x1, x2, x = NA, 
               col1 = 'blue', col2 = 'red', lwd = 3,
               normalization = 'none', norm_independent = FALSE,
               expand = 0.25, fix.ylim = NULL,
               loess = FALSE, sp = 0.5, res.x = 1000, 
               plot.sd = TRUE, 
               xlab = NA, ylab1 = NA, ylab2 = NA, omit.lab2 = FALSE )}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{dataset 1}
  \item{x2}{dataset 2}
  \item{x}{x}
  \item{col1}{color of dataset 1}
  \item{col2}{color of dataset 2}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}% use one of  RShowDoc("KEYWORDS")
