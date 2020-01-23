\name{writeGsea}
\alias{writeGsea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write GMT and CLS files from a matrix of expression}
\description{
Prints a standard script written in bash that can be used for the analysis in HPC environment at CRG.
}
\usage{
writeGsea ( x , cath , filename )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix of expression. Should have rownames and colnames defined.}
  \item{cath}{categorie. In this version only the numeric option is included.}
  \item{filename}{name of the file.}
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
#writeGsea(x,cath,filename)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GSEA}% use one of  RShowDoc('KEYWORDS')
