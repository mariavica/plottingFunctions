\name{printScript}
\alias{printScript}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot two trends}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
printScript ( x = NULL, type , genome , sequencing , loc )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{name of the file}
  \item{type}{one of: \code{'RNAseq'}, \code{'ATACseq'}, \code{'ChipSeq'}, \code{'HiC'}}
  \item{genome}{version of the genoma you are aligning to: \code{'mm10'} or \code{'hg38'}}
  \item{sequencing}{\code{'singe'} (or \code{'singe-end'}) or \code{'paired'} (or \code{'paired-end'})}
  \item{loc}{localization: \code{'geo'} will download data from GEO database, \code{'local'} will download data from the selected folder}
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
\keyword{bash}% use one of  RShowDoc('KEYWORDS')
