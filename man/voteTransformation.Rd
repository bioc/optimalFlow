\name{voteTransformation}
\alias{voteTransformation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{voteTransformation
%%  ~~function to do ... ~~
}
\description{ Transforming votes obtained by using optimalFlowTemplates + OptimalFlowClassification with consenus.method in c("hierarchical", "k-barycenter") and classif.method = "matching" and cost.function = "ellipses" to an appropriate format for using f1ScoreVoting.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
voteTransformation(vote.0, vote.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vote.0}{Values obtained by voteLabelTransfer}
  \item{vote.1}{Original proportions of the clusters after the template obtention.}
%%     ~~Describe \code{x} here~~

}
\value{
A list for the votes on each cell.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
vote.0 = list("1" = data.frame(cell = c(1,2), "compound.proportion" = c(0.7,0.3),
"simple.proportion"= c(0.7,0.3)), "2" = data.frame(cell = c(1,2),
"compound.proportion" = c(0.3,0.7), "simple.proportion"= c(0.3,0.7)))
vote.1.1 = t(c(0.8,0.2))
names(vote.1.1) = c("A","B")
vote.1.2 = t(c(0.2,0.8))
names(vote.1.2) = c("A","B")
vote.1 = list(vote.1.1, vote.1.2)
voteTransformation(vote.0, vote.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
