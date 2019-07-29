\name{f1Score}
\alias{f1Score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{f1Score
%%  ~~function to do ... ~~
}
\description{ Calculates the F1 score fore each group in a partition.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
f1Score(clustering, cytometry, noise.cells)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{The labels of the new classification.}
  \item{cytometry}{Data of the clustering, where the last variable are the original labels.}
  \item{noise.cells}{An array of labels to be considered as noise.}
%%     ~~Describe \code{x} here~~
}
\value{
A matrix where the first row is the F1 score, the seconr row is the Precision and the third row is the Recall.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
E del Barrio, H Inouzhe, JM Loubes, C Matran and A Mayo-Iscar. (2019) optimalFlow: Optimal-transport approach to flow cytometry gating and population matching. arXiv:1907.08006
}
\examples{
f1Score(dplyr::pull(Cytometry3[c(sample(1:250,250),251:(dim(Cytometry3)[1])),],11),
Cytometry3, noise.types)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
