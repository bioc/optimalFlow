\name{f1ScoreVoting}
\alias{f1ScoreVoting}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{f1ScoreVoting
%%  ~~function to do ... ~~
}
\description{ Calculates the F1 score fore each group in a partition, when provided with a fuzzy classification.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
f1ScoreVoting(voting, clustering, cytometry, nivel_sup, noise.cells)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{voting}{A list where each entry is a vote on the respective label}
  \item{clustering}{Labels of the partition}
  \item{cytometry}{Data of the clustering, where the last variable are the original labels.}
  \item{nivel_sup}{level of tolerance for assigning hard clustering. Should be greater or equal than 1. Class A is assigned if class A >nivel_sup*Class B.}
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
database = list(as.data.frame(Cytometry2)[which(match(Cytometry2$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry3)[which(match(Cytometry3$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry4)[which(match(Cytometry4$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry5)[which(match(Cytometry5$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry7)[which(match(Cytometry7$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry8)[which(match(Cytometry8$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry9)[which(match(Cytometry9$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry12)[which(match(Cytometry12$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry13)[which(match(Cytometry13$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry14)[which(match(Cytometry14$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry15)[which(match(Cytometry15$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry16)[which(match(Cytometry16$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry17)[which(match(Cytometry17$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry19)[which(match(Cytometry19$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),],
as.data.frame(Cytometry21)[which(match(Cytometry21$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),])


templates.optimalFlow = optimalFlowTemplates(database = database, templates.number = 5,
cl.paral = 1)


classification.optimalFlow = optimalFlowClassification(as.data.frame(Cytometry1)[
which(match(Cytometry1$`Population ID (name)`,c("Monocytes", "CD4+CD8-",
"Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0), 1:10],database, templates.optimalFlow,
classif.method = "matching", cost.function = "ellipses", cl.paral = 1)


f1ScoreVoting(classification.optimalFlow$cluster.vote, classification.optimalFlow$cluster,
as.data.frame(Cytometry1)[which(match(Cytometry1$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),], 1.01, noise.types)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
