\name{optimalFlowClassification}
\alias{optimalFlowClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{optimalFlowClassification
%%  ~~function to do ... ~~
}
\description{ Performs a supervised classification of the input data where a database and a partition of the database are available
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
optimalFlowClassification(X, database, templates, consensus.method = "pooling",
cov.estimation = "standard", alpha.cov = 0.85,initial.method = "supervized",
alpha.tclust = 0, restr.factor.tclust = 1000,classif.method = "qda",
qda.bar = TRUE, cost.function = "points", cl.paral = 1,
equal.weights.voting = TRUE, equal.weights.template = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Datasample to be classified.}
  \item{database}{A list where each entry is a partition (clustering) represented as dataframe, of the same dimensions, where the last variable represents the labels of the partition.}
  \item{templates}{List of the consensus clusterings for every group in the partition of the database obtained by optimalFlowTemplates}
  \item{consensus.method}{The consensus.method value that was used in optimalFlowTemplates.}
   \item{cov.estimation}{How to estimate covariance matrices in each cluster of a partition. "standard" is for using cov(), while "robust" is for using robustbase::covMcd.}
   \item{alpha.cov}{Only when cov.estimation = "robust". Indicates the value of alpha in robustbase::covMcd.}
   \item{initial.method}{Indicates how to initialize tclust. Currently only supports "supervised".}
   \item{alpha.tclust}{Level of trimming allowed fo tclust.}
   \item{restr.factor.tclust}{Fixes the restr.fact parameter in tclust.}
   \item{classif.method}{Indicates what type of supervised learning we want to do. Takes values on c("matching", "qda", "random forest").}
   \item{qda.bar}{Only if classif.method = "qda". If True then the appropriate consensus clustering (template, prototype) is used for learning. If False, the closest partition in the appropriate group is used.}
   \item{cost.function}{Only if classif.method = "matching". Indicates the cost function, distance between clusters, to be used for label matching.}
   \item{cl.paral}{Number of cores to be used in parallel procedures.}
   \item{equal.weights.voting}{only when classif.method = "qda" and qda.bar =F, or when  classif.method = "random forest". Indicates the weights structure when looking for the most similar partition in a group.}
   \item{equal.weights.template}{If True, weights assigned to every cluster in a partion are uniform (1/number of clusters). If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.}
%%     ~~Describe \code{x} here~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{cluster}{Labels assigned to the input data.}
\item{clusterings}{The different partitions of the data obtained by tclust when initialized with the templates obtained by optimalFlowTemplates.}
\item{assigned.template.index}{Label of the group for which the template is closer to the data.}
\item{cluster.vote}{Only when classif.method = "matching". Vote on the type of every label in the partition of the data.}
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
which(match(Cytometry1$`Population ID (name)`,c("Monocytes", "CD4+CD8-", "Mature SIg Kappa",
"TCRgd-"), nomatch = 0)>0), 1:10], database, templates.optimalFlow, cl.paral = 1)


scoreF1.optimalFlow = optimalFlow::f1Score(classification.optimalFlow$cluster,
as.data.frame(Cytometry1)[which(match(Cytometry1$`Population ID (name)`,
c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0)>0),], noise.types)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
