\name{optimalFlowTemplates}
\alias{optimalFlowTemplates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{optimalFlowTemplates
%%  ~~function to do ... ~~
}
\description{Returns a partition of the input clusterings with a respective consensus cluster for every group in the partition.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
optimalFlowTemplates(database, database.names = NULL, cov.estimation = "standard",
alpha.cov = 0.85,equal.weights.template = TRUE, hclust.method = "complete",
templates.number = NA,minPts = 2, eps = 1, consensus.method = "pooling",
barycenters.number = 37,bar.repetitions = 40, alpha.bar = 0.05,
bar.ini.method = "plus-plus", consensus.minPts = 3, cl.paral = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{database}{A list where each entry is a partition (clustering) represented as dataframe, of the same dimensions, where the last variable represents the labels of the partition.}
  \item{database.names}{Names of the elements in the database.}
  \item{cov.estimation}{How to estimate covariance matrices in each cluster of a partition. "standard" is for using cov(), while "robust" is for using robustbase::covMcd.}
  \item{alpha.cov}{Only when cov.estimation = "robust". Indicates the value of alpha in robustbase::covMcd.}
  \item{equal.weights.template}{If True, weights assigned to every cluster in a partion are uniform (1/number of clusters). If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.}
  \item{hclust.method}{Indicates what kind of hierarchical clustering to  do with the similarity distances matrix of the partitions. Takes values in c("complete", "single", "average", "hdbscan", "dbscan").}
  \item{templates.number}{Only if hclust.method in c("complete", "single", "average"). Indicates the number of clusters to use with cutree. If set to NA (default), plots the hierarchical tree and asks the user to introduce an appropriate number of clusters.}
  \item{minPts}{Only if hclust.method in c("hdbscan", "dbscan"). Indicates the value of argument minPts in dbscan::dbscan and dbscan::hdbscan.}
  \item{eps}{Only if hclust.method = "dbscan". Indicates the value of eps in dbscan::dbscan.}
  \item{consensus.method}{Sets the way of doing consensus clustering when clusters are viewed as Multivariate Distributions. Can take values in c("pooling", "k-barycenter", "hierarchical"). See details.}
  \item{barycenters.number}{Only if consensus.method = "k-barycenter". Sets the number, k, of barycenters when using k-barycenters.}
  \item{bar.repetitions}{Only if consensus.method = "k-barycenter". How many times to repeat the k-barycenters procedure. Equivalent to nstart in kmeans.}
  \item{alpha.bar}{Only if consensus.method = "k-barycenter". The level of trimming allowed during the k-barycenters procedure.}
  \item{bar.ini.method}{Only if consensus.method = "k-barycenter". Takes values in c("rnd", "plus-plus"). See details.}
  \item{consensus.minPts}{Only if consensus.method = "hierarchical". The value of argument minPts for dbscan::hdbscan.}
  \item{cl.paral}{Number of cores to be used in parallel procedures.}
%%     ~~Describe \code{x} here~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{templates}{A list of the consensus clusterings for every group in the partition of the database.}
\item{clustering}{Clustering of the input partitions.}
\item{database.elliptical}{Means, covariances and weights of the clusters in the input partitions.}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
