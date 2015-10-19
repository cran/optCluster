\name{repRankAggreg}
\alias{repRankAggreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Repeat Rank Aggregation}
\description{
  \code{repRankAggreg} repeats rank aggregation of ordered validation measure lists 
  obtained from an \cr object of class \code{"\linkS4class{optCluster}"}. The 
  function returns an object of class \code{"\linkS4class{optCluster}"}.
}
\usage{
  repRankAggreg(optObj, rankMethod = "same", distance = "same", 
  rankVerbose = FALSE, ... ) 
}


\arguments{
  \item{optObj}{An object of class \code{"\linkS4class{optCluster}"}.}
  \item{rankMethod}{A character string
  	providing the method to be used for rank aggregation. As default, the "same" 
  	method as the input \code{"\linkS4class{optCluster}"} object is used. The 
  	cross-entropy Monte Carlo algorithm ("CE") or Genetic algorithm ("GA") can also be 
  	directly specified.  Selection of only one method is allowed.}  
  \item{distance}{A character string providing the type of distance to be used for measuring the similarity between ordered lists  
  	in rank aggregation.  As default, the "same" distance as the input \code{"\linkS4class{optCluster}"} object is used. 
  	The weighted Spearman footrule distance ("Spearman") or the weighted Kendall's tau distance \cr("Kendall") 
  	can also be directly specified.  Selection of only one distance is allowed.}  
  \item{rankVerbose}{If TRUE, current rank aggregation results are displayed at each iteration.}
  \item{\dots}{Additional arguments that can be passed to the internal function \code{RankAggreg}.}
}

\details{
	This function tests the consistency of the rank aggregation results by repeating rank aggregation with the same 
	rank aggregation method, distance measure, clustering algorithm lists, and validation score lists used to create
	the input object of class \code{"\linkS4class{optCluster}"}.  A different rank aggregation algorithm or 
	type of distance measure can also be evaluated using this function, but doing so may affect the final results.
	
}

\value{

  \code{repRankAggreg} returns an object of class \code{"\linkS4class{optCluster}"}.  The class description 
  is provided in the help file.

}
\references{
  
  Pihur, V., Datta, S. and Datta, S. (2007).
  Weighted rank aggregation of cluster validation measures: A Mounte Carlo cross-entropy approach.
  Bioinformatics 23(13): 1607-1615.
     
  Pihur, V., Datta, S. and Datta, S. (2009).
  RankAggreg, an R package for weighted rank aggregation.
  BMC Bioinformatics, 10:62,
  \url{http://www.biomedcentral.com/1471-2105/10/62}.
  
  Sekula, M. (2015). 
  optCluster : An R package for Determining the Optimal Clustering Algorithm and Optimal Number of Clusters.
  Electronic Theses and Dissertations. Paper 2147.
  \url{http://ir.library.louisville.edu/etd/2147}  
}

\author{Michael Sekula, Somnath Datta, Susmita Datta}


\seealso{

  For a description of the \code{RankAggreg} function, including all available arguments that can be
  passed to it, see \code{\link{RankAggreg}} in the \pkg{RankAggreg} package.  
  
  }

\examples{
	
	## These examples may take a few minutes to compute
	\donttest{	
	## Obtain Dataset	
	data(arabid)			
	
	## Normalize Data with Respect to Library Size	
	obj <- t(t(arabid)/colSums(arabid))
	
	## Analysis of Normalized Data using Internal and Stability Validation Measures
	norm1 <- optCluster(obj, 2:4, clMethods = "all")
	print(norm1)
	repCE <- repRankAggreg(norm1)
	print(repCE)
	repGA <- repRankAggreg(norm1, rankMethod = "GA")
	print(repGA)
	}
	
}

\keyword{cluster}
\keyword{optimize}
