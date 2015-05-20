\name{optCluster-class}
\docType{class}
\alias{optClusterObj}
\alias{optCluster-class}

\alias{getDataset,optCluster-method}
\alias{getClValid,optCluster-method}
\alias{methodRanks,optCluster-method}
\alias{scoreRanks,optCluster-method}
\alias{getRankAggreg,optCluster-method}
\alias{topMethod,optCluster-method}
\alias{measureNames,optCluster-method}
\alias{methodNames,optCluster-method}
\alias{clusterResults,optCluster-method}
\alias{valScores,optCluster-method}
\alias{optimalScores,optCluster-method}

\alias{getDataset}
\alias{getClValid}
\alias{methodRanks}
\alias{scoreRanks}
\alias{getRankAggreg}
\alias{topMethod}
\alias{measureNames}
\alias{methodNames}
\alias{clusterResults}
\alias{valScores}
\alias{optimalScores}

\alias{print,optCluster-method}
\alias{show,optCluster-method}
\alias{summary,optCluster-method}

\title{Class "optCluster"}
\description{The class \code{"optCluster"} contains the dataset,
	clustering results, validation measures, ranked lists of clustering algorithms,
	ordered lists of validation scores, and final rank aggregation results from
	the function \code{\link{optCluster}}.
}
\section{Objects from the Class}{
  The function \code{\link{optCluster}} creates objects of the class 
  \code{"optCluster"}.
}

\section{Slots}{
  \describe{
    \item{\code{inputData}:}{Object of class \code{"matrix"} 
    	containing the original dataset.}    
    \item{\code{clVal}:}{Object of class \code{"\linkS4class{clValid}"} 
    	containing the clustering results and validation measures from the
    	internal \code{\link{clValid}} function.}
    \item{\code{ranksWeights}:}{Object of class \code{"list"} 
    	containing the ordered ranks of clustering algorithms and 
    	the ordered validation scores for each measure.}
    \item{\code{rankAgg}:}{Object of class \code{"raggr"}
      containing the rank aggregation results from the internal \cr 
      \code{\link{RankAggreg}} function.}
  }
}

\section{Methods}{
  \describe{
    \item{getDataset}{\code{signature(object = "optCluster")}: Returns the
      original dataset as an object of class \code{"matrix"}.}
    \item{getClValid}{\code{signature(object = "optCluster")}: Returns an
      object of class \code{"\linkS4class{clValid}"}.}
    \item{methodRanks}{\code{signature(object = "optCluster")}: Returns the
    	ranked lists of clustering algorithms for each validation measure.} 
    \item{scoreRanks}{\code{signature(object = "optCluster")}: Returns the
    	ordered lists of scores for each validation measure.} 	     
    \item{getRankAggreg}{\code{signature(object = "optCluster")}: Returns an
      object of class \code{"raggr"}.}
    \item{topMethod}{\code{signature(object = "optCluster")}: Returns the
      name of the optimal clustering algorithm and number of clusters.}
    \item{measureNames}{\code{signature(object = "optCluster")}: Returns the
      names of the validation measures used.}      
    \item{methodNames}{\code{signature(object = "optCluster")}: Returns the
      names of the clustering algorithms used.}      
    \item{clusterResults}{\code{signature(object = "optCluster")}: Returns the
	  clustering results from the selected clustering method.}
	  
	  Additional arguments:
      \describe{
			\item{\code{method = methodNames(object)}}{ The clustering
				algorithm to extract. The selection of only one algorithm 
				is allowed.}
      }      
    \item{valScores}{\code{signature(object = "optCluster")}: Returns the
	  scores from the selected validation measure(s).}
	  
	  Additional arguments:
      \describe{
			\item{\code{measures = measureNames(object)}}{ The validation
				measure(s) to extract.}
      }       
     \item{optimalScores}{\code{signature(object = "optCluster")}: Returns the
      optimal score for each validation measure as well as the 
      corresponding clustering algorithm and number of clusters.}      
            
    \item{print}{\code{signature(x = "optCluster")}: Print method for class
      \code{"optCluster"}.}
    \item{show}{\code{signature(object = "optCluster")}: Same as print.}
    \item{summary}{\code{signature(object = "optCluster")}: Summary method
      for class \code{"optCluster"}.} 
    }  
}

\references{

  Brock, G., Pihur, V., Datta, S. and Datta, S. (2008).
  clValid: An R Package for Cluster Validation.
  Journal of Statistical Software 25(4),
  \url{http://www.jstatsoft.org/v25/i04}.
  
  Datta, S. and Datta, S. (2003).
  Comparisons and validation of statistical clustering techniques for microarray gene expression data.
  Bioinformatics 19(4): 459-466. 
  
  Pihur, V., Datta, S. and Datta, S. (2007).
  Weighted rank aggregation of cluster validation measures: A Mounte Carlo cross-entropy approach.
  Bioinformatics 23(13): 1607-1615.
     
  Pihur, V., Datta, S. and Datta, S. (2009).
  RankAggreg, an R package for weighted rank aggregation.
  BMC Bioinformatics, 10:62,
  \url{http://www.biomedcentral.com/1471-2105/10/62}.
  
  Sekula, M. (2015).
  optCluster: An R Package for Determining the Optimal Clustering Algorithm and Optimal Number of Clusters.
  MS Thesis, University of Louisville.
}
\author{Michael Sekula, Somnath Datta, Susmita Datta}

\seealso{
	
  For details on the function \code{optCluster} see \code{\link{optCluster}}.


  For a description of the \code{clValid} function, including all available arguments that can be
  passed to it, see \code{\link{clValid}} in the \pkg{clValid} package.  
  For a desciption of the class \code{"clValid"} including all available methods see 
  \code{\link{clValid-class}}.
  
  For a description of the \code{RankAggreg} function, including all available arguments that can be 
  passed to it, see \code{\link{RankAggreg}} in the \pkg{RankAggreg} package. 

}


\examples{
	
	## This example may take a few minutes to compute
	\donttest{
	## Obtain Dataset	
	data(arabid)	
	
	## Normalize Data with Respect to Library Size	
	obj <- t(t(arabid)/colSums(arabid))
	
	## Analysis of Normalized Data using Internal and Stability Measures
	norm1 <- optCluster(obj, 2:4, clMethods = "all")

	## View results
	norm1
	topMethod(norm1)
	summary(norm1)
	optimalScores(norm1)
	
	## Extract additional information from slots
	methodNames(norm1)
	clusterResults(norm1, "kmeans")
	measureNames(norm1)
	valScores(norm1)
	methodRanks(norm1)
	scoreRanks(norm1)
	}

}

\keyword{classes}
\keyword{cluster}
\keyword{optimize}