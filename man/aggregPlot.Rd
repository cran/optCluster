\name{aggregPlot}
\alias{aggregPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rank Aggregation Plot}
\description{
  \code{aggregPlot} displays a figure representing the results from rank
  aggregation.
}
\usage{
aggregPlot(x, show.average = TRUE, show.legend = TRUE, colR = "red", ...)
}


\arguments{
  \item{x}{An object of class \code{"\linkS4class{optCluster}"}.}
	\item{show.average}{If \code{TRUE}, the average aggregation will be plotted.}
	\item{show.legend}{If \code{TRUE}, the legend will be displayed.}
	\item{colR}{A character string specifying the line color for the optimal ranked list.}
	\item{\dots}{Additional plotting parameters.}	
}

\details{
	This function will produce a visual representation of rank aggregation results in three plots
	combined into a single figure.  The top left plot displays the minimum values of the objective 
	function as the number of iterations increases.  The top right plot is a histogram of the 
	objective function scores at the last iteration.  The bottom plot shows the individual ranks
	of the clustering algorithms from the data (in grey), the final "optimal" rankings 
	(in red), and the average ranking of each clustering algorithm (in black).
	}

\references{
  
  Pihur, V., Datta, S. and Datta, S. (2009).
  RankAggreg, an R package for weighted rank aggregation.
  BMC Bioinformatics, 10:62,
  \url{http://www.biomedcentral.com/1471-2105/10/62}.
    
}

%\author{Michael Sekula, Somnath Datta, Susmita Datta}


\seealso{

\code{\link{RankAggreg}}, \code{\link{optCluster-class}}
  
}

\examples{
	
	## This example may take a few minutes to compute
	\donttest{
	## Obtain Dataset
	data(arabid)
		
	## Normalize Data with Respect to Library Size	
	obj <- t(t(arabid)/colSums(arabid))
	
	## Analysis of Normalized Data using Internal and Stability Validation Measures
	norm1 <- optCluster(obj, 2:4, clMethods = "all")
	aggregPlot(norm1)
	
	}
}
