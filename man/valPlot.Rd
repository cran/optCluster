\name{valPlot}
\alias{valPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Validation Measure Plots}
\description{
  \code{valPlot} displays a plot of the scores for each selected validation measure.
}
\usage{
valPlot(x, measures = measureNames(x), legend = TRUE, 
	legendLoc = "topright", main = NULL, pch = NULL, type = "b", 
	ask = prod(par("mfcol")) < length(measures) && dev.interactive(), ...)
}


\arguments{
  \item{x}{An object of class \code{"\linkS4class{optCluster}"}.}
	\item{measures}{ Character vector of the names of the validation measures to plot.  Any
		number of choices is allowed.}
	\item{legend}{ If TRUE, provides a legend.}
	\item{legendLoc}{Character string specifying the location of the legend.}
	\item{main}{Character string specifying the title of graph.}
	\item{pch}{ Specifies the plotting characters to use.}
	\item{type}{A character string specifying the type of plot.}
	\item{ask}{ If TRUE, the user is prompted before each plot.}
	\item{\dots}{Additional plotting parameters.}	
}

\details{

	The the biological homogeneity index (BHI), biological stability index (BSI), Dunn index, 
	and silhouette width measures should all be maximized.  
	
	The average proportion of non-overlap (APN), average distance (AD), 
	average distance between means (ADM), figure of merit (FOM), and connectivity measures 
	should all be minimized. }

\references{
  
  Brock, G., Pihur, V., Datta, S. and Datta, S. (2008).
  clValid: An R Package for Cluster Validation.
  Journal of Statistical Software 25(4),
  \url{http://www.jstatsoft.org/v25/i04}.
    
}

\author{Michael Sekula, Somnath Datta, Susmita Datta}


\seealso{

\code{\link{clValid-class}}, \code{\link{optCluster-class}}
  
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
	
	## Plots of Internal and Stability Validation Measures
	valPlot(norm1)
	
	## Plots of Internal Validation Measures in a Single Figure
	par(mfrow = c(2,2))
	valPlot(norm1, measure = c("Dunn", "Silhouette", "Connectivity"), legend = FALSE)
	plot(2:4, valScores(norm1, "Dunn")[,,1], type="n", axes=F, xlab = "", ylab = "")
	legend("center", methodNames(norm1), col=1:9, lty=1:9, pch=paste(c(1:9)), 
			cex=0.8)		
	}
}

\keyword{cluster}
\keyword{optimize}
