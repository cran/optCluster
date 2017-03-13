\name{optHeatmap}
\alias{optHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimal Hierarchical Clustering Heat Map}
\description{
  \code{optHeatmap} creates a heat map from an object of
  class \code{"optCluster"} if the optimal clustering algorithm
  is one of the available hierarchical algorithms: 
  "hierarchical", "agnes", or "diana".
  
}
\usage{
optHeatmap(x, dendroClusters = TRUE, barClusters = FALSE, 
			clusterColors = "rainbow", 
			mapColors = colorRampPalette(c("green", "black", "red"))(256), 
			Colv = FALSE, dendrogram = "row", density.info = "none", ...)
}


\arguments{
  \item{x}{An object of class \code{"\linkS4class{optCluster}"}.}
	\item{dendroClusters}{If \code{TRUE}, the row dendrogram is colored to represent
	the optimal k clusters.}
	\item{barClusters}{If \code{TRUE}, a colored sidebar is produced down the rows 
	to represent the optimal k clusters.}
	\item{clusterColors}{Character vector listing the k colors for the clustering
	paritions.  Default is \cr rainbow(k) where k is the optimal number
	of clusters.}
	\item{mapColors}{Specifies colors to use for heat map image.}
 	\item{Colv}{Determines if and how the column dendrogram should be reordered.
 	If TRUE, a dendrogram is computed and the columns are reordered by means.
	A \code{dendrogram} object or vector of integers can also be used to specify the 
	reordering.}
	\item{dendrogram}{Character string specifying which dendrogram(s) to display. 
		Options include: "none", "row", "column", "both".}
	\item{density.info}{Character string specifying the type of plot to superimpose
		on the color-key.  Options include: "none", "histogram", "density".}
	\item{\dots}{Additional plotting parameters from the \code{heatmap.2} function.}  
}

\details{
	This function utilizes the \code{heatmap.2} function to produce a heat map 
	based on the clustering results from the optimal clustering algorithm and number of
	clusters (as determined by the \code{optCluster} function).  The rows in the heat map 
	are ordered based on the optimal hierarchical clustering algorithm, with the 
	corresponding dendrogram displayed to the left side of the map.  The clustering partition 
	of the rows into the optimal number of clusters can be visualized using 
	the \cr'dendroClusters' argument and/or the 'barClusters' argument. \cr
	  
	Because the \code{optCluster} function performs a cluster analysis on 
	the rows of the dataset, only the rows are reordered as default in the \code{optHeatmap}
	function.  However, columns can be reordered and displayed with a dendrogram using
    the 'Colv' and 'dendrogram' arguments, respectively.  
	The heat map can also be further customized by passing additional agruments to the  
	\code{heatmap.2} function through the \code{optHeatmap} function.
	}

%\author{Michael Sekula, Somnath Datta, Susmita Datta}

\seealso{

\code{\link[gplots]{heatmap.2}}, \code{\link{optCluster-class}}
  
}

\examples{
	## Obtain Dataset	
	data(arabid)			
	
	## Normalize Data with Respect to Library Size	
	obj <- t(t(arabid)/colSums(arabid))	
	
	## Analysis with Only UPGMA using Internal and Stability Validation Measures
	hier1 <- optCluster(obj, 2:10, clMethods = "hierarchical")
	topMethod(hier1)

	## Create Default Heat Map Reordering Rows Only
	optHeatmap(hier1)
	
	## Create a Heat Map Reordering Both Rows and Columns
	optHeatmap(hier1, Colv = TRUE, dendrogram = "both")
	
	## Customized Heat Map Using Several heatmap.2 Arguments
	optHeatmap(hier1, Colv = TRUE, dendrogram = "both", 
	labRow = "", cexCol = 1.0, keysize = 1)


}
