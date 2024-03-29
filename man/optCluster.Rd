\name{optCluster}
\alias{optCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine Optimal Clustering Algorithm and Number of Clusters}
\description{
  \code{optCluster} performs statistical and/or biological validation of
  clustering results and determines the optimal clustering algorithm and
  number of clusters through rank aggreation.  The function returns an 
  object of class \code{"\linkS4class{optCluster}"}.
}

\usage{
  optCluster(obj, nClust, clMethods = c("clara", "diana", "hierarchical", 
  "kmeans", "model", "pam", "som", "sota"), countData = FALSE,
  validation = c("internal", "stability"), hierMethod = "average",
  annotation = NULL, clVerbose = FALSE, rankMethod = "CE", 
  distance = "Spearman", importance = NULL, rankVerbose = FALSE, ...) 
}

\arguments{
  \item{obj}{The dataset to be evaluated as either a data frame, a numeric matrix, or an 
  	\code{\link[Biobase]{ExpressionSet}} object.  Items to be clustered must be the rows of
  	the data.  In the case of data frames, all columns must be numeric.}
  \item{nClust}{A numeric vector providing the range of clusters to be evaluated (e.g. to evaluate 
  	the number of clusters ranging from 2 to 4, input 2:4).
  	A single number can also be provided.}
  \item{clMethods}{A character vector providing the names of the clustering algorithms to be used.
    The available algorithms are: "agnes", "clara", "diana", "fanny", \cr "hierarchical", 
    "kmeans", "model", "pam", "som", "sota", "em.nbinom", \cr "da.nbinom", "sa.nbinom", 
    "em.poisson", "da.poisson", "sa.poisson". Any number of selected methods is allowed.  
    The option "all" may also be used but with some caution, see \bold{Clustering Algorithms} 
    in the `Details' section for more information.}
  \item{countData}{A logical argument, indicating whether the data is count based or not.  Can also be used 
  	in conjuction with the "all" option for the 'clMethods' argument.  If TRUE and 'clMethods' = "all", 
  	all of the clustering algorithms for count data are selected: "em.nbinom", "da.nbinom", "sa.nbinom", 
  	"em.poisson", "da.poisson", "sa.poisson".  If FALSE and 'clMethods' = "all", all of the relevant clustering 
  	algorithms used with continuous data are selected: "agnes", "clara", "diana", "fanny", 
  	"hierarchical", "kmeans", "model", "pam", "som", "sota".}  
   \item{validation}{A character vector providing the names of the types of validation measures to be used. The 
  	options of "internal", "stability", "biological", and "all" are available. Any number or combination 
  	of choices is allowed.}
   \item{hierMethod}{A character string,
   	providing the agglomeration method to be used by the \cr hierarchical clustering options (\code{hclust} and \code{agnes}).  
   	Available choices are \cr "average", "complete", "single", and "ward".}
   \item{annotation}{Used in biological validation. Either a character string providing the name of the
   	Bioconductor annotation package for mapping genes to GO categories, or the names of
   	each functional class and the observations that belong to them in either a list or logical
   	matrix format.}   
   \item{clVerbose}{If TRUE, the progress of cluster validation will be produced as output.}  
  \item{rankMethod}{A character string
  	providing the method to be used for rank aggregation. The two options are the 
  	cross-entropy Monte Carlo algorithm ("CE") or Genetic algorithm ("GA").  Selection of only one method is allowed.}
  \item{distance}{A character string providing the type of distance to be used for measuring the similarity 
  	between ordered lists in rank aggregation.  The two available methods are the 
  	weighted Spearman footrule distance ("Spearman") or the weighted Kendall's tau distance ("Kendall").     
  	Selection of only one distance is allowed.}     
  \item{importance}{Vector of weights indicating the importance of each validation measure list.  
    Default of NULL represents equal weights to each validation measure.  See \bold{Weighted Rank Aggregation} 
    in the `Details' section for more information.}
  \item{rankVerbose}{If TRUE, current rank aggregation results are displayed at each iteration.}
  \item{\dots}{Additional arguments that can be passed to internal functions of \code{clValid} or \code{RankAggreg}:\cr
    Additional \code{clValid} arguments:
    \itemize{
    \item{\code{metric} - Metric used to determine distance matrix in validation measures.  Possible choices are:
    "eucliean" (default), "correlation", and "manhattan".}
    \item{\code{neighbSize} - Integer giving neighborhood size used in "connectivity" validation measure.}
    \item{\code{GOcategory} - For biological valdation, a character string providing which GO category to use.  Options include:
    "BP", "MF", "CC", or "all" (default).}
    \item{\code{goTermFreq} - For BSI validation, the threshold frequency of GO terms to used for functional annotation.}
    \item{\code{dropEvidence} - For biological validation, either NULL or a character vector of GO evidence codes to omit.}
    }
    Additional \code{RankAggreg} arguments:
    \itemize{
    \item{\code{maxIter} - The maximum number of iterations allowed. Default = 1000}
    \item{\code{k} - Size of top-k list in aggregation.}
    \item{\code{convIN} - Stopping criteria for CE and GA algorithms.  The algorithm converges once the "best" solution does not
    change after convIN iterations.  Default: 7 for CE and 30 for GA.}
    \item{\code{N} - Number of samples generated by MCMC in the CE algorithm.  Default = 10*k^2}
    \item{\code{rho} - For CE algorithm, (rho*N) is the qunatile of candidate list sorted by function values.}
    \item{\code{weight} - For CE algorithm, the learning factor used in the probability update feature.  Default = 0.25}
    \item{\code{popSize} - For GA algorithm population size in each generation. Default = 100}
    \item{\code{CP} - For GA algorithm, the crossover probability.  Default = 0.4}
    \item{\code{MP} - For GA algorithm, the mutation probability.  Default = 0.01}
    }
    }
}

\details{
  This function has been created as an extension of the \code{clValid} function.  In addition to the validation
  measures and clustering algorithms available in the \code{clValid} function, six clustering algorithms
  for count data are included in the \code{optCluster} function.  This function also determines a 
  unique solution for the optimal clustering algorithm and number of clusters through rank aggregation of 
  validation measure lists.  A brief description of the available clustering algorithms, validation measures, 
  and rank aggregation algorithms is provided below.  For more details, please refer to the references.
     
  \describe{
  	  	
  	  	\item{\bold{Clustering Algorithms:}}{A total of sixteen clustering algorithms are available for cluster analysis.
  	  	
  	  	\itemize{
  	  	\item Ten clustering algorithms for continuous data are available through the internal function \code{clValid}: 
  	  	"agnes", "clara", "diana", "fanny", "hierarchical", "kmeans", "model", "pam", "som", and "sota". \emph{ \bold{NOTE:} Some
  	  	algorithms (especially Fanny) may have difficulty finding certain numbers of clusters. If warnings or errors are
  	  	produced, the offending algorithm(s) should be removed from the \code{clMethods} argument.}

  		\item Six clustering algorithms for count data are available 
  		through the \pkg{MBCluster.Seq} \cr package: "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", 
  		and \cr "sa.poisson".  The expectation maximization (EM) algorithm, and two of its variations, 
  		the deterministic annealing (DA) algorithm and the simulated annealing (SA) algorithm, have been proposed for 
  		model-based clustering of RNA-Seq count data.  These three methods can be based on a mixture of 
  		either Poisson distributions or negative binomial distributions.  The clustering 
  		options for count data reflect both the algorithm and the distribution being used.   
  		For example, "da.nbinom" represents the deterministic annealing algorithm based on the negative 
  		binomial distribution.  
  		
  		The \pkg{MBCluster.Seq} package uses an adjustment by a normalization factor for these \cr algorithms, 
  		with the default being log(Q3) where Q3 is 75th percentile. A different 
  		normalization factor can be passed through the \code{optCluster} function by using the argument \cr
  		'Normalizer'. \cr  }}
  		
  		\item{\bold{Stability Validation Measures:}}{ Four stability validation are provided: average proportion of 
  		non-overlap (APN), average distance (AD), average distance between means (ADM), and figure of merit (FOM).  
  		These measures compare the clustering partitions established with the full data to the clustering partitions 
  		established while removing each column, one at a time.  For each measure, an average is taken over all of the 
  		removed columns, which should be minimized.
		
		  \itemize{
		  \item The APN determines the average proportion of observations placed in different clusters for
		  both cases.  The APN measure can range from 0 to 1.
		
		  \item The AD calculates the average distance between the observations assigned to the same cluster for 
		  both cases.  The AD measure can range between zero and infinity.
		
		  \item The ADM computes the average distance between the centers of clusters for observations put into 
		  the same cluster for both cases.  ADM values can range between zero and infinity. 
		
		  \item The FOM measures the average intra-cluster variance for the observations in the removed column, 
		  using clustering partitions from the remaining columns. The FOM values can range between zero and infinity. \cr 
		  }}
  		  		
  		\item{\bold{Internal Validation Measures:}}{The three internal validation measures included are:
  		connectivity, Dunn index, and silhouette width. 
  		
  		\itemize{
  		\item Connectivity measures the extent at which neighboring observations are clustered together.  
  		With a value ranging between zero and infinity, this validation measure should be minimized.  
  		
  		\item The Dunn index is the ratio of the minimum distance between observations in different clusters to the 
  		maximum cluster diameter.  With a value between zero and infinity, this measurement should be maximized. 
  		
  		\item Silhouette width is defined as the average of each observation's silhouette value.  The silhouette value 
  		is a measurement of the degree of confidence in an observation's clustering assignment.  Values near 
  		1 mean that the observation is clustered well, while values near -1 mean the observation is poorly clustered.\cr
  		 }}
 		
 		
  		\item{\bold{Biological Validation Measures:}}{The biological homogeneity index (BHI) and 
		  the biological stability index (BSI) are the two biological validation measures.  
		  They were originally proposed to provide guidance in choosing a clustering 
		  technique for microarray data, but can also be used for any other molecular 
		  expression data as well.  Both measures have a range of [0,1] and should be
		  maximized.
		
		  \itemize{
		  \item The BHI evaluates how biologically similar 
		  defined clusters are by calculating the average proportion of paired genes 
		  that are statistically clustered together and have the same functional class.  
		
		  \item The BSI examines the consistency of 
		  clustering similar biologically functioning genes together.  Observations 
		  are removed from the dataset one column at a time and the statistical cluster 
		  assignments of genes with the same functional class are 
		  compared to the cluster assignments based on the full dataset.  \cr
		  }}
		
		  \item{\bold{Rank Aggregation Algorithms:}}{ The cross-entropy Monte Carlo algorithm and the
		  Genetic algorithm are the two approaches available for rank aggregation and come from
		  the \pkg{RankAggreg} package.  Both rank aggregation algorithms can use either the weighted 
		  Spearman footrule distance or the weighted Kendall's tau to measure the "distance" between 
		  any two ordered lists.\cr}
		  
		  \item{\bold{Weighted Rank Aggregation:}}{ A list of weights for each validation measure list
	    can be included using the \code{importance} argument.  The default value of equal weights (NULL) is 
		  represented by rep(1, length(x)), where x is the character vector of validation measure names.  This
		  means each validation measure list has a weight of 1/length(x).
		  To manually change the weights, the order of the validation measures selected needs to be known.
		  The order of validation measures used in \code{optCluster} is provided below:
		  
		  \itemize{
		  \item{When selected, stability measures will ALWAYS be listed first and in the following order: "APN", "AD", "ADM", "FOM".}
		  
		  \item{When selected, internal measures will only precede biological measures.  The order of these measures is:
		  "Connectivity", "Dunn", "Silhouette".}
		  
		  \item{When selected, biological measures will always be listed last and in the following order: "BHI", "BSI".}
		  
		  }
		  }
}} 	


\value{

  \code{optCluster} returns an object of class \code{"\linkS4class{optCluster}"}.  The class description 
  is provided in the help file.

}

\references{

  Sekula, M., Datta, S., and Datta, S. (2017). 
  optCluster: An R package for determining the optimal clustering algorithm. 
  Bioinformation, 13(3), 101.
  \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5450252}

  Brock, G., Pihur, V., Datta, S. and Datta, S. (2008).
  clValid: An R Package for Cluster Validation.
  Journal of Statistical Software 25(4),
  \url{https://www.jstatsoft.org/v25/i04}.
  
  Datta, S. and Datta, S. (2003).
  Comparisons and validation of statistical clustering techniques for microarray gene expression data.
  Bioinformatics 19(4): 459-466. 
  
  Pihur, V., Datta, S. and Datta, S. (2007).
  Weighted rank aggregation of cluster validation measures: A Mounte Carlo cross-entropy approach.
  Bioinformatics 23(13): 1607-1615.
     
  Pihur, V., Datta, S. and Datta, S. (2009).
  RankAggreg, an R package for weighted rank aggregation.
  BMC Bioinformatics, 10:62,
  \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-62}.

  Si, Y., Liu, P., Li, P., & Brutnell, T. (2014).
  Model-based clustering for RNA-seq data. 
  Bioinformatics 30(2): 197-205.
    
}

%\author{Michael Sekula, Somnath Datta, Susmita Datta}
\note{

Prespecifying a list or a logical matrix of genes corresponding to functional classes for biological
validation does not require any additional packages.  If entering a list, each item in the list must 
be a vector providing the genes belonging to a specific biological class.  If entering a logical 
matrix, each column should be a logical vector indicating the genes that belong to a biological class.

If neither a list nor a logical matrix is provided, biological validation will require the \pkg{Biobase}, \cr
\pkg{annotate}, and \pkg{GO.db} packages from Bioconductor in addition to an annotation package for 
your particular data type.
See  \url{http://www.bioconductor.org} for instructions on installing these.

}

\seealso{
  For a description of the \code{clValid} function, including all available arguments that can be
  passed to it, see \code{clValid} in the \pkg{clValid} package.  

  For a description of the \code{RankAggreg} function, including all available arguments that can be 
  passed to it, see \code{\link{RankAggreg}} in the \pkg{RankAggreg} package. 
  
  For details on the clustering algorithm functions for continuous data see 
  \code{\link[cluster]{agnes}}, \code{\link[cluster]{clara}}, \code{\link[cluster]{diana}},
  \code{\link[cluster]{fanny}}, and \code{\link[cluster]{pam}} in package \pkg{cluster},
  \code{\link{hclust}} and \code{\link{kmeans}} in package \pkg{stats}, 
  \code{\link[kohonen]{som}} in package \pkg{kohonen}, 
  \code{\link[mclust]{Mclust}} in package \pkg{mclust}, 
  and \code{sota} in package \pkg{clValid}.
  
  For details the on the clustering algorithm functions for count data see
  \code{\link[MBCluster.Seq]{Cluster.RNASeq}} in package \pkg{MBCluster.Seq}.

  For details on the validation measure functions see
  \code{BHI}, \code{BSI},
  \code{stability}, \code{connectivity} and \code{dunn} 
  in package \pkg{clValid}
  and \code{\link{silhouette}} in package \pkg{cluster}. 
  
}

\examples{
	
	## These examples may each take a few minutes to compute
	\donttest{
	## Obtain Dataset	
	data(arabid)	
		
	## Analysis of Count Data using Internal and Stability Validation Measures
	count1 <- optCluster(arabid, 2:4, clMethods = "all", countData = TRUE)
	summary(count1)
	# Obtain optimal clustering assignment
	optAssign(count1)
	
	
	
	## Normalize Data with Respect to Library Size	
	obj <- t(t(arabid)/colSums(arabid))
		
	## Analysis of Normalized Data using Internal and Stability Validation Measures
	norm1 <- optCluster(obj, 2:4, clMethods = "all")
	summary(norm1)
	# Obtain optimal clustering assignment
	optAssign(norm1)
	#Obtain clustering assignment for diana with 2 clusters
	clusterResults(norm1, "diana", k = 2)$cluster
	
	## Analysis with Only UPGMA using Internal and Stability Validation Measures
	hier1 <- optCluster(obj, 2:10, clMethods = "hierarchical")
	summary(hier1)
	
	 ## Analysis of Normalized Data using All Validation Measures
	 ## Note: These lines of code require the following Bioconductor 
	 ## packages for the biological validation measures:
	 ## "Biobase", "annotate", "GO.db", and "org.At.tair.db".
	 ## If all of these packages are installed, then set
	 ## allBioconductorPackagesInstalled = TRUE
	 
	 allBioconductorPackagesInstalled = FALSE
	 if(allBioconductorPackagesInstalled){
	 require("Biobase")
	 require("annotate")
	 require("GO.db")
	 require("org.At.tair.db")
	 norm2 <- optCluster(obj, 2:4, clMethods = "all", validation = "all", 
					annotation = "org.At.tair.db")
	 summary(norm2)
	 }
	
	}
}

\keyword{cluster}
\keyword{optimize}
