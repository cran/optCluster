\name{clValid-class}
\docType{class}
\alias{clValidObj}
\alias{clValid-class}

\alias{clusterMethods,clValid-method}
\alias{clusters,clValid-method}
\alias{measNames,clValid-method}
\alias{measures,clValid-method}
\alias{nClusters,clValid-method}

\alias{clusterMethods}
\alias{clusters}
\alias{measNames}
\alias{measures}
\alias{nClusters}

\alias{plot,clValid,missing-method}
\alias{print,clValid-method}
\alias{show,clValid-method}
\alias{summary,clValid-method}

\title{Class "clValid"}
\description{The class \code{"clValid"} is taken from the ORPHANED \code{clValid}
  package (Brock et al., 2008). The description below is taken verbatum from the reference manual of \code{clValid}.
}
\section{Objects from the Class}{
  Objects can be created using the function \code{clValid}.
}

\section{Slots}{
  \describe{
    \item{\code{clusterObjs}:}{Object of class \code{"list"}. A list
      containing the results from the clustering methods.}
    \item{\code{measures}:}{Object of class \code{"array"}. A
      3-dimensional array which contains the
      validation measures for the clustering results.  The first
      dimension indicates the validation measures, the second the number of
      clusters, and the third the clustering methods.}
    \item{\code{measNames}:}{Object of class \code{"character"}. The
      names of the validation measures.}
    \item{\code{clMethods}:}{Object of class \code{"character"}. A
      character vector giving the clustering methods.}
    \item{\code{labels}:}{Object of class \code{"character"}. A
      character vector giving the item (gene) labels.}    
    \item{\code{nClust}:}{Object of class \code{"numeric"}. A numeric
      vector giving the numbers of clusters
      which were evaluated.}
    \item{\code{validation}:}{Object of class \code{"character"}. A character vector giving the type of
      validation measures used, consisting of some combination of
      "internal", "stability", or "biological".}
    \item{\code{metric}:}{Object of class \code{"character"}. The metric used to determine the distance
      matrix.}
    \item{\code{method}:}{Object of class \code{"character"}. For
      hierarchical clustering, the agglomeration method used.}
    \item{\code{neighbSize}:}{Object of class \code{"numeric"}. For internal validation, the neighborhood size used for the
      connectivity measure.}
    \item{\code{annotation}:}{Object of class \code{"character or array
	or list"}.
      Either a character string naming the Bioconductor annotation
      package for mapping genes to GO categories, or a list with the names of the functional classes
      and the observations belonging to each class.}
    \item{\code{GOcategory}:}{Object of class \code{"character"}. For biological validation, gives which GO
      categories to use for biological validation.  Can be one of "BP",
      "MF", "CC", or "all"}
    \item{\code{goTermFreq}:}{Object of class \code{"numeric"}. For the
      BSI, what threshold frequency of GO terms to use for functional annotation.}
    \item{\code{call}:}{Object of class \code{"call"}.  Gives the call
      to \code{clValid} used to create the \code{clValid} object.}
  }
}

\section{Methods}{
  \describe{
    \item{clusterMethods}{\code{signature(object = "clValid")}: Returns the
      names of the clustering methods. }
    \item{clusters}{\code{signature(object = "clValid")}: Returns the
      results from the clustering methods.
      
      Additional arguments:
      \describe{
	\item{\code{method = clMethods(object)}}{The clustering method(s)
	to extract.}
      }
    }
    \item{measNames}{\code{signature(object = "clValid")}: Returns the
      names of the validation measures.}
    \item{measures}{\code{signature(object = "clValid")}: Returns
      the validation measures.
      
      Additional arguments:
      \describe{
	\item{\code{measures = measNames(object)}}{ The validation measure(s)
	to extract.}
      }
    }
    \item{nClusters}{\code{signature(object = "clValid")}: Returns the
      numbers of clusters evaluated.}
    \item{plot}{\code{signature(x = "clValid", y = "missing")}: Plots
      the validation measures.
      
      Additional arguments:
      \describe{
	\item{\code{measures=measNames(x)}}{ The validation measures to
	plot.}
	\item{\code{legend=TRUE}}{ If TRUE provides a legend.}
	\item{\code{legendLoc="topright"}}{ The location of the legend.}
	\item{\code{main=NULL}}{ Title of graph.}
	\item{\code{pch=NULL}}{ Plotting characters to use.}
	\item{\code{type="b"}}{ Type of plot.}
	\item{\code{ask=prod(par("mfcol")) < length(measures) &&
	    dev.interactive()}}{ Logical.  If \code{TRUE} the user is
	prompted before each plot.}
      }
    }
    \item{print}{\code{signature(x = "clValid")}: Print method for class
      \code{clValid}.}
    \item{show}{\code{signature(object = "clValid")}: Same as print.}
    \item{summary}{\code{signature(object = "clValid")}: Summary method
      for class \code{clValid}.
      
      Additional arguments:
      \describe{
	\item{\code{digits = max(3,getOption("digits")-3)}}{ The number
	of significant digits to use.}
      }
    }
}}