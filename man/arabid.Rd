\name{arabid}
\alias{arabid}
\docType{data}
\title{Subset of Arabidopsis RNA-Seq Data}
\description{
A subset of RNA-Seq data studying the defense response of Arabidopsis
thaliana to a bacterial infection.  Two treatment groups were analyzed in the study: 
a mock inoculation ("mock") and a bacteria ("hrcc").  Three independent samples 
were used for each treatment.  The 200 genes in this subset were randomly selected 
without replacement from the 26,222 genes in the original dataset. The rows in this matrix
correspond to genes and the columns correspond to samples.  Each individual
cell of the matrix contains the counts of the RNA-Seq reads that are mapped to
a reference database of known genes.
}
\usage{data(arabid)}
\format{
  A 200 by 6 matrix of RNA-Seq read frequencies.
}

\source{
Di, Y., Schafer, D. W., Cumbie, J. S., & Chang, J. H. (2011). 
The NBP negative binomial model for assessing differential gene expression from RNA-Seq.
Statistical Applications in Genetics and Molecular Biology, 10, 1-28.	
}

\examples{
	
	## These examples may each take a few minutes to compute
	\donttest{
	## Obtain Dataset
	data(arabid)
			
	## Analysis of Count Data using Internal and Stability Validation Measures
	count1 <- optCluster(arabid, 2:4, clMethods = "all", countData = TRUE)
	topMethod(count1)
	
	## Analysis of Normalized Data using Internal and Stability Validation Measures
	obj <- t(t(arabid)/colSums(arabid)) ## Normalized with Respect to Library Size
	norm1 <- optCluster(obj, 2:4, clMethods = "all")
	topMethod(norm1)
	
	}

}
\keyword{datasets}
