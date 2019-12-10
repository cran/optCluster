########## optCluter Internal Function ##############

## mbCountVal
## Determine validation scores for count based clustering algorithms
mbCountVal <- function(obj, Normalizer, uniqueTreatment, validation, measNames, nClust, metric, neighbSize, clVerbose, 
 					mbCountMethods, iter.max, TMP, annotation, GOcategory, goTermFreq, dropEvidence){                
		if(!requireNamespace("MBCluster.Seq")) {
      		stop("package 'MBCluster.Seq' required for clustering count data")
    	}

		if(metric=="correlation"){
    		Dist <- as.dist(1-cor(t(obj), use="pairwise.complete.obs"))
    		}  else {
  			Dist <- dist(obj,method=metric)}
  		
  		if("biological" %in% validation){	
  			if(is.character(annotation) & length(grep(".db", annotation))==0) {
    			annotation <- paste(annotation, ".db", sep="")
  			}
  			## Convert annotation list to annotation table
  			if (is.list(annotation)) {
    			if(is.null(rownames(obj))) {
      				stop("rownames of data must be present to specify biological annotation from file")
    			}
    			annotation <- annotationListToMatrix(annotation, genenames=rownames(obj))
  			}

  			if (is.null(annotation)) {
    			stop("annotation must be specified in order to use biological validation")
  			}
  			if (is.character(annotation)) {
    			if(!requireNamespace("Biobase") | !requireNamespace("GO.db") | !requireNamespace("annotate")) {
      				stop("packages 'Biobase', 'GO.db', and 'annotate' required for 2nd type of biological validation \n
					these can be downloaded from Bioconductor (www.bioconductor.org)")
    			}
  			}
  			if(is.null(rownames(obj))) {
    			stop("rownames of data must be present to use biological validation")
  			}			
		}
		## Compute Validation Measures for RNA-Seq clustering
		mbCountData <- MBCluster.Seq::RNASeq.Data(obj, Normalizer = Normalizer, Treatment = uniqueTreatment, GeneID <- rownames(obj))
  		allMeasures <- array(dim=c(length(measNames),length(nClust),length(mbCountMethods)))
  		dimnames(allMeasures) <- list(measNames,nClust,mbCountMethods)
  		colnames(allMeasures) <- nClust

		allClusterObj <- vector("list",length=length(mbCountMethods))
		names(allClusterObj) <- mbCountMethods 
        
        for(i in 1:length(mbCountMethods)){
        	
        	switch(mbCountMethods[i],
        		em.nbinom = {
        			method = "EM"
        			countModel = "nbinom"
        		},
        		da.nbinom = {
        			method = "DA"
        			countModel = "nbinom"
        		},
        		sa.nbinom = {
        			method = "SA"
        			countModel = "nbinom"
        		},
        		em.poisson = {
        			method = "EM"
        			countModel = "poisson"
        		},
        		da.poisson = {
        			method = "DA"
        			countModel = "poisson"
        		},
        		sa.poisson = {
        			method = "SA"
        			countModel = "poisson"
        		}
        		)
        	
        	MBmeasures <- matrix(0,nrow=length(measNames),ncol=length(nClust))
  			rownames(MBmeasures) <- measNames
  			colnames(MBmeasures) <- nClust

			MBclusterObj <- vector("list",length=length(nClust))
			names(MBclusterObj) <- nClust         	           
			ind <- 1
  			for (nc in nClust) {
				centers <- MBCluster.Seq::KmeansPlus.RNASeq(mbCountData, nK = nc)$centers
				noPrint <- capture.output(rnaClust <- MBCluster.Seq::Cluster.RNASeq(mbCountData, model = countModel, centers = centers, method = method,
				iter.max = iter.max, TMP = TMP))
				MBclusterObj[[ind]] <- rnaClust
        		MBcluster <- MBclusterObj[[ind]]$cluster        
        		## Avoid errors in rank aggregation
    			repCount <- 0
    			while(length(table(MBcluster))!=nc) {
					centers <- MBCluster.Seq::KmeansPlus.RNASeq(mbCountData, nK = nc)$centers
					noPrint <- capture.output(rnaClust <- MBCluster.Seq::Cluster.RNASeq(mbCountData, model = countModel, centers = centers, method = method,
						iter.max = iter.max, TMP = TMP))
					MBclusterObj[[ind]] <- rnaClust
        			MBcluster <- MBclusterObj[[ind]]$cluster
					repCount <- repCount + 1
					if(repCount == 5){					
						stop(mbCountMethods[i]," unable to find ",nc," clusters, rank aggregation cannot be performed")
					}
      			}
      			names(MBcluster) <- rownames(obj)      			

    			## internal validation measures
    			if ("internal"%in%validation) {
      				MBmeasures["Dunn",ind] <- dunn(Dist ,MBcluster)
      				MBmeasures["Silhouette",ind] <- mean(silhouette(MBcluster, dmatrix=as.matrix(Dist))[,3])
      				MBmeasures["Connectivity",ind] <- connectivity(Dist ,MBcluster, neighbSize=neighbSize)
      
      			if(clVerbose) print(paste("Finished internal validation", mbCountMethods[i], nc, "clusters"))
    			}
    
    			if("biological"%in%validation) {
      			MBmeasures["BHI",ind] <- BHI(MBcluster,annotation=annotation, names=rownames(obj),
                                 		category=GOcategory, dropEvidence=dropEvidence)
      
        		if(clVerbose & "biological"%in%validation)
        		print(paste("Finished BHI", mbCountMethods[i], nc, "clusters"))      
    			}

    			## stability validation measures
    			if ("stability"%in%validation | "biological"%in%validation) {
      				co.del <- 0 ## for use in verbose printing of progress
      				for (del in 1:ncol(obj)) {
        				objDel <- obj[,-del]
        				if(metric=="correlation") {
          					DistDel <- as.dist(1-cor(t(objDel), use="pairwise.complete.obs"))
        				} else {
          					DistDel <- dist(objDel,method=metric)
        				}
						## identify remaining columns
						uniqueTreatmentDel <- uniqueTreatment[-del]
						if(is.vector(Normalizer)){
							normalizerDel = Normalizer[-del]
						} else if(is.matrix(Normalizer)){
							normalizerDel = Normalizer[,-del]
						} else {
							normalizerDel = NULL
						}
					
						DataDel <- MBCluster.Seq::RNASeq.Data(objDel, Normalizer = normalizerDel, Treatment = uniqueTreatmentDel, GeneID <- rownames(objDel))		
						centersDel <- MBCluster.Seq::KmeansPlus.RNASeq(DataDel, nK = nc)$centers
						noPrintDel <- capture.output(clusterDel <- MBCluster.Seq::Cluster.RNASeq(DataDel, model = countModel, centers = centersDel, 
						method = method, iter.max = iter.max, TMP = TMP)$cluster)
						## Avoid errors in BSI validation
						repCountDel <- 0
    					while(length(table(clusterDel))!=nc) {
							centersDel <- MBCluster.Seq::KmeansPlus.RNASeq(DataDel, nK = nc)$centers
							noPrintDel <- capture.output(clusterDel <- MBCluster.Seq::Cluster.RNASeq(DataDel, model = countModel, centers = centersDel, 
							method = method, iter.max = iter.max, TMP = TMP)$cluster)
							repCountDel <- repCountDel + 1
							if (repCountDel == 5){
								stop(mbCountMethods[i]," unable to find ",nc," clusters when removing column ",del, 
								". \n  Stability measures and/or BSI cannot be calculated")									
							}
      					} 
        				names(clusterDel) <- rownames(objDel) 

        			if("stability"%in%validation) {
          				stabmeas <- stability(obj, Dist, del, MBcluster, clusterDel)
          				MBmeasures["APN",ind] <- MBmeasures["APN",ind] + stabmeas["APN"]
          				MBmeasures["AD",ind]  <- MBmeasures["AD",ind]  + stabmeas["AD"]
          				MBmeasures["ADM",ind] <- MBmeasures["ADM",ind] + stabmeas["ADM"]
          				MBmeasures["FOM",ind] <- MBmeasures["FOM",ind] + stabmeas["FOM"]
        				}
        			if("biological"%in%validation) {        				
						tmp <- BSI(MBcluster,clusterDel,annotation=annotation,
                     		names=rownames(objDel), category=GOcategory, goTermFreq=goTermFreq,
                     		dropEvidence=dropEvidence)
         	 			MBmeasures["BSI",ind] <- MBmeasures["BSI",ind] + tmp
        			}
        
        			if (del/ncol(obj) > 0.25 & co.del==0) 
          			{
            			if(clVerbose & "stability"%in%validation) 
              				print(paste("Stability validation 25% finished", mbCountMethods[i], nc, "clusters"))
            			if(clVerbose & "biological"%in%validation) 
             		 		print(paste("BSI 25% finished", mbCountMethods[i], nc, "clusters"))            
            			co.del <- co.del+1
          			}
        			else if (del/ncol(obj) > 0.50 & co.del==1) 
          			{
            			if(clVerbose & "stability"%in%validation) 
              				print(paste("Stability validation 50% finished", mbCountMethods[i], nc, "clusters"))
            			if(clVerbose & "biological"%in%validation) 
              				print(paste("BSI 50% finished", mbCountMethods[i], nc, "clusters"))            
            			co.del <- co.del+1
          			}
        			else if (del/ncol(obj) > 0.75 & co.del==2) 
          			{
            			if(clVerbose & "stability"%in%validation) 
              				print(paste("Stability validation 75% finished", mbCountMethods[i], nc, "clusters"))
            			if(clVerbose & "biological"%in%validation) 
              				print(paste("BSI 75% finished", mbCountMethods[i], nc, "clusters"))            
            			co.del <- co.del+1
          			}
      			} #END OF del LOOP
      			if(clVerbose & "stability"%in%validation)
        			print(paste("Finished stability validation", mbCountMethods[i], nc, "clusters"))
      			if(clVerbose & "biological"%in%validation)
      				print(paste("Finished BSI", mbCountMethods[i], nc, "clusters"))
    		} #END of STABILITY measures
    		ind <- ind+1  #ind tracks number clusters
  		} #END OF NC LOOP
  
  		if ("stability"%in%validation) {
    		MBmeasures["APN",] <- MBmeasures["APN",]/ncol(obj)
    		MBmeasures["AD",] <-  MBmeasures["AD",]/ncol(obj)
    		MBmeasures["ADM",] <- MBmeasures["ADM",]/ncol(obj)
    		MBmeasures["FOM",] <- MBmeasures["FOM",]/ncol(obj)
  		}
  		if ("biological"%in%validation) {
    		MBmeasures["BSI",] <- MBmeasures["BSI",]/ncol(obj)
  		}
  		allClusterObj[[i]] <- MBclusterObj
  		allMeasures[,,i] <- MBmeasures
  	} # END OF MBCOUNT METHODS LOOP	
	list(clusterObj=allClusterObj, measures=allMeasures)		
}

## Updated code from original clValid (Brock et al., 2008) package
## Fixed vClusters function to avoid multiple somgrid functions
## causing "Error in !toroidal : invalid argument type".
## The change is in the internal vClusters() functions which
## is called by the clValid() function.
##
## Previous fix was to use older version of kohonen package:
## require(devtools)
## install_version("kohonen", version = "2.0.19", repos = "http://cran.us.r-project.org")
##
## Also include changes for enviornment variable in R 4.0.0
## _R_CLASS_MATRIX_ARRAY

clValid2 <- function(obj, nClust, clMethods="hierarchical", validation="stability", maxitems=600,
                    metric="euclidean", method="average", neighbSize=10,
                    annotation=NULL, GOcategory="all", goTermFreq=0.05,
                    dropEvidence=NULL, verbose=FALSE, ...) {
  
  clMethods <- tolower(clMethods)
  clMethods <- match.arg(clMethods,c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes"), several.ok=TRUE)
  
  validation <- match.arg(validation,c("stability","internal","biological"),several.ok=TRUE)
  metric <- match.arg(metric,c("euclidean", "correlation", "manhattan")) ## used for hierarchical, diana, fanny, agnes, pam
  method <- match.arg(method,c("ward", "single", "complete", "average")) ## for hclust, agnes
  GOcategory <- match.arg(GOcategory, c("all","BP","CC","MF"))
  
  ##########################################
  ## Fix to accomodate _R_CLASS_MATRIX_ARRAY (12/10/19)
  ##########################################
  if(is(obj, "matrix")){
    mat <- obj
  }else{
    switch(class(obj), ExpressionSet = mat <- Biobase::exprs(obj), 
           data.frame = {
             if (any(!sapply(obj, class) %in% c("numeric", "integer"))) stop("data frame 'obj' contains non-numeric data")
             mat <- as.matrix(obj)
           }, 
           stop("argument 'obj' must be a matrix, data.frame, or ExpressionSet object"))	
  }
  
  if (nrow(mat)>maxitems) {
    if (interactive()) {
      cat("\nThe number of items to be clustered is larger than 'maxitems'\n")
      cat("The memory and time required may be excessive, do you wish to continue?\n")
      cat("(y to continue, any other character to exit) ")
      ans <- tolower(substr(readLines(n=1),1,1))
      if (ans!="y") {
        stop("Exiting clValid, number of items exceeds 'maxitems'")
      }
    } else {
      stop("The number of items to be clustered is larger than 'maxitems'\n  Either decrease the number of rows (items) or increase 'maxitems'\n")
    }
  }
  
  
  if ("clara"%in%clMethods & metric=="correlation")
    warning("'clara' currently only works with 'euclidean' or 'manhattan' metrics - metric will be changed to 'euclidean'  ")
  
  if(is.character(annotation) & length(grep(".db", annotation))==0) {
    annotation <- paste(annotation, ".db", sep="")
  }
  ## Convert annotation list to annotation table
  
  if (is.list(annotation)) {
    if(is.null(rownames(mat))) {
      stop("rownames of data must be present to specify biological annotation from file")
    }
    annotation <- annotationListToMatrix(annotation, genenames=rownames(mat))
  }
  
  if ("biological"%in%validation & is.null(annotation)) {
    stop("annotation must be specified in order to use biological validation")
  }
  if ("biological"%in%validation & is.character(annotation)) {
    if(!requireNamespace("Biobase") | !requireNamespace("GO.db") | !requireNamespace("annotate")) {
      stop("packages 'Biobase', 'GO.db', and 'annotate' required for 2nd type of biological validation \n
           these can be downloaded from Bioconductor (www.bioconductor.org)")
    }
    }
  if("biological"%in%validation & is.null(rownames(mat))) {
    stop("rownames of data must be present to use biological validation")
  }
  
  #  if (!is.matrix(mat) | !is.numeric(mat))
  #    stop("argument 'mat' must be a numeric matrix")
  
  nClust <- floor(nClust)
  if (any(nClust<1))
    stop("argument 'nClust' must be a positive integer vector")
  
  if(metric=="correlation")
    Dist <- as.dist(1-cor(t(mat), use="pairwise.complete.obs"))  else
      Dist <- dist(mat,method=metric)
  
  clusterObjs <- vector("list",length(clMethods))
  names(clusterObjs) <- clMethods
  
  measures <- c(if("stability"%in%validation) c("APN","AD","ADM","FOM"),
                if("internal"%in%validation) c("Connectivity","Dunn","Silhouette"),
                if("biological"%in%validation) c("BHI","BSI"))
  validMeasures <- array(dim=c(length(measures),length(nClust),length(clMethods)))
  dimnames(validMeasures) <- list(measures,nClust,clMethods)
  
  for (i in 1:length(clMethods)) {
    
    cvalid <- vClusters2(mat,clMethods[i],nClust, validation=validation,
                        Dist=Dist, method=method, metric=metric, annotation=annotation,
                        GOcategory=GOcategory, goTermFreq=goTermFreq, neighbSize=neighbSize,
                        dropEvidence=dropEvidence, verbose=verbose, ...)
    clusterObjs[[i]] <- cvalid$clusterObj
    validMeasures[,,i] <- cvalid$measures
  }
  
  if(is.null(rownames(mat))) {
    rownames(mat) <- 1:nrow(mat)
    warning("rownames for data not specified, using 1:nrow(data)")
  }
  
  new("clValid", clusterObjs=clusterObjs, measures=validMeasures, measNames=measures,
      clMethods=clMethods, labels=rownames(mat), nClust=nClust, validation=validation,
      metric=metric,method=method, neighbSize=neighbSize,  GOcategory=GOcategory,
      goTermFreq=goTermFreq, annotation=annotation,
      call=match.call())
}


## vClusters used inside clValid
## Updated to work with kohonen_3.0.0 or later

vClusters2 <- function(mat,clMethod,nClust,nclustMax, validation,
                      Dist, method, metric, annotation, GOcategory,
                      goTermFreq, neighbSize, dropEvidence, verbose, ... ) {
  
  measNames <- c(if("stability"%in%validation) c("APN","AD","ADM","FOM"),
                 if("internal"%in%validation) c("Connectivity","Dunn","Silhouette"),
                 if("biological"%in%validation) c("BHI","BSI"))
  measures <- matrix(0,nrow=length(measNames),ncol=length(nClust))
  rownames(measures) <- measNames
  colnames(measures) <- nClust
  
  switch(clMethod,
         hierarchical = {
           clusterObj <- hclust(Dist,method)
         },
         diana = {
           clusterObj <- diana(Dist, ...)
         },
         kmeans = {
           clusterObj <- vector("list",length=length(nClust))
           names(clusterObj) <- nClust
           clusterObjInit <- hclust(Dist,method)
         },
         agnes = {
           clusterObj <- agnes(Dist, method=method, ...)
         },
         ## otherwise - sota, fanny, som, model, pam, clara
         { clusterObj <- vector("list",length=length(nClust))
         names(clusterObj) <- nClust })
  
  ind <- 1
  for (nc in nClust) {
    switch(clMethod,
           kmeans = {
             initial <- tapply(mat, list(rep(cutree(clusterObjInit,nc),ncol(mat)),col(mat)),
                               function(x) mean(x, na.rm=TRUE))
             if(length(dup <- which(duplicated(initial)))>0) {
               for(dupi in dup) 
                 initial[dupi,] <- initial[dupi,] + jitter(initial[dupi,])
             }
             dimnames(initial) <- list(NULL,dimnames(mat)[[2]])
             clusterObj[[ind]] <- kmeans(mat,initial,...)
             cluster <- clusterObj[[ind]]$cluster
           },
           fanny = {
             clusterObj[[ind]] <- fanny(Dist, nc, ...)
             cluster <- clusterObj[[ind]]$clustering
           },
           model = {
             clusterObj[[ind]] <- Mclust(mat,nc, ...)
             cluster <- clusterObj[[ind]]$classification
           },
           som = {
             clusterObj[[ind]] <- som(mat, grid=kohonen::somgrid(1,nc), ...) #added kohonen:: to avoid "Error in !toroidal" (12/10/19)
             cluster <- clusterObj[[ind]]$unit.classif
           },
           pam = {
             clusterObj[[ind]] <- pam(Dist, nc, ...)
             cluster <- clusterObj[[ind]]$clustering
           },
           clara = {
             clusterObj[[ind]] <- clara(mat, nc, metric=ifelse(metric=="correlation","euclidean",metric), ...)
             cluster <- clusterObj[[ind]]$clustering
           },
           sota = {
             clusterObj[[ind]] <- sota(mat,nc-1)
             cluster <- clusterObj[[ind]]$clust
             
           },
           ## otherwise - hierarchical, diana, agnes
           {cluster <- cutree(clusterObj,nc)})
    
    if(length(table(cluster))!=nc) {
      warning(paste(clMethod, "unable to find",nc,"clusters, returning NA for these validation measures"))
      measures[,ind] <- NA
      ind <- ind+1
      next()
    }
    
    ## internal validation measures
    if ("internal"%in%validation) {
      measures["Dunn",ind] <- dunn(Dist ,cluster)
      measures["Silhouette",ind] <- mean(silhouette(cluster, dmatrix=as.matrix(Dist))[,3])
      measures["Connectivity",ind] <- connectivity(Dist ,cluster, neighbSize=neighbSize)
      if(verbose) print(paste("Finished internal validation,", clMethod, nc, "clusters"))
    }
    
    if("biological"%in%validation) {
      measures["BHI",ind] <- BHI(cluster,annotation=annotation, names=rownames(mat),
                                 category=GOcategory, dropEvidence=dropEvidence)
      if(verbose & "biological"%in%validation)
        print(paste("Finished BHI,", clMethod, nc, "clusters"))      
    }
    
    ## stability validation measures
    if ("stability"%in%validation | "biological"%in%validation) {
      co.del <- 0 ## for use in verbose printing of progress
      for (del in 1:ncol(mat)) {
        matDel <- mat[,-del]               ## matDel <- as.matrix(matDel)
        if(metric=="correlation") {
          DistDel <- as.dist(1-cor(t(matDel), use="pairwise.complete.obs"))
        } else {
          DistDel <- dist(matDel,method=metric)
        }
        switch(clMethod,
               hierarchical = clusterObjDel <- hclust(DistDel,method),
               kmeans = clusterObjInitDel <- hclust(DistDel,method),
               diana = clusterObjDel <- diana(DistDel, ...),
               agnes = clusterObjDel <- agnes(DistDel, method=method, ...),
               clara = clusterObjDel <- clara(matDel,nc,metric=ifelse(metric=="correlation","euclidean",metric), ...))
        
        
        switch(clMethod,
               kmeans = {
                 initialDel <- tapply(matDel, list(rep(cutree(clusterObjInitDel,nc),
                                                       ncol(matDel)), col(matDel)),
                                      function(x) mean(x, na.rm=TRUE))
                 if(length(dup <- which(duplicated(initialDel)))>0) {
                   for(dupi in dup) 
                     initialDel[dupi,] <- initialDel[dupi,] + jitter(initialDel[dupi,])
                 }
                 dimnames(initialDel) <- list(NULL,dimnames(matDel)[[2]])
                 kmdel <- kmeans(matDel,initialDel, ...)
                 clusterDel <- kmdel$cluster
               },
               fanny = {
                 hfdel <- fanny(DistDel, nc, ...)
                 clusterDel <- hfdel$clustering
               },
               model = {
                 clusterDel <- Mclust(matDel,nc, ...)$classification
               },
               som = {
                 hsdel <- try(som(matDel, grid=kohonen::somgrid(1,nc), ...)) #added kohonen:: to avoid "Error in !toroidal" (12/10/19)
                 clusterDel <- hsdel$unit.classif
               },
               pam = {
                 clusterDel <- pam(DistDel, nc, cluster.only=TRUE, ...)
               },
               clara = {
                 clusterDel <- clusterObjDel$clustering
               },
               sota = {
                 clusterDel <- sota(matDel,nc-1)$clust
               },
               ## otherwise - hierarchical, diana, agnes
               {clusterDel <- cutree(clusterObjDel,nc)})
        
        if("stability"%in%validation) {
          stabmeas <- stability(mat, Dist, del, cluster, clusterDel)
          measures["APN",ind] <- measures["APN",ind] + stabmeas["APN"]
          measures["AD",ind]  <- measures["AD",ind]  + stabmeas["AD"]
          measures["ADM",ind] <- measures["ADM",ind] + stabmeas["ADM"]
          measures["FOM",ind] <- measures["FOM",ind] + stabmeas["FOM"]
        }
        if("biological"%in%validation) {
          tmp <- BSI(cluster,clusterDel,annotation=annotation,
                     names=rownames(mat), category=GOcategory, goTermFreq=goTermFreq,
                     dropEvidence=dropEvidence)
          measures["BSI",ind] <- measures["BSI",ind] + tmp
        }
        ## VERBOSE printing
        if (del/ncol(mat) > 0.25 & co.del==0) 
        {
          if(verbose & "stability"%in%validation) 
            print(paste("Stability validation 25% finished,", clMethod, nc, "clusters"))
          if(verbose & "biological"%in%validation) 
            print(paste("BSI 25% finished,", clMethod, nc, "clusters"))            
          co.del <- co.del+1
        }
        else if (del/ncol(mat) > 0.50 & co.del==1) 
        {
          if(verbose & "stability"%in%validation) 
            print(paste("Stability validation 50% finished,", clMethod, nc, "clusters"))
          if(verbose & "biological"%in%validation) 
            print(paste("BSI 50% finished,", clMethod, nc, "clusters"))            
          co.del <- co.del+1
        }
        else if (del/ncol(mat) > 0.75 & co.del==2) 
        {
          if(verbose & "stability"%in%validation) 
            print(paste("Stability validation 75% finished,", clMethod, nc, "clusters"))
          if(verbose & "biological"%in%validation) 
            print(paste("BSI 75% finished,", clMethod, nc, "clusters"))            
          co.del <- co.del+1
        }
      } #END OF del LOOP
      if(verbose & "stability"%in%validation)
        print(paste("Finished stability validation,", clMethod, nc, "clusters"))
      if(verbose & "biological"%in%validation)
        print(paste("Finished BSI,", clMethod, nc, "clusters"))
    } #END of STABILITY measures
    ind <- ind+1  #ind tracks number clusters
    ## if(verbose) print(paste("Finished with", nc, "clusters"))
  } #END OF NC LOOP
  
  if ("stability"%in%validation) {
    measures["APN",] <- measures["APN",]/ncol(mat)
    measures["AD",] <-  measures["AD",]/ncol(mat)
    measures["ADM",] <- measures["ADM",]/ncol(mat)
    measures["FOM",] <- measures["FOM",]/ncol(mat)  ## little different from Yeung paper (doesn't do this)
  }
  if ("biological"%in%validation) {
    measures["BSI",] <- measures["BSI",]/ncol(mat)
  }
  
  list(clusterObj=clusterObj, measures=measures)
}