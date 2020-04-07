
aTest <- function( fitWithoutAnn, fitWithAnn, vDigit=1000 ) {
	# check correctness of arguments

	if ( !is( fitWithoutAnn, "GPA" ) ) {
		stop( " Input for 'fitWithoutAnn' argument is not 'GPA' class object. Please check the input." )
	}

	if ( !is( fitWithAnn, "GPA" ) ) {
		stop( " Input for 'fitWithAnn' argument is not 'GPA' class object. Please check the input." )
	}

	if ( vDigit %% 10 != 0 | vDigit <= 0 ) {
		stop( "Inappropriate value for 'vDigit' argument. It should be multiples of 10, e.g., 10, 100, ..." )
	}

	# load fits

	empiricalNull <- get_setting(fitWithoutAnn)$empiricalNull

	pis <- get_fit(fitWithoutAnn)$pis
	betaAlpha <- get_fit(fitWithoutAnn)$betaAlpha
	if ( empiricalNull ) {
		betaAlphaNull <- get_fit(fitWithoutAnn)$betaAlphaNull
	}

	# constant

	nBin <- nrow(get_gwasPval(fitWithAnn))
	nGWAS <- ncol(get_gwasPval(fitWithAnn))

	binaryList <- vector( "list", nGWAS )
	for ( k in seq_len(nGWAS) ) {
		binaryList[[k]] <- c( 0, 1 )
	}
	binaryMat <- expand.grid( binaryList )

	nComp <- nrow(binaryMat)
	combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )

	#annMat <- as.matrix(get_annMat(fitWithAnn))
	nAnn <- ncol(get_annMat(fitWithAnn))

	nPis <- length(pis)
	nAlpha <- length(betaAlpha)

	# LRT: numerator

	q1 <- matrix( colSums(get_annMat(fitWithAnn)) / nrow(get_annMat(fitWithAnn)),
								nAnn, nComp )

	betaDist <- matrix( NA, nBin, nComp )
	for ( k in seq_len(nGWAS) ) {
		betaDist[,k] <- dbeta( get_gwasPval(fitWithAnn)[,k], betaAlpha[k], 1 )
	}
	if ( empiricalNull ) {
		betaNullDist <- matrix( NA, nBin, nComp )

		for ( k in seq_len(nGWAS) ) {
			betaNullDist[,k] <- dbeta( get_gwasPval(fitWithAnn)[,k], betaAlphaNull[k], 1 )
		}
	}

    llmat <- matrix( NA, nBin, nComp )

    for ( g in seq_len(nComp) ) {
        # mixing proportion

	    llmat[,g] <- pis[g]

        # emission for GWAS

	    for ( k in seq_len(nGWAS) ) {
	        if( binaryMat[g,k] == 1 ) {
	        	# signal SNP

	        	llmat[,g] <- llmat[,g] * betaDist[,k]
        	} else if ( empiricalNull ) {
	        	# null SNP, if null is empirically estimated

	        	llmat[,g] <- llmat[,g] * betaNullDist[,k]
        	}
	    }

    	# emission for annotation data

    	for ( d in seq_len(nAnn) ) {
        	llmat[,g] <- llmat[,g] *
        		c( ( 1 - q1[d,g] ), q1[d,g] )[ ( get_annMat(fitWithAnn)[,d] + 1 ) ]
    	}
    }

	ll0 <- sum( log( rowSums( llmat ) ) )

	# LRT: denominator

	ll1 <- get_fit(fitWithAnn)$loglik[ length(get_fit(fitWithAnn)$loglik) ]

	# LRT

	LRT <- -2 * ( ll0 - ll1 )

	# calculate p-value

	degree.of.freedom <- nAnn * ( nComp - 1 )
	pvalue <- pchisq( LRT, degree.of.freedom, lower.tail=FALSE )

	# SE for q1

	cov.GPA.wAnn <- cov( fitWithAnn, silent=TRUE )
	seVec <- sqrt(diag(cov.GPA.wAnn))

	# summary

	cat( "Hypothesis testing for annotation enrichment\n" )
	cat( "( Note: This version of test is designed for single annotation data )\n" )
	cat( "--------------------------------------------------\n" )

	cat( "q:\n" )
	cat( "GWAS combination: ", paste( combVec, collapse=" " ), "\n" )
	for ( d in seq_len(nAnn) ) {
		cat( "Annotation #",d,":\n")
    	cat( "    ", paste( round(get_fit(fitWithAnn)$q1[d,]*vDigit)/vDigit,
			collapse=" " ), "\n" )

    	if ( empiricalNull ) {
    		locQ1 <- (nPis-1+2*nAlpha+nPis*(d-1)+1):(nPis-1+2*nAlpha+nPis*d)
		} else {
			locQ1 <- (nPis-1+nAlpha+nPis*(d-1)+1):(nPis-1+nAlpha+nPis*d)
		}
    	cat( "  ( ", paste( round(seVec[locQ1]*vDigit)/vDigit, collapse=" " ), " )\n" )
	}

	cat( "\n" )
	cat( "Ratio of q over baseline (",combVec[1],"):\n" )
	cat( "GWAS combination: ", paste( combVec[-1], collapse=" " ), "\n" )
	q1ratio <- q1ratioSE <- matrix( NA, nrow(get_fit(fitWithAnn)$q1),
																			(ncol(get_fit(fitWithAnn)$q1)-1) )
	for ( d in seq_len(nAnn) ) {
		# estimates

		q1ratio[d,] <- get_fit(fitWithAnn)$q1[d,-1] / get_fit(fitWithAnn)$q1[d,1]

		# SE

		for ( j in seq(2,ncol(q1)) ) {
			qderiv <- rep( 0, nrow(get_fit(fitWithAnn)$q1)*ncol(get_fit(fitWithAnn)$q1) )
			qderiv[ ncol(get_fit(fitWithAnn)$q1) * (d-1) + 1 ] <- - get_fit(fitWithAnn)$q1[d,j] /
							get_fit(fitWithAnn)$q1[d,1]^2
			qderiv[ ncol(get_fit(fitWithAnn)$q1) * (d-1) + j ] <- 1 / get_fit(fitWithAnn)$q1[d,1]

			gderiv <- as.matrix( c( rep( 0, nrow(cov.GPA.wAnn) - nrow(get_fit(fitWithAnn)$q1)*
													ncol(get_fit(fitWithAnn)$q1) ), qderiv ) )
			q1ratioSE[d,(j-1)] <- sqrt( t(gderiv) %*% cov.GPA.wAnn %*% gderiv )
		}
	}
	for ( d in seq_len(nAnn) ) {
		cat( "Annotation #",d,":\n")
		cat( "    ", paste( round(q1ratio[d,]*vDigit)/vDigit, collapse=" " ), "\n" )
		cat( "  ( ", paste( round(q1ratioSE[d,]*vDigit)/vDigit, collapse=" " ), " )\n" )
	}
	cat( "\n" )
    cat( "test statistics: ", paste( round(LRT*vDigit)/vDigit, collapse=" " ), "\n" )
    cat( "p-value: ", pvalue, "\n", collapse=" " )
	cat( "--------------------------------------------------\n" )

	return( list( q=get_fit(fitWithAnn)$q1, statistics=LRT, pvalue=pvalue ) )
}
