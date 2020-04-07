\name{aTest}
\alias{aTest}
\title{
Hypothesis testing for annotation enrichment
}
\description{
Hypothesis testing for annotation enrichment.
}
\usage{
aTest( fitWithoutAnn, fitWithAnn, vDigit=1000 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{fitWithoutAnn}{GPA model fit without using annotation data. }
    \item{fitWithAnn}{GPA model fit with using annotation data. }
	\item{vDigit}{Number of digits for reporting parameter estimates and standard errors.
		For example, setting it to 1000 means printing out values up to three digits below zero.} 
}
\details{
\code{aTest} implements the hypothesis testing for annotation enrichment.
It requires two GPA model fits, 
one fitted with using annotation data and one fitted without using annotation data,
and evaluates annotation enrichment for risk-associated SNPs using the likelihood ratio test.
}
\value{
Returns a list with components:
    \item{q}{ q estimates. }
    \item{statistics}{ Statistics of the test for annotation enrichment. }
    \item{pvalue}{ p-value of the test for annotation enrichment. }
}
\references{
Chung D*, Yang C*, Li C, Gelernter J, and Zhao H (2014), 
"GPA: A statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data,"
PLoS Genetics, 10: e1004787. (* joint first authors)
}
\author{ Dongjun Chung }
\seealso{
\code{\link{pTest}}, \code{\link{GPA}}, 
\code{\linkS4class{GPA}}.
}
\examples{

# simulator function

simulator <- function( risk.ind, nsnp=20000, alpha=0.6 ) {
  
  m <- length(risk.ind)
  
  p.sig <- rbeta( m, alpha, 1 )
  pvec <- runif(nsnp)
  pvec[ risk.ind ] <- p.sig
  
  return(pvec)
}

# run simulation

set.seed(12345)
nsnp <- 1000
alpha <- 0.3
pmat <- matrix( NA, nsnp, 5 )

pmat[,1] <- simulator( c(1:200), nsnp=nsnp, alpha=alpha )
pmat[,2] <- simulator( c(51:250), nsnp=nsnp, alpha=alpha )
pmat[,3] <- simulator( c(401:600), nsnp=nsnp, alpha=alpha )
pmat[,4] <- simulator( c(451:750), nsnp=nsnp, alpha=alpha )
pmat[,5] <- simulator( c(801:1000), nsnp=nsnp, alpha=alpha )

ann <- rbinom(n = nrow(pmat), size = 1, prob = 0.15)
ann <- as.matrix(ann,ncol = 1)

# GPA without annotation data

fit.GPA.noAnn <- GPA( pmat, NULL, maxIter = 100 )
		
# GPA with annotation data

fit.GPA.wAnn <- GPA( pmat, ann, maxIter = 100 )

# hypothesis testing for annotation enrichment

test.annotation <- aTest( fit.GPA.noAnn, fit.GPA.wAnn )

}
\keyword{models}
\keyword{methods}
