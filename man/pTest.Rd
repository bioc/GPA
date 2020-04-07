\name{pTest}
\alias{pTest}
\title{
Hypothesis testing for pleiotropy
}
\description{
Hypothesis testing for pleiotropy.
}
\usage{
pTest( fit, fitH0, vDigit=1000 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{fit}{Fit of the GPA model of interest. }
    \item{fitH0}{GPA model fit under the null hypothesis of pleiotropy test. }
	\item{vDigit}{Number of digits for reporting parameter estimates and standard errors.
		For example, setting it to 1000 means printing out values up to three digits below zero.} 
}
\details{
\code{pTest} implements the hypothesis testing for pleiotropy.
It requires two GPA model fits, one of interest and one under the null hypothesis
(obtained by setting \code{pleiotropyH0=TRUE} when running \code{GPA} function),
and evaluates genetical correlation among multiple phenotypes using the likelihood ratio test.
}
\value{
Returns a list with components:
    \item{pi}{ pi estimates. }
    \item{piSE}{ Standard errors for pi estimates. }
    \item{statistics}{ Statistics of the pleiotropy test. }
    \item{pvalue}{ p-value of the pleiotropy test. }
}
\references{
Chung D*, Yang C*, Li C, Gelernter J, and Zhao H (2014), 
"GPA: A statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data,"
PLoS Genetics, 10: e1004787. (* joint first authors)
}
\author{ Dongjun Chung }
\seealso{
\code{\link{aTest}}, \code{\link{GPA}}, 
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

# GPA without annotation data

fit.GPA.noAnn <- GPA( pmat, NULL, maxIter = 100 )

# GPA under the null hypothesis of pleiotropy test

fit.GPA.pleiotropy.H0 <- GPA( pmat, NULL, pleiotropyH0=TRUE, maxIter = 100 )

# hypothesis testing for pleiotropy

test.pleiotropy <- pTest( fit.GPA.noAnn, fit.GPA.pleiotropy.H0 )

}
\keyword{models}
\keyword{methods}
