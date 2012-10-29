\name{write.impute}
\alias{write.impute}
\title{Write a snpStats object in IMPUTE format}
\usage{
  write.impute(X, a1, a2, bp, pedfile, snp.id = NULL)
}
\arguments{
  \item{pedfile}{Output file name.}

  \item{snp.id}{vector of snp ids}

  \item{X}{SnpMatrix object}

  \item{a1}{vector of first allele at each SNP}

  \item{a2}{vector of second allele at each SNP}

  \item{bp}{vector of base pair positions for each SNP}
}
\value{
  No return value, but has the side effect of writing
  specified output files.
}
\description{
  see \code{\link{write.simple}} for general information
}
\examples{
data(testdata,package="snpStats")
A.small <- Autosomes[1:6,1:10]
pf <- tempfile()

## write in suitable format for IMPUTE
nsnps <- ncol(A.small)
write.impute(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), bp=1:nsnps, pedfile=pf)
unlink(pf)
}
\author{
  Chris Wallace
}
\keyword{manip}

