\name{write.phase}
\alias{write.phase}
\title{Write a snpStats object in PHASE/FastPHASE format}
\usage{
  write.phase(X, a1 = rep(1, ncol(X)),
    a2 = rep(2, ncol(X)), bp = NULL, file)
}
\arguments{
  \item{file}{Output file name.}

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
f <- tempfile()

## write in suitable format for PHASE
nsnps <- ncol(A.small)
write.phase(A.small, file=f)
unlink(f)
}
\author{
  Chris Wallace
}
\keyword{manip}

