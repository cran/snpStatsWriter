\name{write.beagle}
\alias{write.beagle}
\title{Write a snpStats object in beagle format}
\usage{
  write.beagle(X, a1, a2, bp, trait = NULL, gfile, mfile)
}
\arguments{
  \item{trait}{disease trait (0=missing, 1=control,
  2=case)}

  \item{gfile,mfile}{\code{gfile}=genotype file,
  \code{pedfile}=pedigree file}

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
gf <- tempfile()
mf <- tempfile()

## write in suitable format for beagle
nsnps <- ncol(A.small)
write.beagle(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), bp=1:nsnps, gfile=gf, mfile=mf)
unlink(gf)
unlink(mf)
}
\author{
  Chris Wallace
}
\keyword{manip}

