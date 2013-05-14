\name{write.snphap}
\alias{write.snphap}
\title{Write SNPHAP files}
\usage{
  write.snphap(X, a1 = NULL, a2 = NULL, file)
}
\arguments{
  \item{X}{SnpMatrix object}

  \item{a1}{vector of first allele at each SNP}

  \item{a2}{vector of second allele at each SNP}

  \item{file}{Output file name.}
}
\value{
  No return value, but has the side effect of writing
  specified output file.
}
\description{
  Simple wrapper to write.simple to write files in SNPHAP
  format
}
\details{
  If not allele codes are given, a1 and a2 will be set to 1
  and 2 for all SNPs
}
\examples{
data(testdata,package="snpStats")
A.small <- Autosomes[1:6,1:10]
f <- tempfile()
## write in suitable format for snphap
write.snphap(A.small, file=f)
unlink(f)
}
\author{
  Chris Wallace
}

