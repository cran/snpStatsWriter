\name{write.mach}
\alias{write.mach}
\title{Write a snpStats object in mach format}
\usage{
  write.mach(X, a1, a2, pedfile, mfile,
    pedigree = rownames(X), member = rep(1, nrow(X)),
    father = rep(0, nrow(X)), mother = rep(0, nrow(X)),
    sex = rep("M", nrow(X)), snp.names = colnames(X))
}
\arguments{
  \item{pedfile}{Output pedigree file name.}

  \item{mfile}{Output marker file name.}

  \item{pedigree}{Optional
  pedigree/member/father/mother/sex indentifier vectors,
  same order as rows in snpStats object.  If missing,
  pedigree is set to rownames(X) and the others default to
  unrelated males}

  \item{member}{See pedigree}

  \item{father}{See pedigree}

  \item{mother}{See pedigree}

  \item{sex}{See pedigree}

  \item{snp.names}{optional SNP names to include in the
  marker map file.  Defaults to colnames(X).}

  \item{X}{SnpMatrix object}

  \item{a1}{vector of first allele at each SNP}

  \item{a2}{vector of second allele at each SNP}
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
mf <- tempfile()

## write in suitable format for MACH
nsnps <- ncol(A.small)
write.mach(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), pedfile=pf, mfile=mf)
unlink(pf)
unlink(mf)
}
\author{
  Chris Wallace
}
\keyword{manip}

