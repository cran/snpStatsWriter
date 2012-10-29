\name{write.mach}
\alias{write.mach}
\title{Write a snpStats object in mach format}
\usage{
  write.mach(X, a1, a2, pedfile, mfile,
    pedigree = rownames(X), member = rep(1, nrow(X)),
    father = rep(0, nrow(X)), mother = rep(0, nrow(X)),
    sex = rep("M", nrow(X)))
}
\arguments{
  \item{pedfile,mfile}{Output file names.
  \code{pedfile}=pedigree file, \code{mfile}=marker file.}

  \item{pedigree,member,father,mother,sex}{Optional
  pedigree/member/father/mother/sex indentifier vectors,
  same order as rows in snpStats object.  If missing,
  pedigree is set to rownames(X) and the others default to
  unrelated males}

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

