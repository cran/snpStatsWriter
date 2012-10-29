\docType{package}
\name{snpStatsWriter-package}
\alias{snpStatsWriter}
\alias{snpStatsWriter-package}
\title{snpStatsWriter}
\description{
  Write snpMatrix objects to file in other formats.
}
\details{
  \tabular{ll}{ Package: \tab snpStatsWriter\cr Type: \tab
  Package\cr Version: \tab 1.0\cr Date: \tab 2012-10-29\cr
  License: \tab GPL\cr LazyLoad: \tab yes\cr }
}
\examples{
data(testdata,package="snpStats")
A.small <- Autosomes[1:6,1:10]
f <- tempfile()
## write in suitable format for snphap
nsnps <- ncol(A.small)
write.simple(A.small, a1=rep("1",nsnps), a2=rep("2",nsnps), gsep=" ",
             nullallele='0', file=f,
                write.sampleid=FALSE)
unlink(f)
}
\author{
  Maintainer: Chris Wallace <chris.wallace@cimr.cam.ac.uk>
}
\references{
  David Clayton and Hin-Tak Leung (2007). An R package for
  analysis of whole-genome association studies. Hum Hered
  64:45-51.
}
\seealso{
  \code{\link[snpStats:snpStats-package]{snpStats}}
}
\keyword{package}

