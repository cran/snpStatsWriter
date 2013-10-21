\name{write.sbams}
\alias{write.sbams}
\title{write.sbams}
\usage{
  write.sbams(X, response, file)
}
\arguments{
  \item{response}{vector or matrix of response variables.
  rows index subjects, columns index variables}

  \item{X}{SnpMatrix object}

  \item{file}{Output file name.}
}
\value{
  No return value, but has the side effect of writing
  specified output file.
}
\description{
  write an sbams format file
}
\details{
  sbams is software from Xiaoquan Wen at
  https://github.com/xqwen/sbams
}
\examples{
data(testdata,package="snpStats")
A.small <- Autosomes[1:6,1:10]
R <- matrix(rnorm(12),ncol=2)
colnames(R) <- c("var1","var2")
f <- tempfile()

## write in suitable format for sbams
write.sbams(X=A.small, response=R, file=f)
unlink(f)
}
\author{
  Chris Wallace
}
\keyword{manip}

