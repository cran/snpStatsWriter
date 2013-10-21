\name{write.simple}
\alias{write.simple}
\title{Fast and flexible writing of snpStats objects to flat files}
\usage{
  write.simple(X, a1, a2, file, fsep = "\t", gsep = "",
    nullallele = "N", write.header = TRUE,
    transpose = FALSE, write.sampleid = TRUE, bp = NULL,
    num.coding = FALSE)
}
\arguments{
  \item{X}{SnpMatrix object}

  \item{a1}{vector of first allele at each SNP}

  \item{a2}{vector of second allele at each SNP}

  \item{bp}{vector of base pair positions for each SNP}

  \item{fsep,gsep}{Field and genotype separators.}

  \item{nullallele}{Character to use for missing alleles}

  \item{file}{Output file name.}

  \item{write.header}{Write a header line}

  \item{transpose}{Output SNPs as rows, samples as columns
  if \code{TRUE}.  The default is samples as rows, SNPs as
  columns, as represented internally by
  snpStats/SnpMatrix.}

  \item{write.sampleid}{Output sample ids}

  \item{num.coding}{Use alleles 1 and 2 instead of
  supplying allele vectors.}
}
\value{
  No return value, but has the side effect of writing
  specified output files.
}
\description{
  Different genetics phasing and analysis programs (beagle,
  mach, impute, snptest, phase/fastPhase, snphap, etc) have
  different requirements for input files.  These functions
  aim to make creating these files from a SnpMatrix object
  straightfoward.
}
\details{
  It's written in C, so should be reasonably fast even for
  large datasets.

  \code{write.simple} is the most flexible function.  It
  should be able to write most rectangular based formats.

  Additional functions are available tailored to software
  that require a bit more than a rectangular format:
  \code{\link{write.beagle}}, \code{\link{write.impute}},
  \code{\link{write.mach}}, \code{\link{write.phase}}.
}
\note{
  This has been tested with \code{SnpMatrix} objects from
  the package \code{snpStats} but should also work with
  \code{snp.matrix} objects from the package
  \code{snpMatrix}.
}
\section{Warning}{
  Any uncertain genotypes (stored by snpStats as raw codes
  4 to 253) are output as missing.

  The functions use "\\n" as an end of line character,
  unless \code{.Platform$OS.type == "windows"}, when eol is
  "\\r\\n".  I only have access to linux machines for
  testing.

  I have tested these functions with my own data, but it is
  always possible that your data may contain quirks mine
  don't, or that input formats could change for any program
  mentioned here.  Please do have a quick check on a small
  subset of data (eg, as in the example below), that the
  output for your exact combination of options looks
  sensible and matches the specified input format.
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
  Chris Wallace
}
\references{
  David Clayton (2012). snpStats: SnpMatrix and XSnpMatrix
  classes and methods. R package version 1.6.0.
  http://www-gene.cimr.cam.ac.uk/clayton

  phase/fastPhase:
  \url{http://stephenslab.uchicago.edu/software.html}

  beagle:
  \url{http://faculty.washington.edu/browning/beagle/beagle.html}

  IMPUTE:
  \url{http://mathgen.stats.ox.ac.uk/impute/impute_v2.html}

  MACH: \url{http://www.sph.umich.edu/csg/abecasis/MACH}

  snphap:
  \url{https://www-gene.cimr.cam.ac.uk/staff/clayton/software/snphap.txt}
}
\keyword{manip}

