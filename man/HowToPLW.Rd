\name{HowToPLW}
\alias{HowToPLW}

\title{View HowToPLW}
\description{Finds the location of the vignette HowToPLW and optionally opens it.}

\usage{
HowToPLW(view=TRUE)
}

\arguments{
  \item{view}{logical, should the document be opened using the default PDF document reader?}
}

\value{Character string giving the file location.}

\details{
If the operating system is other than Windows, then the PDF viewer used is that given by \code{Sys.getenv("R_PDFVIEWER")}.
The PDF viewer can be changed using \code{Sys.putenv(R_PDFVIEWER=)}.

This function is used by drop-down Vignettes menu when the Rgui interface for Windows is used.
}

\seealso{
\code{\link[utils]{vignette}}, \code{\link[Biobase]{openPDF}}, \code{\link[Biobase]{openVignette}}, \code{\link[base]{Sys.getenv}}, \code{\link[base]{Sys.putenv}}
}

\author{Magnus Astrand}
\examples{
HowToPLW(view=FALSE)
}
\keyword{documentation}
