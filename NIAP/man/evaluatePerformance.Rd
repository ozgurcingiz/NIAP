\name{literatureOverlapAnalysis}
\alias{literatureOverlapAnalysis}

\title{Validating the inferred network using literature data<PPI data>}
\usage{literatureOverlapAnalysis(ppi,anno,mypath="")}
\arguments{
\item{ppi}{Protein- protein interaction data , literature data, which is utilized as validation data in overlap analysis  }

\item{anno}{Annotation data for getting all possible interaction space. }

\item{mypath}{The path of folder which contain all gene coexpression networks whose extension is text file. It has two columns which refers gene- gene interactions.}

}
\value{
  \code{literatureOverlapAnalysis} Return performance evaluation metrics of gene coexpression networks.
}
\description{
  \code{literatureOverlapAnalysis} Overlap analysis of gene networks via literature data.
}



\examples{

data(ppi)
data(annotation)

overlapAnalysisResult <- literatureOverlapAnalysis(ppi,annotation,mypath="")

}
\keyword{misc}

