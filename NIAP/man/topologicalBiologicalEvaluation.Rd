\name{topologicalBiologicalEvaluation}
\alias{topologicalBiologicalEvaluation}

\title{Assessment of gene networks via toplogocial and biological features}
\usage{topologicalBiologicalEvaluation(mypath="",omim)}
\arguments{

\item{mypath}{The path of folder which contain all gene coexpression networks whose extension is text file. It has two columns which refers gene- gene interactions.}

\item{omim}{OMIM disease gene list }

}
\value{
  \code{topologicalBiologicalEvaluation} return environment res that returns biological analysis results using GO- terms as res$biolEval and  topologicaal features as res$topological.
}
\description{
  \code{topologicalBiologicalEvaluation} Returns topological and biological analyses results of gene coexpression networks whose are at folder of given path.
}



\examples{

data(omim)

res <- topologicalBiologicalEvaluation(mypath="",omim)

topologicalResults= res$topological
biologicalAnalyisResults=res$biolEval

}
\keyword{misc}

