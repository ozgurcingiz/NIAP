\name{calculateTopologicalFitness}
\alias{calculateTopologicalFitness}

\title{This function measures how well gene networks follows scale free network topology}
\usage{calculateTopologicalFitness(mypath="")}
\arguments{

\item{mypath}{The path of folder which contain all gene coexpression networks whose extension is text file. It has two columns which refers gene- gene interactions.}

}
\value{
  \code{calculateTopologicalFitness} return  Kolmogorov Smirnov Goodness of fitness and R squared values.
}
\description{
  \code{calculateTopologicalFitness} measures how well gene networks follows scale free network topology.
}


\examples{


topologyFitnessValues <- calculateTopologicalFitness(mypath="")

}
\keyword{misc}

