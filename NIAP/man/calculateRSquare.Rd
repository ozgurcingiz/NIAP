\name{calculateRSquare}
\alias{calculateRSquare}

\title{Return R Squared value }
\usage{calculateRSquare(graph)}
\arguments{

\item{graph}{Igraph object refers gene coexpression network }

}
\value{
  \code{calculateRSquare} Return R-Squared value and goodness of fitness value
}
\description{
  \code{calculateRSquare} Return R Squared value. This function is utilized in calculateTopologicalFitness.R function.
}
\examples{

g = erdos.renyi.game(500, 350, type = "gnm")
rsquareValue= calculateRSquare(g)


}
\keyword{misc}

