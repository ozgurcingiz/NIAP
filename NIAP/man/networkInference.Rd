\name{networkInference}
\alias{networkInference}

\title{Gene coexpression network inference }
\usage{networkInference(dataset,anno,gnimethod)}
\arguments{
\item{dataset}{mRNA,RNA-Seq or any other biological expression dataset file}

\item{anno}{Annotation data for getting all possible interaction space. }

\item{gnimethod}{We used information theory based algorithms. It string value must be any of these values which show popular inference algorthims: c3net, aracne, clr, mrnet or wgcna}

}
\value{
  \code{networkInference} Return gene- gene interactions which forms gene coexpression network
}
\description{
  \code{networkInference} Gene coexpression network inference step
}



\examples{

data(dataset)
data(annotation)

geneNetwork <- networkInference(dataset,annotation,gnimethod="c3net")

}
\keyword{misc}

