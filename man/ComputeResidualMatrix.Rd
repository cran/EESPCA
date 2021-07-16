\name{computeResidualMatrix}
\alias{computeResidualMatrix}
\title{Calculates the residual matrix from the reduced rank reconstruction}
\description{
Utility function for computing the residual matrix formed by 
subtracting from \code{X} a reduced rank approximation of matrix \code{X} generated from 
the top k principal components contained in matrix \code{V}.
}
\usage{
    computeResidualMatrix(X,V)
}
\arguments{
\item{X}{An n-by-p data matrix whose top k principal components are contained in the p-by-k matrix \code{V}.}
\item{V}{A p-by-k matrix containing the loadings for the top k principal components of \code{X}.}
}
\value{
Residual matrix.
}
\examples{
    set.seed(1)
    # Simulate 10x5 MVN data matrix
    X=matrix(rnorm(50), nrow=10)
    # Center columns
    X = scale(X, center=TRUE, scale=FALSE)
    # Perform PCA
    prcomp.out = prcomp(X)
    # Get rank 2 residual matrix
    computeResidualMatrix(X=X, V=prcomp.out$rotation[,1:2])
}