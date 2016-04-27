### Copyright 2016-04 Ghislain DURIF
###
### This file is part of the `hClustering' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


#' @title Hierarchical Clustering based on log-likelihood ratio
#'
#' @description
#' \code{hClust_llr} returns an object of type \code{hclust},
#' results of the clustering on individuals, using the dissimilarity
#' based on log-likelihood ratio test.
#'
#' @details
#' Algorithm to process hierarchical clustering on individuals
#' (from a matrix individuals x variables). The dissimilarity used
#' in the observation space is the log-likelihood ratio, i.e.
#' D(i,j) = log f(Y_i) + log f(Y_j) - log f(Y_i U Y_j)
#' where Y_i is the vector of observations from individual i.
#' The statistical model defining the considered likelihood is
#' standard Gaussian or Zero-Inflated Gaussian (i.e. a Gaussian-Bernoulli
#' mixture). The agglomeration method is also the log-likelihood ratio,
#' however between clusters in this case.
#'
#' Wrapper for Cpp function
#'
#'@author
#'Ghislain Durif, \email{ghislain.durif@univ-lyon1.fr}
#'Franck Picard, \email{franck.picard@univ-lyon1.fr}
#'
#'@seealso \code{\link{plot_heatmap}}
#'
#'@import Rcpp
#'@import RcppEigen
#'@useDynLib ziLRhClust
#'
#'@param X matrix n x p with observations (i.e. individuals) in rows and variables in columns
#'@param ZI boolean, indicating whether using zero-inflated Gaussian model (if TRUE) or
#'just gaussian model (if FALSE), default is TRUE
#'
#'@return an object of type hclust, see \code{\link[stats]{hclust}}
#'\item{diss}{the n x n dissimilarity matrix between individuals}
#'\item{}{(not included in \code{hclust} object)}
#'\item{height}{\code{\link[stats]{hclust}}}
#'\item{merge}{\code{\link[stats]{hclust}}}
#'\item{order}{\code{\link[stats]{hclust}}}
#'\item{labels}{vector of n individual labels}

#### R wrapper to return object understandable by dendrogram function

#' @export
hClust_llr <- function(X, ZI=TRUE) {

    X = as.matrix(X)
    X = apply(X, c(1,2), as.double)
    n = nrow(X)
    p = ncol(X)

    Xvect = as.vector(t(X))

    res = Cpp_hclust(Xvect, n, p, ZI)

    ### return
    rlabels = NULL
    if(!is.null(row.names(X))) {
        rlabels = row.names(X)
    } else {
        rlabels = paste0(1:n)
    }

    row.names(res$diss) = rlabels
    colnames(res$diss) = rlabels

    res$labels=rlabels

    class(res) = "hclust"

    return(res)

}