/// Copyright 2016-04 Ghislain Durif
///
/// This file is part of the `hClustering' library for R and related languages.
/// It is made available under the terms of the GNU General Public
/// License, version 2, or at your option, any later version,
/// incorporated herein by reference.
///
/// This program is distributed in the hope that it will be
/// useful, but WITHOUT ANY WARRANTY; without even the implied
/// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
/// PURPOSE.  See the GNU General Public License for more
/// details.
///
/// You should have received a copy of the GNU General Public
/// License along with this program; if not, write to the Free
/// Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
/// MA 02111-1307, USA

// functions to get hclust output (and eventually reorder it)

#ifndef order_H
#define order_H

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;

// produce vector of individual order (in hclust ouput)
void order(const MatrixXi &merge, VectorXi &iorder, int N);

// reorder merge and iorder for better heatmap aesthetic
// c.f. reorder.hclust from gclus package
void reOrder(const MatrixXd &dist, MatrixXi &merge, VectorXi &iorder, int N);

#endif