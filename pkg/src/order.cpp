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

#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include "order.h"

#define mergeA(i) merge(i-1,0)
#define mergeB(i) merge(i-1,1)
#define _iorder(i) iorder(i-1)
#define _merges(i,j) merges(i-1,j-1)
#define _endpoints(i,j) endpoints(i-1,j-1)
#define _dist(i,j) dist(i-1,j-1)
#define _dir(i,j) dir(i-1,j-1)
#define _clusters(i) clusters[i-1]

using std::vector;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;

// produce vector of individual order (in hclust ouput)
// c.f. hclust.f in stats pkg
void order(const MatrixXi &merge, VectorXi &iorder, int N) {

    _iorder(1) = mergeA(N-1);
    _iorder(2) = mergeB(N-1);
    int loc = 2;
    for(int i=N-2; i>0; i--) {
        for(int j=1; j<loc+1; j++) {
            if(_iorder(j) == i) {
                _iorder(j) = mergeA(i);
                if(j == loc) {
                    loc=loc+1;
                    _iorder(loc) = mergeB(i);
                } else {
                    loc=loc+1;
                    for(int k=loc; k>j+1; k--) {
                        _iorder(k) = _iorder(k-1);
                    }
                    _iorder(j+1) = mergeB(i);
                }
            }
        }
    }

    iorder.array() = (-1) * iorder.array();
}

// reorder merge and iorder for better heatmap aesthetic
// c.f. reorder.hclust from gclus package
void reOrder(const MatrixXd &dist, MatrixXi &merge, VectorXi &iorder, int N) {

    int n = N-1;

    MatrixXi endpoints = MatrixXi::Zero(n,2);
    MatrixXi dir = MatrixXi::Ones(n,2);

    MatrixXi merges(merge);

    for(int i=1; i<n+1; i++) {
        int j = _merges(i,1);
        int k = _merges(i,2);
        if((j < 0) && (k<0)) {
            _endpoints(i,1) = -j;
            _endpoints(i,2) = -k;
        }
        else if(j < 0) {
            j = -j;
            _endpoints(i,1) = j;
            if(_dist(j, _endpoints(k,1)) < _dist(j, _endpoints(k,2))) {
                _endpoints(i,2) = _endpoints(k,2);
            } else {
                _endpoints(i,2) = _endpoints(k,1);
                _dir(i,2) = -1;
            }
        }
        else if(k <0) {
            k = -k;
            _endpoints(i,2) = k;
            if(_dist(k, _endpoints(j,1)) < _dist(k, _endpoints(j,2))) {
                _endpoints(i,1) = _endpoints(j,2);
                _dir(i,1) = -1;
            } else {
                _endpoints(i,1) = _endpoints(j,1);
            }
        } else {
            double d11 = _dist(_endpoints(j,1), _endpoints(k,1));
            double d12 = _dist(_endpoints(j,1), _endpoints(k,2));
            double d21 = _dist(_endpoints(j,2), _endpoints(k,1));
            double d22 = _dist(_endpoints(j,2), _endpoints(k,2));

            double tmp[] = {d11, d12, d21, d22};
            double dmin = *std::min_element(tmp,tmp+4);

            if(dmin == d21) {
                _endpoints(i,1) = _endpoints(j,1);
                _endpoints(i,2) = _endpoints(k,2);
            } else if(dmin == d11) {
                _endpoints(i,1) = _endpoints(j,2);
                _endpoints(i,2) = _endpoints(k,2);
                _dir(i,1) = -1;
            } else if(dmin == d12) {
                _endpoints(i,1) = _endpoints(j,2);
                _endpoints(i,2) = _endpoints(k,1);
                _dir(i,1) = -1;
                _dir(i,2) = -1;
            } else {
                _endpoints(i,1) = _endpoints(j,1);
                _endpoints(i,2) = _endpoints(k,1);
                _dir(i,2) = -1;
            }
        }
    }

    for(int i=n; i>1; i--) {
        if(_dir(i,1) == -1) {
            int m = _merges(i,1);
            if(m > 0) {
                int m1 = _merges(m,1);
                _merges(m,1) = _merges(m,2);
                _merges(m,2) = m1;
                if(_dir(m,1) == _dir(m,2)) {
                    _dir(m,1) = -(_dir(m,1));
                    _dir(m,2) = -(_dir(m,2));
                }
            }
        }
        if(_dir(i,2) == -1) {
            int m = _merges(i,2);
            if(m > 0) {
                int m1 = _merges(m,1);
                _merges(m,1) = _merges(m,2);
                _merges(m,2) = m1;
                if(_dir(m,1) == _dir(m,2)) {
                    _dir(m,1) = -(_dir(m,1));
                    _dir(m,2) = -(_dir(m,2));
                }
            }
        }
    }

//     vector<vector<int> > clusters(n);
//     for(int i=1; i<n+1; i++) {
//         _clusters(i).push_back(i);
//     }

//     for(int i=1; i<n+1; i++) {
//         int j = _merges(i,1);
//         int k = _merges(i,2);
//         if((j<0) && (k<0)) {
//             _clusters(i).clear();
//             _clusters(i).push_back(-j);
//             _clusters(i).push_back(-k);
//         } else if(j < 0) {
//             _clusters(i) = _clusters(k);
//             _clusters(i).insert(_clusters(i).begin(), -j);
//         } else if(k < 0) {
//             _clusters(i) = _clusters(j);
//             _clusters(i).push_back(-k);
//         } else {
//             _clusters(i) = _clusters(j);
//             _clusters(i).insert(_clusters(i).end(), _clusters(k).begin(), _clusters(k).end());
//         }
//     }

    merge = merges;
    //VectorXi tmp(_clusters(n).data(), n);
    //Map<VectorXi> tmp(_clusters(n).data(), n);
    order(merge, iorder, N);



}
