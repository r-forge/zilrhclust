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

// wrapper to construct hclust class object

#include <Rcpp.h>
#include <RcppEigen.h>
#include "hclustZI.h"


using namespace Rcpp;


// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;

// INPUT
//  X: data vector (of length n)
//  K: nb of groups (n if no group and clustering of individuals)

#define mergeA(i) merge(i-1,0)
#define mergeB(i) merge(i-1,1)
#define _iorder(i) iorder(i-1)

// [[Rcpp::export]]
Rcpp::List Cpp_hclust(Map<VectorXd> X, int n, int p, bool ZI) {

    int nbCluster=1; //will be one in the future

    //Rcout << X << std::endl;

    // data for test
    double *data;
    data = X.data();

    int rupture[n];
    for(int i=0; i<n+1; i++) {
        rupture[i] = i*p;
    }

    // construction
    hclustZI myClust(n*p, n, nbCluster, ZI);

    // init
    myClust.Init(data, rupture);


    // distance between individuals
    double **arrayDist = myClust.getDist();
    MatrixXd indivDist( myClust.getDim(), myClust.getDim() );
    indivDist.setZero();

    for(int i=0; i<myClust.getDim()-1; i++) {
        for(int j=0; j<i+1; j++) {
            indivDist(i+1,j) = arrayDist[i][j];
            indivDist(j,i+1) = arrayDist[i][j];
        }
    }

    // clustering
    myClust.CAH();


    // results of clustering
    double **arrayMerge = myClust.getMerge();
    double *arrayHeight = myClust.getHeight();

    VectorXd height = Map<VectorXd>( arrayHeight, myClust.getDim() - 1);

    MatrixXd merge( myClust.getDim()-1, 2 );
    merge.setZero();

    for(int i=0; i<myClust.getDim()-1; i++) {
        merge(i,0) = arrayMerge[i][0];
        merge(i,1) = arrayMerge[i][1];
    }


    //----- order -----//
    VectorXi iorder( myClust.getDim() );
    iorder.setZero();


    // cf hclust.f in stats pkg
    int N = myClust.getDim();
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



    // return
    return Rcpp::List::create(Rcpp::Named("diss") = indivDist,
                              Rcpp::Named("height") = height,
                              Rcpp::Named("merge") = merge,
                              Rcpp::Named("order") = iorder);

}