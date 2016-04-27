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

#include "hclust.h"
#include <fstream>
#include <list>
#include <cmath>
#include <cstring>
#include <math.h>
#include <cstdlib>
#include <memory>
#include <vector>
#include <Rcpp.h>


#define  Dist(i,j) _D[i-2][j-1]
#define  nk(i)     _nk[i-1]
#define  mk(i)     _mk[i-1]
#define  vk(i)     _vk[i-1]
#define  Dtmp(i)   _Dtmp[i-1]


using namespace std;
using namespace Rcpp;


void hclust::CAH() {
    for (int k=_K;k>=_P+1;k--) {

        //Rcout << "##########" << std::endl;
        //Rcout << "k = " << k << std::endl;

        int imin    = 2;
        int jmin    = 1;
        double dmin = Dist(2,1);
        for (int i=2; i<k+1 ; i++) {
            for (int j=1; j<i ; j++) {
                if (Dist(i,j)<=dmin) {
                    dmin = Dist(i,j);
                    imin = i;
                    jmin = j;
                }
            }
        }

        //Rcout << "imin = " << imin << std::endl;
        //Rcout << "jmin = " << jmin << std::endl;

        int round = _K-k;

        _height[round] = dmin;

        // true index of individuals at imin and jmin
        int index_imin = _indiv_index[imin-1]-1;
        int index_jmin = _indiv_index[jmin-1]-1;

        //Rcout << "true imin = " << index_imin << std::endl;
        //Rcout << "true jmin = " << index_jmin << std::endl;

        // if individuals not chosen yet, give the true index, else give the last change
        int merge0 = 0;
        int merge1 = 0;
        if(_indiv_status[index_jmin] == -1) {
            merge0 = -(index_jmin+1);
        } else {
            merge0 = _indiv_status[index_jmin];
        }
        if(_indiv_status[index_imin] == -1) {
            merge1 = -(index_imin+1);
        } else {
            merge1 = _indiv_status[index_imin];
        }

        // order to follow hclust order
        if(merge0 < 0 && merge1 < 0) {
            _merge[round][0] = std::max(merge0, merge1);
            _merge[round][1] = std::min(merge0, merge1);
        } else {
            _merge[round][0] = std::min(merge0, merge1);
            _merge[round][1] = std::max(merge0, merge1);
        }


        // change status of all individuals in the cluster chosen (if not singleton)
        _indiv_status[index_imin] = round+1;
        _indiv_status[index_jmin] = round+1;



        int    ntmp = nk(imin)+nk(jmin);
        double mtmp = (nk(imin)*mk(imin)+nk(jmin)*mk(jmin))/ntmp;
        double vtmp = ( nk(imin)*vk(imin) + nk(jmin)*vk(jmin) + nk(imin)*(mk(imin)-mtmp)*(mk(imin)-mtmp)+nk(jmin)*(mk(jmin)-mtmp)*(mk(jmin)-mtmp) ) / ntmp;
        mk(jmin) = mtmp;
        vk(jmin) = vtmp;
        nk(jmin) = ntmp;

        int i1 = 0;
        int i2 = 0;
        int i3 = 0;

        vector<int> index;

        if (jmin>1) { // 1:(jmin-1)
            i1     = jmin-1;
            for (int h = 1; h<i1+1 ;h++) {
                index.push_back(h);
            }
        } else {
            i1 = 1;
            index.push_back(1);
        }
        if (jmin+1<=imin-1) { // (jmin+1):(imin-1)
            i2     = imin-1-jmin-1+1;
            for (int h=1;h<(i2+1);h++) {
                index.push_back(h+jmin);
            }
        }
        if (imin+1<=(k)) {// (imin+1):k
            i3     = k-imin-1+1;
            for (int h=1;h<(i3+1);h++) {
                index.push_back(h+imin);
            }
        }

        int size_index= index.size();

        for (int h=1;h<size_index+1;h++) {
            int i          = index[h-1];
            /**
            double ybar    =  (nk(i)*mk(i)+ntmp*mtmp)/(nk(i)+ntmp);
            double varpool =  ( ntmp*vtmp + nk(i)*vk(i) + ntmp*(mtmp-ybar)*(mtmp-ybar) + nk(i)*(mk(i)-ybar)*(mk(i)-ybar) ) / (ntmp+nk(i));
            Dtmp(i)        = (nk(i)+ntmp)*log(sqrt(varpool)) -nk(i)*log(sqrt(vk(i)))-ntmp*log(sqrt(vtmp));
             **/
            Dtmp(i) = this->distance(nk(i), mk(i), vk(i), ntmp, mtmp, vtmp);

        }
        double auxm = mk(imin); mk(imin) = mk(k);  mk(k) = auxm;
        double auxv = vk(imin); vk(imin) = vk(k);  vk(k) = auxm;
        int    auxn = nk(imin); nk(imin) = nk(k);  nk(k) = auxn;
        double auxk = Dtmp(k);
        Dtmp(k)     = Dtmp(imin);
        Dtmp(imin)  = auxk;

        int auxindex = _indiv_index[imin-1];
        _indiv_index[imin-1] = _indiv_index[k-1];
        _indiv_index[k-1] = auxindex;

        for(int i=0; i<_K; i++) {
            //Rcout << _indiv_index[i] << " ";
        }

        //Rcout << std::endl;
        for (int i=1; i<(jmin-1+1);i++) {
            Dist(jmin,i) = Dtmp(i);
        }
        for (int i=jmin+1;i<(k-1+1);i++) {
            Dist(i,jmin) = Dtmp(i);
        }
        for (int j=1;j<(jmin-1+1);j++) {
            Dist(imin,j) = Dist(k,j);
        }
        for (int j=jmin+1;j<(imin-1+1);j++) {
            Dist(imin,j) = Dist(k,j);
        }
        for (int i=imin+1;i<(k-1+1);i++) {
            Dist(i,imin) = Dist(k,i);
        }

    } // end k
} //end CAH


void hclust::Init(double *Data, int *RuptVect) {

    memcpy(_x,Data,sizeof(double)*_lengthx);

    for (int k = 0 ;  k < _K+1; k++)
        _Breaks[k] = RuptVect[k];

    for (int k = 0 ;  k < _K; k++) {
        int start  = _Breaks[k];
        int end    = _Breaks[k+1]-1;
        _nk[k]     = end-start+1  ;
        double *xt = &_x[start];
        for (int t = start ; t<end+1 ; t++, xt++) {
            _mk[k]  += *xt;
            _vk[k]  += (*xt)*(*xt);
        }
        xt      = &_x[start];
        _mk[k] /= end-start+1 ;
        _vk[k] /= end-start+1 ;
        _vk[k] -= _mk[k]*_mk[k];
    }

    for (int k = 0 ;  k < _K; k++) {
        _mk0[k] = _mk[k];
        _vk0[k] = _vk[k];
        _nk0[k] = _nk[k];
    }

    for (int k = 0 ;  k < _K-1; k++) {
        for (int r = 0 ;  r < k+1; r++) {
            _D[k][r] = this->distance(_nk[k+1], _mk[k+1], _vk[k+1], _nk[r], _mk[r], _vk[r]);
        }
    }
    for (int k = 0 ;  k < _K-1; k++) {
        Rcout << _D[k][5] << " ";
    }

} //end Init


hclust::hclust(int n, int nbsegments, int nbclusters) : hclustOriginal(n, nbsegments, nbclusters) {

    _N = _lengthx/_K;

    _height = new double[(_K-1)];
    _merge = new double *[(_K-1)];
    _indiv_status = new int[_K];

    //_indiv_index.resize(_N);

    for (int k =0; k<_K-1; k++) {
        _height[k] = 0;
        _merge[k] = new double[2];
        _merge[k][0] = 0;
        _merge[k][1] = 0;
    }

    for (int i=0; i<_K; i++) {
        _indiv_index.push_back(i+1);
        _indiv_status[i] = -1;
    }

} // end constructor

// getter
double** hclust::getDist() {
    return _D;
}

int hclust::getDim() {
    return _K;
}

double** hclust::getMerge() {
    return _merge;
}

double* hclust::getHeight() {
    return _height;
}

// destructor

hclust::~hclust() {

    delete[] _height;
    delete[] _indiv_status;

    for (int k =0; k<_K-1; k++) {
        delete[] _merge[k];
    }

    delete[] _merge;

} // end destructor



// distance
double hclust::distance(int nA, double mA, double vA, int nB, double mB, double vB) {
    double ybar    = (nA*mA+nB*mB) / (nA+nB);
    double varpool =  (  nB*vB + nA*vA + nB*(mB-ybar)*(mB-ybar) + nA*(mA-ybar)*(mA-ybar) ) / (nB+nA);
    double res = (nA+nB) * log(sqrt(varpool)) - nA*log(sqrt(vA)) - nB*log(sqrt(vB));
    if(abs(res) < 1E-16) {
        res = 0;
    }
    return(res);
}
