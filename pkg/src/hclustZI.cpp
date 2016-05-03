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

#include "hclustZI.h"
#include <fstream>
#include <list>
#include <cmath>
#include <cstring>
#include <math.h>
#include <cstdlib>
#include <memory>
#include <vector>
#include <Rcpp.h>


#define  Dist(i,j)  _D[i-2][j-1]
#define  nk(i)      _nk[i-1]
#define  n0k(i)     _n0k[i-1]
#define  n1k(i)     _n1k[i-1]
#define  mk(i)      _mk[i-1]
#define  vk(i)      _vk[i-1]
#define  alphak(i)  _alphak[i-1]
#define  Dtmp(i)    _Dtmp[i-1]


using namespace std;
using namespace Rcpp;

const double pi = 3.141592653589793;


void hclustZI::CAH() {
    for (int k=_K;k>=_P+1;k--) {

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
        int    n0tmp = n0k(imin)+n0k(jmin);
        int    n1tmp = n1k(imin)+n1k(jmin);
        double mtmp = (n1k(imin)*mk(imin)+n1k(jmin)*mk(jmin))/n1tmp;
        double vtmp = ( n1k(imin)*vk(imin) + n1k(jmin)*vk(jmin) + n1k(imin)*(mk(imin)-mtmp)*(mk(imin)-mtmp)+n1k(jmin)*(mk(jmin)-mtmp)*(mk(jmin)-mtmp) ) / n1tmp;
        double alphatmp = (double) n0tmp / (double) ntmp;
        mk(jmin) = mtmp;
        vk(jmin) = vtmp;
        nk(jmin) = ntmp;
        n0k(jmin) = n0tmp;
        n1k(jmin) = n1tmp;
        alphak(jmin) = alphatmp;

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
            int i = index[h-1];
            Dtmp(i) = this->distance(nk(i), n0k(i), n1k(i), mk(i), vk(i), alphak(i), ntmp, n0tmp, n1tmp, mtmp, vtmp, alphatmp, _ZI);
        }
        double auxm = mk(imin); mk(imin) = mk(k);  mk(k) = auxm;
        double auxv = vk(imin); vk(imin) = vk(k);  vk(k) = auxm;
        int    auxn = nk(imin); nk(imin) = nk(k);  nk(k) = auxn;
        int    auxn0 = n0k(imin); n0k(imin) = n0k(k);  n0k(k) = auxn0;
        int    auxn1 = n1k(imin); n1k(imin) = n1k(k);  n1k(k) = auxn1;
        double auxalpha = alphak(imin); alphak(imin) = alphak(k);  alphak(k) = auxalpha;
        double auxk = Dtmp(k);
        Dtmp(k)     = Dtmp(imin);
        Dtmp(imin)  = auxk;

        int auxindex = _indiv_index[imin-1];
        _indiv_index[imin-1] = _indiv_index[k-1];
        _indiv_index[k-1] = auxindex;

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


void hclustZI::Init(double *Data, int *RuptVect) {
    memcpy(_x,Data,sizeof(double)*_lengthx);

    for (int k = 0 ;  k < _K+1; k++)
        _Breaks[k] = RuptVect[k];

    for (int k = 0 ;  k < _K; k++) {
        int start  = _Breaks[k];
        int end    = _Breaks[k+1]-1;
        _nk[k]     = end-start+1;
        double *xt = &_x[start];
        for (int t = start ; t<end+1 ; t++, xt++) {
            _mk[k]  += *xt;
            _vk[k]  += (*xt)*(*xt);
            if( (*xt) == 0) {
                _n0k[k]++;
            } else {
                _n1k[k]++;
            }
        }
        xt      = &_x[start];
        _mk[k] /= _n1k[k];
        _vk[k] /= _n1k[k];
        _vk[k] -= _mk[k]*_mk[k];
        _alphak[k] = (double) _n0k[k] / (double) _nk[k];
    }

    for (int k = 0 ;  k < _K; k++) {
        _mk0[k] = _mk[k];
        _vk0[k] = _vk[k];
        _nk0[k] = _nk[k];

        _alphak0[k] = _alphak[k];
        _n0k0[k] = _n0k[k];
        _n1k0[k] = _n1k[k];
    }

    for (int k = 0 ;  k < _K-1; k++) {
        for (int r = 0 ;  r < k+1; r++) {
            _D[k][r] = this->distance(_nk[k+1], _n0k[k+1], _n1k[k+1], _mk[k+1], _vk[k+1], _alphak[k+1], _nk[r], _n0k[r], _n1k[r], _mk[r], _vk[r], _alphak[r], _ZI);
        }
    }
} //end Init


hclustZI::hclustZI(int n, int nbsegments, int nbclusters, bool ZI) : hclust(n, nbsegments, nbclusters) {

    _ZI = ZI;
    _alphak = new double[_K];
    _n0k = new int[_K];
    _n1k = new int[_K];

    _alphak0 = new double[_K];
    _n0k0 = new int[_K];
    _n1k0 = new int[_K];

    //_indiv_index.resize(_N);

    for (int k =0; k<_K; k++) {
        _alphak[k] = 0;
        _n0k[k] = 0;
        _n1k[k] = 0;

        _alphak0[k] = 0;
        _n0k0[k] = 0;
        _n1k0[k] = 0;
    }

} // end constructor

// getter

// destructor

hclustZI::~hclustZI() {

    delete[] _alphak;
    delete[] _n0k;
    delete[] _n1k;

    delete[] _alphak0;
    delete[] _n0k0;
    delete[] _n1k0;

} // end destructor



// distance
double hclustZI::distance(int nA, int n0A, int n1A, double mA, double vA, double alphaA, int nB, int n0B, int n1B, double mB, double vB, double alphaB, bool ZI) {
    double ybar    = (n1A*mA+n1B*mB) / (n1A+n1B);
    double varpool =  (  n1B*vB + n1A*vA + n1B*(mB-ybar)*(mB-ybar) + n1A*(mA-ybar)*(mA-ybar) ) / (n1B+n1A);
    double alpha   = (double) (n0A + n0B) / (double) (nA + nB);
    double res = 0;

    if(ZI) {
        res = this->loglike(nA, n0A, n1A, vA, alphaA) + this->loglike(nB, n0B, n1B, vB, alphaB) - this->loglike(nA+nB, n0A+n0B, n1A+n1B, varpool, alpha);
    } else {
        res = hclust::distance(nA, mA, vA, nB, mB, vB);
    }
    if(abs(res) < 1E-16) {
        res = 0;
    }
    return(res);
}


double hclustZI::loglike(int n, int n0, int n1, double v, double alpha) {
    double res = 0;
    if(alpha == 0) {
        res = (-n1) * mylog(sqrt(v)) - n1 * log(sqrt(2*pi)) - ( (double) n1 / (double) 2);
    } else if(alpha == 1) {
        res = 0;
    } else {
        res = n1 * mylog(1-alpha) + n0 * mylog(alpha) - n1 * mylog(sqrt(v)) - n1 * log(sqrt(2*pi)) - ( (double) n1 / (double) 2);
    }
    return(res);
}

