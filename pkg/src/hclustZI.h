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

#ifndef hclustZI_H
#define hclustZI_H

#include <iostream>
#include <fstream>
#include <vector>

#include "hclust.h"



class hclustZI : public hclust {

    public:
        bool _ZI;           // indicator of ZI model
        double *_alphak;    // proportion of zeros
        int *_n0k;          // number of zeros
        int *_n1k;          // number of non-zeros

        double *_alphak0;   // proportion of zeros
        int *_n0k0;         // number of zeros
        int *_n1k0;         // number of non-zeros

        hclustZI(int lengthdata, int nbseg, int nbclust, bool ZI);
        void CAH();
        void Init(double *Data, int *RuptVect);
        friend std::ostream & operator << (std::ostream &s, const hclustZI & EMi);
        double distance(int nA, int n0A, int n1A, double mA, double vA, double alphaA, int nB, int n0B, int n1B, double mB, double vB, double alphaB, bool ZI);
        double loglike(int n, int n0, int n1, double v, double alpha);
        ~hclustZI();
};
#endif


