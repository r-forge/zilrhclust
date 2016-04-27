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

#ifndef hclust_H
#define hclust_H

#include <iostream>
#include <fstream>
#include <vector>

#include "hclustOriginal.h"



class hclust : public hclustOriginal {

    public:
        int _N;                     // number of observations (_lengthx/_K)
        double *_height;            // height of merging
        double **_merge;            // merge array as in hclust from R
        std::vector<int> _indiv_index;   // individual index at pos i in dist matrix
        int *_indiv_status;         // individual status (not yet chosen (-1) or last time being chosen)

        hclust(int lengthdata, int nbseg, int nbclust);
        void CAH();
        void Init(double *Data, int *RuptVect);
        friend std::ostream & operator << (std::ostream &s, const hclust & EMi);
        double** getDist();
        int getDim();
        double** getMerge();
        double* getHeight();
        double distance(int nA, double mA, double vA, int nB, double mB, double vB);
        ~hclust();
};
#endif


