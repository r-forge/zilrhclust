/// Copyright 2016-04 Franck Picard
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

#ifndef hclustOriginal_H
#define hclustOriginal_H

#include <iostream>
#include <fstream>


//namespace EMinit{
class hclustOriginal {

    public:
        int       _lengthx;        // size of the data
        int       _K;              // Number of segments
        int       _P;              // Number of clusters
        double   *_phi;            // parameters
        double   *_x;              // data
        double   *_mk;             // empirical means
        double   *_vk;             // empirical var
        int      *_nk;             // length of segments
        double   *_mk0;             // empirical means
        double   *_vk0;             // empirical var
        int      *_nk0;             // length of segments
        int      *_Breaks ;        // breakpoints
        double   **_D;             // distance matrix
        double   *_Dtmp;           // distance vector of the merged groups
        double   **_tau;           // posterior

        hclustOriginal(int lengthdata, int nbseg, int nbclust);
        void CAH();
        void compute_phi();
        void Init(double *Data, int *RuptVect);
        friend std::ostream & operator << (std::ostream &s, const hclustOriginal & EMi);
        double distance(double nA, double mA, double vA, double nB, double mB, double vB);
        ~hclustOriginal();
};
#endif
