#include <fstream>
#include <list>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <memory>
#include <vector>


#include <Rcpp.h>
#include <RcppEigen.h>


using namespace std;
using namespace Rcpp;

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;