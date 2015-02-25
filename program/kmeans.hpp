#ifndef _UBLAS_
#define _UBLAS_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#endif

#ifndef _KMEANS_
#define _KMEANS_

using namespace boost::numeric::ublas;
using namespace boost::numeric;

double dist_xc /* データ点xとクラスタ中心c_jとの距離 */
(const ublas::matrix_row<const ublas::matrix<double> > x,
 const ublas::matrix<double> c,
 const unsigned int j)
{
  double sum = 0;
  for(unsigned int d = 0; d < x.size(); ++d){
    sum += (x[d] - c(j,d)) * (x[d] - c(j,d));    
  }
  return sqrt(sum);
}

double dist_cc /* クラスタ中心間の距離を計算 */
(const ublas::matrix<double> mat1, const unsigned int r1,
 const ublas::matrix<double> mat2, const unsigned int r2) 
{
  if(mat1.size2() != mat2.size2()){
    return -1;
  }else{
    double sum = 0;
    for(unsigned int i = 0; i < mat1.size2(); ++i){
      sum += (mat1(r1, i) - mat2(r2, i)) * (mat1(r1, i) - mat2(r2, i));
    }
    return sqrt(sum);
  }
}

#endif
