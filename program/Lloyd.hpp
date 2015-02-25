#include "my_ublas.hpp"
#include "kmeans.hpp"

#ifndef _LLOYD_
#define _LLOYD_

using namespace boost::numeric::ublas;
using namespace boost::numeric;

int Lloyd_repeat
(const unsigned int n,
 const unsigned int d,
 const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::vector<unsigned int> &a_prev,
 ublas::matrix<double> &c,
 ublas::vector<double> &dist)
{
  ublas::matrix<double> c_sum(k, d, 0);
  ublas::vector<unsigned int> q(k, 0);
  a_prev = a;
  
  for(unsigned int i = 0; i < n; ++i){
    for(unsigned int j = 0; j < k; ++j){
      dist[j] = dist_xc(row(x, i), c, j);
    }
    a[i] = argmin(dist);
    row(c_sum, a[i]) += row(x, i);
    q[a[i]]++;
  }

  for(unsigned int j = 0; j < k; ++j){
    row(c, j) = row(c_sum, j) / q[j];
  }

  int diff = 0; /* クラスタ割り当てが変化した点の数 */
  for(unsigned int i = 0; i < n; ++i){
    if(a[i] != a_prev[i]){ diff++; }
  }
  
  return diff;
}


int Lloyd_main
(const unsigned int n,
 const unsigned int d,
 const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::matrix<double> &c)
{
  /* 作業用メモリを確保 */  
  ublas::vector<double> temp(k,0);
  ublas::vector<unsigned int> a_prev(n, 0);
  
  /* 初期クラスタ中心はすでに得ているので何もしない */
  
  /* 繰り返しステップ */
  int t = 1;
  while(Lloyd_repeat(n, d, k, x, a, a_prev, c, temp) > 0){
    t++; /* 繰り返しステップの実行回数を数える */
  }  
  return t; /* 繰り返しステップの実行回数を返す */
}

#endif
