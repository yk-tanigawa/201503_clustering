#include <algorithm>
#include "my_ublas.hpp"
#include "kmeans.hpp"

#ifndef _ELKAN_
#define _ELKAN_

using namespace boost::numeric::ublas;
using namespace boost::numeric;

int Elkan_initialize
(const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::matrix<double> &c,
 ublas::vector<double> &u,
 ublas::matrix<double> &l,
 ublas::vector<bool> &r,
 ublas::vector<unsigned int> &q,
 ublas::matrix<double> &c_sum,
 ublas::matrix<double> &c_prev,
 ublas::vector<double> &temp)
{
  /* l, q, c_sum はすでに0で初期化されているので何もしない */
  
  /* 各データ点の初期クラスタとu[]を計算 */
  for(unsigned int i = 0; i < a.size(); ++i){
    for(unsigned int j = 0; j < k; ++j){
      temp[j] = l(i, j) = dist_xc(row(x, i), c, j);
    }
    a[i] = argmin(temp); u[i] = temp[a[i]];
    row(c_sum,a[i]) += row(x,i); q[a[i]]++;
  }

  c_prev = c;   /* クラスタ中心の更新 -- O(dk) time  */
  for(unsigned int j = 0; j < k; ++j){
    row(c, j) = row(c_sum, j) / q[j];
  }

  for(unsigned int i = 0; i < a.size(); ++i){ /* u[i]の更新 -- O(nd) time */
    u[i] += dist_cc(c, a[i], c_prev, a[i]);
    r[i] = true;
  }
  
  return 0;
}


int Elkan_repeat
(const unsigned int n,
 const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::matrix<double> &c,
 ublas::vector<double> &u,
 ublas::matrix<double> &l,
 ublas::vector<bool> &r,
 ublas::vector<unsigned int> &q,
 ublas::matrix<double> &c_sum,
 ublas::vector<double> &s,
 ublas::matrix<double> &c_prev,
 ublas::vector<double> &temp)
{
  int count = 0; /* 最近接クラスタが変化したデータ点の数*/
   
  /*  s[]の更新 -- O(dk^2) time
   *  s[j] = min_{jj != j} dist(c[jj], c[j])
   *      = min2_jj dist(c[j], c[jj]) 
   * Elkanの論文では1/2倍しているがHamerlyのアルゴリズムにあわせて
   * 最近接クラスタの距離とする */
  for(unsigned int j = 0; j < k; ++j){
    for(unsigned int jj = 0; jj < k; ++jj){
      temp[jj] = dist_cc(c, j, c, jj);
    }
    s[j] = min2(temp); 
  }

  /* 各点x[i]の最近接クラスタ割り当て -- O(ndk) time
   * 途中，上界と下界を用いて枝刈りを実行 */
  for(unsigned int i = 0; i < a.size(); ++i){
    if(u[i] > (s[a[i]] / 2.0)){
      for(unsigned int j = 0; j < k; ++j){
	if((j != a[i]) &&
	   (u[i] > l(i, j)) &&
	   (2 * u[i] > dist_cc(c, a[i], c, j))){
	  /* x[i]とc[a[i]]との距離を求める */
	  double dist_to_closest = 0; 
	  if(r[i]){ /* u[x] is out of date */
	    dist_to_closest = l(i, a[i]) = dist_xc(row(x, i), c, a[i]);
	    r[i] = false;
	  }else{    /* u[x] is still the tight bound */
	    dist_to_closest = u[i]; /* そのまま代入できる */
	  }
	  
	  if((dist_to_closest > l(i, j)) ||
	     (2 * dist_to_closest > dist_cc(c, a[i], c, j))){
	    /* 実際にx[i]とc[j]の距離を計算して最近接クラスタが変化したか調べる */
	    if((l(i, j) = dist_xc(row(x, i), c, j)) < dist_to_closest){
	      /* x[i] の最近接クラスタは a[i]からjに変化した */
	      row(c_sum,a[i]) -= row(x,i); q[a[i]]--;
	      row(c_sum,j) += row(x,i);    q[j]++;
	      a[i] = j;  count++;
	    }
	  }
	}
      }
    }
  }


  c_prev = c;   /* クラスタ中心の更新 -- O(dk) time  */
  for(unsigned int j = 0; j < k; ++j){
    row(c, j) = row(c_sum, j) / q[j];
  }

  for(unsigned int i = 0; i < n; ++i){ /* l[i, j]の更新  -- O(ndk) time */
    for(unsigned int j = 0; j < k; ++j){
      l(i, j) = std::max((l(i, j) - dist_cc(c_prev, j, c, j)), 0.0);
    }
  }

  for(unsigned int i = 0; i < n; ++i){ /* u[i]の更新 -- O(nd) time */
    u[i] += dist_cc(c, a[i], c_prev, a[i]);
    r[i] = true;
  }
  
  return count; /* 最近接クラスタが変化したデータ点の数 */
}

int Elkan_main
(const unsigned int n,
 const unsigned int d,
 const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::matrix<double> &c)
{

  /* まず，メモリを確保 */  
  /* i番目のデータ点を主語とする記号 */
  ublas::vector<double> u(n, 0);        /* upper bound */
  ublas::matrix<double> l(n, k, 0);     /* lower bound */   
  ublas::vector<bool> r(n, true);
  
  /* j番目のクラスタを主語とする記号  */
  ublas::vector<unsigned int> q(k,0);  /* クラスタ中の点の数 */
  ublas::matrix<double> c_sum(k,d,0);  /* クラスタ中のベクトルの総和 */
  ublas::vector<double> s(k,0);        /* 最近接クラスタの重心との距離 */

  /* 下記は，イテレーションステップで用いる変数のメモリ */
  ublas::matrix<double> c_prev(k,d,0); /* 前のステップでのクラスタ中心 */
  ublas::vector<double> temp(k,0);     /* 計算過程で用いる */

  /* 初期ステップ */
  Elkan_initialize(k, x, a, c, u, l, r, q, c_sum, c_prev, temp);

  /* 繰り返しステップ */
  int t = 1; 
  while(Elkan_repeat(n, k, x, a, c, u, l, r, q, c_sum, s, c_prev, temp) > 0){
    t++; /* 繰り返しステップの実行回数を数える */
  }

  return t; /* 繰り返しステップの実行回数を返す */
}

#endif
