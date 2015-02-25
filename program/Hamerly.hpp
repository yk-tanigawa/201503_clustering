#include <algorithm>
#include "my_ublas.hpp"
#include "kmeans.hpp"

#ifndef _HAMERLY_
#define _HAMERLY_

using namespace boost::numeric::ublas;
using namespace boost::numeric;

int Hamerly_point_all_ctrs
(const unsigned int i,
 const unsigned int k,
 const ublas::matrix_row<const ublas::matrix<double> > xi,
 const ublas::matrix<double> c,
 ublas::vector<unsigned int> &a,
 ublas::vector<double> &u,
 ublas::vector<double> &l,
 ublas::vector<double> &temp)
{
  for(unsigned int j = 0; j < k; ++j){
    temp[j] = dist_xc(xi, c, j);
    /* ここは内積を使って書き換え可能 */
  }
  a[i] = argmin(temp);
  u[i] = dist_xc(xi, c, a[i]);
  l[i] = dist_xc(xi, c, argmin2(temp));
  return 0;
}

int Hamerly_move_centers
(const ublas::matrix<double> c_sum,
 const ublas::vector<unsigned int> q,
 ublas::matrix<double> &c,
 ublas::vector<double> &p,
 ublas::matrix<double> &c_prev)
{
  c_prev = c;
  for(unsigned int j = 0; j < q.size(); ++j){
    row(c, j) = row(c_sum, j) / q[j];
    p[j] = dist_cc(c_prev, j, c, j);
  }
  return 0;
}

int Hamerly_update_bounds
(const ublas::vector<double> p,
 const ublas::vector<unsigned int> a,
 ublas::vector<double> &u,
 ublas::vector<double> &l)
{
  for(unsigned int i = 0; i < u.size(); ++i){
    u[i] += p[a[i]];
    l[i] -= p[argmax(p)];
  }
  return 0;
}

int Hamerly_initialize
(const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::matrix<double> &c,
 ublas::vector<double> &u,
 ublas::vector<double> &l,
 ublas::vector<unsigned int> &q,
 ublas::matrix<double> &c_sum,
 ublas::matrix<double> &c_prev,
 ublas::vector<double> &p)
{
  /* q, c_sum はすでに0で初期化されているので何もしない */

  /* 各データ点の初期クラスタとu[], l[]を計算 */
  for(unsigned int i = 0; i < a.size(); ++i){
    Hamerly_point_all_ctrs(i, k, row(x, i), c, a, u, l, p);    
    row(c_sum,a[i]) += row(x,i); q[a[i]]++;
  }

  /* 初期クラスタ中心に対するデータ割り当てに対応して，
   * クラスター中心の更新や，各データ点のu[i],l[i]の更新が必要 */
  Hamerly_move_centers(c_sum, q, c, p, c_prev);
  Hamerly_update_bounds(p, a, u, l);
  
  return 0;
}

int Hamerly_repeat
(const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::matrix<double> &c,
 ublas::vector<double> &u,
 ublas::vector<double> &l,
 ublas::vector<unsigned int> &q,
 ublas::matrix<double> &c_sum,
 ublas::vector<double> &s,
 ublas::matrix<double> &c_prev,
 ublas::vector<double> &p /*tempのメモリ領域としても使う*/)
{
  int count = 0; /* Hamerlyの命題を満たさないデータ点の数*/
   
  /*  s[]の更新 -- O(dk^2) time
   *  s[j] = min_{jj != j} dist(c[jj], c[j])
   *      = min2_jj dist(c[j], c[jj]) */   
  for(unsigned int j = 0; j < k; ++j){
    for(unsigned int jj = 0; jj < k; ++jj){
      p[jj] = dist_cc(c, j, c, jj);
    }
    s[j] = min2(p);
  }

  /* 各点x[i]の最近接クラスタ割り当て -- O(ndk) time
   * 途中，Hamerlyの命題の条件をみたすか確認 (枝刈りを実行) */
  for(unsigned int i = 0; i < a.size(); ++i){
    double m = std::max((s[a[i]] / 2.0), l[i]);
    if(u[i] > m){     /* Hamerlyの命題の条件を満たさない */
      u[i] = dist_xc(row(x,i),c,a[i]);  /* upper boundをtight boundに */
      if(u[i] > m){   /* それでもHamerlyの命題の条件を満たさない */
	count++; /* Hamerlyの命題を満たさなかった点を数える */
	unsigned int aa = a[i]; /* 最近接クラスタが変化したか調べる */
	Hamerly_point_all_ctrs(i, k, row(x, i), c, a, u, l, p);
	if(aa != a[i]){ /* 最近接クラスタの変化によるパラメタの更新*/
	  row(c_sum,aa)   -= row(x,i); q[aa]--;
	  row(c_sum,a[i]) += row(x,i); q[a[i]]++;
	}
      }
    }
  }

  /* クラスタ中心の更新  -- O(dk) time */
  Hamerly_move_centers(c_sum, q, c, p, c_prev);

  /* u[], l[]の更新  -- O(n + k) time */
  Hamerly_update_bounds(p, a, u, l);

  return count; /* 最近接クラスタを再計算したデータ点の数 */
}
    
int Hamerly_main
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
  ublas::vector<double> l(n, 0);        /* lower bound */   

  /* j番目のクラスタを主語とする記号  */
  ublas::vector<unsigned int> q(k,0);  /* クラスタ中の点の数 */
  ublas::matrix<double> c_sum(k,d,0);  /* クラスタ中のベクトルの総和 */
  ublas::vector<double> s(k,0);        /* 最近接クラスタの重心との距離 */

  /* 下記は，イテレーションステップで用いる変数のメモリ */
  ublas::matrix<double> c_prev(k,d,0); /* 前のステップでのクラスタ中心 */
  ublas::vector<double> p(k,0);        /* クラスター中心の移動距離 */

  /* 初期ステップ */
  Hamerly_initialize(k, x, a, c, u, l, q, c_sum, c_prev, p);

  /* 繰り返しステップ */
  int t = 1;
  while(Hamerly_repeat(k, x, a, c, u, l, q, c_sum, s, c_prev, p) > 0){
    t++; /* 繰り返しステップの実行回数を数える */
  }  
  return t; /* 繰り返しステップの実行回数を返す */
}

#endif
