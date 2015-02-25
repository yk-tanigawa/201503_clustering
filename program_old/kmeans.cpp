#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <cfloat>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "time_bench.hpp"

using namespace std;
using namespace boost::numeric::ublas;
using namespace boost::numeric;

#define BENCH 1
/* k-meansの各種のアルゴリズムのベンチマークを取るときに1*/

enum algorithm {Lloyd, Hamerly}; /* 実行するアルゴリズム */


template <class T>
unsigned int argmin(ublas::vector<T> v)
{ /* vの最小要素の添字を返す */
  if(v.empty()){ return -1; /* エラー */    
  }else{
    unsigned int argmin = 0;
    T minval = v[argmin];
    for(unsigned int i = 0; i < v.size(); ++i){
      if(v[i] < minval){
	argmin = i; minval = v[i];
      }
    }
    return argmin;
  }
}

template <class T>
unsigned int argmin2 (ublas::vector<T> v)
{ /* vの2番目に小さい要素の添字を返す */
  if(v.size() < 2){ return -1; /* エラー */
  }else{
    unsigned int argmin1, argmin2;
    T min1, min2;
    if(v[0] <= v[1]){          
      min1 = v[0]; argmin1 = 0;
      min2 = v[1]; argmin2 = 1;
    }else{           
      min1 = v[1]; argmin1 = 1;
      min2 = v[0]; argmin2 = 0;
    }
    for(unsigned int i = 2; i < v.size(); ++i){
      if(v[i] <= min1){       /* 最小値が発見された */
	min2 = min1; argmin2 = argmin1;       
	min1 = v[i]; argmin1 = i;
      }else if(v[i] <= min2){ /* 準最小値が発見された */
	min2 = v[i]; argmin2 = i;
      }
    }
    return argmin2;
  }
}

template <class T>
unsigned int argmax(ublas::vector<T> v){
  if(v.empty()){ return -1; /* エラー */    
  }else{
    unsigned int argmax = 0; T maxval = v[argmax];
    for(unsigned int i = 0; i < v.size(); ++i){
      if(v[i] > maxval){
	argmax = i; maxval = v[i];
      }
    }
    return argmax;
  }
}

template <class T>
inline T min(ublas::vector<T> v){
  if(v.empty()){ return -1; /* エラー */    
  }else{
    T minval = v[0];
    for(unsigned int i = 0; i < v.size(); ++i){
      if(v.at(i) < minval){	
	minval = v.at(i);
      }
    }
    return minval;
  }
}

template <class T>
inline T min2(ublas::vector<T> v){
  /* vector vの要素で2番目に小さい要素を見つける */
  if(v.size() < 2){ return -1; /* エラー */
  }else{
    T min1 = v[0], min2 = v[1];
    if(min1 > min2){ /* swap する */
      T temp = min1; min1 = min2; min2 = temp;
    }
    for(unsigned int i = 2; i < v.size(); ++i){
      if(v[i] <= min1){       /* 最小値が発見された */
	min2 = min1; min1 = v[i];	
      }else if(v[i] <= min2){ /* 準最小値が発見された */
	min2 = v[i];
      }
    }
    return min2;
  }
}


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


#if 0
ublas::vector<double> map /* ベクトルに関するmap */
(ublas::vector<double> source, function<double(double)> func)
{
  ublas::vector<double> result(source.size());
  for (unsigned int i = 0; i < source.size(); ++i){
    result[i]= func(source[i]);
  }
  return result;
}

ublas::matrix<double> map /* 行列に関するmap */
(ublas::matrix<double> source,
 function<ublas::vector<double>(ublas::vector<double>)> func)
{
  ublas::matrix<double> result(source.size1(), source.size2());
  for (unsigned int i = 0; i < source.size1(); ++i){
    row(result, i)= func(row(source,i));
  }
  return result;
}


template <class T>
inline double norm(ublas::vector<T> v){
  /* ベクトルのノルムを計算する */
  double sum = 0;
  for(unsigned int i = 0; i < v.size(); ++i){
    sum += v[i] * v[i];
  }
  return sqrt(sum);
}

template <class T>
inline double dist(ublas::vector<T> x, ublas::vector<T> y){
  /* ベクトルのノルムから距離を定義 */
  ublas::vector<T> v = x - y;
  return norm(v); 
}

inline double dist(ublas::matrix<double> mat1, unsigned int r1,
		   ublas::matrix<double> mat2, unsigned int r2){
  /* mat1のr1行目のベクトルと, mat2のr2行目のベクトルの距離を計算 */
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

template <class T>
inline unsigned int argmin(ublas::vector<T> v){
  if(v.empty()){ return -1; /* エラー */    
  }else{
    int argmin = 0; T minval = v[argmin];
    for(unsigned int i = 0; i < v.size(); ++i){
      if(v[i] < minval){
	argmin = i; minval = v[i];
      }
    }
    return argmin;
  }
}

template <class T>
inline unsigned int argmax(ublas::vector<T> v){
  if(v.empty()){ return -1; /* エラー */    
  }else{
    int argmax = 0; T maxval = v[argmax];
    for(unsigned int i = 0; i < v.size(); ++i){
      if(v[i] > maxval){
	argmax = i; maxval = v[i];
      }
    }
    return argmax;
  }
}

inline void Hamerly_update_ul(const unsigned int n, const unsigned int k,
			      const ublas::matrix <double> data,
			      const ublas::matrix <double> c,       /* クラスタの重心 */
			      const ublas::vector <unsigned int> a, /* 各点x_iが属するクラスタ */
			      ublas::vector<double> &u,            /* upper bound */
			      ublas::vector<double> &l             /* lower bound */){
  ublas::vector<double> dist_xc(k, 0); /* x[i] とクラスター中心との距離*/
  for(unsigned int i = 0; i < n; ++i){
    for(unsigned int j = 0; j < k; ++j){
      dist_xc[j] = dist(data, i, c, j);
    }
    u[i] = dist_xc[a[i]]; /* u[i] = d(x[i], c[a[i]])*/
    l[i] = min2(dist_xc); /* l[i] = min2_j d(x[i], c[j]) */
    /* lはx[i]から2番目に近いクラスタ重心までの距離 */
    //cerr << i << ": " << u[i] << "\t" << l[i] << "\t" << dist_xc << endl;
  }
  //cerr << u << endl;
  //cerr << l << endl;
  return;
}

bool Hamerly_repeat(const unsigned int n, const unsigned int d, const unsigned int k,
		    const ublas::matrix <double> data,
		    ublas::vector<unsigned int> &q, /* クラスタ中の点の数 */
		    ublas::matrix<double> &c,       /* クラスタの重心 */
		    ublas::matrix<double> &c_sum,   /* クラスタのベクトル和*/
		    ublas::vector<double> &s,       /* 最近接クラスタの重心との距離 */
		    ublas::vector<unsigned int> &a, /* 各点x_iが属するクラスタ */
		    ublas::vector<double> &u,       /* upper bound */
		    ublas::vector<double> &l        /* lower bound */){
  /* Hamerlyの繰り返しステップ クラスタ重心の更新の有無をbooleanで返す */
  /* この繰り返しステップで何らかの更新があったか */
  volatile bool updated = false;

  /* クラスター中心の変化があったか */
  ublas::vector<bool> cluster_unchanged(k, true);

  /* まずs[j]を更新 
   * s[j] = min_{jj != j} dist(c[jj], c[j])
   *      = min2_jj dist(c[j], c[jj]) */   
  for(unsigned int j = 0; j < k; ++j){
    ublas::vector<double> distance(k, 0);
    for(unsigned int jj = 0; jj < k; ++jj){
      distance[jj] = dist(c, j, c, jj);
    }
    s[j] = min2(distance);
  }

  //cerr << "s: " << s << endl;
  
  /* 各点x[i]について，Hamerlyの命題の条件をみたすか確認しながら
   * (つまり，適宜枝刈りを実行しながら)
   * 最近接クラスタの再割当てを行う 
   * 同時に，のちのクラスタ重心の更新のために
   * c_sum[], q[] を正しい値に保つ */
  for(unsigned int i = 0; i < n; ++i){
    double m = max((s[a[i]] / 2.0), l[i]);
    //cerr << "i = " << i << ", m = " << m << ", u[i] = " << u[i] << endl;
    if(u[i] > m){     /* Hamerlyの命題の条件を満たさない */
      //cerr << i << u[i] << " -> ";
      u[i] = dist(data,i,c,a[i]);  /* upper boundを厳密な値に更新 */
      //cerr << u[i] << endl;
      if(u[i] > m){   /* Hamerlyの命題の条件を満たさない */
	unsigned int aa = a[i]; /* 最近接クラスタが変化したか調べる */
	{ /* 最近接クラスタを再計算 */
	  ublas::vector<double> distance(k);
	  for(unsigned int j = 0; j < k; ++j){
	    distance[j] = dist(data,i,c,j);
	  }
	  a[i] = argmin(distance);
	}
	if(aa != a[i]){ /* 最近接クラスタの変化によるパラメタの更新*/
	  /* 点x[i]の最近接クラスタは，aa から a[i] に変化した */
	  updated = true; /* この繰り返しステップ更新があった */
	  /* 後のクラスター重心の更新のために q, c_sum の更新が必要 
	   * 同時に、更新があったクラスターのindexをメモしておく */
	  row(c_sum,aa)   -= row(data,i); q[aa]--;
	  row(c_sum,a[i]) += row(data,i); q[a[i]]++;
	  cluster_unchanged[aa] = false;
	  cluster_unchanged[a[i]] = false;
	  {
	    ublas::vector<double> dist_xc(k, 0);
	    /* x[i] とクラスター中心との距離*/
	    for(unsigned int j = 0; j < k; ++j){
	      dist_xc[j] = dist(data, i, c, j);
	    }
	    u[i] = dist_xc[a[i]]; /* u[i] = d(x[i], c[a[i]])*/
	    l[i] = min2(dist_xc); /* l[i] = min2_j d(x[i], c[j]) */
	  }
	}
      }
    }
  }

  ublas::matrix<double> c_prev = c; /* 変化前のクラスタ重心 */
  if(updated){ /* いずれかのデータ点の最近接クラスタが変化した場合 */
    ublas::vector<double> p(k, 0); /* 各クラスタの重心の移動距離 */
    for(unsigned int j = 0; j < k; ++j){
      if(!cluster_unchanged[j]){ /* j番目のクラスタ重心が移動していた場合 */
	row(c,j) = row(c_sum,j) / q[j]; /* クラスタ重心を更新 */
	p[j] = dist(c_prev, j, c, j); /* クラスター重心の移動距離 */
      }
    }
    unsigned int r = argmax(p);
    for(unsigned int i = 0; i < n; ++i){
      u[i] += p[a[i]];  l[i] -= p[r];
    }
  }

  return updated; /* データ点の最近接クラスタに変化があったかを返す */
}

void Hamerly_init(const unsigned int n, const unsigned int d, const unsigned int k,
		  const ublas::matrix <double> data,
		  ublas::vector<unsigned int> &q,   /* クラスタ中の点の数 */
		  ublas::matrix<double> &c,         /* クラスタの重心 */
		  ublas::matrix<double> &c_sum,     /* クラスタのベクトル和*/
		  ublas::vector<unsigned int> &a,   /* 各点x_iが属するクラスタ */
		  ublas::vector<double> &u,         /* upper bound */
		  ublas::vector<double> &l          /* lower bound */){
  /* 各クラスタの代表点をランダムに選択 */
  for(unsigned int j = 0; j < k; ++j){
    ublas::zero_vector<double> zero(d);
    row(c,j) = row(data, j);
    /* j番目のクラスタ中心をj番目のデータ点と一致させて初期化した */
  }

  /* 全てのデータ点について，初期クラスタa[i]を計算 */
  for(unsigned int i = 0; i < n; ++i){
    { /* argmin_j dist(x(i), c(j))の計算 */
      ublas::vector<double> distance(k);
      for(unsigned int j = 0; j < k; ++j){
	distance[j] = dist(data,i,c,j);
      }
      a[i] = argmin(distance);
    }
    /* クラスタ内の点の数，クラスタの点のベクトル和も更新 */
    q[a[i]]++;  row(c_sum,a[i]) += row(data,i);    
  }
  
  /* upper bound, lower bound の更新 */
  Hamerly_update_ul(n, k, data, c, a, u, l);

  /* 初期クラスタ中心に対するデータ割り当てに対応して，
   * クラスター中心の更新や，各データ点のu[i],l[i]の更新が必要 */
  
  ublas::matrix<double> c_prev = c; /* 変化前のクラスタ重心 */

  ublas::vector<double> p(k, 0); /* 各クラスタの重心の移動距離 */
  for(unsigned int j = 0; j < k; ++j){  
    row(c,j) = row(c_sum,j) / q[j]; /* クラスタ重心を更新 */
    p[j] = dist(c_prev, j, c, j); /* クラスター重心の移動距離 */
  }
  unsigned int r = argmax(p);

  for(unsigned int i = 0; i < n; ++i){
    u[i] += p[a[i]];  l[i] -= p[r];
  }

  
  return;
}

double Hamerly_sqerr(const unsigned int n, const unsigned int k,
		     const ublas::matrix<double>data,
		     const ublas::matrix<double> c,       /* クラスタの重心 */
		     const ublas::vector<unsigned int> a  /* 各点x_iが属するクラスタ */){
  /* クラスタリング結果に基づき二乗平均二乗誤差を計算して返す */
  double sum = 0;
  for(unsigned int i = 0; i < n; ++i){
    sum += dist(data,i,c,a[i]) * dist(data,i,c,a[i]);
  }
  return (sum / n);
}


void Hamerly_main(const unsigned int n, const unsigned int d, const unsigned int k,
	     const ublas::matrix<double> data){
  if(n < k){
    cerr << "data size n is smaller than cluster size d" << endl;
    exit(1);
  }else{
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL); /* 時間計測開始 */

    /* j番目のクラスタを主語とする記号  */
    ublas::vector<unsigned int> q(k,0); /* クラスタ中の点の数 */
    ublas::matrix<double> c(k,d,0);     /* クラスタの重心 */
    ublas::matrix<double> c_sum(k,d,0); /* クラスタ中のベクトルの総和 */
    ublas::vector<double> s(k,0);       /* 最近接クラスタの重心との距離 */
    
    /* i番目のデータ点を主語とする記号 */
    ublas::vector<unsigned int> a(n, 0);  /* 各点x_iが属するクラスタ */
    ublas::vector<double> u(n, 0);        /* upper bound */
    ublas::vector<double> l(n, 0);        /* lower bound */   

    Hamerly_init(n, d, k, data, q, c, c_sum, a, u, l);

    int t = 0; /* 繰り返しステップを何回繰り返したか */
    
    //cerr << "a: " << a << ": t = " << t << endl;
    
    while(Hamerly_repeat(n, d, k, data, q, c, c_sum, s, a, u, l)){
      t++;
      //cerr << "a: " << a << ": t = " << t << endl;
    }
    t++; /* 更新されなかった最後の1回も実行されている */

    gettimeofday(&t_end, NULL); /* 時間計測終了 */

    /* 結果を表示 */
    cerr << "a: " << a << endl;
    cerr << c << endl;
    
    cout << n << "\t"                               /* サンプル数 */
	 << d << "\t"                               /* 次元 */
      	 << k << "\t"                               /* クラスタ数 */
	 << diff_timeval(t_start, t_end) << "\t"    /* 実行時間(us) */
	 << Hamerly_sqerr(n, k, data, c, a) << "\t" /* 平均二乗誤差 */
	 << t << endl;                              /* 繰り返し回数*/

    return;
  }
}

#endif


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
(const unsigned int n,
 const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::matrix<double> &c,
 ublas::vector<double> &u,
 ublas::vector<double> &l,
 ublas::vector<unsigned int> &q,
 ublas::matrix<double> &c_sum,
 ublas::matrix<double> &c_prev,
 ublas::vector<double> &p
 )
{
  /* q, c_sum は0で初期化されていると仮定 */
  for(unsigned int i = 0; i < n; ++i){
    Hamerly_point_all_ctrs(i, k, row(x, i), c, a, u, l, p);

    q[a[i]]++;
    row(c_sum,a[i]) += row(x,i);        
  }

  /* 初期クラスタ中心に対するデータ割り当てに対応して，
   * クラスター中心の更新や，各データ点のu[i],l[i]の更新が必要 */

  Hamerly_move_centers(c_sum, q, c, p, c_prev);

  Hamerly_update_bounds(p, a, u, l);
  
  return 0;
}

int Hamerly_repeat
(const unsigned int n,
 const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::matrix<double> &c,
 ublas::vector<double> &u,
 ublas::vector<double> &l,
 ublas::vector<unsigned int> &q,
 ublas::matrix<double> &c_sum,
 ublas::vector<double> &s,
 ublas::matrix<double> &c_prev,
 ublas::vector<double> &p /*tempとしても使う*/
 )
{
  int count = 0; /* Hamerlyの命題を満たさないデータ点の数*/
   
  /* まずs[j]を更新 
   * s[j] = min_{jj != j} dist(c[jj], c[j])
   *      = min2_jj dist(c[j], c[jj]) */   
  for(unsigned int j = 0; j < k; ++j){
    for(unsigned int jj = 0; jj < k; ++jj){
      p[jj] = dist_cc(c, j, c, jj);
    }
    s[j] = min2(p);
  }

  /* 各点x[i]について，Hamerlyの命題の条件をみたすか確認しながら
   * (つまり，適宜枝刈りを実行しながら)
   * 最近接クラスタの再割当てを行う 
   * 同時に，のちのクラスタ重心の更新のために
   * c_sum[], q[] を正しい値に保つ */
  for(unsigned int i = 0; i < n; ++i){
    double m = max((s[a[i]] / 2.0), l[i]);
    if(u[i] > m){     /* Hamerlyの命題の条件を満たさない */
      u[i] = dist_xc(row(x,i),c,a[i]);  /* upper boundを厳密な値に更新 */
      if(u[i] > m){   /* Hamerlyの命題の条件を満たさない */
	count++;
	unsigned int aa = a[i]; /* 最近接クラスタが変化したか調べる */
	Hamerly_point_all_ctrs(i, k, row(x, i), c, a, u, l, p);
	if(aa != a[i]){ /* 最近接クラスタの変化によるパラメタの更新*/
	  row(c_sum,aa)   -= row(x,i); q[aa]--;
	  row(c_sum,a[i]) += row(x,i); q[a[i]]++;
	}
      }
    }
  }
  Hamerly_move_centers(c_sum, q, c, p, c_prev);
  Hamerly_update_bounds(p, a, u, l);

  return count;
}
    

#if 1
int Hamerly_main
(const unsigned int n,
 const unsigned int d,
 const unsigned int k,
 const ublas::matrix<double> x,
 ublas::vector<unsigned int> &a,
 ublas::matrix<double> &c)
{
  #if BENCH
  struct timeval t_start, t_end;
  gettimeofday(&t_start, NULL); /* 時間計測開始 */
  #endif


  /* メモリ確保 */
  
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

  Hamerly_initialize(n, k, x, a, c, u, l, q, c_sum, c_prev, p);
  int t = 1;
  while(Hamerly_repeat(n, k, x, a, c, u, l, q, c_sum, s, c_prev, p) > 0){
    t++;
  }

  
  #if BENCH
  gettimeofday(&t_end, NULL); /* 時間計測終了 */
  #endif
  
  return 0;
}
#endif

#if 0
int maptest()
{
  ublas::vector<double> vec(5, 0);
  vec[0] = 0.0;
  vec[1] = 1.1;
  vec[2] = -2.4;
  vec[3] = 3.5;
  vec[4] = 4.6;

  cout << vec << endl;
  
  cout <<  map(vec, [](double x){ return x + 1.0; }) << endl;
  return 0;
}
#endif

/* 全てのアルゴリズムに共通 */
int initialize_centers
(const unsigned int k,
 const ublas::matrix<double> x,
 ublas ::matrix<double> &c)
{
  /* とりあえず適当に最初のデータ点を充てる */
  for(unsigned int j = 0; j < k; ++j){
    row(c, j) = row(x, j);
  }
  return 0;
}


int main(const int argc, char *argv[])
{   if(argc < 6){ /* 実行時オプションが少ないときのエラー処理 */
    /* 使い方を表示して異常終了 */
    cerr << "usage: " << argv[0]
	 << " <method> <n> <d> <k> <data> [results] [center]\n"
	 << " <method> = Lloyd | Hamerly\n"
	 << " n : number of data points\n"
      	 << " d : dimension\n"
	 << " k : number of clusters\n"
         << " data    : file path of the data file\n"
	 << " results : (optional) index of the cluster\n"
	 << " center  : (optional) coordinate of cluster center"
	 << endl;
    exit(1);
  }else{
    /* 実行時オプションの処理 */    

    enum algorithm method; /* 実行するアルゴリズム */
    if     (strcmp(argv[1], "Lloyd")   == 0){ method = Lloyd; }
    else if(strcmp(argv[1], "Hamerly") == 0){ method = Hamerly; }
    else{ /* 指定されたアルゴリズムが不明 */
      cerr << "unknown method : " << argv[1] << endl;
      exit(1);
    }

    unsigned int n = atoi(argv[2]), d = atoi(argv[3]), k = atoi(argv[4]);
    char *input_file = argv[5];

    if(n < k){ /* n < d のときのエラー処理 */
      cerr << "data size n is smaller than cluster size d\n"
	   << "n = " << n << ", d = " << d
	   << endl;
      exit(1);
    }
    
    ublas::matrix<double> data(n, d); 
    { /* n個のd次元データをファイルから読み込む */
      ifstream fs(input_file);
      if(fs.fail()){ exit(1); }
      for(unsigned int i = 0; i < n; ++i){
	for(unsigned int j = 0; j < d; ++j){
	  fs >> data(i,j);
	}
      }
      fs.close();
    }

    /* クラスタリングの結果を書き込むメモリ領域の確保 */
    ublas::vector<unsigned int> a(n, 0);  /* 各点x_iが属するクラスタ */
    ublas::matrix<double> c(k,d,0);       /* クラスタの重心 */

    initialize_centers(k, data, c);
	
    if(method == Hamerly){
      Hamerly_main(n, d, k, data, a, c);
    }


    cerr << c << endl;
    cerr << a << endl;
    
    return 0;
  }
}
