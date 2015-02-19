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

#include "time_bench.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::ublas;
using namespace boost::numeric;

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


template <class T>
inline T min(ublas::vector<T> v){
  if(v.empty()){ return -1; /* エラー */    
  }else{
    T minval = v.at(0);
    for(unsigned int i = 0; i < v.size(); ++i){
      if(v.at(i) < minval){	
	minval = v.at(i);
      }
    }
    return argmin;
  }
}

inline void Hamerly_update_ul(const unsigned int n, const unsigned int k,
			      const ublas::matrix <double> data,
			      const ublas::matrix <double> c,       /* クラスタの重心 */
			      const ublas::vector <unsigned int> a, /* 各点x_iが属するクラスタ */
			      ublas::vector<double> &u,            /* upper bound */
			      ublas::vector<double> &l             /* lower bound */){
  for(unsigned int i = 0; i < n; ++i){
    u[i] = dist(data,i,c,a[i]);
    double distance = dist(data,i,c,0);
    double min = distance;
    for(unsigned int j = 1; j < k; ++j){
      if(j != a[i]){
	distance = dist(data, i, c, j);
	if(distance < min){ min = distance; }
      }
    }
    l[i] = min;
  }
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
  bool updated = false;

  /* まずs[j]を更新 
   * s[j] = min_{jj != j} dist(c[jj], c[j]) */
  for(unsigned int j = 0; j < k; ++j){
    double min = DBL_MAX;
    for(unsigned int jj = 0; jj < k; ++jj){
      if(jj != j){
	double distance = dist(c,j,c,jj);
	if(distance < min){ min = distance; }
      }
    }
    s[j] = min;
  }

#if 0
  cerr << "j" << "\t" << "q.at(j)" << "\t" << "c(j)" << endl;
  for(int j = 0; j < k; ++j){
    cerr << j << "\t" << q[j] << "\t" << (ublas::vector<double>)row(c,j) << endl;
  }
#endif

#if 0
  cerr << "!\tq:\t" << q << endl;
  cerr << "\ts:\t" << s << endl;
#endif
  
  ublas::matrix<double> c_prev = c; /* 変化前のクラスタ重心 */
  
  for(unsigned int i = 0; i < n; ++i){
    double m = max((s[a[i]] / 2.0), l[i]);
    if(u[i] > m){     /* Hamerlyの命題の条件を満たさない */
      u[i] = dist(data,i,c,a[i]);  /* upper boundを更新 */           
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
	  updated = true; /* q, c_sum の更新が必要 */
	  row(c_sum,aa)   -= row(data,i);  q[aa]--;
	  // row(c,aa)        = row(c_sum,aa)   / q[aa];
	  row(c_sum,a[i]) += row(data,i);  q[a[i]]++;
	  // row(c,a[i])      = row(c_sum,a[i]) / q[a[i]];	 
	}
      }
    }
  }

  if(updated){ /* いずれかのデータ点の最近接クラスタが変化した場合 */
    ublas::vector<double> p(k, 0); /* 各クラスタの重心の移動距離 */
    for(unsigned int j = 0; j < k; ++j){
      p[j] = dist(c_prev,j, c,j);
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
    q[j] = 0; row(c,j) = row(data, j); row(c_sum,j) = zero;    
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
  
  /* クラスタ重心の計算 */
  for(unsigned int j = 0; j < k; ++j){
    row(c,j) = (1.0 / q[j]) * row(c_sum,j);
  }
  /* upper bound, lower bound の更新 */
  Hamerly_update_ul(n, k, data, c, a, u, l);
  
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


void Hamerly(const unsigned int n, const unsigned int d, const unsigned int k,
	     const ublas::matrix<double> data){
  if(n < k){
    cerr << "data size n is smaller than cluster size d" << endl;
    exit(1);
  }else{
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL); /* 時間計測開始 */

    /* j番目のクラスタを主語とする記号  */
    ublas::vector<unsigned int> q(k);  /* クラスタ中の点の数 */
    ublas::matrix<double> c(k,d);        /* クラスタの重心 */
    ublas::matrix<double> c_sum(k,d);    /* クラスタ中のベクトルの総和 */
    ublas::vector<double> s(k, 0);     /* 最近接クラスタの重心との距離 */
    
    /* i番目のデータ点を主語とする記号 */
    ublas::vector<unsigned int> a(n);  /* 各点x_iが属するクラスタ */
    ublas::vector<double> u(n);        /* upper bound */
    ublas::vector<double> l(n);        /* lower bound */   

    cout << "a: " << a << endl;
    
    Hamerly_init(n, d, k, data, q, c, c_sum, a, u, l);

    cout << "data: " << data << endl;
    
    int t = 1;
    while(Hamerly_repeat(n, d, k, data, q, c, c_sum, s, a, u, l)){ t++; }

    gettimeofday(&t_end, NULL); /* 時間計測終了 */
    
    cout << n << "\t"                               /* サンプル数 */
	 << d << "\t"                               /* 次元 */
      	 << k << "\t"                               /* クラスタ数 */
	 << diff_timeval(t_start, t_end) << "\t"    /* 実行時間(us) */
	 << Hamerly_sqerr(n, k, data, c, a) << "\t" /* 平均二乗誤差 */
	 << t << endl;                              /* 繰り返し回数*/

    cout << "a: " << a << endl;

    return;
  }
}

int main(int argc, char *argv[]){
  if(argc < 5){
    cerr << "usage: " << argv[0]
	 << " <n: data num> <d: dimension> <k: cluster num> <data file>" << endl;
    exit(1);
  }else{
    unsigned int n = atoi(argv[1]), d = atoi(argv[2]), k = atoi(argv[3]);
    char *file = argv[4];
    ublas::matrix<double> data(n, d); /* n個のd次元データを読み込む */
    {
      ifstream fs(file);
      if(fs.fail()){ exit(1); }
      for(unsigned int i = 0; i < n; ++i){
	for(unsigned int j = 0; j < d; ++j){
	  fs >> data(i,j);
	}
      }
      fs.close();
    }
    Hamerly(n, d, k, data);
    return 0;
  }
}
