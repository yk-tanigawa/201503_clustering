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
  volatile bool updated = false; /* この繰り返しステップで何らかの更新があったか */
  ublas::vector<bool> cluster_unchanged(k, true);
#if 0
  cerr << "updated:" << updated 
       << " cluster_unchanged:" << cluster_unchanged
       << " Hamerly_cond:" << Hamerly_cond
       << endl;
#endif
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

  cerr << "s: " << s << endl;
  
  /* 各点x[i]について，Hamerlyの命題の条件をみたすか確認しながら
   * (つまり，適宜枝刈りを実行しながら)
   * 最近接クラスタの再割当てを行う 
   * 同時に，のちのクラスタ重心の更新のために
   * c_sum[], q[] を正しい値に保つ */
  for(unsigned int i = 0; i < n; ++i){
    double m = max((s[a[i]] / 2.0), l[i]);
    cerr << "i = " << i << ", m = " << m << ", u[i] = " << u[i] << endl;
    if(u[i] > m){     /* Hamerlyの命題の条件を満たさない */
      cerr << i << u[i] << " -> ";
      u[i] = dist(data,i,c,a[i]);  /* upper boundを厳密な値に更新 */
      cerr << u[i] << endl;
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

#if 0
  cerr << "updated:" << updated 
       << " cluster_unchanged:" << cluster_unchanged
       << " Hamerly_cond:" << Hamerly_cond
       << endl;
  
  cerr << "u:" << u << endl;
#endif
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
    
    cerr << "a: " << a << ": t = " << t << endl;

    Hamerly_repeat(n, d, k, data, q, c, c_sum, s, a, u, l);
    
    while(Hamerly_repeat(n, d, k, data, q, c, c_sum, s, a, u, l)){
      t++;
      cerr << "a: " << a << ": t = " << t << endl;
    }

    gettimeofday(&t_end, NULL); /* 時間計測終了 */
    
    cout << n << "\t"                               /* サンプル数 */
	 << d << "\t"                               /* 次元 */
      	 << k << "\t"                               /* クラスタ数 */
	 << diff_timeval(t_start, t_end) << "\t"    /* 実行時間(us) */
	 << Hamerly_sqerr(n, k, data, c, a) << "\t" /* 平均二乗誤差 */
	 << t << endl;                              /* 繰り返し回数*/

    cerr << "a: " << a << endl;
    cerr << c << endl;
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
