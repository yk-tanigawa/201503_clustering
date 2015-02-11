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

#include "my_vector.hpp"
#include "time_bench.hpp"

using namespace std;

template <class T>
inline double dist(vector<T> x, vector<T> y){
  return norm(x - y); /* ベクトルのノルムから距離を定義 */
}

template <class T>
inline int argmin(vector<T> v){
  if(v.size() < 2){
    return 0;
  }else{
    int argmin = 0; T minval = v.at(argmin);     
    for(int i = 1; i < v.size(); ++i){
      if(v.at(i) < minval){
	argmin = i; minval = v.at(i);
      }
    }
    return argmin;
  }
}

template <class T>
inline int min(vector<T> v){
  if(v.size() < 2){
    return v.at(0);
  }else{
    T minval = v.at(argmin);     
    for(int i = 1; i < v.size(); ++i){
      if(v.at(i) < minval){
	minval = v.at(i);
      }
    }
    return minval;
  }
}

void Hamerly_update_ul(const int n, const int k,
		       const vector<vector <double> > data,
		       const vector<vector<double> > c, /* クラスタの重心 */
		       const vector<int> a,    /* 各点x_iが属するクラスタ */
		       vector<double> &u,   /* upper bound */
		       vector<double> &l    /* lower bound */);
bool Hamerly_repeat(const int n, const int d, const int k,
		    const vector<vector <double> > data,
		    vector<int> &q,      /* クラスタ中の点の数 */
		    vector<vector<double> > &c,     /* クラスタの重心 */
		    vector<vector<double> > &c_sum, /* クラスタのベクトル和*/
		    vector<double> &s,   /* 最近接クラスタの重心との距離 */
		    vector<int> &a,      /* 各点x_iが属するクラスタ */
		    vector<double> &u,   /* upper bound */
		    vector<double> &l    /* lower bound */);
void Hamerly_init(const int n, const int d, const int k,
		  const vector<vector <double> > data,
		  vector<int> &q,      /* クラスタ中の点の数 */
		  vector<vector<double> > &c,     /* クラスタの重心 */
		  vector<vector<double> > &c_sum, /* クラスタのベクトル和*/
		  vector<int> &a,      /* 各点x_iが属するクラスタ */
		  vector<double> &u,   /* upper bound */
		  vector<double> &l    /* lower bound */);
double Hamerly_sqerr(const int n, const int k,
		     const vector<vector<double> > data,
		     vector<vector<double> > &c,  /* クラスタの重心 */
		     vector<int> &a      /* 各点x_iが属するクラスタ */);
void Hamerly(const int n, const int d, const int k,
	     const vector<vector <double> > data);


/* 下記プログラム本体 */


void Hamerly_update_ul(const int n, const int k,
		       const vector<vector <double> > data,
		       const vector<vector<double> > c, /* クラスタの重心 */
		       const vector<int> a,    /* 各点x_iが属するクラスタ */
		       vector<double> &u,   /* upper bound */
		       vector<double> &l    /* lower bound */){
  for(int i = 0; i < n; ++i){
    u.at(i) = dist(data.at(i), c.at(a.at(i)));
    double distance = dist(data.at(i), c.at(0));
    double min = distance;
    for(int j = 1; j < k; ++j){
      distance = dist(data.at(i), c.at(j));
      if(distance < min){ min = distance; }
    }
    l.at(i) = min;
  }
  return;
}

bool Hamerly_repeat(const int n, const int d, const int k,
		    const vector<vector <double> > data,
		    vector<int> &q,      /* クラスタ中の点の数 */
		    vector<vector<double> > &c,     /* クラスタの重心 */
		    vector<vector<double> > &c_sum, /* クラスタのベクトル和*/
		    vector<double> &s,   /* 最近接クラスタの重心との距離 */
		    vector<int> &a,      /* 各点x_iが属するクラスタ */
		    vector<double> &u,   /* upper bound */
		    vector<double> &l    /* lower bound */){
  /* Hamerlyの繰り返しステップ クラスタ重心の更新の有無をbooleanで返す */
  bool updated = false;

#if 0
  cerr << "***** rep *****" << endl;
  cerr << "j" << "\t" << "q.at(j)" << "\t" << "c(j)" << endl;
  for(int j = 0; j < k; ++j){
    cerr << j << "\t" << q.at(j) << "\t" << c.at(j);
  }

  cerr << "s(j)" << endl;
#endif
  

  /* まずs(j)を更新 
   * s[j] = min_{jj != j} dist(c[jj], c[j]) */
  for(int j = 0; j < k; ++j){
    double min = DBL_MAX;
    for(int jj = 0; jj < k; ++jj){
      if(jj != j){
	double distance = dist(c.at(j), c.at(jj));
	if(distance < min){ min = distance; }
      }
    }
    s.at(j) = min;
#if 0
    cerr << j << "\t" << s.at(j) << endl;
#endif
  }

  for(int i = 0; i < n; ++i){
    double m = max((s.at(a.at(i)) / 2.0), l.at(i));
    if(u.at(i) > m){     /* Hamerlyの命題の条件を満たさない */
      u.at(i) = dist(data.at(i), c.at(a.at(i))); /* upper boundを更新 */
      if(u.at(i) > m){   /* Hamerlyの命題の条件を満たさない */
	int aa = a.at(i); /* 最近接クラスタが変化したか調べる */
	{ /* 最近接クラスタを再計算 */
	  vector<double> distance(k);
	  for(int j = 0; j < k; ++j){
	    distance.at(j) = dist(data.at(i), c.at(j));
	  }
	  a.at(i) = argmin(distance);
	}
	if(aa != a.at(i)){ /* 最近接クラスタの変化によるパラメタの更新*/
	  updated = true; q.at(aa)--; c_sum.at(aa) -= data.at(i);	 
	  q.at(a.at(i))++;  c_sum.at(a.at(i)) += data.at(i);
	  Hamerly_update_ul(n, k, data, c, a, u, l);
	}
      }
    }
  }

  if(updated){ /* いずれかのデータ点の最近接クラスタが変化した場合 */
    /* 各クラスタに関するパラメタの更新 */
    vector<double> p(k, 0); /* 各クラスタの重心の移動距離 */
    double p_max = -1;      /* 番兵で初期化 */
    for(int j = 0; j < k; ++j){
      vector<double> c_prev = c.at(j);
      c.at(j) = (1.0 / q.at(j)) * c_sum.at(j);
      p.at(j) = dist(c_prev, c.at(j));
      if(p_max < p.at(j)){ p_max = p.at(j); }
    }
    for(int i = 0; i < n; ++i){
      u.at(i) += p.at(a.at(i));
      l.at(i) -= p_max;     
    }
  }
  
  return updated; /* データ点の最近接クラスタに変化があったかを返す */
}


void Hamerly_init(const int n, const int d, const int k,
		  const vector<vector <double> > data,
		  vector<int> &q,      /* クラスタ中の点の数 */
		  vector<vector<double> > &c,     /* クラスタの重心 */
		  vector<vector<double> > &c_sum, /* クラスタのベクトル和*/
		  vector<int> &a,      /* 各点x_iが属するクラスタ */
		  vector<double> &u,   /* upper bound */
		  vector<double> &l    /* lower bound */){
  /* 各クラスタの代表点をランダムに選択 */
  for(int j = 0; j < k; ++j){
    vector<double> zero(d, 0);
    q.at(j) = 0;  c.at(j) = data.at(j); c_sum.at(j) = zero;
    /* j番目のクラスタ中心をj番目のデータと一致させて初期化した */
  }
  for(int i = 0; i < n; ++i){
    int argmin_j = 0;
    { /* argmin_j dist(x(i), c(j))の計算 */
      double distance = dist(data.at(i), c.at(0));
      double min = distance;
      for(int j = 1; j < k; ++j){
	distance = dist(data.at(i), c.at(j));
	if(distance < min){
	  min = distance;  argmin_j = j;
	}
      }
    }
    a.at(i) = argmin_j;
    q.at(argmin_j)++;
    c_sum.at(argmin_j) += data.at(i);
  }
  for(int j = 0; j < k; ++j){
    c.at(j) = (1.0 / q.at(j)) * c_sum.at(j);
  }
  /* upper bound, lower bound の更新 */
  Hamerly_update_ul(n, k, data, c, a, u, l);
#if 0
  cerr << "***** init *****" << endl;
  cerr << "j" << "\t" << "q.at(j)" << "\t" << "c(j)" << endl;
  for(int j = 0; j < k; ++j){
    cerr << j << "\t" << q.at(j) << "\t" << c.at(j);
  }
#endif
  return;
}

double Hamerly_sqerr(const int n, const int k,
		     const vector<vector<double> > data,
		     vector<vector<double> > &c,  /* クラスタの重心 */
		     vector<int> &a      /* 各点x_iが属するクラスタ */){
  /* クラスタリング結果に基づき二乗平均二乗誤差を計算して返す */
  double sum = 0;
  for(int i = 0; i < n; ++i){
    sum += dist(data.at(i), c.at(a.at(i))) * dist(data.at(i), c.at(a.at(i)));
  }
  return (sum / n);
}


void Hamerly(const int n, const int d, const int k,
	     const vector<vector <double> > data){
  if(n < k){
    cerr << "data size n is smaller than cluster size d" << endl;
    exit(1);
  }else{
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL); /* 時間計測開始 */

    /* j番目のクラスタを主語とする記号  */
    vector<int> q(k);                 /* クラスタ中の点の数 */
    vector<vector<double> > c(k);     /* クラスタの重心 */
    vector<vector<double> > c_sum(k); /* クラスタ中のベクトルの総和 */
    vector<double> s(k, 0);           /* 最近接クラスタの重心との距離 */
    
    /* i番目のデータ点を主語とする記号 */
    vector<int> a(n);   /* 各点x_iが属するクラスタ */
    vector<double> u(n);   /* upper bound */
    vector<double> l(n);   /* lower bound */   
    
    Hamerly_init(n, d, k, data, q, c, c_sum, a, u, l);

    int t = 0;
    while(Hamerly_repeat(n, d, k, data, q, c, c_sum, s, a, u, l)){ t++; }

    //Hamerly_repeat(n, d, k, data, q, c, c_sum, s, a, u, l);
    //Hamerly_repeat(n, d, k, data, q, c, c_sum, s, a, u, l);

    gettimeofday(&t_end, NULL); /* 時間計測終了 */
    
    cout << n << "\t"                               /* サンプル数 */
	 << d << "\t"                               /* 次元 */
      	 << k << "\t"                               /* クラスタ数 */
	 << diff_timeval(t_start, t_end) << "\t"    /* 実行時間(us) */
	 << Hamerly_sqerr(n, k, data, c, a) << "\t" /* 平均二乗誤差 */
	 << t << endl;                              /* 繰り返し回数*/
    return;
  }
}

int main(int argc, char *argv[]){
  if(argc < 5){
    cerr << "usage: " << argv[0]
	 << " <n: data num> <d: dimension> <k: cluster num> <data file>" << endl;
    exit(1);
  }else{
    int n = atoi(argv[1]), d = atoi(argv[2]), k = atoi(argv[3]);
    char *file = argv[4];
    vector<vector <double> > data(n);
    {
      ifstream fs(file);
      if(fs.fail()){ exit(1); }
      for(int i = 0; i < n; ++i){
	data.at(i).resize(d, 0);
	for(int j = 0; j < d; ++j){
	  fs >> data.at(i).at(j);
	}
      }
      fs.close();
    }
    Hamerly(n, d, k, data);
    return 0;
  }
}
