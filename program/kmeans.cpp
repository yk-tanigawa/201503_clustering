#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#ifndef _UBLAS_
#define _UBLAS_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#endif

#include "time_bench.hpp"
#include "my_ublas.hpp"

using namespace std;
using namespace boost::numeric::ublas;
using namespace boost::numeric;

#define BENCH 1
/* k-meansの各種のアルゴリズムのベンチマークを取るときに1*/

enum algorithm {Lloyd, Hamerly}; /* 実行するアルゴリズム */
static const std::vector<string> algo_name{"Lloyd", "Hamerly"};

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
    double m = max((s[a[i]] / 2.0), l[i]);
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


#if BENCH
int bench_show
(ostream &stream,            /* 出力先 */
 enum algorithm method,      /* 使用したアルゴリズム */
 const unsigned int n,       /* データ点数 */
 const unsigned int d,       /* 次元数 */
 const unsigned int k,       /* クラスタ数 */
 const struct timeval start, /* 開始時刻 */
 const struct timeval end,   /* 終了時刻 */
 const unsigned int rep,     /* 繰り返し回数 */
 const ublas::matrix<double> x,       /* データ点 */ 
 const ublas::vector<unsigned int> a, /* 割り当てクラスタ */
 const ublas::matrix<double> c)       /* クラスタ中心 */
{  
  double err_sum = 0; /* 二乗誤差和を計算する */
  for(unsigned int i = 0; i < n; ++i){
    err_sum += dist_xc(row(x, i), c, a[i]) * dist_xc(row(x, i), c, a[i]);
  }
     
  stream << algo_name[method] << "\t" << n << "\t" << d << "\t" << k << "\t" 	 
	 << diff_timeval(start, end)    << "\t" /* 実行時間(us) */
	 << (err_sum / n)  << "\t" << rep << endl;
  return 0;
}
#endif

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
    ublas::vector<unsigned int> a(n, 0); /* 各点x_iが属するクラスタ */
    ublas::matrix<double> c(k,d,0);      /* クラスタの重心 */
    unsigned int rep = 0;                /* 繰り返し回数 */ 
#if BENCH /* プログラム実行時間のベンチマークを取る */
    struct timeval t_start, t_end;
#endif

    /* 初期クラスタ中心を何らかの方法で得る */
    initialize_centers(k, data, c);

    
#if BENCH
      gettimeofday(&t_start, NULL); /* 時間計測開始 */
#endif
    
      if(method == Hamerly){
	rep = Hamerly_main(n, d, k, data, a, c);
      }

#if BENCH
    gettimeofday(&t_end, NULL); /* 時間計測終了 */
#endif

    { /* 結果を出力 */
      ofstream fs(algo_name[method] + "_"
		  + std::to_string(n) + "_"
		  + std::to_string(d) + "_"
		  + std::to_string(k) + ".txt");
      if(fs.fail()){ exit(1); }
      fs << a << endl;
      fs << c << endl;
      fs.close();
    }

#if BENCH
    bench_show(std::cerr, method, n, d, k, t_start, t_end, rep, data, a, c);
#endif
    
    return 0;
  }
}
