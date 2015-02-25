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
#include "Lloyd.hpp"
#include "Hamerly.hpp"

using namespace std;
using namespace boost::numeric::ublas;
using namespace boost::numeric;

enum algorithm {Lloyd, Hamerly}; /* 実行するアルゴリズム */
static const std::vector<string> algo_name{"Lloyd", "Hamerly"};

#define BENCH 1
/* k-meansの各種のアルゴリズムのベンチマークを取るときに1*/

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
	 << " <method> <n> <d> <k> <data>\n"
	 << " <method> = Lloyd | Hamerly\n"
	 << " n : number of data points\n"
      	 << " d : dimension\n"
	 << " k : number of clusters\n"
         << " data    : file path of the data file\n"
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

    /* 初期クラスタ中心を何らかの方法で得る */
    initialize_centers(k, data, c);

    
#if BENCH /* プログラム実行時間のベンチマークを取る */
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL); /* 時間計測開始 */
#endif
    
    if(method == Lloyd){
      rep = Lloyd_main(n, d, k, data, a, c);
    }else if(method == Hamerly){
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
