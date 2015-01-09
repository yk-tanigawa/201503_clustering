#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#include "my_vector.hpp"
#include "time_bench.hpp"

using namespace std;

void re_clustering(const int n, const int d, const int k,
		   vector<int> &cluster,
		   const vector<vector <double> > representative,
		   const vector<vector <double> > data){
  for(int i = 0; i < n; ++i){
    int closest = 0;
    double dist = norm(data.at(i) - representative.at(closest));
    double min_dist = dist;
    for(int z = 1; z < k; ++z){
      dist = norm(data.at(i) - representative.at(z));
      if(dist < min_dist){
	min_dist = dist;
	closest = z;
      }
    }
    cluster.at(i) = closest;
  }
  return;
}

void calc_representative(const int n, const int d, const int k,
			 const vector<int> cluster,
			 vector<vector <double> > &representative,
			 const vector<vector <double> > data){
  vector<int> count(k, 0); /* 各クラスタに何点あるか数える */
  vector<vector <double> > sum(k); /* 各クラスタの座標の和 */
  {
    vector<double> zero(d, 0);
    for(int z = 0; z < k; ++z){
      sum.at(z) = zero;
    }
  }

  /* それぞれのクラスタについて全データ点の座標の和を計算 */
  for(int i = 0; i < n; ++i){
    count.at(cluster.at(i))++;
    sum.at(cluster.at(i)) += data.at(i);  
  }
  
  /* 代表点を更新する */
  for(int z = 0; z < k; ++z){
    representative.at(z) = (1.0 / count.at(z)) * sum.at(z);
  }
  
  return;
}

double average_square_err(const int n, const int d, const int k,
			 const vector<int> cluster,
			 const vector<vector<double> > representative,
			 const vector<vector<double> > data){
  double sum = 0;
  for(int i = 0; i < n; ++i){
    sum += norm(data.at(i) - representative.at(cluster.at(i))) * norm(data.at(i) - representative.at(cluster.at(i)));
  }
  return (sum / n);
}

double Lloyd_repeat(const int n, const int d, const int k,
		    vector<int> &cluster,
		    vector<vector <double> > &representative,
		    const vector<vector <double> > data){
  re_clustering(n, d, k, cluster, representative, data);
  calc_representative(n, d, k, cluster, representative, data);
  
  return average_square_err(n, d, k, cluster, representative, data);
}

void Lloyd(const int n, const int d, const int k,
	   const vector<vector <double> > data){
  if(n < k){
    cerr << "data size n is smaller than cluster size d" << endl;
    exit(1);
  }else{

    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL); /* 時間計測開始 */

    vector<int> cluster(n);
    /* 各点x_iがどのクラスタに属しているか */
    vector<vector <double> > representative(k);
    /* 各クラスタの代表点 */
    for(int z = 0; z < k; ++z){ /* ランダムに初期化する */
      representative.at(z) = data.at(z);
    }

    double sq_err, sq_err_before = Lloyd_repeat(n, d, k, cluster, representative, data);
    int t = 1; /* 繰り返しを何ステップ行ったか */
    while(1){
      sq_err = Lloyd_repeat(n, d, k, cluster, representative, data);
      if(sq_err == sq_err_before){ break; }
      t++;
      sq_err_before = sq_err;
    }
    
    gettimeofday(&t_end, NULL); /* 時間計測終了 */

    cout << n << "\t"                            /* サンプル数 */
	 << d << "\t"                            /* 次元 */
      	 << k << "\t"                            /* クラスタ数 */
	 << diff_timeval(t_start, t_end) << "\t" /* 実行時間(us) */
	 << sq_err << "\t"                       /* 平均二乗誤差 */
	 << t << endl;                           /* 繰り返し回数*/
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
    Lloyd(n, d, k, data);
    return 0;
  }
}
