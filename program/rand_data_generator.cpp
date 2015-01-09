#include <iostream>
#include <cstdlib>
#include <fstream>
#include <random>

using namespace std;

/* define constant values */
static const double DATA_RANGE_MIN = 0.0;
static const double DATA_RANGE_MAX = 1.0;

#define DEBUG 0

int main(int argc, char *argv[]){
  if(argc < 4){
    cerr << "usage: " << argv[0]
	 << " <n: data num> <d: dimension> <file name> [seed]" << endl;
    exit(1);
  }else{
    int n = atoi(argv[1]), d = atoi(argv[2]), seed;
    char *file = argv[3];
    {
      std::random_device rd;
      if(argc < 5){
#if DEBUG
	seed = 0;
#else
	seed = rd();
#endif
      }else{
	seed = atoi(argv[4]);
      }
    }
    cout << "seed = " << seed << endl;
    {
      ofstream fs(file);
      if(fs.fail()){ exit(1); }
      {
	std::mt19937 mt(seed);
	std::uniform_real_distribution<double> rand(DATA_RANGE_MIN, DATA_RANGE_MAX);
	
	for(int i = 0; i < n; ++i){
	  for(int k = 0; k < d - 1; ++k){
	    fs << rand(mt) << " ";
	  }
	  fs << rand(mt) << endl;
	}
      }
      fs.close();
    }
    return 0;
  }
}
