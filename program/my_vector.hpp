#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

/* ベクトルに対する加算・減算・定数倍などを定義 */
template <class T>
std::vector<T>& operator+=(std::vector<T> &self,
			   const std::vector<T> &other){
  for(int i = 0; i < (int)self.size(); i++)
    self[i] += other[i];
  return self;
}

template <class T>
std::vector<T> operator+(const std::vector<T> &self,
                         const std::vector<T> &other){
  std::vector<T> result = self;
  result += other;
  return result;
}

template <class T>
std::vector<T>& operator-=(std::vector<T> &self,
			   const std::vector<T> &other){
  for(int i = 0; i < (int)self.size(); i++)
    self[i] -= other[i];
  return self;
}

template <class T>
std::vector<T> operator-(const std::vector<T> &self,
                         const std::vector<T> &other){
  std::vector<T> result = self;
  result -= other;
  return result;
}

template <class T>
std::vector<T>& operator*=(std::vector<T> &self, const T &mul){			   
  for(int i = 0; i < (int)self.size(); i++)
    self[i] *= mul;
  return self;
}

template <class T>
std::vector<T> operator*(const std::vector<T> &self,
                         const T &mul){
  std::vector<T> result = self;
  result *= mul;
  return result;
}

template <class T>
std::vector<T> operator*(const T &mul, const std::vector<T> &self){                         
  std::vector<T> result = self;
  result *= mul;
  return result;
}

template <class T>
ostream &operator<<(ostream &stream, vector<T> vector){
  for(int i = 0; i < (int)vector.size() - 1; i++){
    stream << vector.at(i) << ", "; 
  }
  stream << vector.back() << endl;
  return stream;
}

template <class T>
ostream &operator<<(ostream &stream, vector< vector<T> > matrix){
  for(int i = 0; i < (int)matrix.size(); i++){ stream << matrix.at(i); }
  return stream;
}

template <class T>
inline double Euclid_norm(vector<T> vec){
  /* compute Euclid norm of a vector 'vec' */
  double results = 0;
  for(int i = 0; i < (int)vec.size(); ++i){
    results += vec.at(i) * vec.at(i);
  }
  return sqrt(results);
}

template <class T>
inline double sum(vector<T> vec){
  double results = 0;
  for(int i = 0; i < (int)vec.size(); ++i){
    results += vec.at(i);
  }
  return results;
}

template <class T>
double norm(vector<T> vec, int p = 2){
  /* compute norm of a vector 'vec' */
  if(p == 2){
    return Euclid_norm(vec);
  }else if(p == 1){
    return sum(vec);
  }else{
    double results = 0;
    for(int i = 0; i < (int)vec.size(); ++i){
      results += pow(vec.at(i), p);
    }
    return pow(results, 1.0 / p);
  }
}

template <class T>
double inner_product(vector<T> v1, vector<T> v2){
  /* compute inner product of two vectors v1, v2 */
  double results = 0;
  int minsize = (((int)v1.size() < (int)v2.size())
		 ? (int)v1.size() : (int)v2.size());
  for(int i = 0; i < minsize; ++i){
    results += v1.at(i) * v2.at(i);
  }
  return results;
}
