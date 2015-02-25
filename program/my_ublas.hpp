#ifndef _UBLAS_
#define _UBLAS_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#endif

#ifndef _MY_UBLAS_
#define _MY_UBLAS_

using namespace boost::numeric::ublas;
using namespace boost::numeric;

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



#endif
