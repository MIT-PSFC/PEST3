#ifndef __MV_vector__
#define __MV_vector__

#include <math.h>
#include <stdlib.h>
#include "mvvtp.h"

/******************************************************************************

A. Pletzer 

last updated Oct 23 1998

Additional Methods for MV_Vector<template> class:

scalar (double or else)  f
MV_Vector<T> v, w;

max(v)
min(v)
v + w
v + f
f + v
v - w
v - f
f - v
v*w   element by element multiplication
ScalProd(v, w)
sin(v)
cos(v)
tan(v)
asin(v)
acos(v)
atan(v)
exp(v)
log(v)
sqrt(v)
pow(v, f)
pow(v,int)
abs(v)
space(vmin,vmax,dim) generates a vector of length dim made of equally spaced
                     elements between vmin and vmax
range(imin, n)       generates an integer vector (imin, imin+1... imin+n-1)
                     of n elements
max(v, w) is a vector whose elements are the element by element max of v & w
max(v, f) vector whose elements are the max of v & f
max(f, v)
min(v, w)
min(v, f)
min(f, v)
sum(v)    sum of all elements

cat(v1, v2 ) returns a vector which is the concatenation of v1 and v2.


Prototypes:

template<class T> T max(const MV_Vector<T> &a);
template<class T> T min(const MV_Vector<T> &a);
template<class T> MV_Vector<T> operator+=(const T f);
template<class T> MV_Vector<T> operator+=(const MV_Vector<T> &a);
template<class T> MV_Vector<T> operator-=(const T f);
template<class T> MV_Vector<T> operator-=(const MV_Vector<T> &a);
//template<class T> MV_Vector<T> operator*=(const T f);
template<class T> MV_Vector<T> operator*=(const MV_Vector<T> &a);
template<class T> MV_Vector<T> operator/=(const T f);
template<class T> MV_Vector<T> operator/=(const MV_Vector<T> &a);
template<class T> MV_Vector<T> operator+(const MV_Vector<T> &a, const T f);
template<class T> MV_Vector<T> operator+(const T f, const MV_Vector<T> &a);
template<class T> MV_Vector<T> operator-(const MV_Vector<T> &a, const T f);
template<class T> MV_Vector<T> operator-(const T f, const MV_Vector<T> &a);
template<class T> MV_Vector<T> operator-(const MV_Vector<T> &a);
template<class T> T ScalProd(const MV_Vector<T> &a, const MV_Vector<T> &b);
template<class T> MV_Vector<T> operator*(const MV_Vector<T> &a, const MV_Vector<T> &b);
template<class T> MV_Vector<T> operator/(const MV_Vector<T> &a, const T f);
template<class T> MV_Vector<T> operator/(const T f, const MV_Vector<T> &a);
template<class T> MV_Vector<T> operator/(const MV_Vector<T> &a, const MV_Vector<T> &b);
MV_Vector<double> ran(const int n=0);
template<class T> MV_Vector<T> sin(const MV_Vector<T> &a);
template<class T> MV_Vector<T> cos(const MV_Vector<T> &a);
template<class T> MV_Vector<T> tan(const MV_Vector<T> &a);
template<class T> MV_Vector<T> asin(const MV_Vector<T> &a);
template<class T> MV_Vector<T> acos(const MV_Vector<T> &a);
template<class T> MV_Vector<T> atan(const MV_Vector<T> &a);
template<class T> MV_Vector<T> exp(const MV_Vector<T> &a);
template<class T> MV_Vector<T> log(const MV_Vector<T> &a);
template<class T> MV_Vector<T> sqrt(const MV_Vector<T> &a);
template<class T> MV_Vector<T> abs(const MV_Vector<T> &a);
template<class T> MV_Vector<T> pow(const MV_Vector<T> &a, const T exp);
template<class T> MV_Vector<T> pow(const MV_Vector<T> &a, const int exp);
template<class T> MV_Vector<T> space(T xmin=0, T xmax=1, unsigned int n=2);
 MV_Vector<int> range(int imin=0, int nsize=2);
template<class T> MV_Vector<T> max(const MV_Vector<T> &v1, const MV_Vector<T> &v2);
template<class T> MV_Vector<T> max(const MV_Vector<T> &v1, const T f);
template<class T> inline MV_Vector<T> max(const T f, const MV_Vector<T> &v1);
template<class T> MV_Vector<T> min(const MV_Vector<T> &v1, const MV_Vector<T> &v2);
template<class T> MV_Vector<T> min(const MV_Vector<T> &v1, const T f);
template<class T> inline MV_Vector<T> min(const T f, const MV_Vector<T> &v1);
template<class T> T sum(const MV_Vector<T> &v);
template<class T> MV_Vector<T> cat( const MV_Vector<T> &v1, const MV_Vector<T> &v2);


******************************************************************************/

//extern "C" double pow(double, double);
//extern "C" double pow(double, int);


typedef MV_Vector<double> Vec;
typedef MV_Vector<int> Vec_int;

/* 
find max element 
*/

template<class T>  T  max(const MV_Vector<T> &a) {
  T maximum=a(0);
  unsigned int n = a.size();
  for (unsigned int i=0; i<n; ++i) { if( a(i) > maximum ) maximum = a(i); }
  return maximum;
}

/* 
find min element 
*/

template<class T>  T  min(const MV_Vector<T> &a) {
  T minimum=a(0);
  unsigned int n = a.size();
  for (unsigned int i=0; i<n; ++i) { if( a(i) < minimum ) minimum = a(i); }
  return minimum;
}

/*
 +
 */

template<class T> MV_Vector<T> MV_Vector<T>::operator+=(const T f) {
  for (unsigned int j=0; j<size(); ++j) { p_[j] += T(f); }
  return *this;
}
template<class T> MV_Vector<T> MV_Vector<T>::operator+=(const MV_Vector<T> &a) {
  for (unsigned int j=0; j<size(); ++j) { p_[j] += a(j); }
  return *this;
}
  template<class T> MV_Vector<T> operator+(const MV_Vector<T> &a, const T f) {
      MV_Vector<T> c(a);
      return c += f;
    }
   template<class T> MV_Vector<T> operator+(const T f, const MV_Vector<T> &a) {
      MV_Vector<T> c(a);
      return c += f;
    }


/*    template<class T> MV_Vector<T> operator+(const MV_Vector<T> &a, const MV_Vector<T> &b) { */
/*       MV_Vector<T> c(a); */

/*       for (unsigned int j=0; j<a.size(1); ++j) */
/*      { */
/*            c(j) = a(j) + b(j); */
/*      } */
/*       return c; */
/*      } */


/*
 -
 */
     
template<class T> MV_Vector<T> MV_Vector<T>::operator-=(const T f) {
  for (unsigned int j=0; j<size(); ++j) { p_[j] -= T(f); }
  return *this;
}
template<class T> MV_Vector<T> MV_Vector<T>::operator-=(const MV_Vector<T> &a) {
  for (unsigned int j=0; j<size(); ++j) { p_[j] -= a(j); }
  return *this;
}
   template<class T> MV_Vector<T> operator-(const MV_Vector<T> &a, const T f) {
      MV_Vector<T> c(a);
      return c -= f;
    }
   template<class T> MV_Vector<T> operator-(const T f, const MV_Vector<T> &a) {
      MV_Vector<T> c(a.size(), f);
      return c -= a;
    }
   template<class T> MV_Vector<T> operator-(const MV_Vector<T> &a) {
      MV_Vector<T> c(a);
      for (unsigned int j=0; j<a.size(); ++j) { c(j) = - a(j); }
      return c;
    }


/*    template<class T> MV_Vector<T> operator-(const MV_Vector<T> &a, const MV_Vector<T> &b) { */
/*       MV_Vector<T> c(a); */

/*       for (unsigned int j=0; j<a.size(1); ++j) */
/*      { */
/*            c(j) = a(j) - b(j); */
/*      } */
/*       return c; */
/*      } */


/*
 * (multiplication)
 */
     
/*    template<class T> MV_Vector<T> operator*(const MV_Vector<T> &a, const T f) { */
/*       MV_Vector<T> c(a); */

/*       for (unsigned int i=0; i<a.size(); ++i) c(i) = a(i) * f; */
/*       return c; */
/*     } */

/*    template<class T> MV_Vector<T> operator*(const T f, const MV_Vector<T> &a) { */
/*       MV_Vector<T> c(a); */

/*       for (unsigned int i=0; i<a.size(); ++i) c(i) = a(i) * f; */
/*       return c; */
/*     } */

/*
 * Scalar product
 */

   template<class T> T ScalProd(const MV_Vector<T> &a, const MV_Vector<T> &b) {
      T c=0;

      if( a.size() != b.size() ) {
        cerr << "MV_MV_Vectoror::ScalProd incompatible sizes\n";
        return c;
      }
      for (unsigned int j=0; j<a.size(); ++j){ c += ( a(j)*b(j) ); }
      return c;
    }

/*
 element by element multiplication
 */
     
  template<class T> MV_Vector<T> MV_Vector<T>::operator*=(const T f) { 
    for (unsigned int j=0; j<size(); ++j) { p_[j] *= T(f); }  
    return *this;  
  }  
  template<class T> MV_Vector<T> MV_Vector<T>::operator*=(const MV_Vector<T> &a) {  
    for (unsigned int j=0; j<size(); ++j) { p_[j] *= a(j); }  
    return *this;  
  }  

   template<class T> MV_Vector<T> operator*(const MV_Vector<T> &a, const MV_Vector<T> &b) {
      MV_Vector<T> c(a);
      return c*=b;
    }
 

/*
 /  division
 */ 
     
template<class T> MV_Vector<T> MV_Vector<T>::operator/=(const T f) {
  for (unsigned int j=0; j<size(); ++j) { p_[j] /= T(f); }
  return *this;
}
template<class T> MV_Vector<T> MV_Vector<T>::operator/=(const MV_Vector<T> &a) {
  for (unsigned int j=0; j<size(); ++j) { p_[j] /= a(j); }
  return *this;
}
   template<class T> MV_Vector<T> operator/(const MV_Vector<T> &a, const T f) {
     MV_Vector<T> c(a);
      return c/=f;
    }

   template<class T> MV_Vector<T> operator/(const T f, const MV_Vector<T> &a) {
      MV_Vector<T> c(a);
      for (unsigned int i=0; i<a.size(); ++i) c(i) = f / a(i);
      return c;
    }
/*
 / element by element division
 */

   template<class T> MV_Vector<T> operator/(const MV_Vector<T> &a, const MV_Vector<T> &b) {
      MV_Vector<T> c(a);
      return c/=b;
    }



/* Math functions */

     template<class T> MV_Vector<T> sin(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i)  b(i) = sin( a(i) );
       return b;  
     }  

     template<class T> MV_Vector<T> cos(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i) b(i) = cos( a(i) );
       return b;  
     }  

     template<class T> MV_Vector<T> tan(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(0); ++i) b(i) = tan( a(i) );
       return b;  
     }  


     template<class T> MV_Vector<T> asin(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i) b(i) = asin( a(i) );
       return b;  
     }  

     template<class T> MV_Vector<T> acos(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i) b(i) = acos( a(i) );
       return b;  
     }  

     template<class T> MV_Vector<T> atan(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i) b(i) = atan( a(i) );
       return b;  
     }  

     template<class T> MV_Vector<T> exp(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i) b(i) = exp( a(i) );
       return b;  
     }  

     template<class T> MV_Vector<T> log(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i) b(i) = log( a(i) );
       return b;  
     }  

     template<class T> MV_Vector<T> sqrt(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i) b(i) = sqrt( a(i) );
       return b;  
     }
  
     template<class T> MV_Vector<T> abs(const MV_Vector<T> &a) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i) b(i) = fabs( a(i) );
       return b;  
     }  

     template<class T> MV_Vector<T> pow(const MV_Vector<T> &a, const T exp) { 
       MV_Vector<T> b(a); 
 
       for (unsigned int i=0; i<a.size(); ++i) 
	 b(i) = (T)pow( (double)a(i) , (double)exp );
       return b;  
     }  

     template<class T> MV_Vector<T> pow(const MV_Vector<T> &a, const int exp) { 
       MV_Vector<T> b(a);  

       for (unsigned int i=0; i<a.size(); ++i) 
	 b(i) = (T)pow( (double)a(i) , (double)exp );
       return b;  
     }  

/*
grid generator
*/

    
  template<class T> MV_Vector<T> space(T xmin, T xmax, unsigned int n=2) {
      MV_Vector<T> a(n);

      for (unsigned int i=0; i<a.size(); ++i) {
        a(i) = xmin + (xmax-xmin)*T(i)/T(a.size()-1);
      }
      return a;

    }

    
  template<class T> MV_Vector<T> range(T imin, unsigned int nsize=2) {
      MV_Vector<T> a(nsize);
      for (unsigned int i=0; i<nsize; ++i) a(i) = imin + T(i);
      return a;
    }
  template<class T> MV_Vector<T> range(T imin, int nsize=2) {
      MV_Vector<T> a(nsize);
      for (unsigned int i=0; i<nsize; ++i) a(i) = imin + T(i);
      return a;
    }

/* 
max of 2 arguments
*/

template<class T> MV_Vector<T>  max(const MV_Vector<T> &v1, const MV_Vector<T> &v2) {
  unsigned int n = v1.size();
  MV_Vector<T> res(n);
  for (unsigned int i=0; i<n; ++i) { 
    res(i) = v1(i);
    if( v2(i) > v1(i) ) res(i) = v2(i); 
  }
  return res;
}

template<class T> MV_Vector<T> max(const MV_Vector<T> &v1, const T f) {
  unsigned int n = v1.size();
  MV_Vector<T> res(n);
  for (unsigned int i=0; i<n; ++i) { 
    res(i) = v1(i);
    if( f > v1(i) ) res(i) = f; 
  }
  return res;
}

template<class T> inline MV_Vector<T>  max(const T f, const MV_Vector<T> &v1) {
  return max(v1, f);
}

/*
min of 2 arguments
*/

template<class T> MV_Vector<T>  min(const MV_Vector<T> &v1, const MV_Vector<T> &v2) {
  unsigned int n = v1.size();
  if( n > v2.size() ) {
    cerr << 
      "MV_Vector<T>::min(const MV_Vector<T> &, const MV_Vector<T> &) ERROR size() \n";
    exit(1);
  }
  MV_Vector<T> res(n);
  for (unsigned int i=0; i<n; ++i) { 
    res(i) = v1(i);
    if( v2(i) < v1(i) ) res(i) = v2(i); 
  }
  return res;
}

template<class T> MV_Vector<T>  min(const MV_Vector<T> &v1, const T f) {
  unsigned int n = v1.size();
  MV_Vector<T> res(n);
  for (unsigned int i=0; i<n; ++i) { 
    res(i) = v1(i);
    if( f < v1(i) ) res(i) = f; 
  }
  return res;
}

template<class T> inline MV_Vector<T>  min(const T f, const MV_Vector<T> &v1) {
  return min(v1, f);
}

/*
sum of elements
*/

template<class T> T  sum(const MV_Vector<T> &v) {
  unsigned int n = v.size();
  T res = 0;
  for (unsigned int i=0; i<n; ++i) res += v(i);
  return res;
}

/*
concatenate two vectors
*/

template<class T> MV_Vector<T>  cat( const MV_Vector<T> &v1, const MV_Vector<T> &v2) {
  unsigned int n1 = v1.size(); 
  unsigned int n2 = v2.size();
  unsigned int nres = n1 + n2;
  MV_Vector<T> res( nres );

  for (unsigned int i=0; i< n1; ++i) res(i   ) = v1(i);
  for (unsigned int j=0; j< n2; ++j) res(j+n1) = v2(j);

  return res;
}


/*
class methods to fill up elements of an existing vector
*/

    
  template<class TYPE> void MV_Vector<TYPE>::space(TYPE xmin, TYPE xmax) {
      for (unsigned int i=0; i<size(); ++i) {
        p_[i] = xmin + (xmax-xmin)*TYPE(i)/TYPE(size()-1);
      }
    }
    
  template<class TYPE> void MV_Vector<TYPE>::range(TYPE imin) {
      for (unsigned int i=0; i< (*this).size(); ++i) p_[i] = TYPE(i)+imin;
    }

template<class TYPE>
unsigned int index_min(MV_Vector<TYPE> &v){

   /* return index of min(v) */  

  unsigned int ires = 0;
  unsigned int n = v.size();  
  TYPE vmin = v(0);  
  for (unsigned int i=0; i<n; i++) {    
    if (v(i) < vmin) {
      ires = i;
      vmin = v(i);
    }
  }
  return ires;
}

template<class TYPE>
unsigned int index_max(MV_Vector<TYPE> &v){

  /* return index of max(v) */

  unsigned int ires = 0;
  unsigned int n = v.size();
  TYPE vmax = v(0);
  for (unsigned int i=0; i<n; i++) {
    if (v(i) > vmax) {
      ires = i;
      vmax = v(i);
    }
  }
  return ires;
}
 


#endif /* __MV_vector__ */
