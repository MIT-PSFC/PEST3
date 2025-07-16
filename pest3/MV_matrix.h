#ifndef __MV_matrix__
#define __MV_matrix__


#include <math.h>
//#include <stdlib.h>
#include <fstream.h>
#include <iostream.h>
//#include <string>        // <string.h> or <string>
#include "mvvtp.h"
#include "mvmtp.h"
#include "mvvind.h"

typedef MV_ColMat<double> Mat;

/**********************************************************************

 Vector/Matrix class:


MV_ColMat<T> a;

creation

ones(n1,n2)    n1 x n2 matrix of ones
ran(n1,n2)     random matrix with 0 < elements < 1 (default 2x2)
space1(xmin=0, xmax=1, ny=1, nx=11) 
                 create a grid matrix with row points between xmin and xmax.
space1(ymin=0, ymax=1, ny=11, nx=1)
                 create a grid matrix with column points between ymin and ymax.
Eye(n)           create diagonal matrix n x n 
                 

attributes

a.size(0)             number of rows
a.size(1)             number of columns

transformations

a(j1,j2)         element (j1,j2), 0 based indexing

min(a)
max(a)


operations 

a + double
a + a
a - double
double - a
a - a
a * f
f * a
a * a          element by element multiplication
a / f
f / a
a / a          element by element division
MatMult(a,b)   
MatMult(a,b,c) Matrix multiplication

sin(a)
cos(a)
tan(a)
asin(a)
acos(a)
atan(a)
log(a)
exp(a)
sqrt(a)
abs(a)
trans(a)        transpose



a.plotx()        plot against columns
a.ploty()        plot against rows
a.plotxy()       3-d plot of matrix


save/load results

>>
a = load("file")

stretch0( v, n ) converts a vector v to a matrix of dim n times v.size() by
                 copying the row vector n times.
stretch1( v, n ) converts a vector v to a matrix of dim v.size() times n by
                 copying the column vector n times.


Prototypes:
template<class T> MV_ColMat<T> trans(const MV_ColMat<T> &a);
 
template<class T>
 T max(const MV_ColMat<T> &a);
 
template<class T>
 MV_ColMat<T> max(const MV_ColMat<T> &a, const T val);
 
template<class T>
  MV_ColMat<T> max(const MV_ColMat<T> &a, const MV_ColMat<T> &b);

template<class T>
 MV_ColMat<T> max( const T val, const MV_ColMat<T> &a);
 
template<class T>
 T min(const MV_ColMat<T> &a);
 
template<class T>
 MV_ColMat<T> min(const MV_ColMat<T> &a, const T val);
 
template<class T>
 MV_ColMat<T> min(const T val, const MV_ColMat<T> &a);
 
template<class T>
  MV_ColMat<T> min(const MV_ColMat<T> &a, const MV_ColMat<T> &b);

template<class T>
 T average(const MV_ColMat<T> &a);
 
template<class T>
 MV_Vector<T> average0(const MV_ColMat<T> &a);

template<class T>
 MV_Vector<T> average1(const MV_ColMat<T> &a);

template<class T>
 MV_Vector<T> sum0(const MV_ColMat<T> &a);

template<class T>
 MV_Vector<T> sum1(const MV_ColMat<T> &a);

 template<class T> 
 MV_ColMat<T> operator+=(const T f);

 template<class T> 
 MV_ColMat<T> operator+=(const MV_ColMat<T> &a);

 template<class T> 
 MV_ColMat<T> operator-=(const T f);

 template<class T> 
 MV_ColMat<T> operator-=(const MV_ColMat<T> &a);

 template<class T> 
 MV_ColMat<T> operator*=(const T f);

 template<class T> 
 MV_ColMat<T> operator*=(const MV_ColMat<T> &a);

 template<class T> 
 MV_ColMat<T> operator/=(const T f);

 template<class T> 
 MV_ColMat<T> operator/=(const MV_ColMat<T> &a);

 template<class T> 
 MV_ColMat<T> operator+(const MV_ColMat<T> &a, const T f);
 
 template<class T> 
 MV_ColMat<T> operator+(const T f, const MV_ColMat<T> &a);
 
 template<class T> 
 MV_ColMat<T> operator+(const MV_ColMat<T> &a, const MV_ColMat<T> &b);
 
 template<class T> 
 MV_ColMat<T> operator-(const MV_ColMat<T> &a);
 
 template<class T> 
 MV_ColMat<T> operator-(const MV_ColMat<T> &a, const T f);
 
 template<class T> MV_ColMat<T> operator-(const T f, const MV_ColMat<T> &a);
 
 template<class T> 
 MV_ColMat<T> operator-(const MV_ColMat<T> &a, const MV_ColMat<T> &b);
 
 template<class T> 
 MV_ColMat<T> operator*(const MV_ColMat<T> &a, const T f);
 
 template<class T> 
 MV_ColMat<T> operator*(const T f, const MV_ColMat<T> &a);
 
 template<class T> 
 MV_ColMat<T> MatMult(const MV_ColMat<T> &a, const MV_ColMat<T> &b);

 template<class T> 
 MV_Vector<T> MatMult(const MV_ColMat<T> &a, const MV_Vector<T> &b);

 template<class T> 
 MV_Vector<T> ScalProd(const MV_ColMat<T> &a, const MV_Vector<T> &b);

 template<class T> 
 MV_Vector<T> ScalProd(const MV_Vector<T> &b, const MV_ColMat<T> &a);
 
 template<class T> 
 MV_ColMat<T> TripleMatMult(const MV_ColMat<T> &a, const MV_ColMat<T> &b, 
 const MV_ColMat<T> &d);
 

 template<class T> 
 MV_ColMat<T> MatMult(const MV_ColMat<T> &a, const MV_ColMat<T> &b, const MV_ColMat<T> &c);
 
 template<class T> 
 MV_ColMat<T> operator*(const MV_ColMat<T> &a, const MV_ColMat<T> &b);
 
 template<class T> 
 MV_ColMat<T> operator/(const MV_ColMat<T> &a, const T f);
 
 template<class T> 
 MV_ColMat<T> operator/(const T f, const MV_ColMat<T> &a);
 
 template<class T>
 MV_ColMat<T> operator/(const MV_ColMat<T> &a, const MV_ColMat<T> &b);
 
 MV_ColMat<double> 
 Eye(unsigned int n);
 template<class T> 
 MV_ColMat<double> ran(const unsigned int n1, const unsigned int n2);
 
 template<class T> MV_ColMat<T> sin(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> cos(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> tan(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> asin(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> acos(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> atan(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> exp(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> log(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> sqrt(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> abs(const MV_ColMat<T> &a);
  
 template<class T> MV_ColMat<T> pow(const MV_ColMat<T> &a, const T exp);
  
 template<class T> MV_ColMat<T> pow(const MV_ColMat<T> &a, const int exp);
  
 template<class T> 
 MV_ColMat<T> space1(T xmin, T xmax, unsigned int ny=1, unsigned int nx=11);
 
 template<class T> 
 MV_ColMat<T> space0(T ymin, T ymax, unsigned int ny=11, unsigned int nx=1);
 
  
MV_ColMat<double> load(string cfile);
template<class T> MV_ColMat<T>
 stretch0(const MV_Vector<T> &v, const unsigned int n0);

template<class T> MV_ColMat<T>
 stretch1(const MV_Vector<T> &v, const unsigned int n1);

template<class TYPE>
 MV_ColMat<TYPE> TensorProd(const MV_Vector<TYPE> &a, const MV_Vector<TYPE> &b);

template<class TYPE>
 MV_Vector<TYPE> flatten0(const MV_ColMat<TYPE> &a);

template<class TYPE>
 MV_Vector<TYPE> flatten1(const MV_ColMat<TYPE> &a);

**********************************************************************/
//extern "C" double pow(double, double);
//extern "C" double pow(double, int);

template<class T> 
MV_ColMat<T>  trans(const MV_ColMat<T> &a)
    {
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);
      MV_ColMat<T> b(nc,nr);

      for (unsigned int i=0; i<nr; ++i)
	{
	  for (unsigned int j=0; j<nc; ++j) b(j,i) = a(i,j);
	}
      return b;
    }

 
template<class T>
  T      max(const MV_ColMat<T> &a)
    {
      T maximum = a(0,0);
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      if ( a(i,j) > maximum ) maximum = a(i,j);
	    }
	}

      return maximum;
    }

 
template<class T>
  MV_ColMat<T> max(const MV_ColMat<T> &a, const MV_ColMat<T> &b)
    {
      MV_ColMat<T> c(a);
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      if ( b(i,j) > a(i,j) ) c(i,j) = b(i,j);
	    }
	}

      return c;
    }

template<class T>
   MV_ColMat<T> max(const MV_ColMat<T> &a, const T val)
    {
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);
      MV_ColMat<T> b(a);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      if ( a(i,j) < val ) b(i,j) = val;
	    }
	}

      return b;
    }

template<class T>
  MV_ColMat<T> max( const T val, const MV_ColMat<T> &a)
    {
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);
      MV_ColMat<T> b(a);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      if ( a(i,j) < val ) b(i,j) = val;
	    }
	}

      return b;
    }
    

template<class T>
  T      min(const MV_ColMat<T> &a)
    {
      T minimum = a(0,0);
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      if ( a(i,j) < minimum ) minimum = a(i,j);
	    }
	}

      return minimum;
    }


template<class T>
   MV_ColMat<T> min(const MV_ColMat<T> &a, const T val)
    {
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);
      MV_ColMat<T> b(a);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      if ( a(i,j) > val ) b(i,j) = val;
	    }
	}

      return b;
    }


template<class T>
   MV_ColMat<T> min(const T val, const MV_ColMat<T> &a)
    {
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);
      MV_ColMat<T> b(a);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      if ( a(i,j) > val ) b(i,j) = val;
	    }
	}

      return b;
    }

 
template<class T>
  MV_ColMat<T> min(const MV_ColMat<T> &a, const MV_ColMat<T> &b)
    {
      MV_ColMat<T> c(a);
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      if ( b(i,j) < a(i,j) ) c(i,j) = b(i,j);
	    }
	}

      return c;
    }

template<class T>
  T      average(const MV_ColMat<T> &a)
    {
      T av = a(0,0);
      unsigned int nr = a.size(0);
      unsigned int nc = a.size(1);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      av += a(i,j);
	    }
	}

      return av/(nr*nc);
    }
    

template<class T>
  MV_Vector<T>      average0(const MV_ColMat<T> &a)
    {
      /* average along rows */

      unsigned int nr = a.size(0), nc = a.size(1);
      MV_Vector<T> av(nc,0);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      av(j) += a(i,j);
	    }
	}
      return av/T(nr);
    }


template<class T>
  MV_Vector<T>      average1(const MV_ColMat<T> &a)
    {
      /* average along columns */

      unsigned int nr = a.size(0), nc = a.size(1);
      MV_Vector<T> av(nr,0);

      for (unsigned int i=0; i<nr; ++i)
	{
	  for (unsigned int j=0; j<nc; ++j)
	    {
	      av(i) += a(i,j);
	    }
	}
      return av/T(nc);
    }


template<class T>
  MV_Vector<T>      sum0(const MV_ColMat<T> &a)
    {
      /* sum row elements (vertically) */

      unsigned int nr = a.size(0), nc = a.size(1);
      MV_Vector<T> av(nc,0);

      for (unsigned int j=0; j<nc; ++j)
	{
	  for (unsigned int i=0; i<nr; ++i)
	    {
	      av(j) += a(i,j);
	    }
	}
      return av;
    }


template<class T>
  MV_Vector<T>  sum1(const MV_ColMat<T> &a)
    {
      /* sum column elements (horizontally) */

      unsigned int nr = a.size(0), nc = a.size(1);
      MV_Vector<T> av(nr,0);

      for (unsigned int i=0; i<nr; ++i)
	{
	  for (unsigned int j=0; j<nc; ++j)
	    {
	      av(i) += a(i,j);
	    }
	}
      return av;
    }
/*
 +
 */

   template<class T> 
   MV_ColMat<T> MV_ColMat<T>::operator+=(const T f) {
     for (unsigned int j=0; j<size(1)*size(0); ++j) {
       v_[j] += T(f);
     }
   return *this;
   }
   template<class T> 
   MV_ColMat<T> MV_ColMat<T>::operator+=(const MV_ColMat<T> &a) {
     for (unsigned int j=0; j<size(1)*size(0); ++j) {
       v_[j] += a.v_[j];
     }
   return *this;
   }
       
     
   template<class T> 
   MV_ColMat<T> operator+(const MV_ColMat<T> &a, const T f)
    {
      MV_ColMat<T> c(a);
      return c+=f;
    }

   template<class T> 
   MV_ColMat<T> operator+(const T f, const MV_ColMat<T> &a)
    {
      MV_ColMat<T> c(a);
      return c+=f;
    }

   template<class T> 
   MV_ColMat<T> operator+(const MV_ColMat<T> &a, const MV_ColMat<T> &b)
     {
      MV_ColMat<T> c(a);
      return c+=b;
     }

/*
 -
 */
     
   template<class T> 
   MV_ColMat<T> MV_ColMat<T>::operator-=(const T f) {
     for (unsigned int j=0; j<size(1)*size(0); ++j) {
       v_[j] -= T(f);
     }
   return *this;
   }
   template<class T> 
   MV_ColMat<T> MV_ColMat<T>::operator-=(const MV_ColMat<T> &a) {
     for (unsigned int j=0; j<size(1)*size(0); ++j) {
       v_[j] -= a.v_[j];
     }
   return *this;
   }

   template<class T> 
   MV_ColMat<T> operator-(const MV_ColMat<T> &a)
    {
      MV_ColMat<T> c(a);
      for (unsigned int j=0; j<a.size(1); ++j)
	{
	  for (unsigned int i=0; i<a.size(0); ++i)
	    {
	      c(i,j) = -a(i,j);
	    }
	}
      return c;
    }

   template<class T> 
   MV_ColMat<T> operator-(const MV_ColMat<T> &a, const T f)
    {
      MV_ColMat<T> c(a);
      return c-=f;
    }
   template<class T> 
   MV_ColMat<T> operator-(const T f, const MV_ColMat<T> &a)
    {
      MV_ColMat<T> c(a.size(0), a.size(1), f);
      return c-=a;
    }


   template<class T> 
   MV_ColMat<T> operator-(const MV_ColMat<T> &a, const MV_ColMat<T> &b)
     {
      MV_ColMat<T> c(a);
      return c-=b;
     }


/*
 * (multiplication)
 */
     
   template<class T> 
   MV_ColMat<T> MV_ColMat<T>::operator*=(const T f) {
     for (unsigned int j=0; j<size(1)*size(0); ++j) {
       v_[j] *= T(f);
     }
   return *this;
   }
   template<class T> 
   MV_ColMat<T> MV_ColMat<T>::operator*=(const MV_ColMat<T> &a) {
     for (unsigned int j=0; j<size(1)*size(0); ++j) {
       v_[j] *= a.v_[j];
     }
   return *this;
   }

   template<class T> 
   MV_ColMat<T> operator*(const MV_ColMat<T> &a, const T f)
    {
      MV_ColMat<T> c(a);
      return c*=f;
    }
   template<class T> 
   MV_ColMat<T> operator*(const T f, const MV_ColMat<T> &a)
    {
      MV_ColMat<T> c(a);
      return c*=f;
    }

   template<class T> 
   MV_ColMat<T> MatMult(const MV_ColMat<T> &a, const MV_ColMat<T> &b)
    {
      unsigned int i,j,k;
      MV_ColMat<T> c(a.size(0),b.size(1), 0.0);

      for (j=0; j<c.size(1); ++j)
	{
	  for (i=0; i<c.size(0); ++i)
	    {
	      for (k=0; k<a.size(1); k++)
		{
		  c(i,j) += a(i,k) * b(k,j);
		}
	    }
	}
      return c;
    }

   template<class T> 
   MV_Vector<T> MatMult(const MV_ColMat<T> &a, const MV_Vector<T> &b)
    {
      return ScalProd(a,b);
    }

   template<class T> 
   MV_Vector<T> ScalProd(const MV_ColMat<T> &a, const MV_Vector<T> &b)
    {
      unsigned int i,j;
      MV_Vector<T> c(a.size(0),0);

	  for (i=0; i<c.size(); ++i)
	    {
	      for (j=0; j<a.size(1); j++)
		{
		  c(i) += a(i,j) * b(j);
		}
	    }
      return c;
    }

   template<class T> 
   MV_Vector<T> ScalProd(const MV_Vector<T> &b, const MV_ColMat<T> &a)
    {
      unsigned int i,j;
      MV_Vector<T> c(a.size(1),0);

	  for (j=0; j<c.size(); ++j)
	    {
	      for (i=0; i<a.size(0); i++)
		{
		  c(j) += b(i)*a(i,j);
		}
	    }
      return c;
    }


   template<class T> 
   MV_ColMat<T> TripleMatMult(const MV_ColMat<T> &a, const MV_ColMat<T> &b, 
                              const MV_ColMat<T> &d)
    {
      unsigned int i,j,k,m;
      unsigned int n1 = a.size(0), n2 = b.size(0), n3 = d.size(0), n4 = d.size(1);

      MV_ColMat<T> c(n1,n4, 0.0);

      for (j=0; j<n4; ++j)
	{
	  for (i=0; i<n1; ++i)
	    {
	      for (k=0; k<n2; k++)
		{
		  for (m=0; m<n3; m++) c(i,j) = c(i,j) + a(i,k) * b(k,m) * d(m,j);
		}
	    }
	}
      return c;
    }

   template<class T> 
   MV_ColMat<T> MatMult(const MV_ColMat<T> &a, const MV_ColMat<T> &b, const MV_ColMat<T> &c)
    {
      MV_ColMat<T> d = TripleMatMult( a, b, c );
      return d;
    }

/*
 multiplication element by element
 */
     
   template<class T> 
   MV_ColMat<T> operator*(const MV_ColMat<T> &a, const MV_ColMat<T> &b)
    {
      MV_ColMat<T> c(a);
      return c*=b;
    }
 

/*
 /  division
 */
     
   template<class T> 
   MV_ColMat<T> operator/(const MV_ColMat<T> &a, const T f)
    {
      T g=1/f;
      return a * g;
    }

   template<class T> 
   MV_ColMat<T> operator/(const T f, const MV_ColMat<T> &a)
    {
      MV_ColMat<T> c(a.size(0), a.size(1), f);
      return c/=a;
    }
/*
 / element by element division
 */

   template<class T> 
   MV_ColMat<T> MV_ColMat<T>::operator/=(const T f) {
     for (unsigned int j=0; j<size(1)*size(0); ++j) {
       v_[j] /= T(f);
     }
   return *this;
   }
   template<class T> 
   MV_ColMat<T> MV_ColMat<T>::operator/=(const MV_ColMat<T> &a) {
     for (unsigned int j=0; j<size(1)*size(0); ++j) {
       v_[j] /= a.v_[j];
     }
   return *this;
   }

   template<class T>
   MV_ColMat<T> operator/(const MV_ColMat<T> &a, const MV_ColMat<T> &b)
    {
      MV_ColMat<T> c(a);
      return c/=b;
    }

     template<class T> MV_ColMat<T> sin(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)  
         { 
	    for (i=0; i<a.size(0); ++i)  b(i,j) = sin( a(i,j) );
         }  
       return b;  
     }  

     template<class T> MV_ColMat<T> cos(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = cos( a(i,j) );
         }  
       return b;  
     }  

     template<class T> MV_ColMat<T> tan(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = tan( a(i,j) );
         }  
       return b;  
     }  


     template<class T> MV_ColMat<T> asin(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = asin( a(i,j) );
         }  
       return b;  
     }  

     template<class T> MV_ColMat<T> acos(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = acos( a(i,j) );
         }  
       return b;  
     }  

     template<class T> MV_ColMat<T> atan(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = atan( a(i,j) );
         }  
       return b;  
     }  


     template<class T> MV_ColMat<T> exp(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = exp( a(i,j) );
         }  
       return b;  
     }  

     template<class T> MV_ColMat<T> log(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = log( a(i,j) );
         }  
       return b;  
     }  

     template<class T> MV_ColMat<T> sqrt(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = sqrt( a(i,j) );
         }  
       return b;  
     }  
     template<class T> MV_ColMat<T> abs(const MV_ColMat<T> &a)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = fabs( a(i,j) );
         }  
       return b;  
     }  
     template<class T> MV_ColMat<T> pow(const MV_ColMat<T> &a, const T exp)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = (T)pow( (double)a(i,j) , (double)exp );
         }  
       return b;  
     }  

     template<class T> MV_ColMat<T> pow(const MV_ColMat<T> &a, const int exp)  
     { 
       unsigned int i, j;
       MV_ColMat<T> b(a);  

       for (j=0; j<a.size(1); ++j)   
         { 
	   for (i=0; i<a.size(0); ++i) b(i,j) = (T)pow( (double)a(i,j) , (double)exp );
         }  
       return b;  
     }  

    
  template<class T> 
  MV_ColMat<T> space1(T xmin, T xmax, unsigned int ny=1, unsigned int nx=11) 
    {
      unsigned int i,j;
      MV_ColMat<T> a(ny,nx);

      for (i=0; i<a.size(0); ++i) 
	{
	  for (j=0; j<a.size(1); ++j)
	    {
	      a(i,j) = xmin + (xmax-xmin)*T(j)/T(a.size(1)-1);
	    }
	}
      return a;

    }

  template<class T> 
  MV_ColMat<T> space0(T ymin, T ymax, unsigned int ny=11, unsigned int nx=1) 
    {
      unsigned int i,j;
      MV_ColMat<T> a(ny,nx);

      for (i=0; i<a.size(0); ++i) 
	{
	  for (j=0; j<a.size(1); ++j)
	    {
	      a(i,j) = ymin + (ymax-ymin)*T(i)/T(a.size(0)-1);
	    }
	}
      return a;

    }


/* stretch vertically */

template<class T> MV_ColMat<T>
 stretch0(const  MV_Vector<T> &v, const unsigned int n0) {  
    unsigned int nv = v.size();
    MV_ColMat<T> d(n0, nv); 
    for(unsigned int i=0; i<n0; ++i) {
      for (unsigned int j=0; j<nv; ++j) {
	d(i,j) = v(j);
      }
    }
    return d; 
 }  
  
/* stretch horizontally */

  
template<class T> MV_ColMat<T>
  stretch1(const MV_Vector<T> &v, const unsigned int n1) {  
    unsigned int nv = v.size();
    MV_ColMat<T> d(nv, n1); 
    for(unsigned int j=0; j<n1; ++j) {
      for (unsigned int i=0; i<nv; ++i) {
	d(i,j) = v(i);
      }
    }
    return d; 
 }  


/* tensor product of 2 vectors = matrix */

   template<class TYPE>
   MV_ColMat<TYPE> TensorProd(const MV_Vector<TYPE> &a, const MV_Vector<TYPE> &b){
     unsigned int m = a.size(), n = b.size();
     MV_ColMat<TYPE> c(m,n);
     for(unsigned int j=0; j < n; j++){
       for(unsigned int i=0; i < m; i++) c(i,j) = a(i)*b(j);
     }
     return c;
   }

/* convert matrix to vector by flattenning along columns or rows */

template<class TYPE>
MV_Vector<TYPE> flatten0(const MV_ColMat<TYPE> &a){
  unsigned int m = a.size(0), n = a.size(1);
  MV_Vector<TYPE> r(m*n);
  for(unsigned int j=0; j<n; j++){
    for(unsigned int i=0; i<m; i++) r(j*m+i) = a(i,j);
  }
  return r;
}
template<class TYPE>
MV_Vector<TYPE> flatten1(const MV_ColMat<TYPE> &a){
  unsigned int m = a.size(0), n = a.size(1);
  MV_Vector<TYPE> r(m*n);
  for(unsigned int j=0; j<n; j++){
    for(unsigned int i=0; i<m; i++) r(i*n+j) = a(i,j);
  }
  return r;
}
  


#endif /* __MV_matrix__  */

