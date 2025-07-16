
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                   MV++ Numerical Matrix/Vector C++ Library                */
/*                             MV++ Version 1.5                              */
/*                                                                           */
/*                                  R. Pozo                                  */
/*               National Institute of Standards and Technology              */
/*                                                                           */
/*                                  NOTICE                                   */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that this permission notice appear in all copies and             */
/* supporting documentation.                                                 */
/*                                                                           */
/* Neither the Institution (National Institute of Standards and Technology)  */
/* nor the author makes any representations about the suitability of this    */
/* software for any purpose.  This software is provided ``as is''without     */
/* expressed or implied warranty.                                            */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//
//      mvmtp.h  : basic templated numerical matrix class, storage
//                  by columns (Fortran oriented.)
//
//
//

/* 

Additions by A. Pletzer (Aug 98)

The () operator accepts arguments that are integer vectors, as in:
MV_Vector<TYPE> operator()(const MV_Vector<int> &I, const unsigned int j )

slip, a member function that allows one to absorb a vector into a matrix
location given by the integer vectors I and J.
The size of the vector must match those of the I and J integer vectors.
MV_ColMat<TYPE> slip(const MV_Vector<TYPE> &v, const MV_Vector<int> &I, const MV_Vector<int> &J)

*/


#ifndef _MV_MATRIX_H_
#define _MV_MATRIX_H_    

#include "mvvtp.h"

struct Matrix_ 
{
    enum ref_type {  ref = 1 };
};


#include <iostream.h>       // for formatted printing of matrices
#ifdef MV_MATRIX_BOUNDS_CHECK
#   include <assert.h>
#endif


template <class TYPE>
class MV_ColMat
{                                                                      
    private:                                                           
           MV_Vector<TYPE> v_;
           unsigned int dim0_;   // preferred to using dim_[2]. some compilers
           unsigned int dim1_;   // refuse to initalize these in the constructor.
           unsigned int lda_;
           int ref_;   // true if this is declared as a reference vector,
                        // i.e. it does not own the memory space, but 
                        // rather it is a view to another vector or array.
    public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
            MV_ColMat();                             
            MV_ColMat(unsigned int, unsigned int); 

    // some compilers have difficulty with inlined 'for' statements.
    MV_ColMat(unsigned int, unsigned int, const TYPE&);   

    // usual copy by value
    // (can't use default parameter lda=m, because m is not a constant...)
    //
    MV_ColMat(TYPE*, unsigned int m, unsigned int n);
    MV_ColMat(TYPE*, unsigned int m, unsigned int n, unsigned int lda);

    // the "reference" versions
    //
    //
    MV_ColMat(TYPE*, unsigned int m, unsigned int n, Matrix_::ref_type i);
    MV_ColMat(TYPE*, unsigned int m, unsigned int n, unsigned int lda,
                Matrix_::ref_type i);

    MV_ColMat(const MV_ColMat<TYPE>&); 
    ~MV_ColMat();                              
                                                                       

  MV_ColMat<TYPE> operator+=(const TYPE f);
  MV_ColMat<TYPE> operator-=(const TYPE f);
  MV_ColMat<TYPE> operator*=(const TYPE f);
  MV_ColMat<TYPE> operator/=(const TYPE f);
  MV_ColMat<TYPE> operator+=(const MV_ColMat<TYPE> &b);
  MV_ColMat<TYPE> operator-=(const MV_ColMat<TYPE> &b);
  MV_ColMat<TYPE> operator*=(const MV_ColMat<TYPE> &b);
  MV_ColMat<TYPE> operator/=(const MV_ColMat<TYPE> &b);

        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       
    inline TYPE&        operator()(unsigned int, unsigned int); 
    inline const TYPE&  operator()(unsigned int, unsigned int) const; 
    MV_ColMat<TYPE> operator()(const MV_VecIndex &I, const MV_VecIndex &J) ;
    const MV_ColMat<TYPE> operator()(const MV_VecIndex &I, const MV_VecIndex &J) const;
    unsigned int            size(int i) const; 
    MV_ColMat<TYPE>&        newsize(unsigned int, unsigned int);
    int ref() const { return ref_;}
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    MV_ColMat<TYPE> & operator=(const MV_ColMat<TYPE>&);
    MV_ColMat<TYPE> & operator=(const TYPE&);\

    void ran(void);


    //A. Pletzer  friend ostream& operator<<(ostream &s, const MV_ColMat<TYPE> &A);

    MV_ColMat<TYPE> &alloc(const unsigned int nr, const unsigned int nc);

/*     MV_ColMat<TYPE> &MV_ColMat<TYPE>::slip(const MV_Vector<TYPE> &v,   */
/* 				 const MV_Vector<int> &I,  */
/* 				 const MV_Vector<int> &J); */


/* 
Access for integer vector index (need not be contiguous)
A. Pletzer July 1998
*/

MV_ColMat<TYPE> operator()(const MV_Vector<int> &I, const MV_Vector<int> &J ) 
{
       unsigned int nr = I.size();
       unsigned int nc = J.size();
       MV_ColMat<TYPE> y(nr, nc);
       for( unsigned int j=0; j<nc; ++j) {
          for( unsigned int i=0; i<nr; ++i) {
             y(i, j) = v_[J(j)*dim0_ + I(i)];
          }
       }
       return y;
}

MV_ColMat<TYPE> operator()(const MV_Vector<int> &I, const MV_Vector<int> &J ) const
{
       unsigned int nr = I.size();
       unsigned int nc = J.size();
       MV_ColMat<TYPE> y(nr, nc);
       for( unsigned int j=0; j<nc; ++j) {
          for( unsigned int i=0; i<nr; ++i) {
             y(i, j) = v_[J(j)*dim0_ + I(i)];
          }
       }
       return y;
}

MV_Vector<TYPE> operator()(const MV_Vector<int> &I, const unsigned int j ) 
{
       unsigned int n = I.size();
       MV_Vector<TYPE> y(n);
          for( unsigned int i=0; i<n; ++i) {
             y(i) = v_[j*dim0_ + I(i)];
          }
       return y;
}


MV_Vector<TYPE> operator()(const MV_Vector<int> &I, const unsigned int j ) const 
{
       unsigned int n = I.size();
       MV_Vector<TYPE> y(n);
          for( unsigned int i=0; i<n; ++i) {
             y(i) = v_[j*dim0_ + I(i)];
          }
       return y;
}

MV_Vector<TYPE> operator()(const unsigned int i, const MV_Vector<int> &J ) 
{
       unsigned int n = J.size();
       MV_Vector<TYPE> y(n);
          for( unsigned int j=0; j<n; ++j) {
             y(j) = v_[J(j)*dim0_ + i];
          }
       return y;
}


MV_Vector<TYPE> operator()(const unsigned int i, const MV_Vector<int> &J ) const
{
       unsigned int n = J.size();
       MV_Vector<TYPE> y(n);
          for( unsigned int j=0; j<n; ++j) {
             y(j) = v_[J(j)*dim0_ + i];
          }
       return y;
}

/* 
slip vector into matrix. 
absorb vector elements of v into matrix. Every element of v(i) 
takes a matrix slot at I(i) and J(i), i = 1,...v.size().
   */

MV_ColMat<TYPE> slip(const MV_Vector<TYPE> &v, const MV_Vector<int> &I, const MV_Vector<int> &J) {

  unsigned int nv = v.size();

  if( nv != I.size() || nv != J.size() ) {
    cerr << "MV_ColMat<T>::slip(const MV_Vector<T>, const MV_Vector<int> "
         << ", const MV_Vector<int> ) Incompatible sizes \n";
    exit(1);
  }
  for(unsigned int i=0; i<nv; ++i) v_(J(i)*dim0_ + I(i)) = v(i);
  return *this;
}

MV_ColMat<TYPE> slip(const MV_Vector<TYPE> &v, const MV_Vector<int> &I, 
		     const int j) {

  unsigned int nv = v.size();

  if( nv != I.size() ) {
    cerr << "MV_ColMat<T>::slip(const MV_Vector<T>, const MV_Vector<int> "
         << ", const int) Incompatible sizes \n";
    exit(1);
  }
  for(unsigned int i=0; i<nv; ++i) v_(j*dim0_ + I(i)) = v(i);
  return *this;
}

MV_ColMat<TYPE> slip(const MV_Vector<TYPE> &v, const int i, 
		     const MV_Vector<int> &J) {

  unsigned int nv = v.size();

  if( nv != J.size() ) {
    cerr << "MV_ColMat<T>::slip(const int, const MV_Vector<int> "
         << ", const MV_Vector<int> ) Incompatible sizes \n";
    exit(1);
  }
  for(unsigned int j=0; j<nv; ++j) v_(J(j)*dim0_ + i) = v(j);
  return *this;
}

/* stretch: allows to extend a vector into a matrix */
// trying to declare a mixed method
/* friend MV_ColMat<TYPE> MV_Vector<TYPE>::ystretch(const unsigned int ny); */


};                                                                     


template<class TYPE>
unsigned int MV_ColMat<TYPE>::size(int i) const 
{
    if (i==0) return dim0_;
    if (i==1) return dim1_;
    else
    {
     cerr << "Called MV_ColMat::size(" << i << ")  must be 0 or 1 " << endl;
     exit(1);
    }

    // never should be here, but many compilers warn about not
    // returning a value
    return 0;
}

// NOTE: null construct have ref_ flag turned OFF, otherwise, we can
//          never reset the size of matrix....
template <class TYPE>
MV_ColMat<TYPE>::MV_ColMat()  : v_(), dim0_(0), dim1_(0) , lda_(0), ref_(0){}
                                                                

template <class TYPE>
MV_ColMat<TYPE>::MV_ColMat(unsigned int m, unsigned int n) : v_(m*n),
        dim0_(m), dim1_(n), lda_(m), ref_(0) {}

template <class TYPE>
MV_ColMat<TYPE>::MV_ColMat(unsigned int m, unsigned int n, const TYPE &s) : v_(m*n),
        dim0_(m), dim1_(n), lda_(m), ref_(0) 
{
    operator=(s);
}

// operators and member functions



template <class TYPE>
inline TYPE& MV_ColMat<TYPE>::operator()(unsigned int i, unsigned int j)
{
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(0<=i && i<size(0));
    assert(0<=j && j<size(1));
#endif 
    return v_(j*lda_ + i);      // could use indirect addressing
                                // instead...
}

template <class TYPE>
inline const TYPE& MV_ColMat<TYPE>::operator()
                    (unsigned int i, unsigned int j) const
{
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(0<=i && i<size(0));
    assert(0<=j && j<size(1));
#endif
    return v_(j*lda_ + i);
}


template <class TYPE>
MV_ColMat<TYPE>& MV_ColMat<TYPE>::operator=(const TYPE & s) 
{
    unsigned int M = size(0);
    unsigned int N = size(1);

    if (lda_ == M)      // if continuous, then just assign as a ?
        v_ =  s;        // single long vector.

    else
    {
    // this should run much faster than the just accessing each (i,j)
    // element individually 
    //

        MV_VecIndex I(0,M-1);
        for (unsigned int j=0; j<N; j++)
        {
            v_(I) = s;
            I += lda_;
        }
    }

    return *this;
}

template <class TYPE>
MV_ColMat<TYPE>& MV_ColMat<TYPE>::newsize(unsigned int M, unsigned int N)
{
    v_.newsize(M*N);
    dim0_ = M;
    dim1_ = N;
    lda_ = M;

    return *this;
}

template <class TYPE>
MV_ColMat<TYPE>& MV_ColMat<TYPE>::operator=(const MV_ColMat<TYPE> & m) 
{

    unsigned int lM = dim0_;     // left hand arg  (this)
    unsigned int lN = dim1_;

    unsigned int rM = m.dim0_;   // right hand arg (m)
    unsigned int rN = m.dim1_;


    // if the left-hand side is a matrix reference, the we copy the
    // elements of m *into* the region specfied by the reference.
    // i.e. inject().

    if (ref_)
    {
        // check conformance,       
        if (lM != rM  || lN != rN)      
        {
            cerr << "MV_ColMatRef::operator=  non-conformant assignment.\n";
            exit(1);
        }
    }
    else
    {
        newsize(rM,rN);
    }

    // at this point the left hand and right hand sides are conformant

    // this should run much faster than the just accessing each (i,j)
    // element individually 

    // if both sides are contigous, then just copy as one vector
    if ( lM == lda_ && rM == m.lda_)
    {
        MV_VecIndex I(0,rM*rN-1);
        v_(I) = m.v_(I);
    }
    else
    {
        // slower way...

        MV_VecIndex I(0,rM-1);
        MV_VecIndex K(0,rM-1);
        for (unsigned int j=0; j<rN; j++)
        {
            v_(I) = m.v_(K);
            I += lda_;
            K += m.lda_;
        }
    }

    return *this;   
}

template <class TYPE>
MV_ColMat<TYPE>::MV_ColMat(const MV_ColMat<TYPE> & m) : 
        v_(m.dim0_*m.dim1_), dim0_(m.dim0_),
        dim1_(m.dim1_), lda_(m.dim0_), ref_(0)
{

    unsigned int M = m.dim0_;
    unsigned int N = m.dim1_;

    // this should run much faster than the just accessing each (i,j)
    // element individually 

    MV_VecIndex I(0,M-1);
    MV_VecIndex K(0,M-1);
    for (unsigned int j=0; j<N; j++)
    {
        v_(I) = m.v_(K);
        I += lda_;
        K += m.lda_;
    }
}


template <class TYPE>
inline MV_ColMat<TYPE>::MV_ColMat(TYPE* d, unsigned int m, unsigned int n,
        Matrix_::ref_type i ):
            v_(d,m*n, MV_Vector_::ref), dim0_(m), dim1_(n), lda_(m), ref_(i) {}

template <class TYPE>
inline MV_ColMat<TYPE>::MV_ColMat(TYPE* d, unsigned int m, unsigned int n, 
            unsigned int lda, Matrix_::ref_type i) : 
            v_(d, lda*n, MV_Vector_::ref), dim0_(m), dim1_(n), lda_(lda), 
            ref_(i) {}

template <class TYPE>
MV_ColMat<TYPE>::MV_ColMat(TYPE* d, unsigned int m, unsigned int n) :
    v_(m*n), dim0_(m), dim1_(n), lda_(m), ref_(0)
{
    unsigned int mn = m*n;

    // d is contiguous, so just copy 1-d vector
    for (unsigned int i=0; i< mn; i++)
            v_[i] = d[i];
}


template <class TYPE>
MV_ColMat<TYPE>::MV_ColMat(TYPE* d, unsigned int m, unsigned int n, 
        unsigned int lda) :
    v_(m*n), dim0_(m), dim1_(n), lda_(lda), ref_(0)
{
    for (unsigned int j=0; j< n; j++)
        for (unsigned int i=0; i<m; i++)
            operator()(i,j) = d[j*lda + i];   // could be made faster!!
}


template <class TYPE>
MV_ColMat<TYPE> MV_ColMat<TYPE>::operator()(const MV_VecIndex &I, const MV_VecIndex &J)
{
    // check that index is not out of bounds
    //
    if (I.end() >= dim0_  || J.end() >= dim1_)
    {
        cerr << "Matrix index: (" << I.start() << ":" << I.end()  
             << "," << J.start() << ":" << J.end()   
             << ") not a subset of (0:" << dim0_ - 1 << ", 0:" 
             << dim1_-1 << ") " << endl;
        exit(1);
    }

    // this automatically returns a reference
    // 
    return MV_ColMat<TYPE>(&v_[J.start()*lda_ + I.start()], 
            I.end() - I.start() + 1, 
            J.end() - J.start() + 1, lda_, Matrix_::ref);
}

template <class TYPE>
const MV_ColMat<TYPE> MV_ColMat<TYPE>::operator()(const MV_VecIndex &I, 
    const MV_VecIndex &J) const
{

    // cerr << "Const operator()(MV_VecIndex, MV_VecIndex) called " << endl;

    // check that index is not out of bounds
    //
    if (I.end() >= dim0_  || J.end() >= dim1_)
    {
        cerr << "Matrix index: (" << I.start() << ":" << I.end()  
             << "," << J.start() << ":" << J.end()   
             << ") not a subset of (0:" << dim0_ - 1 << ", 0:" 
             << dim1_-1 << ") " << endl;
        exit(1);
    }

    // this automatically returns a reference.  we need to 
    // "cast away" constness here, so the &v_[] arg will
    // not cause a compiler error.
    //
    MV_ColMat<TYPE> *t =  (MV_ColMat<TYPE>*) this;
    return MV_ColMat<TYPE>(&(t->v_[J.start()*lda_ + I.start()]), 
            I.end() - I.start() + 1, 
            J.end() - J.start() + 1, lda_, Matrix_::ref);
}

template <class TYPE>
MV_ColMat<TYPE>::~MV_ColMat() {}

template <class TYPE>
ostream&   operator<<(ostream& s, const MV_ColMat<TYPE>& V)
{
    unsigned int M = V.size(0);
    unsigned int N = V.size(1);

    for (unsigned int i=0; i<M; i++)
    {
        for (unsigned int j=0; j<N; j++)
            s << V(i,j) << " " ;
        s << endl;
    }
    
    return s;
}

  
template <class TYPE>
MV_ColMat<TYPE> &MV_ColMat<TYPE>::alloc(const unsigned int nr, const unsigned int nc) {
  v_.alloc( nr*nc );
  dim0_ = nr;
  dim1_ = nc;
  lda_  = nr;
  return *this;
}


#endif 
// _MV_MATRIX_H_

