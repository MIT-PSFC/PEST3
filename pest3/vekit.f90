!-----------------------------------------------------------------------
         SUBROUTINE pstvekit(val,eps,vepmac,npinv,npoutv)
!
! checked for up-down asym. case ap  17.10_r8 .95
!
!     lcm (vekit)
 USE pstcom

 USE combla

 USE comivi

 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 EPS
      real*8 VEPMAC
      INTEGER ISET
      INTEGER NTYPE
      INTEGER NOUT
      INTEGER NCOMP
      INTEGER IDA
      INTEGER I
      real*8 X1NORM
      real*8 GRAND

!***
!***
!
      data iset/0/
 COMPLEX*16	 	 	 :: zxv
 REAL*8, DIMENSION(*) 	 :: val
 INTEGER, DIMENSION(*) 	 :: npinv
 INTEGER, DIMENSION(*) 	 :: npoutv
      iset = iset + 1
!
!         1._r8    initialize
!
         ntype=npinv(1)
         ng=npinv(2)
         mf=npinv(3)
         ml=npinv(4)
         nitmax=npinv(5)
         nout=npinv(7)
         nda=npinv(8)
         ndb=npinv(9)
         ndlt=npinv(10)
         nds=npinv(14)
!
!
         al0=val(1)
         val(2) =  0.0_r8 
         val(3) =  0.0_r8 
         nconv=0
         nit   = 0
         ndscrc=nds
         neg=0
         m11=mf+1
         ncomp=ng*ml+mf
         m12=mf+ml
         nlong=((mf+ml)*(mf+ml+1))/2
!
         CALL pstfxfu
         if (ntype  /=  2) go to 200
         al0= 0.0_r8 
         ida=nda
         nda=ndb
         ndb=ida
 200     continue
!
!         2._r8    shift
!
         CALL pstbshift(a)
!
!         3._r8    decompose a
!
         CALL pstconalr(vepmac)
         npoutv(1)=neg
         npoutv(4)=nsing
         if (nsing  ==  -1 ) go to 600
         if(ntype /= 0)go to 600
!
!         4._r8    iterations
!
         do 410 nit=1,nitmax
         xnorm= 0.0_r8 
 900  format(/,' eigenvector',/,(10e12.4))
         CALL pstcbxlu(0)
         xnorm=sqrt(xnorm)
!
!     normalize x
!
         CALL pstconorm(eps)
!
!     restore
!
         do  400  i=1,ncomp
         zxv=xt(i)
         xt(i)=vt(i)
         vt(i)=zxv
 400     continue
         if ( nconv  ==  ncomp )  go to 420
  410    continue
  420    continue
         npoutv(2)=nit
         npoutv(3)=nconv
!
!         5._r8    eigenvalues
!
         x1norm=1/xnorm
         CALL psteigeng(eps)
         val(3)=alam
         grand= 1._r8 
         if ((alam-al0)  <   0._r8  )   grand=- 1._r8 
         val(2)=val(1)+x1norm*grand
!!$         CALL pstendxuv
!
!         6._r8    diagnostics
!
  600    continue
      if ( iset  >  1 )  go to 602
         write (nout,700)
         write (nout,710)
         write (nout,720) ng,m12,mf,ncomp,nitmax,nout,nda,ndb,nds,ndlt
  602 continue
         write (nout,730)
         write (nout,740) neg,val(1),nit,nconv
         write (nout,750) val(2),val(3)
      CALL pstempty (nout)
!
!         7._r8    formats
!
!
 700    format (1x,10x,'output from vector iteration',/, &
              10x,'*****************************',//)
!
!
 710     format (//,11x,'input quantities',/, &
              10x,'------------------')
!
!
 720     format (//,11x,'number of blocks        ' ,i5, &
   //,11x,'length of a block        ',i5, &
   //,11x,'number of overlaps       ',i5, &
   //,11x,'number of vectorcomponent',i5, &
   //,11x,'max. number of iterations',i5, &
   //,11x,'output unit              ',i5, &
   //,11x,'unit for matrix a        ',i5, &
   //,11x,'unit for matrix b        ',i5, &
   //,11x,'unit for iteration       ',i5, &
   //,11x,'unit for decomposition   ',i5,//)
!
! c2f90: contracted format
!
 730     format (//,11x,' output quantities',/ &
,              10x,'-------------------')
!
!
 740     format (//,10x,i3,' eigenvalues less than',e13.5,&
           //,10x,i3,' iterations done', &
           //,8x,i5,' vector components have converged')
!
 750     format (//,11x,' eigenvalue from normalization',e13.5,  &
         //,11x,'rayleigh quotient            ',e13.5)
!
         return
         end


