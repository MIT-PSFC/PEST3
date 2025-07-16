!......................................................................
      subroutine pstsaver(nzz)
!
! up-down asym. compatibility check ap  17.10_r8 .95
!
!......................................................................
!
 USE pstcom

 USE combla

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER NZZ
      INTEGER ISET
      INTEGER IS
      INTEGER LNVEC
      INTEGER MAT2
      INTEGER ID1
      INTEGER II1
      INTEGER LQ
      INTEGER LQ1
      INTEGER LL1
      INTEGER LL2
      INTEGER L
      INTEGER NADRES
      INTEGER NADRE2



!
      data iset/0/
 REAL*8, DIMENSION(10) 	 :: eigvag
 REAL*8, DIMENSION(3) 	 :: al
 REAL*8	 	 	 :: epscon
 REAL*8	 	 	 :: epsmac
 INTEGER, DIMENSION(14) 	 :: npin
 INTEGER, DIMENSION(4) 	 :: npout
 COMMON /c2f90eigvag/ eigvag 
 COMMON /c2f90al/ al 
 COMMON /c2f90epscon/ epscon 
 COMMON /c2f90epsmac/ epsmac 
 COMMON /c2f90npin/ npin 
 COMMON /c2f90npout/ npout 
      iset = iset + 1
!
!      first time in create a temporary file called vector to hold
!      the eigenvector for later use.
!
      is = 1
      lnvec = 10000
!***
      mat2 = 2*mat
!***
!
      if ( iset  >  1 )  go to 10
      if ( nzz  ==  1 ) CALL pstzop(outvec,'vector',10000,id1,is, 7000 )
!
   10 continue
!
!      write(outsol,100)
      write(outpst,100)
!
!
!      write(outsol, 102) eigvag(nzz)
!      write(outsol,104)
      write(outpst, 102) eigvag(nzz)
      write(outpst,104)
      do 600 ii1 = 1, mp
      lq = jtot(ii1)
       lq1 = lq + jsub(ii1) - 1
      ll1 = lmin(ii1)
      ll2 = lmax(ii1)
!***      if(ii1 == is)write(outsol,108)(l,l=ll1,ll2)
!***      if(ii1 == is)write(outpst,107)(l,l=ll1,ll2)
!***      write(outsol,109) psinod(ii1), (xt(l),l=lq,lq1)
!***      write(outpst,101)(xt(l),l=lq,lq1)
      ll1 = min0(0,ll1)
      ll2 = max0(5,ll2)
!      if(ii1 == is)write(outsol,108)(l,l=ll1,ll2),(l,l=ll1,ll2)
      if(ii1 == is)write(outpst,107)(l,l=ll1,ll2),(l,l=ll1,ll2)
!      write(outsol,109) psinod(ii1), ( real(xt(l)),l=ll1,ll2),
!     .                               (aimag(xt(l)),l=ll1,ll2)
      write(outpst,101)( real(xt(l)),l=ll1,ll2), &
                 (aimag(xt(l)),l=ll1,ll2)
  600 continue
!
!      now enter into vector
!
!
       nadres = ( nzz-1 )*mat + 1
       nadre2 = ( nzz-1 ) * mat2 + 1
!***      CALL pstzwr ( outvec,xt(1),mat,nadres,is, 7000 )
       CALL pstzwr ( outvec,xt(1),mat2, nadre2,is, 7000 )
!***
!
107    format(8x,5("  re l = ",i3),5("  im l = ",i3))
108    format("psinod  ",5("  re l = ",i3),5("  im l = ",i3))
109    format(f8.6,5(1x,e11.4),5(1x,e11.4))
101    format(8x, 5(1x,e11.4),5(1x,e11.4))
 102  format(5x," eigenvalue =",e14.7/)
 104  format(1x," xipsi  "/)
!
  100 format( 1x," unstable eigenvalues follow : ",1x,/)
 9000 format(i5)
 9100 format(4e20.13)
 9200 format(10i5)
 9300 format(e20.13,10i5)
      return
 7000 CALL psterrmes( outpst,'pstsaver',is )
      end



