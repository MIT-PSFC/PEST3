      SUBROUTINE pstkernel(xobs,zobs,xsce,zsce,grdgre,gren,j1,j2, &
     &  iopw,iops)
!
!..... calculates the kernels of the integral equation for laplace's
!     equation for a torus.
!     j1, j2 correspond to the four blocks in grdgre(j1,j2).  they
!     respectively indicate observer and source.
!     1 and 2 are plasma and wall blocks, respectively.
!     iopw = 1 calculates bval terms on rhs.
!     iops = 1 subtract and add log behavior in bval.  log
!              contribution is done analytically.
!
 USE pstcom

 USE l22com

 USE mtrik1
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER J1
      INTEGER J2
      INTEGER IOPW
      INTEGER IOPS
      INTEGER ISCHK
      INTEGER ISPH
      INTEGER ISHAPE
      real*8 RES
      INTEGER JRES
      INTEGER MTHM
      INTEGER ISGN
      real*8 WSIMP1
      real*8 WSIMP2
      real*8 WSIMP4
      real*8 THIRD
      real*8 ALGDTH
      real*8 SLOG1M
      real*8 SLOG0
      real*8 SLOG1P
      real*8 TDTH
      real*8 ALG
      real*8 ALG0
      real*8 ALG1
      real*8 ALG2
      INTEGER I
      real*8 THETA
      INTEGER J
      real*8 THES
      real*8 AVAL1
      INTEGER IEND
      INTEGER ISTART
      INTEGER MTHS
      INTEGER IC
      real*8 WSIMP
      INTEGER J1J2
      INTEGER JS1
      INTEGER JS2
      INTEGER JS3
      INTEGER JS4
      INTEGER JS5
      INTEGER ILR
      real*8 XL
      real*8 XU
      real*8 AGAUS
      real*8 BGAUS
      real*8 C1GAUS
      real*8 C2GAUS
      real*8 C3GAUS
      real*8 C4GAUS
      INTEGER IG
      real*8 TGAUS0
      real*8 PGAUS
      real*8 PGAUS2
      real*8 WGBG
      real*8 AMM
      real*8 AM
      real*8 A0
      real*8 AP
      real*8 APP
      real*8 RESIDU
      logical clockwise




!
!aplet  18.5_r8 . 94._r8 .. to save values of arrays between SUBROUTINE pstcalls...
!aplet  20.8_r8 . 98._r8 .. removed references to xwal & zwal which caused the code
!                 to crash on alpha's. as a result, the routine may not 
!                 work for a spheromak geometry.
!
      data ischk /0/
!aplet  18.5_r8 .94
      save ischk, isph
!aplet  18.5_r8 .94
!
 REAL*8, DIMENSION(nths)          :: 	 the
 REAL*8, DIMENSION(nths)          :: 	 ww1
 REAL*8, DIMENSION(nths)          :: 	 ww2

! aplet  6.2_r8 .97
!! save the,ww1,ww2
! aplet  6.2_r8 .97

 REAL*8, DIMENSION(nths2,nths2)  :: grdgre
 REAL*8, DIMENSION(nths,nths) 	 :: gren
 REAL*8, DIMENSION(nths) 	 :: xobs
 REAL*8, DIMENSION(nths) 	 :: zobs
 REAL*8, DIMENSION(nths) 	 :: xsce
 REAL*8, DIMENSION(nths) 	 :: zsce
 REAL*8, DIMENSION(8) 	 :: tgaus
 REAL*8, DIMENSION(8) 	 :: wgaus
 REAL*8, DIMENSION(4) 	 :: tlog
 REAL*8, DIMENSION(4) 	 :: wlog
 REAL*8, DIMENSION(nths) 	 :: zxpp
 REAL*8, DIMENSION(nths) 	 :: zzpp
 REAL*8, DIMENSION(3) 	 :: tab
 REAL*8, DIMENSION(nths) 	 :: work
 REAL*8, DIMENSION(nths) 	 :: ww3
 REAL*8, DIMENSION(nths) 	 :: xpr
 REAL*8, DIMENSION(nths) 	 :: zpr
 REAL*8	 	 	 :: aval0
 INTEGER, DIMENSION(2) 	 :: iop
 INTEGER	 	 	 :: mw
 INTEGER	 	 	 :: jtop
 INTEGER	 	 	 :: jbot
 COMMON /c2f90aval0/ aval0 
 COMMON /c2f90mw/ mw 
 COMMON /c2f90jtop/ jtop 
 COMMON /c2f90jbot/ jbot 
!
      ishape=6

      clockwise = .TRUE.
      if( (xinf(2)-xinf(1))*(zinf(2)-zinf(1)) < 0.0_r8 ) clockwise = .FALSE.
!
      zxpp(1) =  0._r8 
      zzpp(1) =  0._r8 
      tab(1) =  0._r8 
!!$      ww1(1) =  0._r8 
!!$      ww2(1) =  0._r8 
!!$      ww3(1) =  0._r8 
!
!........ res will contain integral of sign*k0 at j = jres.
!
      res =  0.0_r8 
      jres = 13
!
      mthm = mth - 1
!
!.....isgn is positive over wall and negative over plasma...
!
      isgn = 2*j2 - 3
!
      gren(1:mth1,1:mth1) = zero
!
!.....weights for simpson quadrature
!
      wsimp1 = dth / three
      wsimp2 = two * dth / three
      wsimp4 = four * dth / three
!
!.....weights for eight point gaussian quadrature.
!
      wgaus(1) =  0.101228536290376_r8 
      wgaus(2) =  0.222381034453374_r8 
      wgaus(3) =  0.313706645877887_r8 
      wgaus(4) =  0.362683783378362_r8 
      wgaus(5) = wgaus(4)
      wgaus(6) = wgaus(3)
      wgaus(7) = wgaus(2)
      wgaus(8) = wgaus(1)
!
!.....calculate some quantities for the singular region
!     for the analytic integral over the log behavior of bval.
!
      third = one/three
      algdth = log(dth)
      slog1m = third * dth * ( algdth - third )
      slog0 = four * third * dth * ( algdth - four*third )
      slog1p = slog1m
      tdth = two * dth
      alg = log(tdth)
      alg0 =  16.0_r8 *dth*(alg -  68.0_r8 /15.0) /  15.0_r8 
      alg1 =  128.0_r8 *dth*(alg -  8.0_r8 /15.0) /  45.0_r8 
      alg2 =  4.0_r8 *dth*( 7.0_r8 *alg -  11.0_r8 /15.0) /  45.0_r8 
!
!.....get x coordinates for wall to check for spher case. put them in
!

the =  0._r8 
the = (/ (  (i-1)* dth, i=1,mth1 ) /)


      if ( ischk  >   0 ) go to 9

      jbot = mth/2 + 1
      jtop = mth/2 + 1
!
! isph = 1 => spheromak?
! isph = 0 => torus
!
      isph = 0
!!
!! aplet 20-Aug-98
!! following lines removed because they use xwal, zwal, which 
!! have not been set during first call. this causes the code to
!! crash on alpha's.
!!
!!$      do 10 i = 1, mth1
!!$      the(i) = (i-1) * dth
!!$!
!!$!......take care of isph case when points on wall are near to the
!!$!     centre line.  that is, lagrange interpolation points are always to
!!$!     right of the center line.
!!$!
!!$      if ( i  ==  mth1 ) go to 10
!!$      if ( xwal(i)*xwal(i+1)  >   zero ) go to 10
!!$      if ( xwal(i)  >  zero ) jbot = i
!!$      if ( xwal(i)  <  zero ) jtop = i+1
!!$      isph = 1
!!$   10 continue
!
    9 continue
      ischk = 1
!
!
!     periodic boundary conditions.
!
      iop(1)  = 4
      iop(2)  = 4

      call pstspl1d1( mth1, the, xsce, zxpp, iop, 1, ww1, ww2, ww3 )
!     ( note.. ww1, ww2, ww3 just       used as temporary working storage here )
      call pstspl1d1( mth1, the, zsce, zzpp, iop, 1, ww1, ww2, ww3 )

!
      do 11 i = 1, mth1
!
      theta = (i-1)*dth
      call pstspl1d2( mth1, the, xsce, zxpp, 1, theta, tab )
      xpr(i) = tab(2)
      call pstspl1d2( mth1, the, zsce, zzpp, 1, theta, tab )
      zpr(i) = tab(2)
!
   11 continue
!
!.....index j runs over the observer points.
!
      do 200 j = 1, mth
!
      xs = xobs(j)
      zs = zobs(j)
      thes = the(j)
     
!.....calculate non-singular part first.
!
!.....  array work will contain grdgre kernel.
!
      work(1:mth1) = zero
!
!.....there is no observer for xobs negative.  put zero in these
!     positions in grdgre.
!
      if ( xs  <  zero ) go to 175
!
      aval1 = zero
!
!.....  in the next do loop, ic indicates the theta values of the
!       source points in the non singular region.   i is the just the
!       integration index.  mths is the number of integration points.
      iend = 2
      if ( .not. (isph  ==  1  .AND.  j2  ==  2) ) go to 15
      if ( jbot-j  ==  1 ) iend = 3
      if ( jbot-j  ==  0 ) iend = 4
      if ( j-jtop  ==  0 ) iend = 0
      if ( j-jtop  ==  1 ) iend = 1
   15 continue
      istart = 4 - iend
!
      mths = mth - (istart+iend-1)
!
      do 25 i = 1, mths
!
      ic = j+i+istart-1
      if ( ic  >=   mth1 ) ic = ic - mth
!     ic = mod ( j+i-1, mth ) +1
      theta = (ic-1) * dth
!
      xt = xsce(ic)
      zt = zsce(ic)
!
!.....there are no sources from negative x. for isph=1, the virtual
!     sources are taken care of through iend above.
!
      if ( xt  <  zero ) go to 25
!
      xtp  = xpr(ic)
      ztp  = zpr(ic)
!
!.....ic should never be equal to j unless ispher=1, in which
!     case we are on the major axis. negligible contribution there.
!
      if ( ic  ==  j ) go to 25
!
      CALL pstgreen
!
!.....weights for simpson integration of grdgre.
!..... note that mths must be odd.
!
!
      wsimp = wsimp2
      if ( (i/2)*2  ==  i ) wsimp = wsimp4
      if ( i  ==  1  .OR.  i  ==  mths ) wsimp = wsimp1
!
      work(ic) = work(ic) +  isgn * aval * wsimp
      gren(j,ic) = gren(j,ic) + bval * wsimp
!
      aval1 = aval1 + aval0 * wsimp
!
   25 continue
!
!.....end of non-singular part.
!
!.....dont compute singular contributions in off-diagonal blocks
!     when xt less than zero.
!
      j1j2 = j1 + j2
      if (j1j2  ==  2 ) go to 28
      if ( (isph == 1)  .AND.  (j >  jbot)  .AND.  (j < jtop) ) go to 175
   28 continue
!
!.....  singular region is done only if iops =  1._r8  an eight point
!.....  gaussian quadrature is used.
!
      thes = the(j)
      js1 = mod ( j-iend+mth-1, mth ) + 1
      js2 = mod ( j-iend+mth, mth ) + 1
      js3 = mod ( j-iend+mth+1, mth ) + 1
      js4 = mod ( j-iend+mth+2, mth ) + 1
      js5 = mod ( j-iend+mth+3, mth ) + 1
!
      do 165 ilr = 1, 2
!
!.....  xl, xu is upper and lower limits of the gaussian integration
!.....  in the singular region. its done in two parts.
!
      xl = thes + (2*ilr - iend - 2)*dth
      xu = xl + tdth
!
      agaus = half*(xu+xl)
      bgaus = half*(xu-xl)
!
      c1gaus =  0.960289856497536_r8  * bgaus
      c2gaus =  0.796666477413627_r8  * bgaus
      c3gaus =  0.525532409916329_r8  * bgaus
      c4gaus =  0.183434642495650_r8  * bgaus
!
      tgaus(1) = agaus - c1gaus
      tgaus(2) = agaus - c2gaus
      tgaus(3) = agaus - c3gaus
      tgaus(4) = agaus - c4gaus
      tgaus(5) = agaus + c4gaus
      tgaus(6) = agaus + c3gaus
      tgaus(7) = agaus + c2gaus
      tgaus(8) = agaus + c1gaus
!
      do 160 ig = 1, 8
!
      tgaus0 = tgaus(ig)
      if ( tgaus0  <  zero ) tgaus0 = twopi + tgaus0
      if ( tgaus0  >=   twopi ) tgaus0 = tgaus0 - twopi
      call pstspl1d2 ( mth1, the, xsce, zxpp, 1, tgaus0, tab )
      xt = tab(1)
      xtp = tab(2)
      call pstspl1d2 ( mth1, the, zsce, zzpp, 1, tgaus0, tab )
      zt = tab(1)
      ztp = tab(2)
!
      CALL pstgreen
!
!..... subtract out log behavior if iops is one.
!
      bval = bval + iops * log( ( thes-tgaus(ig) )**2) / xs
!
!.....chi is approximated in the wall plasma region by
!     a five point lagrange interpolation formula.
!
      pgaus = ( tgaus(ig) - thes-(2-iend)*dth ) / dth
      pgaus2 = pgaus*pgaus
      wgbg = wgaus(ig) * bgaus
!
      amm = (pgaus2-one) * pgaus * ( pgaus - two ) /  24.0_r8 
      amm = amm * wgbg
      work(js1) = work(js1) + isgn*aval*amm
!
      am = -(pgaus-one) * pgaus * (pgaus2-four) /  6.0_r8 
      am = am * wgbg
      work(js2) = work(js2) + isgn*aval*am
!
      a0 = (pgaus2-one) * (pgaus2-four) / four
      work(js3) = work(js3) + isgn*a0*aval*wgbg
!
!
      ap = -(pgaus+one) * pgaus * (pgaus2-four) /  6.0_r8 
      ap = ap * wgbg
      work(js4) = work(js4) + isgn*aval*ap
!
      app = (pgaus2-one) * pgaus * ( pgaus+two) /  24.0_r8 
      app = app * wgbg
      work(js5) = work(js5) + isgn*aval*app
!
!.....add in integral over k0 in singular region.
!
      work(j) = work(j) - isgn*aval0*wgbg
!
      if ( j  ==  jres ) res = res - isgn*aval0*wgbg
!
!..... dont need gren if iopw is zero.
!
      if ( iopw  ==  0 ) go to 160
!
      gren(j,js1) = gren(j,js1) + bval*amm
      gren(j,js2)  = gren(j,js2)  + bval*am
      gren(j,js3)   = gren(j,js3)   + bval*a0*wgbg
      gren(j,js4)  = gren(j,js4)  + bval*ap
      gren(j,js5) = gren(j,js5) + bval*app
!
  160 continue
  165 continue
!
!..... add in integral over k0 in nonsingular region.
!
!     work(j) = work(j) - isgn*( 2.0_r8 *(2-j1) + aval1)
      residu =  0.0_r8 
      if ( j1  ==  j2 ) residu =  2.0_r8 
      if ( ishape  <  10 ) residu = -isgn *  2.0_r8 *(2-j1)
      if(.NOT. clockwise) residu = - residu
      
      work(j) = work(j) - isgn*aval1 + residu
!
      if ( j  ==  jres ) res = res - isgn*aval1
!
      if ( iops  /=  1  .OR.  iopw  ==  0 ) go to 170
!
!.......add in analytical integrals over log behavior.
!
      gren(j,js1) = gren(j,js1) - alg2/xs
      gren(j,js2)  = gren(j,js2) - alg1/xs
      gren(j,js3)   = gren(j,js3) - alg0/xs
      gren(j,js4)  = gren(j,js4) - alg1/xs
      gren(j,js5) = gren(j,js5) - alg2/xs
!
  170 continue
!
  175 continue
!
      if ( (xs  <  zero)  .AND.  (j2  ==  2) ) work(j) =  1.0_r8 
!
!.....now put work vector in appropriate place in grdgre .....
!
      do ic = 1, mth
         grdgre ( (j1-1)*mth+j, (j2-1)*mth+ic ) = work(ic)
      end do
!
  200 continue
!
      if ( checkd ) write ( outmod,330 ) jres, j1,j2, res
  330 format ( /,1x,"jres, j1,j2, res*isgn = ", 3i3, e14.5,/ )
!
      return
      END SUBROUTINE pstkernel

