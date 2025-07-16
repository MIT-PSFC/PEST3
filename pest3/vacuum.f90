!...................................................................
      SUBROUTINE pstvacuum
!...................................................................
!
!      3.2_r8           solution of the vacuum integral equations.
!
!
 
 USE pstcom
 USE newmet

 USE l22com

 USE mtrik1
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IER
      INTEGER ISHAPE
      INTEGER IDGT
      INTEGER JMAX
      INTEGER LMIN1
      real*8 Q
      INTEGER LGIVUP
      INTEGER LENGTH
      INTEGER NADRES
      INTEGER MTHSQ
      INTEGER LMTH
      INTEGER J1V
      INTEGER J2V
      INTEGER L1
      INTEGER LL
      INTEGER J
      real*8 THETA
      real*8 ELTH
      real*8 ELTHNQ
      real*8 SINLTH
      real*8 COSLTH
      INTEGER I
      INTEGER IWP
      INTEGER MTH12
      INTEGER II
      INTEGER I1
      INTEGER I2
      real*8 SUM11
      real*8 SUM12
      INTEGER LMAX2
      real*8 CN0
      INTEGER IOK
      INTEGER IERR
      real*8 ZRCOND
      real*8 ZERR
      INTEGER NFM12
      real*8 AIRS
      real*8 ARIS
      real*8 VACIS
      INTEGER LL1
      real*8 AL1NQ
      INTEGER L2
      real*8 AL2NQ
      INTEGER J1
      real*8 ALNQ1
      INTEGER J2
      real*8 ALNQ2
      real*8 ARRD
      real*8 AIID
      real*8 VACRD
      real*8 VACRSO
      real*8 VACRSD
      real*8 ASYMV1
      real*8 ASYMV2
      INTEGER KK
      INTEGER JJ
      real*8 OUTPEST


!
!
 LOGICAL	 	 	 :: ldelta
 LOGICAL	 	 	 :: lmapck
 LOGICAL	 	 	 :: lvacdsk
 REAL*8	 	 	 :: nq
 REAL*8, DIMENSION(nths2) 	 :: workl
 REAL*8, DIMENSION(nths2,nths2) 	 :: grdgre
!!$ REAL*8, DIMENSION(nths) 	 :: xwal
!!$ REAL*8, DIMENSION(nths) 	 :: zwal
 REAL*8, DIMENSION(nths,nths) 	 :: grwp
 REAL*8, DIMENSION(nths2,nfm21) 	 :: grri
 REAL*8, DIMENSION(nfm,nfm) 	 :: arr
 REAL*8, DIMENSION(nfm,nfm) 	 :: aii
 REAL*8, DIMENSION(nfm,nfm) 	 :: ari
 REAL*8, DIMENSION(nfm,nfm) 	 :: air
 REAL*8, DIMENSION(2) 	 :: summ
 REAL*8, DIMENSION(nfm*nfm) 	 :: work
 REAL*8, DIMENSION(nfm*nfm) 	 :: work1
 REAL*8, DIMENSION(nths) 	 :: delt
 INTEGER	 	 	 :: tmth
 INTEGER, DIMENSION(nths2) 	 :: ipvt
 INTEGER	 	 	 :: mw
 INTEGER	 	 	 :: jtop
 INTEGER	 	 	 :: jbot

! arrays for lapack linear system routine

 REAL*8, dimension(:,:),allocatable :: zamat, ZAF
 REAL*8, dimension(:),allocatable ::  ZBERR,ZFERR
 REAL*8, dimension(:),allocatable  :: zwork
 REAL*8, dimension(:),allocatable ::  zc, zr
 REAL*8, dimension(:,:),allocatable :: zsource, zx , zxerr

 real(r8) clockwise

INTEGER,dimension(:),allocatable :: JPV, JWORK

! arrays for lapack linear system routine

 COMMON /c2f90ldelta/ ldelta 
 COMMON /c2f90lmapck/ lmapck 
 COMMON /c2f90lvacdsk/ lvacdsk 
 COMMON /c2f90mw/ mw 
 COMMON /c2f90jtop/ jtop 
 COMMON /c2f90jbot/ jbot 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( wall )  then
      write(outmod,9300)
      write(itty  ,9300)
      end if
      if ( infwal   .AND.   ( .not. wall) ) then
      write(outmod,9400)
      write(itty  ,9400)
      end if
      if ( .not. wall   .AND.   (.not. infwal) ) then
         if( a <= -10 ) then
            write(outmod,9500) b
            write(itty  ,9500) b
         else 
!
! a and b represent the normalized offset from the magnetic axis 
! both must be < 1 and > -1
!
            if( b > 1   .OR.   b < -1 ) then
               write(outmod,*)' WARNING: b out of range will reset to zero'
               write(  itty,*)' WARNING: b out of range will reset to zero'
               b =  0.0_r8 
            end if
            if( a > 1   .OR.   a < -1 ) then
               write(outmod,*)' WARNING: a out of range will reset to zero'
               write(  itty,*)' WARNING: a out of range will reset to zero'
               a =  0.0_r8 
            end if
            write(outmod,9501) a, b, aw, bw, dw, sw
            write(itty  ,9501) a, b, aw, bw, dw, sw
         end if
      end if
      if ( wall )  write( itty,9300)
      if ( infwal   .AND.   ( .not. wall) ) write( itty,9400)
!
!      check for wall
!
      if(.not. infwal) go to 1300
      a=  22.0_r8 
      b= 22.0_r8 
 1300 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$      xwal(1) =  0._r8 
!!$      zwal(1) =  0._r8 
      ier = 0
      ishape = 6
      idgt = 6
!
      jmax = 2*lmax(1) + 1
!aplet  17.5_r8 .94      jmax1 = lmax(1) - lmin(1) + 1
      lmax1 = lmax(1) + 1
      lmin1 = lmin(1) + 1
      q = qa(nusurf)
      nq = n*q
!     
!     check for vacmat from disk
!     
      if(lvacdsk) then
!aplet  20.5_r8 .94      CALL pstvacdsk
      CALL pstvacds2
      return
      endif
!
!      get delta on interface...
!
      !!lgivup = 1
      !!length = nths * nsf
      !1nadres = 50 + 2*mth2 + nosurf + 10*length + (nosurf-1)*nths
      !!CALL pstzrd ( iomode,delt(1),mth2,nadres,lgivup, 7000 )

      !delt(1:mth2) = delta(nosurf, 1:mth2) !WOW The BUG DPB 7/10
      delt(1:mth2) = delta(1:mth2,nusurf)
!
      tmth = 2*mth
      mthsq = tmth * tmth
      lmth = tmth * 2*jmax1
!
      infwal = .false.
      if ( a  >=    10.0_r8  ) infwal = .true.
!
!.....j1v and j2v are input to vacout and are the size of the matrices.
!
      j1v = nfm
      j2v = nfm
!
      grri = zero
!
      grdgre = zero
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
      CALL pstkernel(xinf,zinf,xinf,zinf,grdgre,grwp,1,1,1,1)
!
      if ( checkd ) &
      CALL pstmatwrt(grwp,nths,nths,9,9,"grwp                " )
!
!.....fourier analyse source points. real and imag. parts.
!     pack into one-d array. ( ... not for leqtif or leqt2f... )
!
!..... the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
!
      do 135 l1 = 1, jmax1
!
      ll = l1 - 1 + lmin(1)
      do 130 j = 1, mth
!
      theta = (j-1) * dth
      elth = ll * theta
      elthnq = elth + nq*delt(j)
      sinlth = sin(elthnq)
      coslth = cos(elthnq)
!
      do 140 i = 1, mth
      grri(i,l1) = grri(i,l1) + coslth*grwp(i,j)
      grri(i,jmax1+l1) = grri(i,jmax1+l1) + sinlth*grwp(i,j)
  140 continue
!
  130 continue
  135 continue
!
      iwp = 1
!
      if ( infwal ) go to 34
!
!.....get wall coordinates....
!
!aplet  16.10_r8 .95      CALL pstvacdat ( infwal, xwal, zwal )
!!$      CALL pstvacda2( infwal, xwal, zwal )
!CALL pstvacda2( infwal )
CALL pstvacda2
!
      CALL pstkernel(xinf,zinf,xwal,zwal,grdgre,grwp,1,2,0,0)
      CALL pstkernel(xwal,zwal,xwal,zwal,grdgre,grwp,2,2,0,0)
      CALL pstkernel(xwal,zwal,xinf,zinf,grdgre,grwp,2,1,1,0)
!
!.....fourier analyse source points. real and imag. parts.
!     pack into one-d array. ( ... not for leqtif or leqt2f... )
!
!..... the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
!
      do 435 l1 = 1, jmax1
!
      ll = l1 - 1 + lmin(1)
      do 430 j = 1, mth
!
      theta = (j-1) * dth
      elth = ll * theta
      elthnq = elth + nq*delt(j)
      sinlth = sin(elthnq)
      coslth = cos(elthnq)
!
      do 440 i = 1, mth
      grri(mth+i,l1) = grri(mth+i,l1) + coslth*grwp(i,j)
      grri(mth+i,jmax1+l1) = grri(mth+i,jmax1+l1) + sinlth*grwp(i,j)
  440 continue
!
  430 continue
  435 continue
!
      iwp = 2
!
   34 continue
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
      mth12 = iwp*mth
!
      if ( n  >    0.1_r8  ) go to 46
!
!.....examine possibility of sum rules......
!
!
      do 45 j = 1, mth12
!
      do 44 ii = 1, 2
!
      summ(ii) = zero
      i1 = (ii-1) * mth + 1
      i2 = ii * mth
!
      do 40 i = i1, i2
!
!     is = 2*(i-1)*mth + j
!     summ(ii) = summ(ii) + grdgre(is)
      summ(ii) = summ(ii) + grdgre ( j, i )
!
   40 continue
!
   44 continue
!
      sum11 = summ(1)
      sum12 = summ(2)
      if ( checkd ) &
      write ( outmod, 9002 ) j, sum11, sum12
 9002 format (  1x, i4, 3x, 5e14.5 )
!
   45 continue
!
   46 continue
!
!.....repack grdgre if wall is far away.
!
!     if ( .not. infwal ) go to 49
!
!     do 48 j = 1, mth-1
!     do 48 i = 1, mth
!
!     ii = mth*j + i
!     jj = 2*mth*j + i
!     grdgre(ii) = grdgre(jj)
!
!  48 continue
!  49 continue
      if ( checkd ) &
      write(outmod,8010) (grdgre(i,1),i=1,9)
 8010 format ( 1x, "grdgre = ",/, 10e11.4 )
!
!.....gaussian elimination.
 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
      lmax2 = 2*jmax1
!
!     if ( abs(n)  >    1.0E-5_r8   ) go to 190
      if ( (abs(n)  >    1.E-5_r8 )  .OR.  infwal  .OR.  (ishape  >  10) ) &
      go to 190
!
!.......   add arbitrary constant to grdgre to make matrix nonsingular.
!
!!$      cn0 =  1.0_r8 
!!$!
!!$      write (outmod,8180) cn0
!!$ 8180 format ( /, 5x, "constant,cn0, added to grdgre = ", f7.3 )
!!$!
!!$      do 180 i = 1, mth12
!!$      do 180 j = 1, mth12
!!$      grdgre(i,j) = grdgre(i,j) + cn0
!!$  180 continue
!!$!
  190 continue
!
      if (checkd ) then
         CALL pstmatwrt( grdgre   ,nths2,nths2,mth12,lmax2  , "grdgre            " )
         CALL pstmatwrt( grri   ,nths2,nfm21,mth12,lmax2  , "grri              " )
      end if

!aplet  28.1_r8 .97     CALL pstleqt1f(grdgre,lmax2,mth12,nths2,grri,idgt,workl,ier)
! aplet 6/2/97
! replaced imsl routine by lapack routine dgesvx (sgesvx for crays)

allocate( zamat(mth12,mth12), zsource(mth12,lmax2), zx(mth12,lmax2) )
allocate ( ZAF(mth12,mth12),jpv(mth12), zr(mth12), zc(mth12),ZFERR(lmax2), &
 ZBERR(lmax2), ZWORK(4*mth12),JWORK(mth12), stat=iok    )


if (iok /= 0 ) print *,'error allocation in vacuum routine'

zamat(1:mth12,1:mth12) = grdgre(1:mth12,1:mth12)
zaf(1:mth12,1:mth12) = zamat(1:mth12,1:mth12)
zsource(1:mth12,1:lmax2) = grri(1:mth12,1:lmax2)
zx(1:mth12,1:lmax2) = grri(1:mth12,1:lmax2)


CALL dgesv( MTH12, LMAX2, ZAMAT, mth12, jpv,  ZX, mth12, ierr)

!!$if ( checkd ) then
!!$print *,'zrcond = ',zrcond
!!$print *, 'zsource = ', zsource
!!$print *, 'zx = ', zx
!!$zerr = sum( sum( matmul(zaf,zx) - zsource, 2),1 )
!!$print *,'error estimate = ',zerr
!!$endif


grri(1:mth12,1:lmax2) = zx(1:mth12,1:lmax2)

if (checkd ) then
   CALL pstmatwrt(grri,nths2,nfm12,mth,lmax2,"grri after   dgesvx" )
end if
!


!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!!$      write ( outmod,8050 ) ier
!!$      write ( 6, 8050 ) ier
!!$ 8050 format ( 1x, " ier in  = leqt1f/gelg", i5 )
!
      vacmat(1:jmax1,1:jmax1) = zero
      vacmti(1:jmax1,1:jmax1) = zero
!
      arr(1:jmax1,1:jmax1) = zero
      aii(1:jmax1,1:jmax1) = zero
      ari(1:jmax1,1:jmax1) = zero
      air(1:jmax1,1:jmax1) = zero
!
!.....integrate over observer points with trapezoidal rule.
!
      airs =  0.0_r8 
      aris =  0.0_r8 
      vacis =  0.0_r8 
!
!
      do     l1 = 1, jmax1
         !
         ll1 = l1 - 1 + lmin(1)
         al1nq = ll1 - nq
         do     l2 = 1, jmax1
            !
            ll = l2 - 1 + lmin(1)
            al2nq = ll - nq

            do    i = 1, mth
               !
               theta = (i-1)*dth
               elth = ll*theta
               elthnq = elth + nq*delt(i)
               !
               sinlth = sin(elthnq)
               coslth = cos(elthnq)
               !
               !     il11 = (l1-1) * mth12 + i
               !     il12 = jmax1 * mth12 + il11
               !
               arr(l2,l1) = arr(l2,l1) + dth*coslth*grri(i,l1)
               aii(l2,l1) = aii(l2,l1) + dth*sinlth*grri(i,jmax1+l1)
               ari(l2,l1) = ari(l2,l1) + dth*coslth*grri(i,jmax1+l1)
               air(l2,l1) = air(l2,l1) + dth*sinlth*grri(i,l1)
               !
               !
            end do
            !
            airs = airs + abs ( air(l2,l1) )
            aris = aris + abs ( ari(l2,l1) )
            !
            !
         end do
      end do
!
    if ( checkd ) then
      CALL pstmatwrt ( arr, j1v,j2v, jmax1,jmax1, "arr                 " )
      CALL pstmatwrt ( aii, j1v,j2v, jmax1,jmax1, "aii                 " )
      CALL pstmatwrt ( ari, j1v,j2v, jmax1,jmax1, "ari                 " )
      CALL pstmatwrt ( air, j1v,j2v, jmax1,jmax1, "air                 " )
    end if
!
!.....calculate the vacuum matrix.
!
    clockwise = +1.0_r8
    if( (xinf(2)-xinf(1))*(zinf(2)-zinf(1)) < 0.0_r8 ) clockwise = -1.0_r8
!
      do 300 j1 = 1, jmax1
!
      l1 = j1 + lmin(1) - 1
      alnq1 = l1 - nq
!
      do 300 j2 = 1, jmax1
!
      l2 = j2 + lmin(1) - 1
      alnq2 = l2 - nq
!     
!     introduced the (l-nq) terms here and removed them from addvac 3/93
!
      vacmat(j1,j2) = clockwise*  alnq1*alnq2 * ( arr(j1,j2) + aii(j1,j2) )
      vacmti(j1,j2) = clockwise*  alnq1*alnq2 * ( air(j1,j2) - ari(j1,j2) )
!
!
      vacis = vacis + abs ( vacmti(j1,j2) )
!
  300 continue
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
if ( checkd ) then
      CALL pstmatwrt(vacmat,j1v,j2v,jmax1,jmax1,"vacmat         " )
      CALL pstmatwrt(vacmti,j1v,j2v,jmax1,jmax1,"vacmti         " )
end if
!
!
!....... calculate asymmetries...
!
      arrd =  0.0_r8 
      aiid =  0.0_r8 
      vacrd =  0.0_r8 
      vacrso =  0.0_r8 
      vacrsd =  0.0_r8 
!
      do 110 l1 = 1, jmax1
      vacrsd = vacrsd + vacmat(l1,l1)
      do 110 l2 = 1, l1
      arrd = arrd + abs ( arr(l1,l2) - arr(l2,l1) )
      aiid = aiid + abs ( aii(l1,l2) - aii(l2,l1) )
      vacrd = vacrd + abs ( vacmat(l1,l2) - vacmat(l2,l1) )
      if ( l1  /=  l2 ) &
 vacrso = vacrso + abs ( vacmat(l1,l2)+vacmat(l2,l1) )
  110 continue
!
      asymv1 = vacrd / vacrso
      asymv2 = vacrd / ( vacrso + vacrsd )
!
      write ( outmod,112 ) airs,aris,vacis,arrd,aiid, vacrd, &
                     vacrsd, vacrso, asymv1,asymv2
!
!
  112 format( /,1x," airs=",1pe12.4," aris=",1pe12.4," vacis=",1pe12.4, &
   /,1x," arrd=",1pe12.4," aiid=",1pe12.4," vacrd=",1pe12.4,  &
  /,1x," vacrsd=",1pe12.4, " vacrso=",1pe12.4, &
   /,1x," asymv1=",1pe12.4," asymv2=",1pe12.4 )
!
!
      write ( outmod, 500 ) n,q, nj
  500 format (//,1x,'n,q, nj = ',2e12.4,i5,/ )
!..... symmetrize.
      do l1 = 1, jmax1
         do l2 = l1, jmax1
            vacmat(l1,l2) =  0.5_r8  * ( vacmat(l1,l2)+vacmat(l2,l1) )
            vacmti(l1,l2) =  0.5_r8  * ( vacmti(l1,l2)-vacmti(l2,l1) )
         enddo
      enddo
      do l1 = 1, jmax1
         do l2 = l1, jmax1
            vacmat(l2,l1) = vacmat(l1,l2)
            vacmti(l2,l1) =-vacmti(l1,l2)
         enddo
      enddo
!
      if ( checkv ) then
!
!     find eigenvalues of vacuum matrix....
!
      do 9195 i = 1, nfv0
      work1(i) = zero
 9195  continue
!
      kk = 1
      do 9200 i = 1, jmax1
        l1 = i + lmin(1) - 1
        alnq1 = l1 - nq
          do 9200 j = 1, i
            l2 = j + lmin(1) - 1
            alnq2 = l2 - nq
            work1(kk) =   alnq1 * alnq2 *vacmat(i,j)
            kk = kk + 1
 9200    continue
      CALL psteigen ( work1,work, jmax1, 0 )
      kk = 1
      do 9210 i = 1, jmax1
      ii = ( i*  (i+1) ) / 2
      j = (i-1)*jmax1 + 1
      jj = j + jmax1 - 1
!
      write(outmod,9003) i, work1(ii), ( work(j1), j1=j,jj )
!
!
 9003 format(1x,/,1x, "eigenvector no.",i4 &
 , "  eigenvalue = ",e12.5,/,             10(1x,10e13.6,/ ) )
!
 9210 continue
!
      endif
!

deallocate(zamat, zsource, zx)
deallocate(ZAF,jpv, zr, zc,ZFERR, &
 ZBERR, ZWORK,JWORK)

      return
  999 CALL psterrmes ( outpest,'vacuum' )
 9300 format(//," conducting wall on plasma surface ",//)
 9400 format(//," conducting wall at infinity ",//)
 9500 format(//," conformal conducting wall in vacuum at b =",e12.5,//)
 9501 format(//," conducting wall geometry: ",/, &
                " x shift/R0      a = ", f8.4, /, &
                " z shift/R0      b = ", f8.4, /, &
                " half-width/R0  aw = ", f8.4, /, &
                " half-height/R0 bw = ", f8.4, /, &
                " triangularity     = ", f8.4, /, &
                " squareness        = ", f8.4, /)
      end

