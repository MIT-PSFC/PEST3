!
!
!      3.1_r8    conducting wall specification.
!.......................................................................
!!$      SUBROUTINE pstvacda2(infwal,xwal,zwal)
!!$      SUBROUTINE pstvacda2(infwal)
      SUBROUTINE pstvacda2
!.......................................................................
!      conducting wall specification
!
!         a, b, aw, bw, dw and sw from namelist in modin.
!
!        if a greater than equal 10, then cond. wall at infinity
!        if a less than equal - 10._r8  then wall at a fixed distance
!        given by b.
!        if a=- 10._r8  the wall is constructed numerically by forming the
!        local outwards normal at each interface grid point and moving
!        out b*plasma radius (on outboard midplane). the indentation
!        can be  adjusted by varying bw from 0 (straight inboard side)
!        to 1 (similar 'bean' to plasma )
!
!        if a<-10 the wall is formed by joining the points obtained
!        by extending each ray from the magnetic axis to the interface
!        grid points by the amount b*max plasma radius.then the points
!        are moved along the wall so that each point is close to it's
!        corresponding interface grid point.
!
!        if not then set wall by the expn.
!
!         xw = a + xma + aw*(1+sw*sin(theta)**2) *cos(theta + dw*sin(theta))
!         zw = zma + bw*(1+sw*cos(theta)**2)*sin(theta)
!
!        all lengths are scaled to the larger of the max. vertical
!        distance ,max. horizontal distance.
!
!      in the third option 'a' defines the shift of the wall shape 'axis',
!        'dw' the triangularity, 'aw' the minor radius, ' bw/aw' 
!        the elongation and 'sw' the squareness.
!
!      versions:   chance, manickam    1978
!      modified:   grimm   dec, 82._r8 
!
!...........................................................
!
 USE pstcom

 USE l22com

 USE mtrik1
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IFIRST
      real*8 BSAV
      real*8 XSHIFT
      real*8 ZSHIFT
      real*8 XMAX
      real*8 XMIN
      real*8 ZMAX
      real*8 ZMIN
      real*8 PLRAD
      real*8 XMAJ
      real*8 DELTA
      real*8 TOLD
      real*8 THETO
      INTEGER I
      real*8 RR
      real*8 RO
      real*8 THE
      real*8 XMX
      INTEGER INSIDE
      INTEGER J
      INTEGER J1
      real*8 TT
      INTEGER JPLUS
      INTEGER JMNUS
      real*8 ZDTH
      real*8 ZDXINF
      real*8 ZDZINF
      real*8 ZTX
      real*8 ZTZ
      real*8 ZMOD
      real*8 ZFOLD
      real*8 TTOLD
      real*8 ZXTOLD
      real*8 ZZTOLD
      real*8 ZXTMIN
      real*8 ZZTMIN
      real*8 ZCOEF
      INTEGER JITER
      INTEGER K
      real*8 XMX1
      real*8 ZMZ1
      real*8 ZZ
      real*8 ZZP
      real*8 ZZPP
      real*8 ZF
      real*8 ZFP
      real*8 ZFPP
      real*8 DELTZ
      real*8 ZINTER


!
!
      data ifirst/0/
!
!!$ LOGICAL	 	 	 :: infwal
 LOGICAL	 	 	 :: lfix
 LOGICAL	 	 	 :: insect
!!$ REAL*8, DIMENSION(1) 	 :: xwal
!!$ REAL*8, DIMENSION(1) 	 :: zwal
 REAL*8, DIMENSION(nths) 	 :: xpp
 REAL*8, DIMENSION(nths) 	 :: zpp
 REAL*8, DIMENSION(nths) 	 :: ww1
 REAL*8, DIMENSION(nths) 	 :: ww2
 REAL*8, DIMENSION(nths) 	 :: ww3
 REAL*8, DIMENSION(nths) 	 :: thet
 REAL*8, DIMENSION(3) 	 :: tabx
 REAL*8, DIMENSION(3) 	 :: tabz
 REAL*8, DIMENSION(nths) 	 :: xwnew
 REAL*8, DIMENSION(nths) 	 :: zwnew
 INTEGER, DIMENSION(2) 	 :: iop
      if(ifirst  ==  0) bsav = b
!     


      lfix = .false.
      if( a  >   - 10._r8 ) lfix=.true.
      if( a  >=    10._r8 ) then
         infwal = .true.
         return
      end if
!
      xshift = a*xma
      zshift = b*xma
!
!      check b. should be greater than  0.0_r8 
!
      if( b  <    0.0_r8   .AND.  a <= - 10._r8  ) then
         write(6,1100)b
         write(outpst,1100)b
         write(outmod,1100)b
         CALL psterrmes(outpst,'vacdat')
1100     format(" **error** b is less than zero=",e12.5)
         infwal = .true.
         return
      end if
!
      xmax = maxval( xinf(1:mth ) )
      xmin = minval( xinf(1:mth ) )
      zmax = maxval( zinf(1:mth ) )
      zmin = minval( zinf(1:mth ) ) 

      plrad = ( xmax - xmin ) /  2.0_r8 
      xmaj = ( xmax + xmin ) /  2.0_r8 
!
      if ( a  ==  - 10._r8  )  go to 200
!
      write(outmod,1101) xma, zma, xshift, zshift
 1101  format( / ' xma = ',f10.6,' zma = ', f10.6, &
                 ' xshift = ',f10.6,' zshift = ', f10.6/)
!
      aw = aw * xma
      bw = bw * xma
      delta = b * max( xmax-xma, xma-xmin )
!
      told = atan2( zinf(1)-zma, xinf(1)-xma )
      theto=  0._r8 
      do 100 i = 1, mth
!
! thet must increase monotonically from 0 to 2*pi
!
      rr = (xinf(i)-xma)**2 + (zinf(i)-zma)**2
      ro = sqrt(rr)
      the = atan2( zinf(i)-zma, xinf(i)-xma )
      zt = abs(the - told)
      if( zt >    0.8_r8 *twopi ) zt = zt - twopi
!
      thet(i) = theto + abs(zt)
      told = the
      theto = thet(i)
!     
      if(lfix) go to 50
!
      ro = ro + delta
      xwal(i) = xma + ro * cos ( the )
      zwal(i) = zma + ro * sin ( the )
      go to 60
!
!      employ todds equation for a wall
!
   50 continue
!
      xwal(i) = xshift + xma +  &
                aw * (1 + sw*sin(the)**2/2)* cos(the + dw * sin(the))
      zwal(i) = zshift + zma + bw* (1 + sw*cos(the)**2/2) * sin(the)
!
   60 continue
!
 100  continue
!
      go to 300
!
  200 continue
!
!      conforming wall, numerically.
!
      xmx =   0.5_r8  * ( xmin + xmax )
      b = bsav *  0.5_r8  * ( xmax - xmin )
      write(outmod,7359)b
      write(6,7359)b
 7359 format(' conforming wall at ',f5.2,' meters')
!
      told = atan2( zinf(1)-zma, xinf(1)-xma )
      theto=  0._r8 
      do 220 i = 1, mth
!
! thet must increase monotonically from 0 to 2*pi
!
      rr = (xinf(i)-xma)**2 + (zinf(i)-zma)**2
      ro = sqrt(rr)
      the = atan2( zinf(i)-zma, xinf(i)-xma )
      zt = abs(the - told)
      if( zt >    0.8_r8 *twopi ) zt = zt - twopi
!
      thet(i) = theto + abs(zt)
      told = the
      theto = thet(i)
!
      xwal(i) = xinf(i) + b*cos(the)
      zwal(i) = zinf(i) + b*sin(the)
!
  220 continue
!
 300  continue
!
! complete
!
      thet(mth1) = thet(1) + twopi
      thet(mth2) = thet(2) + twopi
      xwal(mth1) = xwal(1)
      xwal(mth2) = xwal(2)
      zwal(mth1) = zwal(1)
      zwal(mth2) = zwal(2)
!
!      now check for wall intersecting the plasma surface
!      and realign wall points to lie along the direction
!      normal to the plasma surface at ( xinf(i), zinf(i) )
!
!
      insect = .false.
      inside = 0
!
! check that thet is monotonically increasing
!
      j = 0
      do i = 1, mth1
      j1 = i+1
      if( thet(j1)  <  thet(i) ) then
      j = j + 1
      write(6,*) &
   ' *** thet nonmonotonic thet(',i,') = ',thet(i), &
   ' thet(',j1,') = ', thet(j1)
      write(outmod,*) &
   ' *** thet nonmonotonic thet(',i,') = ',thet(i), &
   ' thet(',j1,') = ', thet(j1)
       end if
      end do
      if ( j >  0 ) stop 'ERROR in vacuum: non monotonic theta'
!
!     write(outmod,*)' wall point before alignment'
!     write(outmod,1460)
!     do i = 1,mth
!     write(outmod,1480)i,thet(i),xinf(i),zinf(i),
!    .                     xwal(i),zwal(i)
!     end do
!
!    fit spline to get xw,zw at new intervals of theta
!    periodic boundary conditions
!
      iop(1) = 4
      iop(2) = 4
!
xpp =  0._r8 
zpp =  0._r8 

      call pstspl1d1(mth1,thet,xwal,xpp,iop,1,ww1,ww2,ww3)
      call pstspl1d1(mth1,thet,zwal,zpp,iop,1,ww1,ww2,ww3)

!
! realign xwal and zwal perpendicularly to the tangent vector at (xinf,zinf)
!
      do 125 i=2,mth1
      xs = xinf(i)
      zs = zinf(i)
      xt = xwal(i)
      zt = zwal(i)
      tt = thet(i)
!
      jplus = i+1
      jmnus = i-1
!
! compute tangential unit vector [ ztx ztz ]
!
      zdth =   thet(jplus) - thet(jmnus)
      zdxinf = xinf(jplus) - xinf(jmnus)
      zdzinf = zinf(jplus) - zinf(jmnus)
      ztx  = zdxinf/zdth
      ztz  = zdzinf/zdth
      zmod = sqrt( ztx**2 + ztz**2 )
      ztx = ztx/zmod
      ztz = ztz/zmod
!
!      now use a newton-raphson iteration scheme to minimize  0.5_r8 *(t*d)**2
!      where d = [ (xt-xs)   (zt-zs) ]
!
      zfold  =  0.5_r8 * (ztx*(xt-xs) + ztz*(zt-zs))**2
      ttold = tt
      zxtold = xt
      zztold = zt
      zxtmin = xt
      zztmin = zt
      zcoef =  1._r8 
!
      jiter = 100 
      do 120 k=1, jiter
      call pstspl1d2(mth1,thet,xwal,xpp,1,tt,tabx)
      call pstspl1d2(mth1,thet,zwal,zpp,1,tt,tabz)
      xt = tabx(1)
      zt = tabz(1)
      xmx1 = xt - xs
      zmz1 = zt - zs
!
      zz =   ztx*xmx1 + ztz*zmz1
      zzp =  ztx*tabx(2) + ztz*tabz(2)
      zzpp = ztx*tabx(3) + ztz*tabz(3)
      zf =  0.5_r8 * zz**2
      zfp = zz*zzp
      zfpp = zzp**2 + zz*zzpp
!
      deltz = -zfp/zfpp
!
! check if new zf is smaller, if not take previous value
!
      zxtmin = xt
      zztmin = zt
!
      if( zf >  zfold+ 1.E-10_r8  ) then
      zcoef = zcoef/ 2._r8 
      tt = ttold
      zxtmin = zxtold
      zztmin = zztold
      go to 120
      end if
!
! check for presence of a maximum
!
      if( zfpp  <   0.05_r8  ) then
!
! use linear extrapolation to zero
!
      deltz = -zf/zfp
      else
!
! reset coef if previous value is not identical
!
      if( abs((zf-zfold)/(zf+ 1.E-8_r8 )) >    1.E-8_r8  ) zcoef =  1._r8 
      end if
!
      if ( abs(deltz)  >    1.0_r8  ) deltz = deltz /  10.0_r8 
! aplet  7.8_r8 .95      if ( abs(deltz)  >   pye/ 3._r8 ) deltz = pye/ 3._r8 
      if(abs(deltz) <  1.0E-4_r8    .AND.  abs(zfp)  <   1.0E-4_r8  ) go to 124
!
      ttold = tt
      zxtold = xt
      zztold = zt
      zfold = zf
!
      tt = tt + deltz*zcoef
 3     if (tt <  0._r8  ) then
      tt = tt + twopi
      go to 3
      end if
 4     if (tt >  twopi ) then
      tt = tt - twopi
      go to 4
      end if
!
  120 continue
!
      write(6,1430)i,xs,zs,xt,zt,tt,zf,zfp,zfpp
      write(outpst,1430)i,xs,zs,xt,zt,tt,zf,zfp,zfpp
!      CALL psterrmes(outpst,'vacdat')
!
!
 1430 format(" newton scheme for finding the nearest wall point", &
 /,   "  did not converge ", &
 /,   " i=",i4," (x,z)_plasma=(",2e15.8,") (x,z)_wall=(",2e15.8,")"&
 /  " theta=",1e10.3/" zf=",1e15.8," zfp=",1e15.8," zfpp=",1e15.8)
  124 continue
!
      xt = zxtmin
      zt = zztmin
      xwnew(i) = xt
      zwnew(i) = zt
!
! test for intersection
!
      zinter =  (xt-xs)*(xs-xma) &
        +  (zt-zs)*(zs-zma)
      if( zinter  <  + 1.E-6_r8  ) then
      insect = .true.
      inside = inside + 1
      end if
!
  125 continue
!
      do i = 2, mth1
      xwal(i) = xwnew(i)
      zwal(i) = zwnew(i)
      end do
!
      xwal(1)    = xwal(mth1)
      xwal(mth2) = xwal(2)
      zwal(1)    = zwal(mth1)
      zwal(mth2) = zwal(2)
!
      if(ifirst  >   0) return
      ifirst = 1
!
! aplet  8.4_r8 .94
!      CALL pstdrawc2(zwal,xwal,zinf,xinf,1,mth1,"z","x",xmx,xma,zma)
!
      if( imode ==  1 ) then
      write(outmod, 1470)
      do i = 1, mth1
      write(outmod, 1480) i, xinf(i), zinf(i), xwal(i), zwal(i)
      end do
      end if
!
!
 1470 format(//' plasma and wall boundaries:' &
 / 1x,' i ',1x,5x,'xinf',  11x,'zinf',11x,'xwal',11x,'zwal')
!
 1480 format(1x,i3,5(1x,e14.7))
!
! aplet  8.4_r8 .94
!
      if(.not. insect) return
      write(6,1450)inside
      write(outpst,1450)inside
      write(outmod,1450)inside
      CALL psterrmes(outpst,'vacdat')
!
 1460 format(1x," i ",5x,"thet",10x,"xinf",10x,"zinf", &
 10x,"xwal",  10x,"zwal")
 1450 format(" there are at least ",i3," wall points in the plasma")
!
      stop 'ERROR in vacuum'
      end

