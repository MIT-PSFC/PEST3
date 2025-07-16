subroutine eqi_bdfind_circ(ivec,imask,zR,zZ,iseg, &
     ztheta_near,zdist, &
     nsrch1,zRsurf,zZsurf,zTheta,krzrad,rzrad,ierr)

  ! internal routine to compute distance based on circles interpolation
  ! approximation: each sequence of the consecutive points in
  ! {zRsurf,zZsurf} is fit to a circle, which defines a distance to each
  ! target point.  This gives two distances-- one for each end of the
  ! segment to which the target point is closest.  A single distance is
  ! chosen as an interpolation between these two, linearly according to 
  ! the average of 2 locations of intersection of lines from target points
  ! to the fit circles's centers.  I.e. if these intersections are 75% of
  ! the way to segment endpoint #2, then 
  !    d = 75%*[distance to circle fit to endpoint #2] +
  !        25%*[distance to circle fit to endpoint #1]
  ! is used.

  implicit NONE

  integer, intent(in) :: ivec   ! size of vector of target points
  logical, intent(in) :: imask(ivec) ! skip vector elements j if imask(j) TRUE

  real*8, intent(in) :: zR(ivec) ! target points' R coordinates
  real*8, intent(in) :: zZ(ivec) ! target points' Z coordinates
  integer, intent(in) :: iseg(ivec)    ! segment to use (cf eqi_bdfind)

  !  the following are intent(inout) so that outputs corresponding to masked
  !  points are unmodified:

  real*8, intent(inout) :: ztheta_near(ivec)   ! estimated nearest theta
  real*8, intent(inout) :: zdist(ivec)         ! estimated distance

  integer, intent(in) :: nsrch1  ! # of search segment end points

  real*8, intent(in) :: zRsurf(0:nsrch1+1),zZsurf(0:nsrch1+1) ! segment end pts
  !  augmented at start and end so that the triple involving an end point and
  !  its nearest neighbours is always available.

  real*8, intent(in) :: zTheta(0:nsrch1+1) ! poloidal angle at each seg end pt.

  integer, intent(inout) :: krzrad(nsrch1) ! status for each end point:
  !  0 -- circle fit not done
  !  1 -- circle fit done
  !  2 -- circle fit deemed "degenerate"; treat end pt triple as colinear
  !    status may be updated here...

  real*8, intent(inout) :: rzrad(nsrch1,3) ! circle data [Rcen,Zcen,radius]
  !  for each segment end point... may be updated here... 

  integer, intent(out) :: ierr    ! error status code, 0=OK
  !----------------------------------------------------------------
  integer :: i,jsrch1,jsrch2
  real*8 :: d1,d2,thet1,thet2,thetav

  real*8, parameter :: TWO = 2.0d0
  real*8, parameter :: HALF = 0.5d0
  real*8, parameter :: C4TH = 0.25d0
  real*8, parameter :: EPS10= 1.0d-10
  real*8, parameter :: ZERO = 0.0d0
  !----------------------------------------------------------------

  ierr=0

  do i=1,ivec
     if(imask(i)) cycle

     jsrch1=iseg(i)
     jsrch2=jsrch1+1

     if(krzrad(jsrch1).eq.0) call circ_fit(jsrch1)
     if(krzrad(jsrch2).eq.0) call circ_fit(jsrch2)

     call dist(jsrch1,zR(i),zZ(i),d1,thet1)
     call dist(jsrch2,zR(i),zZ(i),d2,thet2)

     thetav=HALF*(thet1+thet2)

     ztheta_near(i)=thetav

     zdist(i)=(d1*(zTheta(jsrch2)-thetav) + d2*(thetav-zTheta(jsrch1)))/ &
          (zTheta(jsrch2)-zTheta(jsrch1))
  enddo

  contains

    subroutine circ_fit(j)

      !  circle fit the three points j-1:j+1 on the search contour
      !  finding the 1 point that is equidistant from 3 consecutive
      !  boundary points is pretty quick-- reduces to a linear equation.

      !  on exit, the sign of the radius indicates if the center of
      !  the circle is displaced towards the inside or the outside of
      !  the plasma.

      !  if the radius is too large a flag is set and the sequence of 3
      !  points is treated as being colinear.

      integer, intent(in) :: j

      !-------------------
      integer :: jp1,jm1
      real*8 :: Rp,Zp,Rm,Zm,Rhalf,Zhalf,Rj,Zj,delR,delZ,del2
      real*8 :: delRj,delZj,delc2
      real*8 :: zrhs,zlamd,zcoeff
      !-------------------

      jp1=j+1
      Rp=zRsurf(jp1); Zp=zZsurf(jp1)

      jm1=j-1
      Rm=zRsurf(jm1); Zm=zZsurf(jm1)

      Rhalf=HALF*(Rm+Rp)
      Zhalf=Half*(Zm+Zp)

      Rj=zRsurf(j)
      Zj=zZsurf(j)

      delR = Rp-Rm
      delZ = Zp-Zm
      del2 = delR*delR + delZ*delZ

      delRj=(Rhalf-Rj)
      delZj=(Zhalf-Zj)
      delc2= delRj*delRj + delZj*delZj

      zrhs = delc2 - del2*C4TH
      zcoeff = TWO*(delZ*delRj-delR*delZj)
      if(abs(zcoeff).lt.1.0d-6*abs(zrhs)) then
         !  treat as colinear segment
         krzrad(j)=2
         
      else
         krzrad(j)=1
         zlamd = zrhs/zcoeff
         rzrad(j,1)=Rhalf - delZ*zlamd
         rzrad(j,2)=Zhalf + delR*zlamd

         ! **testing** only one is needed: all should have the same value
         ! (verified dmc May 2006)

         rzrad(j,3)=sqrt((Rj-rzrad(j,1))**2+(Zj-rzrad(j,2))**2)
         ! rzrad(j,3)=sqrt((Rp-rzrad(j,1))**2+(Zp-rzrad(j,2))**2)
         ! rzrad(j,3)=sqrt((Rm-rzrad(j,1))**2+(Zm-rzrad(j,2))**2)

         if(zlamd.lt.0.0d0) rzrad(j,3)=-rzrad(j,3)
      endif

    end subroutine circ_fit

    subroutine dist(j,Rpt,Zpt,d,th)

      !  distance and theta of nearest approach estimate
      !  relative to point j on the search contour

      integer, intent(in) :: j       ! search contour point
      real*8, intent(in) :: Rpt,Zpt  ! point from which distance is wanted

      real*8, intent(out) :: d       ! distance estimate
      !  signed: +1 for outside & -1 for inside the contour...

      real*8, intent(out) :: th      ! theta (interp. zTheta of search contour)

      !--------------------------
      integer :: j1,j2,jm1,jp1
      real*8 :: dmsq,dpsq,th1,th2
      real*8 :: delR,delZ,znum,zdenom,zlamd,Rx,Zx,d1,d2,R1,Z1,R2,Z2,zdir
      real*8 :: Rc,Zc
      real*8 :: dc,rad,signtest
      !--------------------------
      jm1=j-1
      jp1=j+1

      if(krzrad(j).eq.2) then

         ! pts treated as colinear (very rare in practice)

         dmsq = (Rpt-zRsurf(jm1))**2 + (Zpt-zZsurf(jm1))**2
         dpsq = (Rpt-zRsurf(jp1))**2 + (Zpt-zZsurf(jp1))**2

         ! use segment with nearer end point...
         if(dmsq.lt.dpsq) then
            j1=jm1
            j2=j
         else
            j1=j
            j2=jp1
         endif
         th1=zTheta(j1)
         th2=zTheta(j2)
         R1=zRsurf(j1); Z1=zZsurf(j1)
         R2=zRsurf(j2); Z2=zZsurf(j2)

         delR = R2-R1
         delZ = Z2-Z1

         ! find intersection of line through (Rpt,Zpt) normal to the segment
         ! with the segment itself or a linear extension thereof

         zdenom = delR*delR + delZ*delZ  ! guaranteed non-zero
         znum = (Z1-Zpt)*delR - (R1-Rpt)*delZ

         zlamd = znum/zdenom

         ! this is the intercept

         Rx = Rpt - delZ*zlamd
         Zx = Zpt + delR*zlamd

         ! this is the distance: + for outside, - for inside the plasma surface
         ! (this sign convention applies because the surfaces are drawn
         ! counter-clockwise around the plasma).

         if(zlamd.gt.ZERO) then
            d=sqrt((Rx-Rpt)**2+(Zx-Zpt)**2)
         else
            d=-sqrt((Rx-Rpt)**2+(Zx-Zpt)**2)
         endif

         ! find which end point is closer.  Extrapolate if past an endpoint;
         ! this should be rare/small.  Define output theta by interpolation 
         ! using ratio of distances from the end points to the intercept.

         if((delR*(Rx-R1)+delZ*(Zx-Z1)).ge.ZERO) then
            d1=sqrt((Rx-R1)**2+(Zx-Z1)**2)
         else
            d1=-sqrt((Rx-R1)**2+(Zx-Z1)**2)  ! this should be rare
         endif

         d2=sqrt(delR**2+delZ**2)-d1

         th = (d2*th1 + d1*th2)/(d1+d2)

      else

         ! normal case -- use circle fit data associated with point j

         ! signtest: usually (Rpt,Zpt) is on the same side of the circle
         ! center as the bdy points, but... not always...

         Rc = rzrad(j,1)  ! circle center...
         Zc = rzrad(j,2)

         signtest = (Rpt-Rc)*(zRsurf(j)-Rc) + &
              (Zpt-Zc)*(zZsurf(j)-Zc)
         if(signtest.lt.0.0d0) then
            signtest=-1.0d0
         else
            signtest=1.0d0
         endif

         dc = sqrt((Rpt-Rc)**2+(Zpt-Zc)**2)*signtest

         rad = rzrad(j,3)

         if(rad.lt.ZERO) then
            d = -(rad+dc)  ! - if inside surface, + if outside
         else
            d = dc-rad     ! - if inside surface, + if outside
         endif
         
         if(abs(dc).lt.EPS10*rad) then
            th=zTheta(j)
         else
            Rx=Rc + abs(rad)*(Rpt-Rc)/dc
            Zx=Zc + abs(rad)*(Zpt-Zc)/dc
         endif

         zdir = (zRsurf(jp1)-zRsurf(j))*(Rx-zRsurf(j)) + &
              (zZsurf(jp1)-zZsurf(j))*(Zx-zZsurf(j))

         if(zdir.ge.ZERO) then
            j1=j
            j2=jp1
         else
            j1=jm1
            j2=j
         endif
         th1=zTheta(j1)
         th2=zTheta(j2)
         R1=zRsurf(j1); Z1=zZsurf(j1)
         R2=zRsurf(j2); Z2=zZsurf(j2)

         d1=(Rx-Rc)*(Zx-Z1)-(Zx-Zc)*(Rx-R1)
         d2=(Rx-Rc)*(Z2-Zx)-(Zx-Zc)*(R2-Rx)

         th = (d2*th1 + d1*th2)/(d1+d2)

      endif
    end subroutine dist

end subroutine eqi_bdfind_circ
