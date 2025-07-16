subroutine eqi_tbdy_circs(zdata,inum,rc,zc,rad,zR0,zZ0)

  !  add in circular limiters

  implicit NONE

  !  input:
  integer inum                      ! no. of circular limiters

  !  output:
  real*8 zdata(4,inum)              ! 4 data items / circular limiter

  !  input:
  real*8 rc(inum),zc(inum)          ! (R,Z) of circle centers
  real*8 rad(inum)                  ! circle radii

  real*8 zR0,zZ0                    ! mag. axis of plasma (for sign test)

  !----------------------------
  integer i
  real*8 zradtest
  !----------------------------

  do i=1,inum
     zdata(1,i)=rc(i)
     zdata(2,i)=zc(i)
     zdata(3,i)=rad(i)

     !  set "zsign"

     zradtest=sqrt((zR0-rc(i))**2+(zZ0-zc(i))**2)
     if(zradtest.gt.rad(i)) then
        zdata(4,i)=-1           ! plasma outside the circle
     else
        zdata(4,i)=+1           ! plasma inside the circle
     endif
  enddo

  !  for a given point (R,Z), let r=sqrt((R-Rc)**2+(Z-Zc)**2)
  !  then d = zsign*(r-rad); if d<0, inside plasma-allowed region;
  !  if d>0 then inside plasma-excluded region

end subroutine eqi_tbdy_circs

!-------------------------------------------------------------
subroutine eqi_tbdy_lines(zdata,ilines,rl,zl,thl,zR0,zZ0)

  !  add line limiters, including "sanity" limiters

  implicit NONE

  !  input:
  integer ilines                    ! no. of user specified line limiters

  !  output:
  real*8 zdata(6,ilines)            ! limiter info

  !  input:
  real*8 rl(ilines),zl(ilines)      ! point on line
  real*8 thl(ilines)                ! line orientation, DEGREES

  real*8 zR0,zZ0                    ! plasma mag. axis

  !-----------------------------------------------

  integer i

  real*8 zpi,zconv                      ! radians/degree

  real*8 zrl(4),zzl(4),zthl(4)
  !-----------------------------------------------

  zpi=3.1415926535897931D+00

  zconv=zpi/180

  do i=1,ilines
     call eqi_tbdy_line0(zdata(1,i),rl(i),zl(i),thl(i)*zconv, zR0,zZ0)
  enddo

end subroutine eqi_tbdy_lines

!----------------------------------------------------------------
subroutine eqi_tbdy_line0(zdata,rlin,zlin,ztheta,zR0,zZ0)

  !  set up one line limiter

  implicit NONE

  !  output:
  real*8 zdata(6)                   ! limiter data vector

  !  input:
  real*8 rlin,zlin                  ! point on line
  real*8 ztheta                     ! orientation of line

  real*8 zR0,zZ0                    ! point inside plasma
  real*8 ztest

  !------------------------

  zdata(1)=rlin                     ! point on line...
  zdata(2)=zlin
  zdata(3)=-sin(ztheta)             ! normal vector to line...
  zdata(4)=cos(ztheta)
  zdata(5)=zdata(1)*zdata(3)+zdata(2)*zdata(4)
                                        ! signed distance, (0,0) to line

  ztest = zR0*zdata(3)+zZ0*zdata(4) - zdata(5)
  if(ztest.lt.0.0d0) then
     zdata(6)=1                 ! inside plasma
  else
     zdata(6)=-1
  endif

end subroutine eqi_tbdy_line0
