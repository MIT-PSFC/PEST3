subroutine eqi_cbdy_pts(zdata,ipts,rlim,zlim,ierr)

  !  create the piecewise linear contour limiter data array.
  !  make the following checks:
  !    ipts.gt.2
  !    rlim(1)=rlim(ipts)
  !    zlim(1)=zlim(ipts)
  !    adjacent points are not equal
  !    sum of angles btween pieces = +/- 2*pi

  use xplasma_obj
  use eqi_rzbox_module
  implicit NONE

  !  input:
  integer ipts                      ! 1 + number of linear pieces
  real*8 rlim(ipts),zlim(ipts)      ! the (R,Z) contour point sequence

  !  output:
  integer ierr                      ! completion code, 0=OK
  real*8 zdata(6,ipts)              ! limiter data

  !        zdata(1,i) -- Rlim(i)
  !        zdata(2,i) -- Zlim(i)
  !        zdata(3,i) -- -sin(theta)
  !        zdata(4,i) --  cos(theta)
  !                             theta=atan2(dZ,dR), (dR,dZ)=displacement
  !                                                  from pt i to pt i+1
  !        zdata(5,i) -- -sin(theta)*Rlim(i)+cos(theta)*Zlim(i)
  !        zdata(6,i) -- +/- 1

  !  with elements 3-5 can quickly compute distance from any point to the
  !  infinite line which is a continuation of the given segment, and the
  !  location on that infinite line of the nearest approach, which we can
  !  quickly determine if its in the segment or out.

  !----------------------------------------------
  real*8 zR0,zZ0     ! a test point presumed inside the plasma, computed here

  real*8 zR,zZ,zdr,zdz,zth,zthp,zdth,zthsum,zsdlsum,zdl
  real*8 zsign,period
  integer i
  character*128 zmsg

  !----------------------------------------------
  real*8, parameter :: cone = 1.0d0
  real*8, parameter :: czero = 0.0d0
  real*8, parameter :: ceps4 = 1.0d-4
  real*8, parameter :: cpi = 3.1415926535897931D+00
  real*8, parameter :: c2pi= 6.2831853071795862D+00
  !----------------------------------------------

  ierr=0

  if(ipts.le.2) then
     zmsg=' '
     write(zmsg,*) ' ?eqm_cbdy:  ipts=',ipts,'.  At least 3 pts are required.'
     call xplasma_errmsg_append(sp,zmsg)
     ierr=3000
     return
  endif

  if((rlim(ipts).ne.rlim(1)).or.(zlim(ipts).ne.zlim(1))) then
     zmsg=' '
     write(zmsg,*) ' ?eqm_cbdy:  limiter contour not closed:'
     call xplasma_errmsg_append(sp,zmsg)
     zmsg=' '
     write(zmsg,*) '  R1,Z1 = ',rlim(1),zlim(1), &
          '  Rn,Zn = ',rlim(ipts),zlim(ipts)
     call xplasma_errmsg_append(sp,zmsg)
     ierr=3000
     return
  endif

  do i=1,ipts-1
     if((rlim(i).eq.rlim(i+1)).and.(zlim(i).eq.zlim(i+1))) then
        zmsg=' '
        write(zmsg,*) ' ?eqm_cbdy:  (R,Z) pts identical:  ',i,i+1, &
             '  (R,Z)=',rlim(i),zlim(i)
        call xplasma_errmsg_append(sp,zmsg)
        ierr=3000
        return
     endif
  enddo

  !  OK... loop through the points

  call find_ref_pt  ! find (R0,Z0) estimate first...

  zthsum=czero
  zsdlsum=czero

  zdr=rlim(ipts)-rlim(ipts-1)
  zdz=zlim(ipts)-zlim(ipts-1)
  zth=atan2(zdz,zdr)

  do i=2,ipts
     zthp=zth
     zR=rlim(i-1)
     zZ=zlim(i-1)
     zdr=rlim(i)-rlim(i-1)
     zdz=zlim(i)-zlim(i-1)
     zth=atan2(zdz,zdr)
     zdth=zth-zthp
!  check range
     period=c2pi
     if((zdth.lt.-cpi).or.(zdth.gt.cpi)) then
        zdth=mod(zdth+cpi,period)
        if(zdth.lt.czero) zdth=zdth+period
        zdth=zdth-cpi
     endif
     zdl=sqrt(zdr**2+zdz**2)
     zthsum=zthsum+zdth

     call eqi_tbdy_line0(zdata(1,i-1),zR,zZ,zth,zR0,zZ0)
     zsdlsum=zsdlsum+zdata(6,i-1)*zdl
  enddo

  do i=1,6
     zdata(i,ipts)=zdata(i,1)
  enddo

  !  check sum

  if(min(abs(zthsum-c2pi),abs(zthsum+c2pi)).gt.ceps4) then
     zmsg=' '
     write(zmsg,*) ' ?eqm_cbdy:  contour bend angle check sum failure!'
     call xplasma_errmsg_append(sp,zmsg)
     zmsg=' '
     write(zmsg,*) '  zthsum=',zthsum,' not +/- 2pi'
     ierr=3000
     return
  endif

  !  OK, looks good...
  !  standardize segment plasma side signs...

  if(zsdlsum.lt.czero) then
     zsign=-cone
  else
     zsign=+cone
  endif
  do i=1,ipts
     zdata(6,i)=zsign
  enddo

  contains
    subroutine find_ref_pt

      !  find a reference point that is considered to be "inside" the
      !  limiter contour.  This is defined as follows:
      !     Zmid = (Zmax+Zmin)/2
      !     R = (R2+R1)/2
      !  where R2 is the contour intercept of Zmid with the largest R value
      !  and R1 is the contour intercept of Zmid with the next largest R value.

      integer :: ii
      real*8 :: zmin,zmax,zmid,zintrp,rintrp,r1,r2

      real*8, parameter :: bigneg = -1.0d34
      real*8, parameter :: half = 0.5d0
      real*8, parameter :: one = 1.0d0

      zmin=zlim(1)
      zmax=zlim(1)
      do ii=2,ipts-1
         zmin=min(zmin,zlim(ii))
         zmax=max(zmax,zlim(ii))
      enddo

      zmid = (zmin+zmax)*half

      r1 = bigneg
      r2 = bigneg

      do ii=1,ipts-1
         if ( ((zlim(ii).lt.zmid).and.(zmid.le.zlim(ii+1))) .or. &
              ((zlim(ii).gt.zmid).and.(zmid.ge.zlim(ii+1))) ) then
            if(zlim(ii).eq.zlim(ii+1)) then
               zintrp=half
            else
               zintrp=(zmid-zlim(ii))/(zlim(ii+1)-zlim(ii))
            endif
            rintrp = rlim(ii)*(one-zintrp) + rlim(ii+1)*zintrp

            if(r1<0.d0) then
               r1=rintrp
            else if(r2<0.d0) then
               r2=rintrp
            else
               call errmsg_exit("?find_ref_pt: too many midplane crossings")
            endif
         endif
      enddo

      zR0 = half*(r1+r2)
      zZ0 = zmid

      if (zR0<0.d0) call errmsg_exit("?find_ref_pt: failed to find two limiter crossings??")
    end subroutine find_ref_pt

end subroutine eqi_cbdy_pts
