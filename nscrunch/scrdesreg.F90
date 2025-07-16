SUBROUTINE scrdesreg_clean(rmb_in,zmb_in,nmb_in,mpolin,icon,ioflag, &
     fsq,lasym0,lzsmlonly,nmbout,rmbout,zmbout,ier)
  !
  !  call scrDESCUR -- but first:
  !   a) remove duplicate points
  !   b) using scrdesreg subroutine, regularize the contour using equal
  !  steps in theta as defined by elongation-corrected atan2((Z-Z0),(R-R0))
  !
  !  arguments are as in scrDESCUR.f90
  !  as in scrDESCUR, ier=1 means NORMAL exit
  !
  implicit NONE
  !
  integer, intent(in) :: nmb_in  ! no. of contour (R,Z) data pts
  real*8, intent(in) :: rmb_in(nmb_in),zmb_in(nmb_in) 
  ! closed contour (R,Z) data pts
  !   expecting rmb(1)=rmb(nmb) and zmb(1)=zmb(nmb)

  integer, intent(in) :: mpolin  ! no. of moments (0:mpolin-1)
  integer, intent(in) :: icon    ! debug arg -- .gt.0 to have scrDESCUR
  !   write a debug file on fortran LUN icon

  logical, intent(in) :: ioflag  ! .TRUE. for extra debug output messages

  real*8, intent(out) :: fsq  ! output measure of error (fft contour deviation)

  logical, intent(in) :: lasym0  ! .TRUE. -- updown assymmetry
  logical, intent(in) :: lzsmlonly  ! .TRUE. -- Z=Z0+Z1*sin(theta) only
  !   (all other Z moments to be set to zero)

  integer, intent(in) :: nmbout  ! array dimension for output arrays
  !   note FULL DIMESIONS are (0:NMBOUT,2)

  !  R cos & sin moments
  real*8, intent(out) :: rmbout(0:nmbout,2),zmbout(0:nmbout,2)

  integer, intent(out) :: ier    ! completion code, OK = 1

  !------------------------------------
  !  local:
  real*8 :: rmb(nmb_in),zmb(nmb_in)
  real*8 :: zzr,zzrp,zzz,zzzp,zdl
  integer :: i,nmb,lunmsg
  real*8, parameter :: ZERO = 0.0d0
  !------------------------------------

  lunmsg = 6

  nmb=1
  rmb(1)=rmb_in(1)
  zmb(1)=zmb_in(1)
  zzr=rmb(1)
  zzz=zmb(1)

  do i=2,nmb_in
     zzrp=zzr
     zzzp=zzz
     zzr=rmb_in(i)
     zzz=zmb_in(i)
     zdl=sqrt((zzr-zzrp)**2+(zzz-zzzp)**2)
     if(zdl.eq.ZERO) then
        continue  ! exact overlap of points...
     else
        nmb=nmb+1
        rmb(nmb)=rmb_in(i)
        zmb(nmb)=zmb_in(i)
     endif
  enddo

  if(nmb.lt.nmb_in) then
     write(lunmsg,*) ' %scrdesreg_clean: ',nmb_in-nmb, &
          ' overlapping points removed.'
  endif

  call scrdesreg(rmb,zmb,nmb,mpolin,icon,ioflag, &
       fsq,lasym0,lzsmlonly,nmbout,rmbout,zmbout,ier)

end SUBROUTINE scrdesreg_clean

SUBROUTINE scrdesreg (rmb,zmb,nmb,mpolin,icon,ioflag, &
     fsq,lasym0,lzsmlonly,nmbout,rmbout,zmbout,ier)
!
!  call scrDESCUR -- but first -- regularize the contour using equal
!  steps in theta as defined by elongation-corrected atan2((Z-Z0),(R-R0))
!
!  arguments are as in scrDESCUR.f90
!  as in scrDESCUR, ier=1 means NORMAL exit
!
  implicit NONE
!
  integer, intent(in) :: nmb  ! no. of contour (R,Z) data pts
  real*8, intent(in) :: rmb(nmb),zmb(nmb) ! closed contour (R,Z) data pts
  !   expecting rmb(1)=rmb(nmb) and zmb(1)=zmb(nmb)

  integer, intent(in) :: mpolin  ! no. of moments (0:mpolin-1)
  integer, intent(in) :: icon    ! debug arg -- .gt.0 to have scrDESCUR
  !   write a debug file on fortran LUN icon

  logical, intent(in) :: ioflag  ! .TRUE. for extra debug output messages

  real*8, intent(out) :: fsq  ! output measure of error (fft contour deviation)

  logical, intent(in) :: lasym0  ! .TRUE. -- updown assymmetry
  logical, intent(in) :: lzsmlonly  ! .TRUE. -- Z=Z0+Z1*sin(theta) only
  !   (all other Z moments to be set to zero)

  integer, intent(in) :: nmbout  ! array dimension for output arrays
  !   note FULL DIMESIONS are (0:NMBOUT,2)

  !  R cos & sin moments
  real*8, intent(out) :: rmbout(0:nmbout,2),zmbout(0:nmbout,2)

  integer, intent(out) :: ier    ! completion code, OK = 1

  !--------------------------------------------
  !  local
  integer ith0,ith1
  integer nmb2,i,j,jprev,idir,lunmsg,inter
  real*8, dimension(:), allocatable :: rmb2,zmb2,rmbx,zmbx,zthx,zdx,zdla
  real*8 :: zl,zdl,zzl,zzr,zzrp,zzz,zzzp,R0,Z0,ztmp
  real*8 :: zrmax,zrmin,zzmin,zzmax,zaovb,zaovb2,zintrp
  real*8 :: zdthmin,zdthmax,zrat,zbranch
  real*8 :: zddth1,zddth2,zdth,zdelth,za,zb,zth,zd
  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: HALF = 0.5d0
  real*8, parameter :: ONE = 1.0d0
  real*8, parameter :: PI = 3.1415926535897931D+00
  !--------------------------------------------
#ifdef __DEBUG

  real*8 :: r0atan,z0atan,shif2pi
  real*8 :: thatan(nmb),datan(nmb)
  real*8, dimension(:), allocatable :: thb2,db2

  integer, parameter :: nthspl=301
  real*8 :: thvmec(nthspl),thvma(nthspl),dvmec(nthspl), &
       rvmec(nthspl),zvmec(nthspl)
  real*8 :: thvmcc(nthspl-1),dthth(nthspl-1),ddth(nthspl-1)

  real*8 :: sinm(mpolin-1),cosm(mpolin-1)

#endif
  !--------------------------------------------

  lunmsg = 6

  ier=1
  rmbout=ZERO
  zmbout=ZERO

  !  double length array for regularized polar angle step
  nmb2=nmb*2
  allocate(rmb2(nmb2),zmb2(nmb2))

  !  extended contour -- 1 point duplicated on each end
  allocate(rmbx(0:nmb+1),zmbx(0:nmb+1),zthx(0:nmb+1),zdx(0:nmb+1))
  allocate(zdla(1:nmb+1))
  rmbx(1:nmb)=rmb; zmbx(1:nmb)=zmb
  rmbx(0)=rmb(nmb-1); rmbx(nmb+1)=rmb(2)
  rmbx(1)=HALF*(rmbx(1)+rmbx(nmb)); rmbx(nmb)=rmbx(1)
  zmbx(0)=zmb(nmb-1); zmbx(nmb+1)=zmb(2)
  zmbx(1)=HALF*(zmbx(1)+zmbx(nmb)); zmbx(nmb)=zmbx(1)

  !  define Z0 estimate: line average
  !  also error check -- all pts distinct
  !  also find min & max

  zl=ZERO
  zzl=ZERO
  zzr=rmbx(1); zrmin=zzr; zrmax=zzr
  zzz=zmbx(1); zzmin=zzz; zzmax=zzz
  do i=2,nmb
     zzrp=zzr
     zzzp=zzz
     zzr=rmbx(i); zrmin=min(zrmin,zzr); zrmax=max(zrmax,zzr)
     zzz=zmbx(i); zzmin=min(zzmin,zzz); zzmax=max(zzmax,zzz)
     zdl=sqrt((zzr-zzrp)**2+(zzz-zzzp)**2)
     zdla(i)=zdl
     if(zdl.eq.ZERO) then
        ier=99
        write(lunmsg,*) '?scrdesreg: 2 consecutive points overlap exactly.'
        exit
     endif
     zl=zl+zdl
     zzl=zzl+zdl*(zzz+zzzp)*HALF
  enddo
  zdla(1)=zdla(nmb)
  zdla(nmb+1)=zdla(2)

  if(ier.ne.1) return

  Z0 = zzl/zl
  zaovb = (zrmax-zrmin)/(zzmax-zzmin)  ! inverse ellipticity approx.
  zaovb2 = zaovb**2

  !  define R0 estimate: half way btw midplane intercepts
  !  via piecewise linear estimation

  inter=0
  zzz=zmbx(1)
  zzr=ZERO
  do i=2,nmb
     zzzp=zzz
     zzz=zmbx(i)
     if(zzz.eq.Z0) cycle
     if((zzz-Z0)*(zzzp-Z0).le.ZERO) then
        zintrp=(Z0-zzzp)/(zzz-zzzp)
        inter=inter+1
        zzr=zzr + rmbx(i-1)*(ONE-zintrp)+rmbx(i)*zintrp
     endif
  enddo

  if(inter.ne.2) then
     write(lunmsg,*) '?scrdesreg: 2 surface/midplane intercepts not found:', &
          inter,' were found.'
     ier=100
     return
  endif

  R0=HALF*zzr

  !  define theta sequence... and distance sequence; (a/b) modified space...

  idir=0
  zbranch=ZERO
  zthx(0)=atan2(zaovb*(zmbx(0)-Z0),(rmbx(0)-R0))
  zdx(0)=sqrt(zaovb2*(zmbx(0)-Z0)**2+(rmbx(0)-R0)**2)
  do i=1,nmb+1
     zthx(i)=atan2(zaovb*(zmbx(i)-Z0),(rmbx(i)-R0)) + zbranch
     zdx(i)=sqrt(zaovb2*(zmbx(i)-Z0)**2+(rmbx(i)-R0)**2)
     if((zthx(i)-zthx(i-1)).gt.PI) then
        zbranch=zbranch-2*pi
        zthx(i)=zthx(i)-2*pi
     endif
     if((zthx(i)-zthx(i-1)).lt.-PI) then
        zbranch=zbranch+2*pi
        zthx(i)=zthx(i)+2*pi
     endif
     if(idir.eq.0) then
        if(zthx(i).gt.zthx(i-1)) then
           idir=1
        else
           idir=-1
        endif
        zdthmin=(zthx(i)-zthx(i-1))/zdla(i)
        zdthmax=zdthmin
     else
        zdthmin=min(zdthmin,(zthx(i)-zthx(i-1))/zdla(i))
        zdthmax=max(zdthmax,(zthx(i)-zthx(i-1))/zdla(i))
     endif
  enddo

  if((zdthmin*zdthmax).le.ZERO) then
     write(lunmsg,*) '?scrdesreg: d(theta)/dl changes direction.'
     ier=101
     return
  endif

  zrat=max(abs(zdthmin),abs(zdthmax))/min(abs(zdthmin),abs(zdthmax))
  if(zrat.gt.nmb) then
     write(lunmsg,*) '?scrdesreg: too much variation in |d(theta)/dl|'
     ier=102
     return
  endif

  !  reverse theta & D(theta) order if necessary

  if(zthx(0).gt.zthx(nmb+1)) then
     ith0=0
     ith1=nmb+1
     do
        ztmp=zthx(ith0)
        zthx(ith0)=zthx(ith1)
        zthx(ith1)=ztmp
        ztmp=zdx(ith0)
        zdx(ith0)=zdx(ith1)
        zdx(ith1)=ztmp
        ith0=ith0+1
        ith1=ith1-1
        if(ith0.ge.ith1) exit
     enddo
  endif

  !  OK  have D(theta) based on original contour; new contour is formed
  !      using equal steps in theta; end points are shared with the old

  rmb2(1)=rmbx(1)
  rmb2(nmb2)=rmbx(nmb)
  zmb2(1)=zmbx(1)
  zmb2(nmb2)=zmbx(nmb)

  !  Hermite interpolation will be used, based on numerical derivatives
  !  dD/dtheta at original contour (R,Z) pts

  j=1
  zddth1=(zdx(j+1)-zdx(j-1))/(zthx(j+1)-zthx(j-1))
  zddth2=(zdx(j+2)-zdx(j))/(zthx(j+2)-zthx(j))
  call deriv(j,zddth1)
  call deriv(j+1,zddth2)
  zdelth=zthx(j+1)-zthx(j)
  za=(2*(zdx(j)-zdx(j+1))/zdelth + (zddth1+zddth2))/(zdelth*zdelth)
  zb=(3*(zdx(j+1)-zdx(j))/zdelth - (zddth2+2*zddth1))/zdelth
  do i=2,nmb2-1
     zth=((i-1)*zthx(nmb)+(nmb2-i)*zthx(1))/(nmb2-1)
     jprev=j
     do
        if(zth.le.zthx(j+1)) exit
        j=j+1
     enddo
     if(j.ne.jprev) then
        zddth1=(zdx(j+1)-zdx(j-1))/(zthx(j+1)-zthx(j-1))
        zddth2=(zdx(j+2)-zdx(j))/(zthx(j+2)-zthx(j))
        call deriv(j,zddth1)
        call deriv(j+1,zddth2)
        zdelth=zthx(j+1)-zthx(j)
        za=(2*(zdx(j)-zdx(j+1))/zdelth + (zddth1+zddth2))/(zdelth*zdelth)
        zb=(3*(zdx(j+1)-zdx(j))/zdelth - (zddth2+2*zddth1))/zdelth
     endif

     !  in interval j

     zdth=zth-zthx(j)

     zd=((za*zdth+zb)*zdth+zddth1)*zdth+zdx(j)
     zzr = zd*cos(zth)
     zzz = zd*sin(zth)/zaovb

     rmb2(i)=zzr+R0
     zmb2(i)=zzz+Z0
  enddo

#ifdef __DEBUG
  !  call r8_grafc2(rmb,zmb,nmb, rmb2,zmb2,nmb2, 'm', 'm', &
  !       'scrDESREG debug plot','zmb(rmb), zmb2(rmb2)','should be similar')
  z0atan=zmb(1)
  r0atan=sum(rmb(1:nmb))/nmb
  call bdtan(nmb,rmb,zmb,thatan,datan)
  !  allocate(thb2(nmb2),db2(nmb2))
  !  call bdtan(nmb2,rmb2,zmb2,thb2,db2)
  !  call r8_grafc2(thatan,datan,nmb, thb2,db2,nmb2, 'rad', 'm', &
  !       'scrDESREG debug plot','d(th)orig, d(th)b2','should be similar')
  !  deallocate(thb2,db2)
#endif

  call scrDESCUR (rmb2,zmb2,nmb2, mpolin, &
       icon,ioflag, &
       fsq,lasym0,lzsmlonly,nmbout,rmbout,zmbout,ier)

#ifdef __DEBUG
  do i=1,nthspl
     thvmec(i)=(i-1)*2*PI/(nthspl-1)
     call r8sincos(thvmec(i),mpolin-1,sinm,cosm)
     rvmec(i)=rmbout(0,1)
     zvmec(i)=zmbout(0,1)
     do j=1,mpolin-1
        rvmec(i)=rvmec(i)+cosm(j)*rmbout(j,1)+sinm(j)*rmbout(j,2)
        zvmec(i)=zvmec(i)+cosm(j)*zmbout(j,1)+sinm(j)*zmbout(j,2)
     enddo
  enddo
  call r8_grafc2(rmb,zmb,nmb, rvmec,zvmec,nthspl, 'm', 'm', &
       'scrDESREG debug plot','zmb(rmb), zvmec(rvmec)','should be similar')
  call bdtan(nthspl,rvmec,zvmec,thvma,dvmec)
  call r8_grafc2(thatan,datan,nmb, thvma,dvmec,nthspl, 'rad', 'm', &
       'scrDESREG debug plot','d(th)orig, d(th)vmec','should be similar')
  do i=1,nthspl-1
     thvmcc(i)=(thvmec(i)+thvmec(i+1))/2
     dthth(i)=(thvma(i+1)-thvma(i))/(thvmec(i+1)-thvmec(i))
     ddth(i)=(dvmec(i+1)-dvmec(i))/(thvmec(i+1)-thvmec(i))
  enddo
  call r8_grafx2(thvmcc,ddth,dthth,nthspl-1,'rad','mixed', &
       'scrDESREG debug plot','d(dvmec)/dth & d(thvmec)/dth',' ')
#endif
  deallocate(rmb2,zmb2,rmbx,zmbx,zthx,zdx)

  contains
    subroutine deriv(ij,zddth)
      integer, intent(in) :: ij
      real*8,intent(out) :: zddth
      real*8 :: zw1,zw2
      real*8, parameter :: Z08 = 0.8d0
      real*8, parameter :: Z09 = 0.9d0

      if(zdla(ij).le.zdla(ij+1)) then
         zw1=(ONE/zdla(ij))/(ONE/zdla(ij) + ONE/zdla(ij+1))
         if(zw1.gt.Z09) then
            zw1=ONE
         else if(zw1.gt.Z08) then
            zw1=min(ONE,(zw1+2*(zw1-Z08)))
         endif
         zw2=ONE-zw1
      else
         zw2=(ONE/zdla(ij+1))/(ONE/zdla(ij) + ONE/zdla(ij+1))
         if(zw2.gt.Z09) then
            zw2=ONE
         else if(zw2.gt.Z08) then
            zw2=min(ONE,(zw2+2*(zw2-Z08)))
         endif
         zw1=ONE-zw2
      endif
      
      zddth=( zw2*(zdx(ij+1)-zdx(ij))/(zthx(ij+1)-zthx(ij)) + &
              zw1*(zdx(ij)-zdx(ij-1))/(zthx(ij)-zthx(ij-1)) )

    end subroutine deriv

#ifdef __DEBUG    
    subroutine bdtan(nmb,rmb,zmb,thatan,datan)
      integer :: nmb
      real*8 :: rmb(nmb),zmb(nmb),datan(nmb),thatan(nmb)

      shif2pi=0
      do i=1,nmb
         datan(i)=sqrt((rmb(i)-r0atan)**2+(zmb(i)-z0atan)**2)
         thatan(i)=atan2((zmb(i)-z0atan),(rmb(i)-r0atan)) + shif2pi
         if(i.gt.1) then
            if(thatan(i).lt.thatan(i-1)-PI) then
               shif2pi=shif2pi+2*PI
               thatan(i)=thatan(i)+shif2pi
            endif
         endif
      enddo
    end subroutine bdtan
#endif

end SUBROUTINE scrdesreg
