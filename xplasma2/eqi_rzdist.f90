subroutine eqi_rzdist(ivec,zR,zZ,zphi,zdist, &
     ilim_type,iwork,r8work, &
     ilpt,rlim,zlim,philim, ierr)

  !  **vectorized**
  !  compute distance from mechanical limiter: .gt.0 means outside 
  !  .lt.0 means inside, # gives actual distance in (R,Z) units
  !  relative to nearest limiting object.

  !  if switch is set, also return coordinates of nearest contact
  !  point on the limiting object.

  implicit NONE

  !  input:

  integer ivec                      ! vector dimension
  real*8 zR(ivec),zZ(ivec),zphi(ivec) ! location

  integer :: ilim_type              ! type of limiter
  !  100= contour limiter; 101= "circles and lines" limiter

  integer :: iwork(6)               ! limiter data pointers
  real*8 :: r8work(*)               ! limiter data

  integer ilpt                      ! =1 to calculate nearest limiter pt

  !  output:
  real*8 zdist(ivec)                ! signed distance from bdy/limiter
  integer ierr                      ! completion code, 0=OK

  real*8 rlim(*),zlim(*),philim(*)  ! (1:ivec) if used:
                 ! nearest limiter contact pt

  !---------------------------------------------------
  real*8 zdist0(ivec)
  real*8 zrlim(ivec),zzlim(ivec),zphilim(ivec)
  integer i,iv,j,iadr,iadz,imin
  real*8 zrnex,zznex,zrcur,zzcur
  real*8 ztest(ivec),zdsum,zdis2(ivec),zdl2nex
  real*8 zrsave(ivec),zzsave(ivec)

  integer :: ilines,icircs,iadrl,iadrc

  real*8, parameter :: czero = 0.0d0
  !---------------------------------------------------

  ierr = 0

  if((ilim_type.lt.100).or.(ilim_type.gt.101)) then
     ierr=9999
     return
  endif

  iadrl =iwork(1)
  ilines=iwork(4)

  iadrc =iwork(2)
  icircs=iwork(5)

  if(ilim_type.eq.101) then

     !  circles and lines...

     if(ilpt.eq.1) then
        philim(1:ivec)=zphi
     endif

     zdist=-1.0d30

     iadr=iadrc
     do i=1,icircs
        call eqi_rzdist_c(ivec,r8work(iadr),zR,zZ,zdist0,ilpt,zrlim,zzlim)
        do iv=1,ivec
           if(zdist0(iv).gt.zdist(iv)) then
              zdist(iv)=zdist0(iv)
              if(ilpt.eq.1) then
                 rlim(iv)=zrlim(iv)
                 zlim(iv)=zzlim(iv)
              endif
           endif
        enddo
        iadr=iadr+4
     enddo

     iadr=iadrl
     do i=1,ilines
        call eqi_rzdist_l(ivec,r8work(iadr),zR,zZ,zdist0,ilpt,zrlim,zzlim)
        do iv=1,ivec
           if(zdist0(iv).gt.zdist(iv)) then
              zdist(iv)=zdist0(iv)
              if(ilpt.eq.1) then
                 rlim(iv)=zrlim(iv)
                 zlim(iv)=zzlim(iv)
              endif
           endif
        enddo
        iadr=iadr+6
     enddo

  else if(ilim_type.eq.100) then

     !  axisymmetric piecewise linear limiter description

     if(ilpt.eq.1) then
        philim(1:ivec)=zphi
     endif

     !  find nearest segment where closest approach of point to infinite
     !  extension of segment is on the segment

     zdis2=1.0d30
     do i=1,ilines-1
        iadr=iadrl+6*(i-1)
        zrcur=r8work(iadr)
        zzcur=r8work(iadr+1)
        zrnex=r8work(iadr+6)
        zznex=r8work(iadr+7)
        zdl2nex=(zrnex-zrcur)**2+(zznex-zzcur)**2

        call eqi_rzdist_l(ivec,r8work(iadr),zR,zZ,zdist,1,zrlim,zzlim)

        do iv=1,ivec
           ztest(iv)=max( &
                ((zRlim(iv)-zRcur)**2+(zZlim(iv)-zZcur)**2), &
                ((zRlim(iv)-zRnex)**2+(zZlim(iv)-zZnex)**2))

           if(ztest(iv).le.zdl2nex) then

              !  qualifying segment

              if(abs(zdist(iv)).lt.abs(zdis2(iv))) then
                 zdis2(iv)=zdist(iv)
                 zrsave(iv)=zrlim(iv)
                 zzsave(iv)=zzlim(iv)
              endif
           endif
        enddo

     enddo

     do iv=1,ivec
        zdist(iv)=zdis2(iv)
        zdis2(iv)=zdist(iv)*zdist(iv)

        !  find nearest endpoint

        imin=0
        do j=2,ilines
           iadr=iadrl+6*(j-1)
           iadz=iadr+1
           zdist0(iv)=(zR(iv)-r8work(iadr))**2+(zZ(iv)-r8work(iadz))**2
           if(zdist0(iv).lt.zdis2(iv)) then
              imin=j
              zdis2(iv)=zdist0(iv)
           endif
        enddo

        if(imin.eq.0) then

           !  segment was closest

           zrlim(iv)=zrsave(iv)
           zzlim(iv)=zzsave(iv)

        else

           !  endpoint was closest

           zdsum=czero
           iadr=iadrl+6*(imin-2)
           call eqi_rzdist_l(1,r8work(iadr),zR(iv),zZ(iv),zdist(iv), &
                0,zrlim(iv),zzlim(iv))
           zdsum=zdsum+zdist(iv)
           iadr=iadr+6
           call eqi_rzdist_l(1,r8work(iadr),zR(iv),zZ(iv),zdist(iv), &
                0,zrlim(iv),zzlim(iv))
           zdsum=zdsum+zdist(iv)

           if(zdsum.ge.czero) then
              zdist(iv)=sqrt(zdis2(iv))
           else
              zdist(iv)=-sqrt(zdis2(iv))
           endif

           zrlim(iv)=r8work(iadr)
           zzlim(iv)=r8work(iadr+1)

        endif
        if(ilpt.eq.1) then
           rlim(iv)=zrlim(iv)
           zlim(iv)=zzlim(iv)
        endif

     enddo

  endif

end subroutine eqi_rzdist

!------------------------
!  xplasma internal routines for eqi_rzdist
!------------------------

subroutine eqi_rzdist_c(ivec,zdata,zR,zZ,zdist,ilpt,zrlim,zzlim)

  !  distance parameter from a circular limiter

  implicit NONE

  !  input:
  integer ivec                      !vector dimension
  real*8 zdata(4)                   ! Rc,Zc,rad,sign
  real*8 zR(ivec),zZ(ivec)          ! points from which to eval distance

  integer ilpt                      ! =1 to calculate nearest limiter pt

  !  output:
  real*8 zdist(ivec)                ! (zR,zZ) distance from limiter

  real*8 zrlim(ivec),zzlim(ivec)    ! nearest limiter contact points

  !  if zdist.lt.0.0 the point is inside
  !  if zdist.gt.0.0 the point is outside
  !-------------------------

  integer iv
  real*8 zrad,zdr,zdz

  !-------------------------

  do iv=1,ivec
     zdr=zR(iv)-zdata(1)
     zdz=zZ(iv)-zdata(2)
     zrad=sqrt(zdr**2+zdz**2)
     zdist(iv)=zdata(4)*(zrad-zdata(3))

     if(ilpt.eq.1) then
        zrlim(iv)=zdata(1)+zdata(3)*(zdr/zrad)
        zzlim(iv)=zdata(2)+zdata(3)*(zdz/zrad)
     endif
  enddo

end subroutine eqi_rzdist_c

!------------------------
subroutine eqi_rzdist_l(ivec,zdata,zR,zZ,zdist,ilpt,zrlim,zzlim)

  !  distance parameter from a line limiter

  implicit NONE

  !  input:
  integer ivec                      ! vector dimension
  real*8 zdata(6)                   ! Rl,Zl,<normal-vec>,param,sign
  real*8 zR(ivec),zZ(ivec)          ! points in question

  integer ilpt                      ! =1 to calculate nearest limiter pt

  !  output:
  real*8 zdist(ivec)                ! (zR,zZ) distance from limiter
  
  real*8 zrlim(ivec),zzlim(ivec)    ! nearest limiter contact point

  !  if zdist.lt.0.0 the point is inside
  !  if zdist.gt.0.0 the point is outside
  !-------------------------

  integer iv
  real*8 ztest

  !-------------------------

  do iv=1,ivec
     ztest = zdata(3)*zR(iv)+zdata(4)*zZ(iv) - zdata(5)
     zdist(iv) = zdata(6)*ztest

     if(ilpt.eq.1) then
        zrlim(iv)=zR(iv)-zdata(3)*ztest
        zzlim(iv)=zz(iv)-zdata(4)*ztest
     endif
  enddo

end subroutine eqi_rzdist_l
