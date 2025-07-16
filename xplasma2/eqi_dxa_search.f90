subroutine eqi_dxa_search(iin,iout,nvec,distv,atanv, &
     inth,th,inrho,rho,distf,atanf,atan0, &
     inthb,thb,distb,atanb, &
     rhov,thv,nregion)

  implicit NONE

  !  **utility routine** for fast inverse map based on distance from axis...

  !  GIVEN profiles
  !     distf(ith,irho) = sqrt((R(ith,irho)-Raxis)**2+(Z(ith,irho)-Zaxis)**2)
  !     atanf(ith,irho) = atan2((Z(ith,irho)-Zaxis),(R(ith,irho)-Raxis))
  !  and a vector of points at locations (distv(...),atanv(...)) wrt the 
  !  magnetic axis, return estimates of the corresponding rho and theta
  !  locations.
  
  !  Do not return an answer for points that are closer to the magnetic axis
  !  than the closest points in the mapped region;
  
  !  Do not return an answer for points (j) with nregion(j).gt.0.

  !  If iout>0: set nregion(j)=iout if point (j) is beyond the mapped region;
  !  If iout=0: force point(j) to be considered inside the mapped region.

  !  For points found or considered to be within the mapped region, set
  !  nregion(j)=iin and output the correspond rho & theta values in rhov(j)
  !  and thv(j).

  !  For each target point this is a binary search in both of two dimensions
  !  taking advantage of the fact that d(theta,rho) is monotonic increasing
  !  vs. rho along the ray from the magnetic axis to the target point, and,
  !  the adjusted arctan(theta,rho) is monotonic increasing in theta at fixed
  !  rho.

  !---------------------------
  !  PASSED:

  integer, intent(in) :: iin  ! marker code for items in range
  integer, intent(in) :: iout ! marker code for items beyond range
  ! (force into range if iout=0)

  integer, intent(in) :: nvec ! size of vector of points sought

  real*8, intent(in) :: distv(nvec)  ! distances from mag axis to points sought
  real*8, intent(in) :: atanv(nvec)  ! angle to points sought (wrt mag axis)

  integer, intent(in) :: inth        ! theta (poloidal angle) grid of profiles
  real*8, intent(in) :: th(inth)     ! theta grid values

  integer, intent(in) :: inrho       ! rho (radial flux coordinate) grid
  real*8, intent(in) :: rho(inrho)   ! rho grid values

  real*8, intent(in) :: distf(inth,inrho)  ! distance profile   (from mag axis)
  real*8, intent(in) :: atanf(inth,inrho)  ! arctangent profile (wrt mag axis)
  !  above btw 0 and 2pi, monotonic over 1:inth at each rho, relative to:
  real*8, intent(in) :: atan0(inrho)       ! arctan of th(1) half ray

  !  in some configurations xplasma has higher resolution information at/beyond
  !  a boundary.  So, these inputs are provided s.t. they can be used.
  !  the grid thb(1:inthb) is a superset of th(1:inth).  Any th(j) is matched
  !  by some thb(k).  If inthb=inth, then the grids are identical

  integer, intent(in) :: inthb       ! theta grid @bdy (can be different)
  real*8, intent(in) :: thb(inthb)   ! theta grid values @bdy
  real*8, intent(in) :: distb(inthb) ! distance profile @bdy
  real*8, intent(in) :: atanb(inthb) ! atan info @bdy

  !  arctan((Z(ith,irho)-Zaxis)/(R(ith,irho)-Raxis)) = 
  !                                        atanf(ith,irho) - atan0(irho)

  ! "inout" because some values might not get set on this call...

  real*8, intent(inout) :: rhov(nvec)      ! rho grid values found
  real*8, intent(inout) :: thv(nvec)       ! theta grid values found
  integer, intent(inout) :: nregion(nvec)  ! location code (iin/iout)

  !---------------------------
  !  LOCAL:
  integer :: ivec,nseek
  logical, dimension(:), allocatable :: isdone
  integer, dimension(:), allocatable :: jrho,jrho_min,jrho_max
  real*8, dimension(:), allocatable :: th_minrho,th_maxrho
  real*8, dimension(:,:), allocatable :: dsave,thsave
  real*8, parameter :: C2PI = 6.2831853071795862D+00
  real*8, parameter :: CPI =  3.1415926535897931D+00
  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE = 1.0d0
  real*8 :: d,ata,df,thf,dup,ddown,hd,thup,thdown,dmax,dmin,dfac
  integer :: ith,ithb,krho,krho_new,krho_old,iter
  logical :: step_up
  !---------------------------
  !  work vectors

  allocate(isdone(nvec),jrho(nvec),jrho_min(nvec),jrho_max(nvec))
  allocate(th_minrho(nvec),th_maxrho(nvec))
  allocate(thsave(inrho,nvec),dsave(inrho,nvec))

  nseek=0
  do ivec=1,nvec
     jrho_min(ivec)=0
     jrho_max(ivec)=inrho
     jrho(ivec)=inrho
     if(nregion(ivec).eq.0) then
        isdone(ivec)=.FALSE.
        nseek=nseek+1
     else
        isdone(ivec)=.TRUE.
     endif
  enddo

  iter=0
  do
     if(nseek.eq.0) exit
     iter=iter+1
     if(iter.gt.2*inrho) then
        call errmsg_exit(' fatal error, algorithm failure, eqi_dxa_search.f90')
     endif

     do ivec=1,nvec
        if(isdone(ivec)) cycle

        krho=jrho(ivec)
        d=distv(ivec)

        ata=atanv(ivec)-atan0(krho)
        if(ata.lt.ZERO) ata = ata+C2PI
        if(ata.gt.C2PI) ata = C2PI    ! now guaranteed btw 0 and 2pi

        ith = iget_ith(ata,atanf(1:inth,krho),inth)
        ithb = ith

        if(krho.eq.inrho) then
           if(inthb.gt.inth) then
              ithb = iget_ith(ata,atanb,inthb)
           endif
           df = hfun(ata,0,ithb,distb,atanb,inthb)
           thf = hfun(ata,1,ithb,thb,atanb,inthb)

        else
           df = hfun(ata,0,ith,distf(1:inth,krho),atanf(1:inth,krho),inth)
           thf = hfun(ata,1,ith,th,atanf(1:inth,krho),inth)

        endif

        dsave(krho,ivec)=df
        thsave(krho,ivec)=thf

        if(d.gt.df) then
           if(krho.eq.inrho) then
              ! beyond last surface
              isdone(ivec)=.TRUE.
              nseek=nseek-1
              thv(ivec)=thf
              if(iout.eq.0) then
                 df=d  ! force to boundary
                 rhov(ivec)=rho(inrho)
                 nregion(ivec)=iin
              else
                 rhov(ivec)=rho(inrho)*(d/df)
                 nregion(ivec)=iout
              endif
              cycle
           else
              jrho_min(ivec)=krho
              step_up=.TRUE.
           endif

        else if(d.lt.df) then
           if(krho.eq.1) then
              nregion(ivec)=0
              rhov(ivec)=rho(1)
              thv(ivec)=0
              isdone(ivec)=.TRUE.
              nseek=nseek-1
              cycle
           else
              jrho_max(ivec)=krho-1
              step_up=.FALSE.
           endif

        else
           ! exact match
           isdone(ivec)=.TRUE.
           nseek=nseek-1
           nregion(ivec)=iin
           rhov(ivec)=rho(krho)
           thv(ivec)=thf
           cycle

        endif
           
        if(jrho_min(ivec).eq.jrho_max(ivec)) then
           isdone(ivec)=.TRUE.
           nseek=nseek-1

           krho=jrho_min(ivec)

           dup=dsave(krho+1,ivec)
           ddown=dsave(krho,ivec)

           hd=dup-ddown

           thup=thsave(krho+1,ivec)
           thdown=thsave(krho,ivec)

           !  check for branch cut
           if(thup.gt.thdown+CPI) then
              thup=thup-C2PI
           else if(thup.lt.thdown-CPI) then
              thup=thup+C2PI
           endif

           nregion(ivec)=iin
           rhov(ivec)=(rho(krho)*(dup-d)+rho(krho+1)*(d-ddown))/hd
           thv(ivec)=(thdown*(dup-d)+thup*(d-ddown))/hd
           cycle
        endif

        if(step_up) then
           krho_old=jrho_max(ivec)+1
           dmax=dsave(krho_old,ivec)
           dfac=(d-df)/(dmax-df)
           krho_new = krho+dfac*(krho_old-krho)
           krho_new = max(krho+1,min(krho_old-1,krho_new))

        else
           krho_old=jrho_min(ivec)
           if(krho_old.eq.0) then
              krho_new=max(1,min(krho-1,krho/4))
           else
              dmin=dsave(krho_old,ivec)
              dfac=(d-df)/(dmin-df)
              krho_new=krho+dfac*(krho_old-krho)
              krho_new=max(krho_old+1,min(krho-1,krho_new))
           endif

        endif

        jrho(ivec) = krho_new

     enddo
  enddo

  contains

    integer function iget_ith(ata,atanp,inthp)

      real*8, intent(in) :: ata  ! arctan parameter
      real*8, intent(in), dimension(:) :: atanp  ! arctan grid
      integer, intent(in) :: inthp  ! grid size

      integer :: ith,ithmin,ithmax

      ! binary search in theta
      ithmin=1
      ithmax=inthp

      do 
         if(ata.ge.atanp(ithmax-1)) then
            ith=ithmax-1
            exit
         else if(ata.le.atanp(ithmin+1)) then
            ith=ithmin
            exit
         else
            ithmin=ithmin+1
            ithmax=ithmax-1

            if(ithmax.eq.ithmin+1) then
               ith=ithmin
               exit
            endif

            ith = ithmin + (ithmax-ithmin)* &
                 (ata-atanp(ithmin))/ &
                 (atanp(ithmax)-atanp(ithmin))
            if(ith.ge.ithmax) ith=ithmax-1   ! btw ithmin and (ithmax-1) now.

            if(ata.gt.atanp(ith+1)) then
               ithmin=ith+1

            else if(ata.lt.atanp(ith)) then
               ithmax=ith

            else
               exit
            endif
         endif
      enddo

      iget_ith = ith

    end function iget_ith

    real*8 function hfun(ata,ibranch,ith,funp,atanp,inthp)
      !  Hermite interpolation within a theta zone on a surface
      real*8, intent(in) :: ata  ! arctan of target
      integer, intent(in) :: ibranch  ! 1 if funp has a branch cut
      integer, intent(in) :: ith  ! theta zone (must contain ata)
      real*8, dimension(:), intent(in) :: funp  ! function vs. atan param.
      real*8, dimension(:), intent(in) :: atanp  ! atan param. grid
      integer, intent(in) :: inthp  ! size of funp & atanp

      real*8 :: d0,d1,d0p,d1p,h

      d0=funp(ith)
      d1=funp(ith+1)
      if(ith.gt.1) then
         d0p=(funp(ith+1)-funp(ith-1))/ &
              (atanp(ith+1)-atanp(ith-1))
      else
         d0p=(funp(ith+1)+ibranch*C2PI-funp(inthp-1))/ &
              (atanp(ith+1)+C2PI-atanp(inthp-1))
      endif
      if(ith+1.lt.inthp) then
         d1p=(funp(ith+2)-funp(ith))/ &
              (atanp(ith+2)-atanp(ith))
      else
         d1p=(funp(2)+ibranch*C2PI-funp(ith))/ &
              (atanp(2)+C2PI-atanp(ith))
      endif

      h = (atanp(ith+1)-atanp(ith))

      hfun = hermfa(ata-atanp(ith),h,d0,d0p,d1,d1p)

    end function hfun

    real*8 function hermfa(avar,h,f0,f0p,f1,f1p)

      !  Hermite cubic evaluation...
      real*8, intent(in) :: avar,f0,f0p,f1,f1p,h

      real*8 :: hinv,hi2,f00,f01,f10,f11,avari

      hinv = ONE/h
      hi2 = hinv*hinv

      avari=h-avar

      f00 = hi2*avari*avari*(3-2*avari*hinv)*f0
      f01 = hi2*avar*avar*(3-2*avar*hinv)*f1
      f10 = hinv*avari*avari*(1-avari*hinv)*f0p
      f11 = hinv*avar*avar*(avar*hinv-1)*f1p

      hermfa = f00+f01+f10+f11

    end function hermfa

end subroutine eqi_dxa_search
