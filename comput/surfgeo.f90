!
! compute various flux surface parameters derived from the position of 
! the flux surface on the R,Z grid
!
! the midplane intercepts are defined as the major radiaus at the height of the
! centroid of the flux surface.
!
! notes from gelong:
!c
!c> Date: Wed, 22 Sep 1999 11:59:40 -0400
!c> From: Glenn Bateman <bateman@fusion.physics.lehigh.edu>
!c> 	Consider a toroidal magnetic surface (also called a flux surface).
!c> Use the following variables
!c>
!c> R_out   = major radius to the outboard edge of the flux surface
!c>           that is, the largest value of major radius anywhere on the surface
!c> R_in    = major radius to the inboard edge (ie, smallest value of major
!c> radius)
!c> R_top   = major radius to the top edge of the flux surface
!c>           that is, to the highest point on the flux surface cross section
!c> R_bottom = major radius to the lowest point on the flux surface
!c> height  = distance between the elevations of the highest point
!c>           and the lowest point
!c> R_dent_in  = major radius to the most indented inboard part of the surface
!c>           in cases where the flux surface is bean shaped
!c>           with a dent on the inboard edge
!c>           Note: R_dent_in = R_in whenever the flux surface is not indented
!c> R_dent_out = major radius to the most indented outboard part of the surface
!c>
!c> r_minor = half width = (R_out - R_in ) / 2     (usually called minor radius)
!c> r_major = ( R_out + R_in ) / 2                 (usually called major radius)
!c>
!c> elongation:	kappa = height / ( R_out - R_in )
!c>
!c> tiangularity:	delta = ( r_major - min ( R_top, R_bottom ) ) / r_minor
!c>
!c> indentation:	indent = max ( ( R_dent_in - R_in ) / rminor ,
!c> 			       ( R_dent_out - R_out ) / rminor )
!c>
!c> 	Note: for normal convex surfaces, indent = 0.0
!c> For bean-shaped surfaces with indentation on the inboard edge, indent > 0
!c> For bean-shaped surfaces with indentation on the outboard edge, indent < 0
!c>
!c> 	For surfaces with the point of the triangle facing outward,
!c> delta > 0.0, while if the triangle points inward, delta < 0.0
!c>
!c> 	For vertically elongated surfaces, kappa > 1.0
!c> For horizontally elongated surfaces, kappa < 1.0.
!c>
!c> squareness based on Holcomb Phys. Plasmas 16, 056116 (2009), see fig. 1
!c>
!c> Z_rmax = elevation where R==R_out
!c>
!c> Gyro squareness attempts to match the GYRO moment representation
!c>
!c>    R = Rmajor + Rminor*cos(theta + asin(triang)sin(theta))
!c>    Z = Zrmax  + elong*Rminor*sin(theta + square*sin(2*theta))
!c>
!c> by finding the point Rmk = R(pi/4.), interpolating from the boundary
!c> data to Zmk and then using the moment representation at Z(pi/4.) to
!c> derive the squareness.
!c>
!c>
!c>
!
subroutine surfgeo(itheta, ivec, ar, az, &
     elong, triang, indent, zmidp, rmnmidp, rmjmidp, nlim, limits, &
     idebug, mline, ier)

  implicit none

  integer, intent(in) :: itheta           ! number of poloidal points, first and last R,Z points in theta must be the same
  integer, intent(in) :: ivec             ! number of radial points
  real*8,  intent(in) :: ar(itheta, ivec) ! R(theta,rho)
  real*8,  intent(in) :: az(itheta, ivec) ! Z(theta,rho)

  real*8,  intent(out) :: elong(ivec)    ! surface elongation
  real*8,  intent(out) :: triang(ivec)   ! triangularity
  real*8,  intent(out) :: indent(ivec)   ! indentation

  real*8,  intent(out) :: zmidp (ivec)   ! height of flux surface centroid
  real*8,  intent(out) :: rmnmidp(ivec)  ! inner midplane intercept at centroid height
  real*8,  intent(out) :: rmjmidp(ivec)  ! outer midplane intercept at centroid height

  integer, intent(in)  :: nlim              ! number of limits to return
  real*8,  intent(out) :: limits(nlim,ivec) ! surface limits used for elong,triang,indent
                                          ! 1  -> R_in
                                          ! 2  -> R_out
                                          ! 3  -> R_top
                                          ! 4  -> Z_top
                                          ! 5  -> R_bottom
                                          ! 6  -> Z_bottom
                                          ! 7  -> R_dent_in
                                          ! 8  -> R_dent_out
                                          ! 9  -> R_centroid
                                          ! 10 -> Z_centroid
                                          ! 11 -> lower Holcomb squareness
                                          ! 12 -> upper Holcomb squareness
                                          ! 13 -> Z_rmax
                                          ! 14 -> Rmk_low
                                          ! 15 -> Zmk_low
                                          ! 16 -> lower Gyro squareness
                                          ! 17 -> Rmk_upp
                                          ! 18 -> Zmk_upp
                                          ! 19 -> upper Gyro squareness

  logical,      intent(in)  :: idebug     ! .true. for some debug output
  character*(*),intent(out) :: mline      ! error message
  integer,      intent(out) :: ier        ! nonzero on error
    

  integer, parameter :: R8=SELECTED_REAL_KIND(12,100)
  real*8,  parameter :: TWOPI      = 6.2831853071795862d0
  real*8,  parameter :: PI4        = TWOPI/8.d0
  real*8,  parameter :: RMINAXIS   = 1.D-6       ! threshold of zrminor/zrmajor for an axis surface

  real*8,  allocatable :: zr(:),zz(:)

  integer :: kaxis       ! index of last surface which is axis like
  integer :: i,ipt,k,iv  ! loop variable
  integer :: mlim        ! number of limits to return

  real*8  :: zrcentr, zzcentr     ! centroid
  real*8  :: tr,tz,trp,tzp        ! points in centroid loop
  real*8  :: zytest               ! >0 when segment spans midplane
  real*8  :: zrans                ! radius at midplane crossing
  real*8  :: zrmin,zrmax          ! limits of R around surface
  real*8  :: zymin,zymax          ! limits of Z around surface
  real*8  :: zrymin, zrymax       ! R values at Z limits
  real*8  :: zrminor, zrmajor     ! estimate from R major limits
  real*8  :: zyrmax               ! Y value at Rmax limit
  real*8  :: zra(3),zya(3)        ! three points around target
  real*8  :: xlimits(20)          ! actual limits on a single surface
  real*8  :: zzR0,zzAR,zzBR,zzY0  ! temps
  real*8  :: zzAY,zzBY,zx,zrx,zyx ! temps
  real*8  :: zryx,zyrx,z          ! temps
  
  logical :: o_bound                   ! temporary
  logical :: o_inside, o_outside       ! .true. if tm,tp bound the inboard and outboard sides
  real*8  :: tm,tp                     ! theta indices bounding top and bottom 
  real*8  :: tymin,tymax,trmin,trmax   ! normalized theta indices where limits are achieved
  integer :: iext(2)                   ! number of mimnimum and maximum extremum
  real*8,allocatable :: text(:,:)      ! normalized theta indices of min,max extremums  (2,itheta)
  real*8,allocatable :: rext(:,:)      ! min,max extremums (2,itheta)
  real*8  :: zrindent_in, zrindent_out ! indetation on inboard and outboard side
  real*8  :: zt                        ! normalized theta index

  integer :: ith                ! number of previous points encountered

  real*8 :: sqlow, squpp        ! lower and upper squareness
  real*8 :: th0, th, thp        ! target angle, and angle of current and prvious points
  real*8 :: a,b,ak,bk,ap,bp     ! R,Z of enclosing rectangle and current and previous points relative to ellipse center
  real*8 :: asq, bsq            ! a,b interpolated to th0
  real*8 :: ab,ad,ac,cd         ! see Hocomb, fig 1.

  integer :: imk                ! index closest to Gyro R squareness
  integer :: i0,i1,i2           ! closest indices to Millelsen R
  integer :: ii,jj              ! temp
  integer :: i47,i48            ! I/O units, if idebug=.TRUE.

  integer, parameter :: nsamp = 10 ! for debugging, surface sample #

  integer :: indxn(nsamp)         ! Surface indices for debug output
  real*8 ::  rhoxn(nsamp)         ! Corresponding rho values

  real*8 :: sqmk_low, sqmk_upp  ! Gyro lower and upper squareness
  real*8 :: klow, kupp          ! lower and upper elongation
  real*8 :: dlow, dupp          ! lower and upper triangularity
  real*8 :: rmklow, rmkupp      ! value of R at pi/4 of moment representation
  real*8 :: zmklow, zmkupp      ! value of Z corresponding to rmklow, rmkupp

  real*8 :: zsq

  !---------------------------------------------
  ! allocate temporary storage
  !
  ier   = 0
  mline = " "

  if (itheta<3) then
     mline = "?surfgeo: not enough theta points"
     ier=1 ; return
  end if

  ! debug output setup
  if(idebug) then
     call find_io_unit_list(2,indxn(1:2))
     i47=indxn(1)
     i48=indxn(2)
     if(ivec.gt.nsamp) then
        open(unit=i48,file='eq_RZout.dat',status='replace')
        write(i48,*) '#  Ntheta, Nsurf'
        write(i48,'(1x,2i6)') itheta,nsamp
        write(i48,*) '#  No., Orig. Index, rho value (est.)'
        do ii=1,nsamp
           indxn(ii) = ((ivec-1)*ii)/nsamp + 1
           rhoxn(ii) = ((indxn(ii)-1)*1.0d0)/(ivec-1)
           write(i48,'(1x,2i6,1x,f6.3)') ii,indxn(ii),rhoxn(ii)
        enddo
        write(i48,*) '#  iTheta, iSurf, Rgyro, Zgyro, Rorig, Zorig'
     endif
  endif

  allocate(zr(itheta), zz(itheta))
  allocate(text(2,itheta),rext(2,itheta))

  mlim   = min(nlim, size(xlimits))
  limits = 0._r8     ! make sure all elements are set
  kaxis  = 0

  !
  ! all exits through 999
  !

  jj=1  ! for debug sampling

  do iv = 1, ivec    ! loop over surfaces
     zr = ar(:,iv)
     zz = az(:,iv)

     if (max(abs(zr(itheta)-zr(1)),abs(zz(itheta)-zz(1)))>1.e-10_r8*max(1.e-6_r8,maxval(zr))) then
        mline = "?surfgeo: R,Z surfaces are not closed"
        ier=1 ; goto 999
     end if

     zr(itheta) = zr(1)   ! i=itheta is the same point as i=1
     zz(itheta) = zz(1)
     
     !
     ! --- centroid ---
     !
     call r8_plcentr(zr, zz, itheta, zrcentr, zzcentr)

     zmidp(iv) = zzcentr

     !
     ! snatched from r8_getmpa
     ! find approximate midplane intercept locations
     !
     rmnmidp(iv)=zrcentr  ! in case at magnetic axis
     rmjmidp(iv)=zrcentr
     tr=zr(1)
     tz=zz(1)
     zx=1.D-8*abs(zrcentr) ! really small

     do i=2,itheta
        trp=tr
        tzp=tz
        tr=zr(i)
        tz=zz(i)
        
        zytest=(tz-zzcentr)*(zzcentr-tzp)
        if(zytest.ge.0.0D0) then
           if (abs(tz-tzp)>zx) then
              zrans=(tr*(zzcentr-tzp)+trp*(tz-zzcentr))/(tz-tzp)  ! linear interpolation
              if(zrans.gt.zrcentr) then
                 rmjmidp(iv)=zrans
              else
                 rmnmidp(iv)=zrans
              endif
           end if
        endif
     enddo

     !
     ! --- elongation,triangulation,indentation,limits ---
     ! snatched from gelong
     !
     zrmin=zr(1)
     zrmax=zr(1)
     zymin=zz(1)
     zymax=zz(1)
     zrymin=zr(1)
     zrymax=zr(1)
     zyrmax=zz(1)

     !
     ! the 't' indices are the floating point equivalent of 'i' going from 0 to itheta-1
     ! on the surface.  The 't' domain is [0.,itheta-1) and theta=2.*pi*t/(itheta-1).
     !
     tymin=0.D0  
     tymax=0.D0
     trmin=0.D0
     trmax=0.D0
     iext=0
     ipt=1
     do i=-1,itheta-1                     ! i=0,i=itheta-1  corresponds to theta=0.
        k = mod(i+itheta-1,itheta-1)      ! i restricted to the range [0,itheta-2]
        tr = zr(k+1)
        tz = zz(k+1)
        
        if (tr<zrmin) then
           zrmin=tr
           trmin=k
        end if
        if (tr>zrmax) then
           zrmax=tr
           zyrmax=tz
           trmax=k
        end if
        if(tz.lt.zymin) then
           zymin=tz
           zrymin=tr
           tymin=k
        endif
        if(tz.gt.zymax) then
           zymax=tz
           zrymax=tr
           tymax=k
        endif
        
        if(i.le.0) then
           ipt=ipt+1
           zra(ipt)=tr                 ! ipt=2, then 3, ...
           zya(ipt)=tz
        else
           zra(1)=zra(2)               ! i=1 is first real point
           zra(2)=zra(3)
           zra(3)=tr
           zya(1)=zya(2)
           zya(2)=zya(3)
           zya(3)=tz
           !        
           !  parabolic fitting of 3 pt sequence to improve estimate of min,max.
           !
           zzR0=zra(2)
           zzAR=(zra(3)-zra(1))*0.5_R8
           zzBR=(zra(3)+zra(1))*0.5_R8- zzR0
           
           zzY0=zya(2)
           zzAY=(zya(3)-zya(1))*0.5_R8
           zzBY=(zya(3)+zya(1))*0.5_R8- zzY0
           
           if((zra(2).ge.max(zra(1),zra(3))).and. &
                (zra(2).gt.1.000001_R8*min(zra(1),zra(3)))) then
              zx=-zzAR/(2.0_R8*zzBR)
              zrx=zzR0+zx*(zzAr+ zx*zzBR)
              zyrx=zzY0+zx*(zzAy+ zx*zzBY)
              zt=k+zx-1
              if (zrx>ZRMAX) then
                 ZRMAX=zrx
                 ZYRMAX=zyrx
                 trmax=zt
              end if
              
              iext(2) = iext(2)+1    ! store this R maximum
              text(2,iext(2))=zt
              rext(2,iext(2))=zrx
           endif
           
           if((zra(2).le.min(zra(1),zra(3))).and. &
                (zra(2).lt.0.999999_R8*max(zra(1),zra(3)))) then
              zx=-zzAR/(2.0_R8*zzBR)
              zrx=zzR0+zx*(zzAr+ zx*zzBR)
              zt=k+zx-1
              if (zrx<ZRMIN) then
                 ZRMIN=zrx
                 trmin=zt
              end if
              
              iext(1) = iext(1)+1    ! store this R minimum
              text(1,iext(1))=zt
              rext(1,iext(1))=zrx
           endif
           
           if((zya(2).ge.max(zya(1),zya(3))).and. &
                (zya(2).gt.1.000001_R8*min(zya(1),zya(3)))) then
              zx=-zzAY/(2.0_R8*zzBY)
              zyx=zzY0+zx*(zzAy+ zx*zzBY)
              zryx=zzR0+zx*(zzAr+ zx*zzBR)
              zt=k+zx-1
              if(zyx.gt.ZYMAX) then
                 ZYMAX=zyx
                 ZRYMAX=zryx
                 tymax=zt
              endif
           endif
           
           if((zya(2).le.min(zya(1),zya(3))).and. &
                (zya(2).lt.0.999999_R8*max(zya(1),zya(3)))) then
              zx=-zzAY/(2.0_R8*zzBY)
              zyx=zzY0+zx*(zzAy+ zx*zzBY)
              zryx=zzR0+zx*(zzAr+ zx*zzBR)
              zt=k+zx-1
              if(zyx.lt.ZYMIN) then
                 ZYMIN=zyx
                 ZRYMIN=zryx
                 tymin=zt
              endif
           endif
           
        endif        
     enddo
     !
     !  elongation
     !
     elong(iv) = (zymax-zymin)/max(zrmax-zrmin, 1.D-30)
     
     !
     ! triangularity
     !
     zrmajor=0.5_R8*(zrmin+zrmax)
     zrminor=max(0.5_R8*(zrmax-zrmin), 1.D-30)
     
     if (zrminor<RMINAXIS*zrmajor) kaxis=iv

     triang(iv) = (zrmajor-min(zrymin,zrymax))/zrminor
     !ztriangU=(zrmajor-zrymax)/zrminor
     !ztriangL=(zrmajor-zrymin)/zrminor
     
     !
     ! indentation
     !
     zrindent_in  = zrmin   ! defaults will give zero indentation
     zrindent_out = zrmax
     trmin=tnorm(trmin)
     trmax=tnorm(trmax)
     tymin=tnorm(tymin)
     tymax=tnorm(tymax)
     if (tymin/=tymax) then
        tm = min(tymin,tymax)  ! define interval between top and bottom
        tp = max(tymin,tymax)
        
        o_bound   = (tm==trmin .or. tp==trmin .or. tm==trmax .or. tp==trmax)
        o_inside  = (tm<trmin .and. trmin<tp)  ! true if tm,tp bounds inboard side
        o_outside = (tm<trmax .and. trmax<tp)  ! true if tm,tp bounds outboard side
        
        if ((o_inside .neqv. o_outside) .and. (.not. o_bound)) then
           ! there is sufficient separation of points to evaluate indentation
           do i = 1, iext(1)                  ! loop over minimum extremums
              zt = tnorm(text(1,i))           ! index at extremum
              if (zt==tm .or. zt==tp) cycle   ! skip if near the top or bottom
              o_bound = (tm<zt .and. zt<tp)   ! .true. if within tm,tp limits
              if (o_bound .eqv. o_outside) then
                 zrindent_out = min(zrindent_out, rext(1,i))
              end if
           end do
           do i = 1, iext(2)                  ! loop over maximum extremums
              zt = tnorm(text(2,i))           ! index at extremum
              if (zt==tm .or. zt==tp) cycle   ! skip if near the top or bottom
              o_bound = (tm<zt .and. zt<tp)   ! .true. if within tm,tp limits
              if (o_bound .eqv. o_inside) then
                 zrindent_in = max(zrindent_in, rext(2,i))
              end if
           end do
        end if
     end if

     indent(iv) = max ( (zrindent_in-zrmin)/zrminor, (zrindent_out-zrmax)/zrminor )
     
     !
     ! -- upper squareness --
     !
     squpp    = 0.d0
     sqmk_upp = 0.d0

     a = zrmax - zrymax
     b = zymax - zyrmax

     kupp = b/max(zrminor, 1.D-30)
     dupp = min(1.d0, max(-1.d0,(zrmajor-zrymax)/max(zrminor, 1.D-30)))

     th0 = atan2(b,a)

     ith = 0
     do k=-1, itheta
        i = k + int(trmax+0.5d0) - 2
        i = mod(i+itheta-1,itheta-1) + 1   ! [1,itheta-1]

        if (ith>0) then
           ap  = ak
           bp  = bk
           thp = th
        end if

        ak = zr(i)-zrymax
        bk = zz(i)-zyrmax

        if (ak<=0.d0) then
           ith=0 ; cycle
        end if

        th  = atan2(bk,ak)
        ith = ith+1

        if (ith>1 .and. th0>=thp .and. th0<=th .and. th>thp) then
           z = (th0-thp)/(th-thp)

           asq = ap + z*(ak-ap)
           bsq = bp + z*(bk-bp)

           ab = sqrt(asq**2+bsq**2)
           ad = sqrt(a**2+b**2)
           ac = ad/sqrt(2.0d0)
           cd = ad-ac
           squpp = (ab-ac)/cd
           exit
        end if
     end do

     rmkupp = zrmajor + zrminor*cos(pi4 + asin(dupp)/sqrt(2.d0))

     imk = 1
     do i= 1, itheta-1
        if (zz(i)<zyrmax .or. zr(i)<zrymax) cycle        ! outside ellipse quadrant
        if (abs(zr(i)-rmkupp)<abs(zr(imk)-rmkupp)) imk=i
     end do

     i0 = mod(imk-1 +itheta-2,itheta-1) + 1
     i1 = mod(imk   +itheta-2,itheta-1) + 1
     i2 = mod(imk+1 +itheta-2,itheta-1) + 1

     !if (iv>1 .and. iv==ivec) then
     !   print *, 'at upp ivec'
     !end if

     if (rmkupp>=zr(i2) .and. rmkupp<=zr(i0)) then  ! only do this if in bounds
        call zsolve(rmkupp, zmkupp, ii)
        if (ii==0) then
           a = (zmkupp-zyrmax)/max(kupp*zrminor,1.d-30)
           if (abs(a)<=1.d0) then
              sqmk_upp = asin(a)-pi4
           end if
        end if
     end if
     
     !
     ! -- lower squareness --
     !
     sqlow    = 0.d0
     sqmk_low = 0.d0

     a = zrmax - zrymin
     b = zymin - zyrmax  ! <0

     klow = -b/max(zrminor, 1.D-30)
     dlow = min(1.d0, max(-1.d0,(zrmajor-zrymin)/max(zrminor, 1.D-30)))

     th0 = atan2(b,a)

     ith = 0
     do k=-1, itheta
        i = -k + int(trmax+0.5d0)
        i = mod(i+itheta-1,itheta-1) + 1   ! [1,itheta-1]

        if (ith>0) then
           ap  = ak
           bp  = bk
           thp = th
        end if

        ak = zr(i)-zrymin
        bk = zz(i)-zyrmax

        if (ak<=0.d0) then
           ith=0 ; cycle
        end if

        th  = atan2(bk,ak)
        ith = ith+1

        if (ith>1 .and. th0<=thp .and. th0>=th .and. th<thp) then
           z = (th0-thp)/(th-thp)

           asq = ap + z*(ak-ap)
           bsq = bp + z*(bk-bp)

           ab = sqrt(asq**2+bsq*bsq)
           ad = sqrt(a**2+b*b)
           ac = ad/sqrt(2.0d0)
           cd = ad-ac
           sqlow = (ab-ac)/cd
           exit
        end if
     end do
     
     rmklow = zrmajor + zrminor*cos(-pi4 - asin(dlow)/sqrt(2.d0))

     imk = 1
     do i= 1, itheta-1
        if (zz(i)>zyrmax .or. zr(i)<zrymin) cycle        ! outside ellipse quadrant
        if (abs(zr(i)-rmklow)<abs(zr(imk)-rmklow)) imk=i
     end do

     i0 = mod(imk-1 +itheta-2,itheta-1) + 1
     i1 = mod(imk   +itheta-2,itheta-1) + 1
     i2 = mod(imk+1 +itheta-2,itheta-1) + 1

     !if (iv>1 .and. iv==ivec) then
     !   print *, 'at low ivec'
     !end if

     if (rmklow>=zr(i0) .and. rmklow<=zr(i2)) then  ! only do this if in bounds
        call zsolve(rmklow, zmklow, ii)
        if (ii==0) then
           a = (zmklow-zyrmax)/max(klow*zrminor,1.d-30)
           if (abs(a)<=1.d0) then
              sqmk_low = -asin(a)-pi4
           end if
        end if
     end if

     !
     ! limits
     !
     xlimits = (/zrmin,    zrmax,       zrymax,       zymax,    zrymin,  &
                 zymin,    zrindent_in, zrindent_out, zrcentr,  zzcentr, &
                 sqlow,    squpp,       zyrmax,       rmklow,   zmklow,  &
                 sqmk_low, rmkupp,      zmkupp,       sqmk_upp, 0.d0    /)

     limits(1:mlim,iv) = xlimits(1:mlim)

     if (idebug .and. iv>1 .and. iv==ivec) then
        !
        !plot 'eq_bdy.dat' using 3:4 with lines title 'bdy points', 'eq_mom.dat' using 3:4 with lines title 'moments', \
        !'eq_pts.dat' index 0 with points title 'R0,Z0', '' index 1 with points title 'upper pi/4', '' index 2 with points title 'lower pi/4'
        !
        write(6,*) ' **surgeo.f90** debug output files:'
        write(6,*) '   eq_bdy.dat, eq_mom.dat, eq_pts.dat, eq_RZout.dat'
        call find_io_unit(i47)
        open(unit=i47, file="eq_bdy.dat", status="replace")

        do i=1, itheta
           if (zz(i)>=zyrmax) then
              ak = zr(i)-zrymax
              bk = zz(i)-zyrmax
           else
              ak = zr(i)-zrymin
              bk = zz(i)-zyrmax
           end if

           th = atan2(bk,ak)

           write(i47,'(i5,3(1x,es14.6))') i, th, zr(i), zz(i)
        end do
        close(i47)

        open(unit=i47, file="eq_mom.dat", status="replace")

        do i=1, itheta
           th0 = twopi*(i-1)/real(itheta-1,r8)

           if (th0<=twopi/2.d0) then
              ak = zrmajor + zrminor*cos(th0 + asin(dupp)*sin(th0))
              bk = zyrmax  + kupp*zrminor*sin(th0 + sqmk_upp*sin(2.d0*th0))
           else
              ak = zrmajor + zrminor*cos(th0 + asin(dlow)*sin(th0))
              bk = zyrmax  + klow*zrminor*sin(th0 + sqmk_low*sin(2.d0*th0))
           end if

           th = atan2(bk,ak)

           write(i47,'(i5,4(1x,es14.6))') i, th, ak, bk, th0
        end do
        close(i47)

        open(unit=i47, file="eq_pts.dat", status="replace")
        write(i47,'(2(1x,es14.6)//)') zrmajor, zyrmax
        write(i47,'(2(1x,es14.6)//)') rmkupp, zmkupp
        write(i47,'(2(1x,es14.6)//)') rmklow, zmklow
        close(i47)
     end if

     if(idebug.AND.(indxn(jj).eq.iv)) then

        do i=1, itheta
           th0 = twopi*(i-1)/real(itheta-1,r8)

           zsq = 0.5d0*(sqmk_upp+sqmk_low)

           ak = zrmajor + zrminor*cos(th0 + asin(triang(iv))*sin(th0))
           bk = zyrmax  + elong(iv)*zrminor*sin(th0 + zsq*sin(2.0d0*th0))

           write(i48,'(1x,2i6,4(1x,1pe14.7))') i,jj,ak,bk,zr(i),zz(i)
        enddo

        jj=min(nsamp,jj+1)
     endif
  end do

  if (kaxis>0 .and. kaxis<ivec) then   ! replace metrics which have no meaning at the axis
     elong(1:kaxis)  = elong(kaxis+1)
     triang(1:kaxis) = triang(kaxis+1)
     indent(1:kaxis) = indent(kaxis+1)

     if (nlim>=11) limits(11,1:kaxis) = limits(11,kaxis+1)
     if (nlim>=12) limits(12,1:kaxis) = limits(12,kaxis+1)
     if (nlim>=16) limits(16,1:kaxis) = limits(16,kaxis+1)
     if (nlim>=19) limits(19,1:kaxis) = limits(19,kaxis+1)
  end if

  !
  ! exit after allocation
  !
999 continue
  deallocate(zr,zz)
  deallocate(text,rext)

  if(idebug) then
     close(unit=i48)
  endif

contains
  !
  ! restrict the argument to the range [0,itheta-1)
  !
  function tnorm(x) result(t)
    real*8, intent(in) :: x    ! normalized theta index
    real*8 :: t                ! restricted domain
    real*8 :: p                ! itheta-1
    p = itheta-1
    t = x
    do while (t<0.D0)
       t = t+p
    end do

    do while (t>=p)
       t = t-p
    end do
  end function tnorm

  !
  ! given indices i0,i1,i2 and target rtar, return the corrsponding ztar
  ! which is on the quadratic curve through i0,i1,i2 of zr(),zz().  Based
  ! on interpolation over angle from (Rmajor, Zrmax).
  !
  subroutine zsolve(rtar, ztar, ii)
    real*8,  intent(in)  :: rtar     ! target R
    real*8,  intent(out) :: ztar     ! returned Z
    integer, intent(out) :: ii       ! nonzero on error

    real*8 :: r0,r1,r2   ! R at i0,i1,i2 relative to Rmajor
    real*8 :: z0,z1,z2   ! Z at i0,i1,i2 relative to Z at Rmax
    real*8 :: d0,d1,d2   ! theta at i0,i1,i2
    real*8 :: rt         ! Rtar relative to Rmajor

    real*8 :: del1,del2         ! (theta1-theta0), (theta2-theta0)
    real*8 :: d                 ! temp
    real*8 :: m11,m12,m21, m22  ! inverse matrix

    real*8 :: ar,br,az,bz       ! quadratic coefficients for R,Z equation
    real*8 :: da,db             ! del solutions
    
    ii=0

    r0 = zr(i0)-zrmajor ; z0 = zz(i0)-zyrmax ; d0 = btan2(z0,r0)
    r1 = zr(i1)-zrmajor ; z1 = zz(i1)-zyrmax ; d1 = btan2(z1,r1)
    r2 = zr(i2)-zrmajor ; z2 = zz(i2)-zyrmax ; d2 = btan2(z2,r2)
    rt = rtar-zrmajor

    del1 = d1-d0
    del2 = d2-d0

    d = del1*del2*(del1-del2)

    if (d==0.d0) goto 100

    m11 =  del2/d
    m12 = -del1/d
    m21 = -del2*del2/d
    m22 =  del1*del1/d

    ar = m11*(r1-r0) + m12*(r2-r0)
    br = m21*(r1-r0) + m22*(r2-r0)

    az = m11*(z1-z0) + m12*(z2-z0)
    bz = m21*(z1-z0) + m22*(z2-z0)

    d = br*br-4.d0*ar*(r0-rt)

    if (d<0.d0 .or. ar==0.d0) goto 100

    d = sqrt(d)

    da = (-br + d)/(2.d0*ar)
    db = (-br - d)/(2.d0*ar)

    if (abs(db-del1)<abs(da-del1)) da=db

    ztar = da*(da*az+bz) + z0 + zyrmax
    
    return

100 continue
    ii=1
  end subroutine zsolve

  function btan2(z,r) result(b)
    real*8, intent(in) :: z,r
    real*8 :: b
    
    if (z==0.d0 .and. r==0.d0) then
       b = 0.d0
    else
       b = atan2(z,r)
    end if
  end function btan2
end subroutine surfgeo
