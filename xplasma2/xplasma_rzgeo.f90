module xplasma_rzgeo

  !  methods for defining flux surface geometry -- 2d axisymmetric tokamak

  use xplasma_obj
  use xplasma_ctran
  use xplasma_flxint, only: xplasma_rho_zonint, xplasma_2d_zonint
  use xplasma_sol, only: xplasma_limcheck,xplasma_lim_RZminmax
  implicit NONE

  private

  public :: xplasma_rzmagf,xplasma_rzmag
  public :: xplasma_jacheck,xplasma_eqClear,xplasma_eqcheck
  public :: xplasma_gen_q,xplasma_gen_p,xplasma_gen_torflx
  public :: xplasma_gen_boozer
  public :: xplasma_volume,xplasma_vol1,xplasma_voln
  public :: xplasma_area,xplasma_area1,xplasma_arean
  public :: xplasma_create_RZgrid
  public :: xplasma_max_maptype
  public :: xplasma_surfgeo

  interface xplasma_volume
     module procedure xplasma_vol1
     module procedure xplasma_voln
  end interface

  interface xplasma_area
     module procedure xplasma_area1
     module procedure xplasma_arean
  end interface

  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: THREE= 3.0d0
  real*8, parameter :: C2PI = 6.2831853071795862D+00

  contains

    !=======================================================================
    subroutine xplasma_rzmagf(s,rzmapper,id_rho,id_theta,id_R,id_Z,ierr, &
         zdiff)

      !  ** public **

      !  create 2d functions R(rho,chi), Z(rho,chi)
      !  using passed external routine "rzmapper"
      !  rzmapper is a dummy argument for a subroutine which MUST have the
      !  following interface:

      !      subroutine rzmapper(rho,chi,R,Z,ierr)

      !      real*8 rho,chi      ! INPUT (rho,chi) - radial, poloid coordinates
      !      real*8 R,Z          ! OUTPUT R(rho,chi), Z(rho,chi)
      !      integer IERR        ! OUTPUT completion code, 0=OK
      !
      !-------------------------------
      !  arguments

      type (xplasma), pointer :: s      ! object, state modified by this call

      ! input:
      external rzmapper                 ! passed subroutine, as described

      integer, intent(in) :: id_rho     ! rho grid id
      integer, intent(in) :: id_theta   ! theta grid id

      ! output:
      integer, intent(out) :: id_R      ! id for R spline
      integer, intent(out) :: id_Z      ! id for Z spline
      integer, intent(out) :: IERR      ! completion code, 0=OK
  
      ! OPTIONAL input (default=1.0d-3)

      real*8, intent(in), optional :: zdiff  ! finite difference factor for 
      ! boundary condition calculation-- use fraction zdiff of the 1st/last
      ! rho grid zone in evaluating dR/drho and dZ/drho finite difference
      ! formulae for first derivative boundary conditions

      !--------------------------
      !  local

      integer :: irho,ith,nrho,nth,ibc

      real*8, dimension(:), allocatable :: zrho,zth
      real*8, dimension(:,:), allocatable :: R,Z
      real*8, dimension(:), allocatable :: Rbc0,Rbc1,Zbc0,Zbc1

      real*8 R0sum,Z0sum,zdiffi,zrho0x,zrho1x,zdrho0,zdrho1

      real*8, parameter :: ONE = 1.0d0
      real*8, parameter :: EPS4 = 1.0d-4
      !--------------------------

      id_R=0
      id_Z=0
      ierr=0

      if(present(zdiff)) then
         zdiffi=min(ONE,max(EPS4,zdiff))
      else
         zdiffi=1.0d-3
      endif

      R0sum=ZERO
      Z0sum=ZERO

      call xplasma_grid_size(s,id_rho,nrho,ierr)
      if(ierr.ne.0) return

      call xplasma_grid_size(s,id_theta,nth,ierr)
      if(ierr.ne.0) return

      allocate(zrho(nrho),zth(nth),R(nrho,nth),Z(nrho,nth))
      allocate(Rbc0(nth),Rbc1(nth),Zbc0(nth),Zbc1(nth))

      do
         call xplasma_grid(s,id_rho,zrho,ierr)
         if(ierr.ne.0) exit

         zdrho0=zdiffi*(zrho(2)-zrho(1))
         zrho0x=zrho(1)+zdrho0

         zdrho1=zdiffi*(zrho(nrho)-zrho(nrho-1))
         zrho1x=zrho(nrho)-zdrho1

         call xplasma_grid(s,id_theta,zth,ierr)
         if(ierr.ne.0) exit

         do ith=1,nth

            call rzmapper(zrho0x,zth(ith),Rbc0(ith),Zbc0(ith),ierr)
            if(ierr.ne.0) exit

            call rzmapper(zrho1x,zth(ith),Rbc1(ith),Zbc1(ith),ierr)
            if(ierr.ne.0) exit

            do irho=1,nrho

               call rzmapper(zrho(irho),zth(ith),R(irho,ith),Z(irho,ith),ierr)
               if(ierr.ne.0) exit

               if(irho.eq.1) then
                  R0sum=R0sum+R(irho,ith)
                  Z0sum=Z0sum+Z(irho,ith)
               endif
            enddo
            if(ierr.ne.0) exit
         enddo
         if(ierr.ne.0) exit

         !  eval finite difference BC derivative estimates

         do ith=1,nth
            Rbc0(ith)=(Rbc0(ith)-R(1,ith))/zdrho0
            Zbc0(ith)=(Zbc0(ith)-Z(1,ith))/zdrho0

            Rbc1(ith)=(R(nrho,ith)-Rbc1(ith))/zdrho1
            Zbc1(ith)=(Z(nrho,ith)-Zbc1(ith))/zdrho1
         enddo

         !  force axis to be exactly fixed vs. chi

         R0sum=R0sum/nth
         Z0sum=Z0sum/nth

         R(1,1:nth)=R0sum
         Z(1,1:nth)=Z0sum

         ibc=1
         call xplasma_rzmag(s,id_rho,id_theta,R,Z,id_R,id_Z,ierr, &
              ibc,Rbc0,Zbc0, ibc,Rbc1,Zbc1)

         exit
      enddo

      if(ierr.ne.0) then
         call xplasma_errmsg_append(s,' error in xplasma_rzmagf')
         deallocate(zrho,zth,R,Z,Rbc0,Rbc1,Zbc0,Zbc1)
      endif

    end subroutine xplasma_rzmagf

    !=======================================================================
    subroutine xplasma_rzmag(s,idx1,idx2,R,Z,id_R,id_Z,ierr, &
              ibc0,Rbc0,Zbc0, ibc1,Rbc1,Zbc1, lhermite, &
              Rderivs, Zderivs)

      !  ** public **

      !  create 2d functions R(x,theta) and Z(x,theta) from data arrays
      !  provided, with optional boundary condition specifications

      !  NOTE: xplasma stores {R,Z} data internally against a poloidal angle
      !  coordinate which increases in the counter-clockwise direction.  IF
      !  {R,Z} data is given with a clockwise orientation, the theta grid
      !  (identified with idx1 or idx2) has to have been defined to xplasma
      !  as having a clockwise orientation as well, for results to be correct.

      !  This code checks the data orientation, but it cannot re-check the 
      !  original theta grid.  NB if theta is evenly spaced, there will
      !  be no error.

      !-------------------------------
      !  arguments

      type (xplasma), pointer :: s      ! object, state modified by this call

      integer, intent(in) :: idx1       ! grid id of 1st dimension of data
      integer, intent(in) :: idx2       ! grid id of 2nd dimension of data

      real*8, intent(in), dimension(:,:) :: R,Z ! R(x,theta) & Z(x,theta)
      !  OR R(theta,x), z(theta,x) -- sized to match specified grids

      ! output:
      integer, intent(out) :: id_R      ! id for R spline
      integer, intent(out) :: id_Z      ! id for Z spline
      integer, intent(out) :: IERR      ! completion code, 0=OK
  
      ! OPTIONAL input

      integer, intent(in), optional :: ibc0  ! boundary condition control
      !  for R & Z at x=rho=0
      !  ibc0=0: default "not a knot" (Rbc0 & Zbc0 not used)
      !  ibc0=1: 1st derivative BC: Rbc0(j) = dR/dx @ x=0,theta=theta(j)
      !          j refers to jth point in theta grid (id_theta)
      !  ibc0=2: 2nd derivative BC: Zbc0(k) = d2Z/dx2 @ x=0, theta=theta(k)
      !          k refers to kth point in theta grid (id_theta)

      real*8, intent(in), dimension(:), optional :: Rbc0,Zbc0  ! BC data
      
      integer, intent(in), optional :: ibc1  ! boundary condition control
      !  for R & Z at x=rho=1
      !  ibc1=0: default "not a knot" (Rbc1 & Zbc1 not used)
      !  ibc1=1: 1st derivative BC: Rbc1(j) = dR/dx @ x=1,theta=theta(j)
      !          j refers to jth point in theta grid (id_theta)
      !  ibc1=2: 2nd derivative BC: Zbc1(k) = d2Z/dx2 @ x=1, theta=theta(k)
      !          k refers to kth point in theta grid (id_theta)

      real*8, intent(in), dimension(:), optional :: Rbc1,Zbc1  ! BC data

      ! OPTION: convert representation to Hermite when done...

      logical, intent(in), optional :: lhermite

      ! OPTIONAL input -- Hermite derivatives at nodal points
      !   Rderivs(1,:,:) d/dx1 (ref idx1)
      !   Rderivs(2,:,:) d/dx2 (ref idx2)
      !   Rderivs(3,:,:) d2/dx1dx2
      !     ...similarly for Zderivs

      ! if one of these is present BOTH must be present
      ! if these are present, BC data must NOT be present

      real*8, intent(in), dimension(:,:,:), optional :: Rderivs,Zderivs

      !---------------------------------------------------
      !  local: 
      integer :: jbc0,jbc1
      real*8, dimension(:), allocatable :: zRbc0,zRbc1,zZbc0,zZbc1
      real*8, dimension(:), allocatable :: rbdy,zbdy,atan00
      real*8, dimension(:,:), allocatable :: ruse,zuse,dist0,atan0
      real*8, dimension(:,:,:), allocatable :: rdduse,zdduse
      integer :: i,inrho,inchi,iertmp,ichi,irho,icount,isign,idiff
      integer :: id_rho,id_chi,icoord,isw_rho,idtmp,ishift
      integer :: ichk_derivs,ichk_bc
      real*8 :: ccwsum,dra,drb,dza,dzb,R0,Z0,raxis,zaxis,zdra,zdza,zatan
      real*8 :: thmin
      logical :: ccwflag,khermite
      real*8, parameter :: CZERO = 0.0d0
      real*8, parameter :: C2PI  = 6.2831853071795862D+00
      character*120 atan00_label,zmsg
      character*32 zname
      !---------------------------------------------------

      ierr=0
      id_R=0
      id_Z=0

      ichk_bc=0
      ichk_derivs=0

      if(present(ibc0)) ichk_bc=ichk_bc + 1
      if(present(ibc1)) ichk_bc=ichk_bc + 1

      if(present(rbc0)) ichk_bc=ichk_bc + 1
      if(present(zbc0)) ichk_bc=ichk_bc + 1

      if(present(rbc1)) ichk_bc=ichk_bc + 1
      if(present(zbc1)) ichk_bc=ichk_bc + 1

      if(present(Rderivs)) ichk_derivs=ichk_derivs + 1
      if(present(Zderivs)) ichk_derivs=ichk_derivs + 1

      ierr = 0
      if((ichk_bc.gt.0).and.(ichk_derivs.gt.0)) then
         ierr = 2001
         call xplasma_errmsg_append(s, &
              ' ?(R,Z) boundary condition data and nodal derivative data cannot BOTH be present!')
      endif

      if(ichk_derivs.eq.1) then
         ierr = 2001
         call xplasma_errmsg_append(s, &
              ' ?If nodal derivative data is provided for either R or Z it must be provide for BOTH!')
      endif

      khermite=.FALSE.
      if(present(lhermite)) khermite=lhermite

      if(ierr.ne.0) return

      do
         id_rho=0
         id_chi=0

         !  check grids; find which one is rho

         call xplasma_grid_info(s,idx1,ierr, coord=icoord); if(ierr.ne.0) exit

         if(icoord.eq.xplasma_rho_coord) then
            id_rho=idx1
            isw_rho=1
         else if(icoord.eq.xplasma_theta_coord) then
            id_chi=idx1
            isw_rho=2
         endif

         call xplasma_grid_info(s,idx2,ierr, coord=icoord); if(ierr.ne.0) exit

         if(icoord.eq.xplasma_rho_coord) then
            id_rho=idx2
         else if(icoord.eq.xplasma_theta_coord) then
            id_chi=idx2
         endif

         if((id_rho.eq.0).or.(id_chi.eq.0)) then
            ierr=2001
            call xplasma_errmsg_append(s, &
                 ' ?both theta and rho grid identifiers are required.')
            exit
         endif

         call xplasma_grid_info(s,id_chi,ierr, xmin=thmin)

         jbc0=0
         if(present(ibc0)) jbc0=ibc0

         jbc1=0
         if(present(ibc1)) jbc1=ibc1

         call xplasma_grid_size(s,id_rho,inrho,ierr)
         if(ierr.ne.0) exit

         call xplasma_grid_size(s,id_chi,inchi,ierr)
         if(ierr.ne.0) exit

         !-------------------------------
         !  for backwards compatibility:  if grid names are not __RHO,__CHI
         !  create copies with these names...

         call xplasma_grid_info(s,id_rho,ierr, name=zname)
         if(ierr.ne.0) exit
         if(zname.ne.'__RHO') call gcopy(id_rho,inrho,'__RHO')
         if(ierr.ne.0) exit

         call xplasma_grid_info(s,id_chi,ierr, name=zname)
         if(ierr.ne.0) exit
         if(zname.ne.'__CHI') call gcopy(id_chi,inchi,'__CHI')
         if(ierr.ne.0) exit

         !-------------------------------
         if(isw_rho.eq.1) then
            if((size(R,1).ne.inrho).or.(size(R,2).ne.inchi)) then
               ierr=ierr+1
               call xplasma_errmsg_append(s, &
                    ' ?R data dimensions inconsistent with grids given.')
            endif
            if((size(Z,1).ne.inrho).or.(size(Z,2).ne.inchi)) then
               ierr=ierr+1
               call xplasma_errmsg_append(s, &
                    ' ?Z data dimensions inconsistent with grids given.')
            endif
            if(ichk_derivs.gt.0) then
               if((size(Rderivs,2).ne.inrho).or.(size(Rderivs,3).ne.inchi)) then
                  ierr=ierr+1
                  call xplasma_errmsg_append(s, &
                       ' ?Rderivs data dimensions inconsistent with grids given.')
               endif
               if((size(Zderivs,2).ne.inrho).or.(size(Zderivs,3).ne.inchi)) then
                  ierr=ierr+1
                  call xplasma_errmsg_append(s, &
                       ' ?Zderivs data dimensions inconsistent with grids given.')
               endif
               if((size(Rderivs,1).ne.3).or.(size(Zderivs,1).ne.3)) then
                  ierr=ierr+1
                  call xplasma_errmsg_append(s, &
                       ' ?1st dimension size of Rderivs & Zderivs should = 3.')
               endif
            endif
         else
            if((size(R,2).ne.inrho).or.(size(R,1).ne.inchi)) then
               ierr=ierr+1
               call xplasma_errmsg_append(s, &
                    ' ?R data dimensions inconsistent with grids given.')
            endif
            if((size(Z,2).ne.inrho).or.(size(Z,1).ne.inchi)) then
               ierr=ierr+1
               call xplasma_errmsg_append(s, &
                    ' ?Z data dimensions inconsistent with grids given.')
            endif
            if(ichk_derivs.gt.0) then
               if((size(Rderivs,3).ne.inrho).or.(size(Rderivs,2).ne.inchi)) then
                  ierr=ierr+1
                  call xplasma_errmsg_append(s, &
                       ' ?Rderivs data dimensions inconsistent with grids given.')
               endif
               if((size(Zderivs,3).ne.inrho).or.(size(Zderivs,2).ne.inchi)) then
                  ierr=ierr+1
                  call xplasma_errmsg_append(s, &
                       ' ?Zderivs data dimensions inconsistent with grids given.')
               endif
               if((size(Rderivs,1).ne.3).or.(size(Zderivs,1).ne.3)) then
                  ierr=ierr+1
                  call xplasma_errmsg_append(s, &
                       ' ?1st dimension size of Rderivs & Zderivs should = 3.')
               endif
            endif
         endif
         if(ierr.ne.0) then
            ierr=2001
            exit
         endif

         if(ichk_derivs.eq.0) then
            allocate(zRbc0(inchi),zZbc0(inchi),zRbc1(inchi),zZbc1(inchi))
            zRbc0=0; zRbc1=0; zZbc0=0; zZbc1=0

            iertmp=0

            if(present(Rbc0)) call bchk('R @ x=0',Rbc0,zRbc0,ierr)
            iertmp=max(iertmp,ierr)

            if(present(Rbc1)) call bchk('R @ x=1',Rbc1,zRbc1,ierr)
            iertmp=max(iertmp,ierr)

            if(present(Zbc0)) call bchk('Z @ x=0',Zbc0,zZbc0,ierr)
            iertmp=max(iertmp,ierr)

            if(present(Zbc1)) call bchk('Z @ x=1',Zbc1,zZbc1,ierr)
            iertmp=max(iertmp,ierr)

            ierr=iertmp
            if(ierr.ne.0) exit
         endif

         !---------------------------------------
         !  OK, determine theta orientation by examining outer bdy.  Does 
         !   increasing theta trace out the boundary counter-clockwise or 
         !   clockwise?

         allocate(rbdy(inchi),zbdy(inchi))
         if(isw_rho.eq.1) then
            rbdy=r(inrho,1:inchi)
            zbdy=z(inrho,1:inchi)
         else
            rbdy=r(1:inchi,inrho)
            zbdy=z(1:inchi,inrho)
         endif

         R0 = sum(rbdy)/inchi
         Z0 = sum(zbdy)/inchi

         ccwsum = 0
         do i = 1,inchi-1
            dra=Rbdy(i)-R0
            dza=Zbdy(i)-Z0
            drb=Rbdy(i+1)-Rbdy(i)
            dzb=Zbdy(i+1)-Zbdy(i)
            ccwsum = ccwsum - dza*drb + dra*dzb
         enddo

         deallocate(rbdy,zbdy)

         if(ccwsum.gt.CZERO) then

            !  normal orientation:  CCW
            ccwflag = .TRUE.
            isign=1

         else if(ccwsum.lt.CZERO) then

            !  reversed orientation:  CW
            ccwflag = .FALSE.
            isign=-1

         else
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_rzgeo: CW/CCW orientation determination failed.')
            call xplasma_errmsg_append(s, &
                 '  is the data all ZERO?')
            ierr=2001
         endif
         if(ierr.ne.0) exit

         ! OK -- proceed...
         ! delete {R,Z} and derived items...

         call xplasma_eqClear(s,.TRUE.,iertmp)

         !  create the {R,Z} splines... set BC on right dimension...

         iertmp = 0

         allocate(ruse(inchi,inrho),zuse(inchi,inrho))
         allocate(dist0(inchi,inrho),atan0(inchi,inrho),atan00(inrho))
         if(ichk_derivs.gt.0) then
            allocate(rdduse(3,inchi,inrho),zdduse(3,inchi,inrho))
         endif

         if(isw_rho.eq.1) then

            do ichi=1,inchi
               do irho=1,inrho
                  ruse(ichi,irho)=r(irho,ichi)
                  zuse(ichi,irho)=z(irho,ichi)
                  if(ichk_derivs.gt.0) then
                     rdduse(2,ichi,irho)=Rderivs(1,irho,ichi) ! note order
                     rdduse(1,ichi,irho)=Rderivs(2,irho,ichi) ! is swapped
                     rdduse(3,ichi,irho)=Rderivs(3,irho,ichi)

                     zdduse(2,ichi,irho)=Zderivs(1,irho,ichi) ! note order
                     zdduse(1,ichi,irho)=Zderivs(2,irho,ichi) ! is swapped
                     zdduse(3,ichi,irho)=Zderivs(3,irho,ichi)
                  endif
               enddo
            enddo

         else

            ruse = R
            zuse = Z
            if(ichk_derivs.gt.0) then
               rdduse = Rderivs
               zdduse = Zderivs
            endif

         endif

         !-----------------------------------
         ! dmc-- added inequality test; only reset axis if any points
         ! are unequal.  This improves reproducibility of representations
         ! as the averaging formula can move the axis ever so slightly even
         ! if all points are equal.

         idiff=0
         do i=2,inchi
            if(ruse(i-1,1).ne.ruse(i,1)) then
               idiff=1
               exit
            endif
         enddo
         if(idiff.gt.0) then
            ruse(1:inchi,1)=sum(ruse(1:inchi,1))/inchi
         endif

         idiff=0
         do i=2,inchi
            if(zuse(i-1,1).ne.zuse(i,1)) then
               idiff=1
               exit
            endif
         enddo
         if(idiff.gt.0) then
            zuse(1:inchi,1)=sum(zuse(1:inchi,1))/inchi
         endif
         !-----------------------------------

         raxis=ruse(1,1)
         zaxis=zuse(1,1)

         ruse(1,1:inrho)=(ruse(1,1:inrho)+ruse(inchi,1:inrho))/2
         ruse(inchi,1:inrho)=ruse(1,1:inrho)
         zuse(1,1:inrho)=(zuse(1,1:inrho)+zuse(inchi,1:inrho))/2
         zuse(inchi,1:inrho)=zuse(1,1:inrho)

         if(ichk_derivs.eq.0) then
            zrbc0(1)=(zrbc0(1)+zrbc0(inchi))/2
            zrbc0(inchi)=zrbc0(1)
            zrbc1(1)=(zrbc1(1)+zrbc1(inchi))/2
            zrbc1(inchi)=zrbc1(1)

            zzbc0(1)=(zzbc0(1)+zzbc0(inchi))/2
            zzbc0(inchi)=zzbc0(1)
            zzbc1(1)=(zzbc1(1)+zzbc1(inchi))/2
            zzbc1(inchi)=zzbc1(1)

         else
            rdduse(1,1:inchi,1)=0  ! d(R,Z)/dtheta = 0 on axis
            zdduse(1,1:inchi,1)=0
         endif

         do irho=inrho,1,-1
            do ichi=1,inchi-1
               if(irho.eq.1) exit
               zdra=ruse(ichi,irho)-raxis
               zdza=zuse(ichi,irho)-zaxis
               dist0(ichi,irho)=sqrt(zdra*zdra+zdza*zdza)
               atan0(ichi,irho)=atan2(zdza,zdra)
               if(irho.eq.2) then
                  dist0(ichi,1)=0.0d0
                  atan0(ichi,1)=atan0(ichi,2)
               endif
            enddo
            dist0(inchi,irho)=dist0(1,irho)
            !  now make atan0(theta,rho) be relative to the 1st theta atan.
            zatan=atan0(1,irho)
            atan00(irho)=zatan  ! arctan at 1st theta
            atan0(1,irho)=0.0d0
            atan0(inchi,irho)=isign*C2PI
            ishift=0
            do ichi=2,inchi-1
               if(ishift.eq.0) then
                  if(isign*(atan0(ichi,irho)-zatan).lt.CZERO) then
                     ishift=1
                  else
                     zatan=atan0(ichi,irho)
                  endif
               endif
               atan0(ichi,irho)=atan0(ichi,irho)-atan00(irho)+isign*ishift*C2PI
            enddo
         enddo

         !  sanity check

         do ichi=2,inchi
            if(isign*(atan0(ichi,inrho)-atan0(ichi-1,inrho)).lt. &
                 C2PI/(100*inchi)) then
               ierr=2003
               zmsg=' '
               write(zmsg,*) ' at boundary interval #',ichi-1
               call xplasma_errmsg_append(s,zmsg)
               zmsg=' '
               write(zmsg,*) ' from (R,Z) = ',ruse(ichi-1,inrho),zuse(ichi-1,inrho)
               call xplasma_errmsg_append(s,zmsg)
               zmsg=' '
               write(zmsg,*) '   to (R,Z) = ',ruse(ichi,inrho),zuse(ichi,inrho)
               call xplasma_errmsg_append(s,zmsg)
               zmsg=' '
               write(zmsg,*) ' distances to mag. axis: ', &
                    dist0(ichi-1,inrho),dist0(ichi,inrho)
               call xplasma_errmsg_append(s,zmsg)
               zmsg=' '
               write(zmsg,*) ' arc-tangent of displacement from axis: ', &
                    atan0(ichi-1,inrho)-atan00(inrho), &
                    atan0(ichi,inrho)-atan00(inrho)
               call xplasma_errmsg_append(s,zmsg)
               exit
            endif
         enddo
         iertmp=max(iertmp,ierr)

         if(ichk_derivs.eq.0) then
            !  these will be splines of order s%rzOrder...

            call xplasma_create_prof(s,'R',id_chi,id_rho, &
                 Ruse,id_R,ierr, &
                 ibcx2a=ibc0,zbcx2a=zRbc0,ibcx2b=ibc1,zbcx2b=zRbc1, &
                 ccwflag1=ccwflag)

            if((ierr.eq.0).and.khermite) call xplasma_hermitize(s,id_R,ierr)

            iertmp=max(iertmp,ierr)

            call xplasma_create_prof(s,'Z',id_chi,id_rho, &
                 Zuse,id_Z,ierr, &
                 ibcx2a=ibc0,zbcx2a=zZbc0,ibcx2b=ibc1,zbcx2b=zZbc1, &
                 ccwflag1=ccwflag)

            if((ierr.eq.0).and.khermite) call xplasma_hermitize(s,id_Z,ierr)

            iertmp=max(iertmp,ierr)
            ierr=iertmp

         else
            !  these will be Hermites with the specified nodal derivatives

            call xplasma_create_prof(s,'R',id_chi,id_rho, &
                 Ruse,id_R,ierr, &
                 ccwflag1=ccwflag, ispline=1, coeffs=Rdduse)
            iertmp=max(iertmp,ierr)
            call xplasma_create_prof(s,'Z',id_chi,id_rho, &
                 Zuse,id_Z,ierr, &
                 ccwflag1=ccwflag, ispline=1, coeffs=Zdduse)

            iertmp=max(iertmp,ierr)
            ierr=iertmp
         endif

         if(ierr.eq.0) then
            call xplasma_rzcheck(s,id_rho,inrho,ierr)
         endif

         if(ierr.eq.0) then
            call xplasma_author_set(s,xplasma_xmhd,iertmp)

            !  these will be piecewise linear functions...

            iertmp=0

            call xplasma_create_prof(s,'__dist_axis',id_chi,id_rho, &
                 dist0,idtmp,ierr, ccwflag1=ccwflag, &
                 label='distance to axis',units='m')
            iertmp=max(iertmp,ierr)

            call xplasma_create_prof(s,'__atan_axis',id_chi,id_rho, &
                 atan0,idtmp,ierr, ccwflag1=ccwflag, &
                 label='adjusted arctan((Z-Zaxis)/(R-Raxis))',units='rad')
            iertmp=max(iertmp,ierr)

            atan00_label=' '
            write(atan00_label, &
                 '(" arctan((Z-Zaxis)/(R-Raxis)) vs. rho @theta=",1pd11.4)') &
                 thmin

            call xplasma_create_prof(s,'__atan00_axis',id_rho, &
                 atan00,idtmp,ierr, &
                 label=trim(atan00_label),units='rad')

            ierr=iertmp

            call xplasma_author_clear(s,xplasma_xmhd,iertmp)

         endif

         deallocate(ruse,zuse,dist0,atan0,atan00)
         exit
      enddo

      if(allocated(zRbc0)) deallocate(zRbc0,zZbc0,zRbc1,zZbc1)
      if(allocated(rdduse)) deallocate(rdduse,zdduse)

      contains
        subroutine bchk(str,bcin,bcuse,ierr)
          character*(*), intent(in) :: str  ! for error msg
          real*8, dimension(:) :: bcin,bcuse
          integer, intent(out) :: ierr
          character*128 msgbuf

          ierr = 0

          if(size(bcin).ne.size(bcuse)) then
             ierr=2001
             call xplasma_errmsg_append(s, ' ?xplasm_rzgeo: '// &
                  str//' boundary condition array data size mismatch.')
             write(msgbuf,*) '   size expected: ',size(bcuse), &
                  '; size received: ',size(bcin)
             call xplasma_errmsg_append(s,msgbuf)
             call xplasma_errmsg_append(s,'  expected:'// &
                  ' boundary condition array dimension = theta dimension.')

          else

             bcuse = bcin

          endif
        end subroutine bchk

        subroutine gcopy(idg,ing,zname)
          ! copy a grid
          integer, intent(in) :: idg,ing  ! id & size
          character*(*), intent(in) :: zname  ! name

          !-----------------------------
          real*8 :: zgrid(ing),zgrid2(ing)
          integer :: id2,ing2
          logical :: imatch
          !-----------------------------

          call xplasma_grid(s,idg,zgrid,ierr)
          if(ierr.ne.0) return

          call xplasma_gridID(s,zname,id2)
          if(id2.eq.0) then
             imatch=.TRUE.
             call xplasma_create_grid(s,zname,idg,zgrid,idtmp,ierr, &
                  label=trim(zname)//': copy of grid with standard name.')
          else
             imatch=.TRUE.
             call xplasma_grid_info(s,id2,ierr,size=ing2)
             if(ierr.ne.0) return
             if(ing.ne.ing2) then
                imatch=.FALSE.
             else
                call xplasma_grid(s,id2,zgrid2,ierr)
                if(ierr.ne.0) return
                if(maxval(abs(zgrid-zgrid2)).gt.1.0d-12) imatch=.FALSE.
             endif
          endif

          if(.not.imatch) then
             ierr=9999
             call xplasma_errmsg_append(s, &
                  'grid supplied to xplasma_rzmag does not match existing grid: '//trim(zname))
          endif
        end subroutine gcopy

    end subroutine xplasma_rzmag

    !=======================================================================
    subroutine xplasma_eqclear(s,rzonly,ier)

      !  ** public **

      ! delete MHD equilibrium from xplasma (e.g. before inserting a new one).
      ! if (rzonly) delete {R,Z} & derived quantities
      ! if (.not.rzonly) delete {g,Psi,P,R,Z} & derived quantities

      type (xplasma), pointer :: s
      logical, intent(in) :: rzonly  ! R,Z only flag
      integer, intent(out) :: ier    ! status code =0 on exit
      !  ier only set if s is uninitialized...

      !------------------
      integer :: id_R,id_Z,id_g,id_psi,id_P
      !------------------

      call xplasma_common_ids(s,ier, &
           id_g=id_g,id_psi=id_psi,id_P=id_P,id_R=id_R,id_Z=id_Z)

      if(ier.ne.0) return

      call xplasma_author_set(s,xplasma_root,ier)
      if(ier.ne.0) return

      if(.not.rzonly) then
         if(id_g.ne.0) call xplasma_remove_item(s,id_g,ier)
         if(id_p.ne.0) call xplasma_remove_item(s,id_p,ier)
         if(id_psi.ne.0) call xplasma_remove_item(s,id_psi,ier)
      endif

      if(id_R.ne.0) call xplasma_remove_item(s,id_R,ier)
      if(id_Z.ne.0) call xplasma_remove_item(s,id_Z,ier)

      call xplasma_author_clear(s,xplasma_root,ier)

      call xplasma_remove_author(s,xplasma_xmhd,ier)

    end subroutine xplasma_eqclear

    !=======================================================================
    subroutine xplasma_rzcheck(s,idrho,inrho,ier)

      !  ** PRIVATE **

      use eqi_rzbox_module

      !  call this when {R,Z} are updated

      ! check update of R(rho,theta) and Z(rho,theta)
      !   check Jacobian
      !   get magnetic axis
      !   get RZ min and max of boundary surface
      !   set up integrator for profiles of enclosed volume and cross 
      !     sectional area
      !   check equilibrium against limiter, if one is present.
      !   if scrapeoff region flag is set, compute extrapolated surfaces.

      type (xplasma), pointer :: s
      integer, intent(in) :: idrho     ! flux surface grid id
      integer, intent(in) :: inrho     ! no. of radial flux surfaces
      integer, intent(out) :: ier      ! completion status 0=OK

      !-------------------------------
      integer :: id_R,id_Z
      integer :: idx1,idc1,idx2,idc2,id_rho,id_integ
      integer :: i,id,iertmp,idlim

      real*8 :: raxis,zaxis
      real*8, dimension(:), allocatable :: zrho

      real*8, parameter :: zero = 0.0d0
      real*8, parameter :: rhobdy = 1.0d0

      logical :: sol
      !-------------------------------

      ier=0

      call xplasma_common_ids(s,ier, id_R=id_R,id_Z=id_Z)
      if(ier.ne.0) return
 
      call xplasma_jacheck(s,idrho,ier)
      if(ier.ne.0) return

      !--------------------------------------------------

      !  get magnetic axis... (stored as xplasma list)

      call xplasma_mag_axis(s,ier, raxis=raxis, zaxis=zaxis)
      if(ier.ne.0) return

      !  get R and Z min and max (stored as list)

      call xplasma_RZminmax(s,rhobdy,ier)
      if(ier.ne.0) return

      !--------------------------------------------------
      !  create area and volume
      !  --converted to surface integral by Gauss scheme...

      call xplasma_author_set(s,xplasma_xmhd,iertmp)

      call xplasma_prof_info(s,id_R,iertmp, gridId1=idx1, gridId2=idx2)

      call xplasma_grid_info(s,idx1,iertmp, coord=idc1)

      if(idc1.eq.xplasma_rho_coord) then
         id_rho=idx1
      else
         id_rho=idx2
      endif

      allocate(zrho(inrho))

      call xplasma_grid(s,id_rho,zrho,iertmp)

      call xplasma_create_integ(s,'__RZ_RHO_INTEGRATOR', &
           zrho,id_integ,ier, cache_enable=.TRUE.)

      if(ier.eq.0) then
         call xplasma_label_item(s,id_integ,iertmp, &
              label='flux surface integrals, rho grid of {R,Z} surfaces.')
         ier = max(ier,iertmp)
      endif

      deallocate(zrho)

      call xplasma_author_clear(s,xplasma_xmhd,iertmp)

      if(ier.eq.0) then
         call xplasma_find_item(s,'__LIMITER',idlim,iertmp,nf_noerr=.TRUE.)
         if(idlim.gt.0) then
            call xplasma_limcheck(s,ier)
         endif
      endif

      if(ier.eq.0) then
         call xplasma_global_info(s,ier, scrapeoff=sol)
         if(sol) then

            sp => s
            call eqi_xtrz(ier)
            
         endif
      endif

    end subroutine xplasma_rzcheck

    !=======================================================================
    subroutine xplasma_vol1(s,ans,ier)
      !  ** public **

      !  return total plasma volume

      type (xplasma), pointer :: s
      real*8, intent(out) :: ans
      integer, intent(out) :: ier

      !------------
      real*8 xv(1),ansv(1)
      !------------

      xv=1

      call xplasma_voln(s,xv,ansv,ier)

      ans=ansv(1)

    end subroutine xplasma_vol1
    
    subroutine xplasma_voln(s,xs,ans,ier, ideriv)

      !  ** public **

      !  return Volume(x) or dVol/dx -- calculate if needed.

      type (xplasma), pointer :: s
      real*8, dimension(:), intent(in) :: xs      ! evaluation vector
      real*8, dimension(:), intent(out) :: ans    ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      integer, intent(in), optional :: ideriv     ! derivative option 0,1,2

      !-------------
      integer :: id_integ,id_vol,iertmp,id_R,id_Z,id_rho,inrho
      integer :: icounter, icount_rz
      real*8, dimension(:), allocatable :: zvol
      !-------------

      ans = 0.0d0

      call xplasma_find_item(s,"__RZ_RHO_INTEGRATOR", id_integ, ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_volume: numeric integration dataset unavailable.')
         return
      endif

      call xplasma_find_item(s,"XPLASMA_VOLUME",id_vol, ier, nf_noerr=.TRUE.)
      !  do not check ier -- id_vol=0, not found, possibility is expected here.

      !  id_R and id_Z must be defined, since __RZ_RHO_INTEGRATOR exists...

      call xplasma_common_ids(s,ier, id_R=id_R, id_Z=id_Z)
      if(ier.ne.0) return

      if(id_vol.gt.0) then

         !  make sure that volume info is not stale

         call xplasma_prof_info(s,id_R,ier, counter=icount_rz)
         call xplasma_prof_info(s,id_Z,ier, counter=icounter)
         icount_rz=max(icount_rz,icounter)
         call xplasma_prof_info(s,id_vol,ier, counter=icounter)

         if(icounter.lt.icount_rz) id_vol=0
      endif

      if(id_vol.eq.0) then

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         do

            call xplasma_prof_info(s,id_R,ier, gridId2=id_rho)
            if(ier.ne.0) exit

            call xplasma_grid_size(s,id_rho,inrho,ier)
            if(ier.ne.0) exit

            allocate(zvol(inrho))

            call xplasma_rho_zonint(s,id_integ,'VOL',zvol,ier)
            if(ier.ne.0) exit

            call xplasma_create_1dprof(s,'XPLASMA_VOLUME', &
                 id_rho,zvol,id_vol,iertmp, &
                 ispline=2,ibca=1,zbca=zero, &
                 label='Plasma Volume',units='m^3')
            if(ier.ne.0) id_vol=0

            deallocate(zvol)

            exit
         enddo

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)

      endif

      if(id_vol.eq.0) return

      call xplasma_eval_prof(s,id_vol,xs,ans,ier, ideriv1=ideriv)

    end subroutine xplasma_voln

    !=======================================================================
    subroutine xplasma_area1(s,ans,ier)
      !  ** public **

      !  return total plasma poloidal plane cross sectional area

      type (xplasma), pointer :: s
      real*8, intent(out) :: ans
      integer, intent(out) :: ier

      !------------
      real*8 xv(1),ansv(1)
      !------------

      xv=1

      call xplasma_arean(s,xv,ansv,ier)

      ans=ansv(1)

    end subroutine xplasma_area1
    
    subroutine xplasma_arean(s,xs,ans,ier, ideriv)

      !  ** public **

      !  return plasma poloidal plane cross sectional area
      !  return Area(x) or dArea/dx -- calculate if needed.

      type (xplasma), pointer :: s
      real*8, dimension(:), intent(in) :: xs      ! evaluation vector
      real*8, dimension(:), intent(out) :: ans    ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      integer, intent(in), optional :: ideriv     ! derivative option 0,1,2

      !-------------
      integer :: id_integ,id_area,iertmp,id_R,id_Z,id_rho,inrho
      integer :: icounter, icount_rz
      real*8, dimension(:), allocatable :: zarea
      !-------------

      ans = 0.0d0

      call xplasma_find_item(s,"__RZ_RHO_INTEGRATOR", id_integ, ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_area: numeric integration dataset unavailable.')
         return
      endif

      call xplasma_find_item(s,"XPLASMA_AREA",id_area, ier, nf_noerr=.TRUE.)
      !  do not check ier -- id_area=0, not found, is expected here.

      !  id_R and id_Z must be defined, since __RZ_RHO_INTEGRATOR exists...

      call xplasma_common_ids(s,ier, id_R=id_R, id_Z=id_Z)
      if(ier.ne.0) return

      if(id_area.gt.0) then

         !  make sure that volume info is not stale

         call xplasma_prof_info(s,id_R,ier, counter=icount_rz)
         call xplasma_prof_info(s,id_Z,ier, counter=icounter)
         icount_rz=max(icount_rz,icounter)
         call xplasma_prof_info(s,id_area,ier, counter=icounter)

         if(icounter.lt.icount_rz) id_area=0
      endif

      if(id_area.eq.0) then

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         do

            call xplasma_prof_info(s,id_R,ier, gridId2=id_rho)
            if(ier.ne.0) exit

            call xplasma_grid_size(s,id_rho,inrho,ier)
            if(ier.ne.0) exit

            allocate(zarea(inrho))

            call xplasma_rho_zonint(s,id_integ,'AREA',zarea,ier)
            if(ier.ne.0) exit

            call xplasma_create_1dprof(s,'XPLASMA_AREA', &
                 id_rho,zarea,id_area,iertmp, &
                 ispline=2,ibca=1,zbca=zero, &
                 label='Plasma Area',units='m^2')
            if(ier.ne.0) id_area=0
            exit
         enddo

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)
      endif

      if(id_area.eq.0) return

      call xplasma_eval_prof(s,id_area,xs,ans,ier, ideriv1=ideriv)

    end subroutine xplasma_arean

    !=======================================================================
    subroutine xplasma_gen_boozer(s,id_iboozer,id_zetacor,ierr)

      !  **public**

      !  compute quantities related to Boozer coordinates

      !  compute two 2d spline profiles:
      !     I(theta,rho) (T*m)
      !     nu(theta,rho) (rad) = zeta - phi; Boozer toroidal coord - phi

      !  The formulas are:
      !
      !     (dnu/dtheta) = q - (d(Lp)/dtheta)*g/(R**2*Bpol)
      !     I = Bpol*(d(Lp)/dtheta) - g*(dnu/dtheta)
      !
      !     g = R*Bphi, T*m; q = d[toroidal_flux]/d[poloidal_flux]
      !     Psi = [poloidal_flux]/2pi, Wb/rad.
      !     Lp = poloidal path length around flux surface, m.
      !     (rho,theta) are the radial and poloidal flux coordinates
      !                 used by xplasma.
      !     Bphi and Bpol the (axisymmetric) toroidal and poloidal fields (T).
      !
      !  The above formulas will satisfy covariant and contravariant
      !  formulations of B:
      !    verifying...
      ! 
      !  B = g*grad(zeta) + I*grad(theta) + delta*grad(Psi)
      !         delta*grad(Psi) cancels the components of B*grad(zeta)
      !         and I*grad(theta) in the grad(Psi) direction;
      !             grad(zeta) = grad(phi) + grad(nu)
      !             component of grad(nu) = (dnu/dtheta)*(dtheta/d(Lp))
      !                 (tangent to the flux surface in the poloidal direction)
      !                 (the normal component cancelled by delta*grad(Psi))
      !             component of grad(theta): (dtheta/d(Lp)).
      !                 (tangent to the flux surface in the poloidal direction)
      !  Bphi = g*grad(phi)  -- OK
      !  Bpol = g*(dnu/dtheta)*(dtheta/d(Lp)) + I*(dtheta/d(Lp))
      !       = Bpol when formula for I is substituted -- OK
      !
      !  B = q*grad(Psi)^grad(theta) + grad(zeta)^grad(Psi)
      !  B = q*grad(Psi)^grad(theta) + grad(phi)^grad(Psi) + grad(nu)^grad(Psi)
      !  Bpol = grad(phi)^grad(Psi)  -- OK
      !  Bphi = q*grad(Psi)*(dtheta/d(Lp)) - 
      !                         grad(Psi)*(dnu/dtheta)*(dtheta/d(Lp))
      !
      !  Bphi = (q - (dnu/dtheta))*grad(Psi)*(dtheta/d(Lp))
      !       = (d(Lp)/dtheta)*(g/(R**2*Bpol))*(R*Bpol)*(dtheta/d(Lp))
      !       = g/R  -- OK

      !--------------------------------------
      !  arguments:

      type (xplasma), pointer :: s
      integer, intent(out) :: id_iboozer   ! I spline id
      integer, intent(out) :: id_zetacor   ! nu spline (zetacor, zeta-phi)
      integer, intent(out) :: ierr   ! completion status (0=OK)
      !--------------------------------------
      !  local:
      integer :: id_integ,idi2,id_R,id_g,id_th,id_rho,iertmp
      integer :: ids(5),iderivs(5),id_Z,id_BR,id_BZ,id_nu
      integer :: iderivs_0x(5),iderivs_0t(5)
      integer :: inth,inrho,ith,irho,ith2
      real*8, dimension(:), allocatable :: thgrid,rhogrid
      real*8, dimension(:), allocatable :: qovg    ! q/g from integration
      real*8, dimension(:), allocatable :: gfun    ! g
      real*8, dimension(:), allocatable :: wk1,rhowk
      real*8, dimension(:,:), allocatable :: qovg_2d,fiwk
      real*8, dimension(:,:), allocatable :: ii,nu,wk5
      real*8, dimension(:), pointer :: thvec,wth
      real*8 :: rhomin
      !--------------------------------------
      !  find 1d integrator

      ierr=0

      call xplasma_profId(s,'__I_BOOZER',id_iboozer)
      call xplasma_profId(s,'__NU_BOOZER',id_zetacor)
      if(min(id_iboozer,id_zetacor).gt.0) return

      !  Need to evaluate these profiles and set up their splines...

      call xplasma_find_item(s,"__RZ_RHO_INTEGRATOR", id_integ, ierr)
      if(ierr.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_gen_boozer: numeric integration dataset unavailable.')
         return
      endif

      !  find R and Z and g and BR and BZ
      call xplasma_common_ids(s,iertmp,id_R=id_R,id_g=id_g,id_Z=id_Z, &
           id_BR=id_BR,id_BZ=id_BZ)

      ids(1)=id_BR
      ids(2)=id_BZ
      ids(3)=id_R
      ids(4)=id_Z
      ids(5)=id_R  ! for now

      ! ids(5) will be "nu" after it is computed

      iderivs(1:2)=0    ! no derivative
      iderivs(3:5)=1    ! d/dtheta

      ! for axis limit:  d/dtheta:
      iderivs_0t(1:2)=0; iderivs_0t(3:4)=1; iderivs_0t(5)=0

      ! for axis limit:  d/drho:
      iderivs_0x(1:2)=1; iderivs_0x(3:4)=1; iderivs_0x(5)=0

      !  find theta and rho -- grids used for equilibrium

      call xplasma_prof_info(s,id_R,ierr, gridId1=id_th, gridId2=id_rho)
      if(ierr.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_gen_boozer: could not find grids from "R"')
         return
      endif

      !  grid sizes

      call xplasma_grid_size(s,id_th,inth,iertmp)
      call xplasma_grid_size(s,id_rho,inrho,iertmp)

      !  copy of theta grid

      allocate(thgrid(inth),rhogrid(inrho))

      call xplasma_grid(s,id_rho,rhogrid,ierr)
      if(ierr.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_gen_boozer: could not fetch rho grid from "R"')
         deallocate(rhogrid,thgrid)
         return
      endif

      call xplasma_grid(s,id_th,thgrid,ierr)
      if(ierr.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_gen_boozer: could not fetch theta grid from "R"')
         deallocate(rhogrid,thgrid)
         return
      endif

      !  extend to 2d integrator

      call xplasma_author_set(s,xplasma_xmhd,iertmp)

      call xplasma_augment_integ(s,"__RZ_RHOTH_INTEGRATOR", id_integ, &
           thgrid,idi2,ierr)

      call xplasma_author_clear(s,xplasma_xmhd,iertmp)

      if(ierr.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_gen_boozer: failed to setup (rho,th) integrator.')
         deallocate(rhogrid,thgrid)
         return
      endif

      !  allocate work arrays

      allocate(qovg(inrho),gfun(inrho))
      allocate(qovg_2d(inth-1,inrho),wk1(inth),wk5(inth,5))
      allocate(ii(inth,inrho),nu(inth,inrho))

      !  explicit integration for axis ...L'Hopital rule expression, avoid 0/0:
      !  (dl/dtheta)/(R**2*Bpol) -> (d2l/dtheta*drho)/(R**2*d(Bpol)/dtheta)

      call xplasma_integ_info(s,idi2,iertmp, &
           thvec_ptr=thvec, thwts_ptr=wth, rhomin=rhomin)
      
      allocate(fiwk(size(thvec),5),rhowk(size(thvec))); rhowk=ZERO

      do
         !  compute integral contributions over theta grid

         ierr=0
         call xplasma_2d_zonint(s,idi2,'BOOZER_QOVG',qovg_2d,iertmp)
         if(iertmp.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_gen_boozer: BOOZER_QOVG integration failed.')
            ierr=max(ierr,iertmp)
         endif

         if(ierr.ne.0) exit

         ! fill in axial data

         call xplasma_eval_prof(s,ids, &
              xplasma_theta_coord,thvec, xplasma_rho_coord, rhowk, &
              fiwk, ierr, ideriv1s=iderivs_0t, ideriv2s=iderivs_0x)
         if(ierr.ne.0) exit

         ith=1
         do ith2=1,size(thvec)
            if(ith.lt.inth) then
               if(thvec(ith2).gt.thgrid(ith+1)) ith=ith+1
            endif
            qovg_2d(ith,1)=qovg_2d(ith,1) + wth(ith2) * &
                 sqrt(fiwk(ith2,3)**2+fiwk(ith2,4)**2)/ &
                 (sqrt(fiwk(ith2,1)**2+fiwk(ith2,2)**2)* &
                 C2PI*fiwk(ith2,5)**2)
         enddo

         call xplasma_eval_prof(s,id_G,rhogrid,gfun,iertmp)

         qovg=ZERO
         do ith=1,inth-1
            qovg = qovg + qovg_2d(ith,1:inrho)     ! (q/g)
         enddo

         nu(1,1:inrho)=ZERO
         nu(inth,1:inrho)=ZERO

         !  ii = I*thB    -- periodic derivative
         !  nu = zeta-phi -- strictly periodic; nu(theta(1))=nu(theta(inth))=0

         do irho=1,inrho
            do ith=2,inth-1
               nu(ith,irho)=nu(ith-1,irho) + &
                    qovg(irho)*gfun(irho)*(thgrid(ith)-thgrid(ith-1)) - &
                    qovg_2d(ith-1,irho)*gfun(irho)*C2PI
            enddo
         enddo

         ierr=0
         call xplasma_author_set(s,xplasma_xmhd,iertmp)
         ierr=max(ierr,iertmp)

         !  create profiles with default Bdy Conds in rho, periodic in theta

         call xplasma_create_2dprof(s,'__NU_BOOZER', &
              id_th,id_rho,nu,id_zetacor,iertmp, &
              ispline=2, &
              label='Toroidal coordinate adjust: Zeta_Boozer - Phi', &
              units='rad')
         ierr=max(ierr,iertmp)

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)
         ierr=max(ierr,iertmp)

         if(ierr.ne.0) then
            call xplasma_errmsg_append(s, &
                 'xplasma_gen_boozer: 2d nu profile creation failure.')
            exit
         endif

         id_nu=id_zetacor
         ids(5)=id_nu

         !  compute I(theta,rho)

         do irho=1,inrho
            wk1 = rhogrid(irho)
            call xplasma_eval_prof(s,ids, &
                 xplasma_theta_coord,thgrid, xplasma_rho_coord,wk1, &
                 wk5, ierr, ideriv1s=iderivs)
            if(ierr.ne.0) exit
            do ith = 1,inth
               ii(ith,irho) = sqrt((wk5(ith,1)**2+wk5(ith,2)**2)* &
                    (wk5(ith,3)**2+wk5(ith,4)**2)) - gfun(irho)*wk5(ith,5)
            enddo
         enddo

         if(ierr.eq.0) then
            call xplasma_author_set(s,xplasma_xmhd,iertmp)
            ierr=max(ierr,iertmp)

            !  create profiles with default rho boundary conditions;
            !  automatic periodic boundary condition in theta direction.

            call xplasma_create_2dprof(s,'__I_BOOZER', &
                 id_th,id_rho,ii,id_iboozer,iertmp, &
                 ispline=2, &
                 label='Boozer I(theta,rho) profile', &
                 units='Amps*mu0/2pi')
            ierr=max(ierr,iertmp)

            call xplasma_author_clear(s,xplasma_xmhd,iertmp)
            ierr=max(ierr,iertmp)
         endif

         if(ierr.ne.0) then
            call xplasma_errmsg_append(s, &
                 'xplasma_gen_boozer: 2d I profile creation failure.')
            exit
         endif

         exit

      enddo
      deallocate(qovg,gfun,rhogrid,thgrid)
      deallocate(qovg_2d,ii,nu,wk1,wk5)
      deallocate(fiwk,rhowk)

    end subroutine xplasma_gen_boozer

    !=======================================================================
    subroutine xplasma_jacheck(s,idrho,ierr)
 
      !  ** public **

      ! check for singular or near singular Jacobian
      ! axisymmetric equilibrium only
 
      !-------------------------------------------
      !  arguments
 
      type (xplasma), pointer :: s
      integer, intent(in) :: idrho   ! flux surface grid
      integer, intent(out) :: ierr   ! completion status (0=OK)
 
      !-------------------------------------------
      !  local stuff
 
      logical allzones_ok,badzone
      integer inrow
      integer inth,nrho
      integer izone,istat,iertmp
 
      real*8 zrho,zth,detj,detjav,ajac_maxvar

      logical, dimension(:), allocatable :: zone_ok ! (nrho-1)
      real*8, dimension(:), allocatable :: xrho     ! (nrho)
 
      real*8  zscreen1
 
      real*8, parameter :: c2pi = 6.2831853071795862D+00
      real*8, parameter :: czero = 0.0d0

      character*80 zmsg1

      !-------------------------------------------
      !  executable code
 
      !  first screen jacobian values at modest resolution
      !  if this shows a lot of variation, try again at higher resolution
 
      inrow=4
      inth=100
 
      allzones_ok=.TRUE.
      badzone=.FALSE.
 
      call xplasma_grid_size(s,idrho,nrho,ierr)
      if(ierr.ne.0) return

      allocate(zone_ok(nrho-1),xrho(nrho))
      zone_ok=.FALSE.

      call xplasma_grid(s,idrho,xrho,ierr)
      if(ierr.ne.0) return
 
      call xplasma_global_info(s,iertmp, ajac_maxvar=ajac_maxvar)

      zscreen1=4*ajac_maxvar

      do izone=1,nrho-1
 
         call jac_screen(izone,zscreen1,istat,zrho,zth,detj,detjav)
         if(istat.eq.0) then
            zone_ok(izone)=.TRUE.
         else if(istat.eq.1) then
            allzones_ok=.FALSE.   ! det(J) stays positive but varies too much
         else
            allzones_ok=.FALSE.   ! det(J) changes sign
            badzone=.TRUE.
            call jac_error('changes sign')
         endif
 
      enddo
 
      if(allzones_ok) then
         ierr=0
      else
         if(badzone) then
            ierr=2002
         endif
      endif
 
      ! result ambiguous
      ! some zones may not be OK -- look with higher resolution
      
      if((.not.badzone).and.(.not.allzones_ok)) then

         inrow=40
         inth=1000
 
         allzones_ok=.TRUE.

         do izone=1,nrho-1
 
            if(.not.zone_ok(izone)) then
 
               call jac_screen(izone,ajac_maxvar,istat,zrho,zth,detj,detjav)
               if(istat.eq.0) then
                  zone_ok(izone)=.true.
               else if(istat.eq.1) then
                  allzones_ok=.FALSE.
                  call jac_error( &
                       'too small -- may need to adjust tol via xplasma_ajac_maxvar_set')
                  call xplasma_errmsg_append(s, &
                       ' f77 style maxvar adjust call: eqi_jacheck_maxvar_set')
               else
                  allzones_ok=.FALSE.
                  call jac_error('changes sign')
               endif
 
            endif
 
         enddo
 
         if(allzones_ok) then
            ierr=0
         else
            ierr=2002
         endif
      endif

      if(ierr.eq.2002) then
         write(zmsg1,'(a,1x,1pe13.6)') &
              'max variation parameter (ajac_maxvar):',ajac_maxvar
         call xplasma_errmsg_append(s,zmsg1)
      endif

      deallocate(xrho,zone_ok)

    contains
 
      subroutine jac_error(zmsg)
        character*(*) zmsg      ! error message about det(J)

        character*128 zmsg2

        call xplasma_errmsg_append(s, &
             '===========================================')
        call xplasma_errmsg_append(s, &
             ' Xplasma (xplasma_jacheck): det(J) '//trim(zmsg)//'.')
        zmsg2=' '
        write(zmsg2,1001) zrho,zth,detj
        call xplasma_errmsg_append(s,trim(zmsg2))
        zmsg2=' '
        write(zmsg2,1002) detjav
        call xplasma_errmsg_append(s,trim(zmsg2))
        call xplasma_errmsg_append(s, &
             '===========================================')

1001    format(' ==> at rho=',1pd12.5,' theta=',1pd12.5,' det(J)=',1pd12.5)
1002    format('     average det(J) in vicinity: ',1pd12.5)
 
      end subroutine jac_error
 
      subroutine jac_screen(izone,zscreen,istat,zrho,zth,detj,detjav)
        integer,intent(in) :: izone    ! rho(izone:izone+1)
        real*8, intent(in) :: zscreen  ! max allowed rel. variation
 
        integer,intent(out) :: istat  ! 0:ok,1:too much variation,2:sign change
        real*8, intent(out) :: zrho,zth ! location of min det(J)
        real*8, intent(out) :: detj,detjav ! min det(J), avg det(J) in region
 
        !--------------
        real*8, dimension(:), allocatable :: zrhovec,zthvec,zphvec,zdetjvec
        real*8 zrhocur,zrho1,zrho2,zsign,zdetjmax
        integer ierr,ith,ir
        integer ithmax,irmax
        !--------------
 
        allocate(zrhovec(inth),zthvec(inth),zphvec(inth),zdetjvec(inth))
 
        zphvec = 0
        do ith=1,inth
           zthvec(ith)=(ith-1)*c2pi/inth
        enddo
 
        zrho1=xrho(izone)
        zrho2=xrho(izone+1)
 
        istat=0
 
        do ir=1,inrow
           zrhocur=(zrho1*(inrow-ir)+zrho2*ir)/inrow  ! zrho1+drho to zrho2
           zrhovec=zrhocur
 
           call xplasma_rzjac(s,zrhovec,zthvec,ierr, rzdetj=zdetjvec)
 
           if(izone.eq.1) then
              !  1/rho scaling of det(J), izone=1 only.
              zdetjvec=(zrho2/(zrhocur-xrho(izone)))*zdetjvec
           endif
 
           if(ierr.ne.0) then
            
              call xplasma_errmsg_append(s, &
                   ' ?? unexpected xplasma_rzjac error in xplasma_jacheck.')
              zrho=0; zth=0; detj=0; detjav=0; istat=2
              exit
           endif
 
           if(ir.eq.1) then
 
              !  first row; init avgs & min/max scans
              !  determine predominant sign of det(J)
 
              detjav=sum(zdetjvec)/(inth*inrow)
              if(detjav.lt.0.0) then
                 zsign=-1
                 zdetjvec=zdetjvec*zsign
                 detjav=detjav*zsign
              else if(detjav.gt.0.0) then
                 zsign=1
              else
                 call xplasma_errmsg_append(s, &
                      ' ?? xplasma_rzjac returned zeroes in xplasma_jacheck.')
                 zrho=0; zth=0; detj=0; detjav=0; istat=2
                 exit
              endif
 
              zrho=zrhovec(1)
              zth=zthvec(1)
              detj=zdetjvec(1)
 
              ithmax=1
              irmax=1
              zdetjmax=detj
 
           else
              zdetjvec=zdetjvec*zsign
              detjav = detjav + sum(zdetjvec)/(inth*inrow)
           endif
 
           do ith=1,inth
              if(zdetjvec(ith).gt.zdetjmax) then
                 zdetjmax=zdetjvec(ith)
                 ithmax=ith
                 irmax=ir
              endif
              if(zdetjvec(ith).lt.detj) then
                 zrho=zrhovec(ith)
                 zth=zthvec(ith)
                 detj=zdetjvec(ith)
              endif
           enddo
 
        enddo
 
        if(detj.le.czero) then
           istat=2
        else if((detj/zdetjmax).lt.zscreen) then
           istat=1
        endif
 
        deallocate(zrhovec,zthvec,zphvec,zdetjvec)
 
      end subroutine jac_screen
 
    end subroutine xplasma_jacheck

    !=======================================================================
    subroutine xplasma_eqcheck(s,iforce,icount,ier)

      !  **public**

      ! update equilibrium and increment counter IF all input profiles are
      ! available, and if Bphi and Jphi directions are set, and if one or more
      ! of the input profiles are new

      type (xplasma), pointer :: s
      logical, intent(in) :: iforce    ! TRUE to FORCE equilibrium update
      integer, intent(out) :: icount   ! counter set if equilibrium updated
      integer, intent(out) :: ier      ! completion status 0=OK

      ! use iforce=T to force an equilibrium update (i.e. calculation of
      ! field profiles) even if it appears that not all profiles have been
      ! modified.  If necessary data is missing, an error code will be set.

      !-------------------------------
      integer :: bccw,jccw,ieq_prev,inmax_inp
      integer :: iinp_R,iinp_RZ,iinp_Z,iinp_g,iinp_psi
      integer :: id_g,id_psi,id_R,id_Z,idtmp,ix,ith
      integer :: iertmp
      integer :: idx1,idx2,idrho,inrho,idth,inth,ispline,jcoord
      integer :: id5(5),idth5(5),idrho5(5)

      real*8,dimension(:), allocatable :: xrho,xth,xwk,fg,fpsi
      real*8,dimension(:), allocatable :: fr,fdrdrho,fdrdth,fdzdrho,fdzdth
      real*8,dimension(:), allocatable :: fdPsidrho,zdeti
      real*8,dimension(:,:), allocatable :: bmod,br,bz,f5

      real*8, parameter :: rhobdy = 1.0d0

      !-------------------------------

      ier=0
      icount=0

      call xplasma_global_info(s,ier, bphi_ccw=bccw,jphi_ccw=jccw, &
           eq_counter=ieq_prev)
      if(ier.ne.0) return

      if((bccw.eq.0).or.(jccw.eq.0)) then
         if(iforce) then
            ier=2003
            call xplasma_errmsg_append(s, &
                 '   Signs of toroidal field and current unspecified.')
         endif
         return   ! signs not specified
      endif

      !  see if all required profiles are present...

      call xplasma_common_ids(s,iertmp, &
           id_g=id_g,id_psi=id_psi,id_R=id_R,id_Z=id_Z)

      if(min(id_g,id_psi,id_R,id_Z).eq.0) then
         if(iforce) then
            ier=2003
            call xplasma_errmsg_append(s, &
                 '   Not all profiles in set {G,PSI,R,Z} are defined.')
         endif
         return ! don't have all needed profiles.
      endif

      !--------------------------------------------------
      !  see if any of the profiles are newer than the derived field profiles...
      inmax_inp=0

      call xplasma_prof_info(s,id_g,iertmp, counter=iinp_g)
      inmax_inp=max(inmax_inp,iinp_g)

      call xplasma_prof_info(s,id_psi,iertmp, counter=iinp_psi)
      inmax_inp=max(inmax_inp,iinp_psi)

      call xplasma_prof_info(s,id_R,iertmp, counter=iinp_R)

      call xplasma_prof_info(s,id_Z,iertmp, counter=iinp_Z)

      iinp_RZ=min(iinp_R,iinp_Z)

      inmax_inp=max(inmax_inp,iinp_RZ)

      if((.not.iforce).and.(inmax_inp.le.ieq_prev)) return

      !--------------------------------------------------
      !  *** an input profile has changed *** or iforce=T
  
      !  B field components... writing as xplasma_xmhd

      call xplasma_global_info(s,iertmp,rzOrder=ispline)

      call xplasma_author_set(s,xplasma_xmhd,iertmp)

      !  use same x axes as R(rho,theta)

      call xplasma_prof_info(s,id_R,iertmp,gridId1=idx1,gridId2=idx2)

      call xplasma_grid_info(s,idx1,iertmp, coord=jcoord)
      if(jcoord.eq.xplasma_rho_coord) then
         idrho=idx1
         idth=idx2
      else
         idrho=idx2
         idth=idx1
      endif

      !  get grid sizes

      call xplasma_grid_size(s,idrho,inrho,ier)
      if(ier.ne.0) return

      call xplasma_grid_size(s,idth,inth,ier)
      if(ier.ne.0) return

      !  get grids

      allocate(xrho(inrho),xth(inth),xwk(inth))

      call xplasma_grid(s,idrho,xrho,ier)
      if(ier.ne.0) return

      call xplasma_grid(s,idth,xth,ier)
      if(ier.ne.0) return

      do

         !  compute B_phi = g/R; use to start bmod...

         allocate(bmod(inth,inrho),fg(inrho),fr(inth))

         call xplasma_eval_prof(s,id_g,xrho,fg,ier)
         do ix=1,inrho
            xwk = xrho(ix)
            call xplasma_eval_prof(s,id_R,idth,xth,idrho,xwk,fr,ier)
            if(ier.ne.0) exit
            if(ix.eq.1) then
               fr=sum(fr)/inth
            endif
            bmod(1:inth,ix)=bccw*fg(ix)/fr(1:inth)  ! just BPHI at present
         enddo

         deallocate(fg,fr)
         if(ier.ne.0) then
            deallocate(bmod)
            exit
         endif

         call xplasma_create_prof(s,'BPHI',idth,idrho,bmod,idtmp,ier, &
              ispline=ispline,label='B_phi',units='T')

         !  compute [BR,BZ] = (1/R)grad(Psi) ^ e_phi

         allocate(br(inth,inrho),bz(inth,inrho),fr(inth))
         allocate(fdPsidrho(inrho),zdeti(inth),f5(inth,5))
         allocate(fdRdrho(inth),fdRdth(inth),fdZdrho(inth),fdZdth(inth))

         call xplasma_eval_prof(s,id_psi,xrho,fdPsidrho,ier, ideriv1=1)
         do ix=1,inrho

            if(ix.eq.1) then

               bR(1:inth,1)=0
               bZ(1:inth,1)=0
               bmod(1:inth,1)=abs(bmod(1:inth,1))
            else

               id5(1:3)=id_R ! fetch: R, dR/dtheta, dR/drho, dZ/dtheta, dZ/drho
               id5(4:5)=id_Z
               idth5 = (/ 0, 1, 0, 1, 0/)
               idrho5= (/ 0, 0, 1, 0, 1/)

               xwk=xrho(ix)
               call xplasma_eval_prof(s,id5,idth,xth,idrho,xwk,f5,ier, &
                    ideriv1s=idth5, ideriv2s=idrho5)
               if(ier.ne.0) exit

               fr = f5(1:inth,1)
               fdrdth = f5(1:inth,2)
               fdrdrho = f5(1:inth,3)
               fdzdth = f5(1:inth,4)
               fdzdrho = f5(1:inth,5)

               zdeti = -1/(fdrdth*fdzdrho-fdrdrho*fdzdth)

               bR(1:inth,ix) = -jccw*fdPsidrho(ix)*fdrdth*zdeti/fr
               bZ(1:inth,ix) = -jccw*fdPsidrho(ix)*fdzdth*zdeti/fr

               do ith=1,inth

                  bmod(ith,ix)=sqrt(bmod(ith,ix)*bmod(ith,ix) + &
                       br(ith,ix)*br(ith,ix) + bz(ith,ix)*bz(ith,ix))

               enddo

            endif
         enddo

         if(ier.eq.0) then
            call xplasma_create_prof(s,'BR',idth,idrho,bR,idtmp,ier, &
                 ispline=ispline,label='B_R',units='T')
            call xplasma_create_prof(s,'BZ',idth,idrho,bZ,idtmp,ier, &
                 ispline=ispline,label='B_Z',units='T')
            call xplasma_create_prof(s,'BMOD',idth,idrho,bmod,idtmp,ier, &
                 ispline=ispline,label='mod(B)',units='T')
         endif

         deallocate(br,bz,bmod,fr,fdpsidrho,fdrdrho,fdrdth,fdzdrho,fdzdth,zdeti,f5)
         exit

      enddo

      call xplasma_author_clear(s,xplasma_xmhd,iertmp)

      deallocate(xth,xrho,xwk)

      if(ier.eq.0) then
         icount = inmax_inp
         call xplasma_eqCount_incr(s,ier)
      endif
      
    end subroutine xplasma_eqcheck

    subroutine xplasma_gen_q(s,qname,iorder,id_q,ier)

      !  create a q profile from equilibrium information-- this involves a
      !  flux surface integral; the standard equilibrium grid is used for this.

      type (xplasma), pointer :: s

      character*(*), intent(in) :: qname  ! user selected name for profile
      integer, intent(in) :: iorder       ! fit order
      !  0 for piecewise linear, 1 for Hermite, 2 for Spline

      integer, intent(out) :: id_q        ! id of q profile created 
      integer, intent(out) :: ier         ! status code 0=OK

      !  if an error occurs, id_q=0 will be returned.

      !------------------------------
      integer :: id_rho,id_psi,id_g,id_R,id_integ,inrho,i,iertmp
      real*8, dimension(:), allocatable :: zrho,zg,zdvdrho,zr2i,zq
      real*8, dimension(:,:), allocatable :: zpsi
      real*8 :: zdpsmin
      integer :: ids(3),iderivs(3)
      real*8, parameter :: C2PI=6.2831853071795862D+00
      !------------------------------

      id_q=0
      call xplasma_common_ids(s,ier,id_psi=id_psi,id_g=id_g,id_R=id_R)
      if(ier.ne.0) return

      if(id_psi.eq.0) then
         ier=106
         call xplasma_errmsg_append(s,' ?xplasma_gen_q: no Psi(x) profile.')
      endif

      if(id_g.eq.0) then
         ier=106
         call xplasma_errmsg_append(s,' ?xplasma_gen_q: no R*Bphi(x) profile.')
      endif

      if(id_R.eq.0) then
         ier=106
         call xplasma_errmsg_append(s,' ?xplasma_gen_q: no R(theta,x).')
      endif
      if(ier.ne.0) return

      call xplasma_prof_info(s,id_R,ier, gridId2=id_rho)
      if(ier.ne.0) return

      call xplasma_grid_size(s,id_rho,inrho,ier)
      if(ier.ne.0) return

      allocate(zrho(inrho),zpsi(inrho,3),zg(inrho),zdvdrho(inrho),zr2i(inrho))
      allocate(zq(inrho))

      ids(1:2)=id_psi; iderivs(1)=0; iderivs(2)=1
      ids(3)=id_g; iderivs(3)=0

      do
         call xplasma_grid(s,id_rho,zrho, ier)
         if(ier.ne.0) exit

         !  get zpsi(:,1)=Psi(rho); zpsi(:,2)=d(psi)/d(rho).
         !    also putting zpsi(:,3)=g(rho) ***

         call xplasma_eval_prof(s,ids,zrho,zpsi,ier, ideriv1s=iderivs)
         if(ier.ne.0) exit

         !  force dPsi/drho > 0 everywhere...

         do i=2,inrho
            zdpsmin=zpsi(i,1)/zrho(i)/10
            if(zdpsmin.gt.0) exit
         enddo

         zpsi(2:inrho,2)=max(zdpsmin,zpsi(2:inrho,2))

         zg=zpsi(1:inrho,3)  ! and copy g

         !  get (dV/drho); this also sets up the flux surface integrator

         call xplasma_voln(s,zrho,zdvdrho,ier, ideriv=1)
         if(ier.ne.0) exit

         call xplasma_find_item(s,"__RZ_RHO_INTEGRATOR", id_integ, ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_gen_q: xplasma_voln call failed to create integrator.')
            exit
         endif

         call xplasma_rho_zonint(s,id_integ,'<1/R^2>S',zr2i,ier)
         if(ier.ne.0) exit

         zq(2:inrho)=zdvdrho(2:inrho)*zr2i(2:inrho)* &
              zg(2:inrho)/(zpsi(2:inrho,2)*C2PI*C2PI)
         zq(1)=(4*zq(2)-zq(3))/3  ! simple extrapolation to axis

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         call xplasma_create_prof(s,qname,id_rho,zq,id_q,ier, &
              ispline=iorder,label='q profile (xplasma_gen_q)',units=' ')

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)

         exit
      enddo
      deallocate(zrho,zpsi,zg,zdvdrho,zr2i,zq)

    end subroutine xplasma_gen_q

    subroutine xplasma_gen_torflx(s,tfname,iorder,id_tf,ier)

      !  create a profile of enclosed toroidal flux -- this involves a
      !  flux surface integral; the standard equilibrium grid is used for this.

      type (xplasma), pointer :: s

      character*(*), intent(in) :: tfname ! user selected name for profile
      integer, intent(in) :: iorder       ! fit order
      !  0 for piecewise linear, 1 for Hermite, 2 for Spline

      integer, intent(out) :: id_tf       ! id of tf profile created 
      integer, intent(out) :: ier         ! status code 0=OK

      !  if an error occurs, id_tf=0 will be returned.

      !------------------------------
      integer :: id_rho,id_g,id_R,id_integ,inrho,i,iertmp
      real*8, dimension(:), allocatable :: zrho,zg,zdvdrho,zr2i,ztf
      real*8 :: zdrho,zrhoav,zgav,zr2iav,zdvdrhoav
      real*8, parameter :: C2PI=6.2831853071795862D+00
      !------------------------------

      id_tf=0
      call xplasma_common_ids(s,ier,id_g=id_g,id_R=id_R)
      if(ier.ne.0) return

      if(id_g.eq.0) then
         ier=106
         call xplasma_errmsg_append(s,' ?xplasma_gen_q: no R*Bphi(x) profile.')
      endif

      if(id_R.eq.0) then
         ier=106
         call xplasma_errmsg_append(s,' ?xplasma_gen_q: no R(theta,x).')
      endif
      if(ier.ne.0) return

      call xplasma_prof_info(s,id_R,ier, gridId2=id_rho)
      if(ier.ne.0) return

      call xplasma_grid_size(s,id_rho,inrho,ier)
      if(ier.ne.0) return

      allocate(zrho(inrho),zg(inrho),zdvdrho(inrho),zr2i(inrho))
      allocate(ztf(inrho))

      do
         call xplasma_grid(s,id_rho,zrho, ier)
         if(ier.ne.0) exit

         !  get g(rho)

         call xplasma_eval_prof(s,id_g,zrho,zg,ier)
         if(ier.ne.0) exit

         !  get (dV/drho); this also sets up the flux surface integrator

         call xplasma_voln(s,zrho,zdvdrho,ier, ideriv=1)
         if(ier.ne.0) exit

         call xplasma_find_item(s,"__RZ_RHO_INTEGRATOR", id_integ, ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_gen_q: xplasma_voln call failed to create integrator.')
            exit
         endif

         call xplasma_rho_zonint(s,id_integ,'<1/R^2>S',zr2i,ier)
         if(ier.ne.0) exit

         ztf(1)=0
         do i=1,inrho-1
            zdrho=zrho(i+1)-zrho(i)
            zrhoav=(zrho(i)+zrho(i+1))/2
            zgav=(zg(i)+zg(i+1))/2
            zr2iav=(zr2i(i)+zr2i(i+1))/2
            zdvdrhoav=zdvdrho(i)+ &
                 (zrhoav+zrho(i))*(zrhoav-zrho(i))*(zdvdrho(i+1)-zdvdrho(i))/ &
                 (zrho(i+1)+zrho(i))*(zrho(i+1)-zrho(i))
            ztf(i+1)=ztf(i) + zdrho*zgav*zdvdrhoav*zr2iav/C2PI
         enddo

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         call xplasma_create_prof(s,tfname,id_rho,ztf,id_tf,ier, &
              ispline=iorder, &
              label='enclosed toroidal flux (xplasma_gen_torflx)',units='Wb')

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)

         exit
      enddo
      deallocate(zrho,zg,zdvdrho,zr2i,ztf)

    end subroutine xplasma_gen_torflx

    subroutine xplasma_gen_p(s,pname,iorder,id_p,ier)

      !  create a P profile from equilibrium information-- this involves use
      !  of flux surface integrals, to solve the surface averaged 
      !  Grad-Shafronov equation JxB = grad(P)

      type (xplasma), pointer :: s

      character*(*), intent(in) :: pname  ! user selected name for profile
      integer, intent(in) :: iorder       ! fit order
      !  0 for piecewise linear, 1 for Hermite, 2 for Spline

      integer, intent(out) :: id_p        ! id of p profile created 
      integer, intent(out) :: ier         ! status code 0=OK

      !  if an error occurs, id_q=0 will be returned.

      !------------------------------
      integer :: id_rho,id_psi,id_g,id_R,id_integ,inrho,i,iertmp
      real*8, dimension(:), allocatable :: zrho,zg,zr2i,zarea,zp
      real*8, dimension(:), allocatable :: zpsi,zri,zx2r2i,zb2,zitor
      real*8 :: zmu0,zga,zr2ia,zria,zdpsidx,zdgdx,zx2r2ia,zb2a,zja,zjdotb
      real*8 :: zdpdxm
      real*8 :: get_jdotb_from_j   ! source/comput utility function
      real*8, parameter :: C2PI=6.2831853071795862D+00
      !------------------------------

      id_p=0
      call xplasma_common_ids(s,ier,id_psi=id_psi,id_g=id_g,id_R=id_R)
      if(ier.ne.0) return

      if(id_psi.eq.0) then
         ier=106
         call xplasma_errmsg_append(s,' ?xplasma_gen_q: no Psi(x) profile.')
      endif

      if(id_g.eq.0) then
         ier=106
         call xplasma_errmsg_append(s,' ?xplasma_gen_q: no R*Bphi(x) profile.')
      endif

      if(id_R.eq.0) then
         ier=106
         call xplasma_errmsg_append(s,' ?xplasma_gen_q: no R(theta,x).')
      endif
      if(ier.ne.0) return

      call xplasma_prof_info(s,id_R,ier, gridId2=id_rho)
      if(ier.ne.0) return

      call xplasma_grid_size(s,id_rho,inrho,ier)
      if(ier.ne.0) return

      allocate(zrho(inrho),zpsi(inrho),zg(inrho),zr2i(inrho),zp(inrho))
      allocate(zri(inrho),zx2r2i(inrho),zarea(inrho),zitor(inrho))

      do
         call xplasma_grid(s,id_rho,zrho, ier)
         if(ier.ne.0) exit

         !  get zpsi(:)=Psi(rho); get zg(:)=g(rho)

         call xplasma_eval_prof(s,id_psi,zrho,zpsi,ier)
         if(ier.ne.0) exit

         call xplasma_eval_prof(s,id_g,zrho,zg,ier)
         if(ier.ne.0) exit

         !  get (dV/drho); this also sets up the flux surface integrator

         call xplasma_arean(s,zrho,zarea,ier)
         if(ier.ne.0) exit

         call xplasma_find_item(s,"__RZ_RHO_INTEGRATOR", id_integ, ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_gen_q: xplasma_voln call failed to create integrator.')
            exit
         endif

         call xplasma_rho_zonint(s,id_integ,'<1/R^2>S',zr2i,ier)
         if(ier.ne.0) exit

         call xplasma_rho_zonint(s,id_integ,'<1/R>S',zri,ier)
         if(ier.ne.0) exit

         call xplasma_rho_zonint(s,id_integ,'<GRAD(RHO)^2/R^2>S',zx2r2i,ier)
         if(ier.ne.0) exit

         call xplasma_rho_zonint(s,id_integ,'ITOR',zitor,ier)
         zp(inrho)=0

         zmu0 = C2PI*2.0d-7

         do i=inrho-1,1,-1
            !  g, <1/R^2>, <1/R>, dPsi/dx, dg/dx, <grad(x)^2/R^2>, <B^2>, j
            zga=(zg(i)+zg(i+1))/2
            zr2ia=(zr2i(i)+zr2i(i+1))/2
            zria=(zri(i)+zri(i+1))/2
            zdpsidx=(zpsi(i+1)-zpsi(i))/(zrho(i+1)-zrho(i))
            zdgdx=(zg(i+1)-zg(i))/(zrho(i+1)-zrho(i))
            zx2r2ia=(zx2r2i(i)+zx2r2i(i+1))
            zb2a=zga*zga*zr2ia + zdpsidx*zdpsidx*zx2r2ia
            zja=(zitor(i+1)-zitor(i))/(zarea(i+1)-zarea(i))

            !  <j> -> <j.B>
            zjdotb = get_jdotb_from_j(zga,zr2ia,zb2a,zria,zja)

            ! -dP/dx

            zdpdxm = zb2a*zdgdx/(zmu0*zga) + zdpsidx*zjdotb/zga

            zp(i) = zp(i+1) + zdpdxm*(zrho(i+1)-zrho(i))
         enddo

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         call xplasma_create_prof(s,pname,id_rho,zp,id_p,ier, &
              ispline=iorder,label='P profile (xplasma_gen_p)',units='Pascals')

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)

         exit
      enddo
      deallocate(zrho,zpsi,zg,zr2i,zri,zx2r2i,zarea,zp,zitor)

    end subroutine xplasma_gen_p
      
    subroutine xplasma_create_RZgrid(s,Rgrid,Zgrid,id_Rgrid,id_Zgrid,ier)

      !  create (R,Z) rectangular grid
      !  this should at least cover the plasma & limiters.
      !  the R grid should start at a point greater than zero.

      use eqi_rzbox_module

      type (xplasma), pointer :: s
      real*8, dimension(:), intent(in) :: Rgrid,Zgrid  ! R & Z grids
      integer, intent(out) :: id_Rgrid,id_Zgrid        ! ids returned
      integer, intent(out) :: ier

      !---------------------
      integer :: id_lim,itype,iwarn,inr,inz
      real*8 :: zrmin,zrmax,zzmin,zzmax
      character*128 zmsg
      !---------------------
      !  limiter must be defined first...

      call xplasma_find_item(s,'__LIMITER',id_lim,ier, nf_noerr=.TRUE.)
      if(ier.ne.0) return

      if(id_lim.eq.0) then
         ier=3004
         call xplasma_errmsg_append(s, &
              '  need to define limiter prior to specifying R & Z rectangle grid.')
         return
      endif

      inr=size(Rgrid)
      inz=size(Zgrid)
      call xplasma_lim_RZminmax(s,zRmin,zRmax,zZmin,zZmax,ier,itype=itype)
      if(ier.ne.0) return

      if((zRmin.lt.Rgrid(1)).or.(zRmax.gt.Rgrid(inr)).or. &
           (zZmin.lt.Zgrid(1)).or.(zZmax.gt.Zgrid(inz))) then
         iwarn=1
      else
         iwarn=0
      endif

      if(Rgrid(1).le.0.0d0) iwarn=2

      call xplasma_create_grid(s, &
           '__Rgrid',xplasma_R_coord,Rgrid,id_Rgrid,ier, &
           label="R grid")
      if(ier.ne.0) return

      call xplasma_create_grid(s, &
           '__Zgrid',xplasma_Z_coord,Zgrid,id_Zgrid,ier, &
           label="Z grid")
      if(ier.ne.0) return

      if(iwarn.ne.0) then
         ier=3006
         if(iwarn.eq.1) then
            zmsg=' '
            write(zmsg,*) ' Rmin & Rmax of limiter: ',zRmin,zRmax
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) ' Rmin & Rmax of grid: ',Rgrid(1),Rgrid(inr)
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) ' Zmin & Zmax of limiter: ',zZmin,zZmax
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) ' Zmin & Zmax of grid: ',Zgrid(1),Zgrid(inz)
            call xplasma_errmsg_append(s,zmsg)
         endif
      else

         sp => s
         call xplasma_misc_set(s,ier, scrapeoff=.TRUE.)
         call eqi_xtrz(ier)

      endif

    end subroutine xplasma_create_RZgrid

    subroutine xplasma_max_maptype(s,maptype,ier)

      !  return the maximum maptype: 2 if there is no SOL; 3 if R&Z grids
      !  for SOL exist, so that bilinear maps can be constructed.

      type(xplasma), pointer :: s

      integer, intent(out) :: maptype
      integer, intent(out) :: ier

      !--------------------
      logical :: sol,axisymm
      !--------------------

      maptype=0

      call xplasma_global_info(s,ier, scrapeoff=sol,axisymm=axisymm)
      if(ier.ne.0) then
         sol=.FALSE.
         axisymm=.TRUE.
         return
      endif

      if(.not.axisymm) then
         ier=2509
         return
      endif

      if(sol) then
         maptype=3
      else
         maptype=2
      endif

    end subroutine xplasma_max_maptype

    ! DMC Apr 2011: surfgeo wrapper imported into xplasma module...
    ! original comments copied here...

    ! compute various flux surface parameters derived from the position of 
    ! the flux surface on the R,Z grid
    !
    ! the midplane intercepts are defined as the major radiaus at the height
    ! of the centroid of the flux surface.
    !
    ! the rho value will be forced to be larger then  1.e-5*(rho_bdy-rho_axis)
    !
    ! notes from gelong:
    !c
    !c> Date: Wed, 22 Sep 1999 11:59:40 -0400
    !c> From: Glenn Bateman <bateman@fusion.physics.lehigh.edu>
    !c> 	Consider a toroidal magnetic surface (also called a flux 
    !c>              surface).
    !c> Use the following variables
    !c>
    !c> R_out   = major radius to the outboard edge of the flux surface
    !c>           that is, the largest value of major radius anywhere on the 
    !c>           surface
    !c> R_in    = major radius to the inboard edge (ie, smallest value of major
    !c>           radius)
    !c> R_top   = major radius to the top edge of the flux surface
    !c>           that is, to the highest point on the flux surface cross 
    !c>           section
    !c> R_bottom = major radius to the lowest point on the flux surface
    !c> height  = distance between the elevations of the highest point
    !c>           and the lowest point
    !c> R_dent_in  = major radius to the most indented inboard part of the 
    !c>           surface;
    !c>           in cases where the flux surface is bean shaped
    !c>           with a dent on the inboard edge
    !c>           Note: R_dent_in = R_in whenever the flux surface is not 
    !c>                 indented
    !c> R_dent_out = major radius to the most indented outboard part of the 
    !c>           surface
    !c>
    !c> r_minor = half width = (R_out - R_in ) / 2     
    !c>           (usually called minor radius)
    !c> r_major = ( R_out + R_in ) / 2
    !c>           (usually called major radius)
    !c>
    !c> elongation:	kappa = height / ( R_out - R_in )
    !c>
    !c> tiangularity:	delta = ( r_major - min ( R_top, R_bottom ) ) / r_minor
    !c>
    !c> indentation:	indent = max ( ( R_dent_in - R_in ) / rminor ,
    !c> 			       ( R_dent_out - R_out ) / rminor )
    !c>
    !c> 	Note: for normal convex surfaces, indent = 0.0
    !c> For bean-shaped surfaces with indentation on inboard edge, indent > 0
    !c> For bean-shaped surfaces with indentation on outboard edge, indent < 0
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
    !c> Mikkelsen squareness attempts to match moment representation
    !c>
    !c>    R = Rmajor + Rminor*cos(theta + asin(triang)sin(theta))
    !c>    Z = Rmajor + elong*Rminor*sin(theta + square*sin(2*theta))
    !c>
    !c> by finding the point Rmk = R(pi/4.), interpolating from the boundary
    !c> data to Zmk and then using the moment representation at Z(pi/4.) to
    !c> derive the squareness.
    !c>

    subroutine xplasma_surfgeo(s, rho, &
         elong, triang, indent, &
         zmidp, rin_midp, rout_midp, &
         ierr, &
         ntheta, phival, auxdata, idebug)

      !------------------------------------
      type (xplasma), pointer :: s
      !------------------------------------

      ! NOTE: all 1d arrays must be of the same size;
      ! this size must also match size(auxdata,2) if auxdata is present

      real*8, dimension(:), intent(in) :: rho ! radial coord. of flux surfaces


      real*8, dimension(:), intent(out) :: elong    ! surface elongation
      real*8, dimension(:), intent(out) :: triang   ! triangularity
      real*8, dimension(:), intent(out) :: indent   ! indentation
      
      real*8, dimension(:), intent(out) :: zmidp    ! height of centroid
                                                    ! of each flux surface
      real*8, dimension(:), intent(out) :: rin_midp ! inner midplane intercept 
                                                    ! at centroid height
      real*8, dimension(:), intent(out) :: rout_midp ! outer midplane intercept
                                                    ! at centroid height

      !----------------
      integer, intent(out) :: ierr        ! completion status; 0=OK
      !----------------

      ! OPTIONAL arguments...

      integer, intent(in), optional :: ntheta ! #theta pts to evaluate
                                              ! (default 400)

      real*8, dimension(:), intent(in), optional :: phival  ! toroidal coord.
                                              ! of flux surfaces

      real*8, dimension(:,:), intent(out), optional :: auxdata  ! aux. outputs
                                          ! (item index, surface index)
                                          !   decoding of item index:
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
                                          ! 16 -> lower Mikkelsen squareness
                                          ! 17 -> Rmk_upp
                                          ! 18 -> Zmk_upp
                                          ! 19 -> upper Mikkelsen squareness
                                          !   NOTE: if actual index dimension
                                          !         is <19, there is no error:
                                          !         less information is output

      logical, intent(in), optional :: idebug   ! .TRUE. for extra debugging output


      !-------------------------------------------
      !  local declarations

      integer, parameter :: R8=SELECTED_REAL_KIND(12,100)
      integer, parameter :: ntheta_def = 400    ! default ntheta
      real*8,  parameter :: ZERO       = 0.0d0
      real*8,  parameter :: TWOPI      = 6.2831853071795862d0
      real*8,  parameter :: PI4        = TWOPI/8.d0
      real*8,  parameter :: RHOZERO    = 1.D-5  ! minimum allowed rho

      character*128 :: mline    ! for messages

      logical :: jdebug         ! .true. for some debug output

      real*8,  allocatable :: phi(:)
      real*8,  allocatable :: zrho(:),zphi(:),ztheta(:),zr(:),zz(:)
      real*8,  allocatable :: ar(:,:), az(:,:), locdata(:,:)
      integer, allocatable :: nregion(:)       ! xplasma region code

      integer :: i,iv        ! loop variable
      integer :: itheta      ! number of poloidal points

      real*8  :: rrmin       ! minimum allowed rmin
      real*8  :: tr,tz       ! temp

      integer :: ivec        ! 1d array dimension
      integer :: ivec_min,ivec_max  ! for error checking...
      !-------------------------------------------
      !  executable code

      !  clear output variables...

      ierr = 0
      elong = ZERO
      triang = ZERO
      indent = ZERO

      zmidp = ZERO
      rin_midp = ZERO
      rout_midp = ZERO

      if(present(auxdata)) auxdata = ZERO

      !-------------------------------
      !  dimension size check

      ivec_min = size(rho)
      ivec_max = ivec_min

      if(present(phival)) then
         ivec_min = min(ivec_min,size(phival))
         ivec_max = max(ivec_max,size(phival))
      endif

      if(present(auxdata)) then
         ivec_min = min(ivec_min,size(auxdata,2))
         ivec_max = max(ivec_max,size(auxdata,2))
      endif

      ivec_min = min(ivec_min,size(elong))
      ivec_min = min(ivec_min,size(triang))
      ivec_min = min(ivec_min,size(indent))

      ivec_min = min(ivec_min,size(zmidp))
      ivec_min = min(ivec_min,size(rin_midp))
      ivec_min = min(ivec_min,size(rout_midp))

      ivec_max = max(ivec_max,size(elong))
      ivec_max = max(ivec_max,size(triang))
      ivec_max = max(ivec_max,size(indent))

      ivec_max = max(ivec_max,size(zmidp))
      ivec_max = max(ivec_max,size(rin_midp))
      ivec_max = max(ivec_max,size(rout_midp))

      if(ivec_min.ne.ivec_max) then
         ierr=524
         mline = 'xplasma_rzgeo:xplasma_surfgeo argument array sizes do not match.'
         call xplasma_errmsg_append(s,mline)
         return
      endif

      !-------------------------------
      ! OK... array sizes match

      ivec = ivec_min ! = ivec_max

      !-------------------------------
      !  initial allocations...

      mline = ' '
      jdebug = .false.
      if(present(idebug)) jdebug=idebug

      if(.NOT.present(ntheta)) then
         itheta = ntheta_def
      else
         if (ntheta>16) then
            itheta = ntheta
         else
            itheta = ntheta_def
         endif
      endif

      allocate(phi(ivec))
      if(present(phival)) then
         phi = phival
      else
         phi = ZERO
      endif

      allocate(ar(itheta,ivec), az(itheta,ivec))
      allocate(zrho(itheta),ztheta(itheta),zphi(itheta))
      allocate(zr(itheta), zz(itheta), nregion(itheta))

      if(.NOT.present(auxdata)) then
         allocate(locdata(19,ivec))
      endif

      !
      ! i=1,itheta -> theta=0.0
      do i = 1, itheta-1
         ztheta(i) = TWOPI*(i-1)/(itheta-1)  ! equal steps in theta
      end do
      ztheta(itheta)=0.D0

      rrmin = RHOZERO

      !
      ! all exits through 999
      !
      do iv = 1, ivec
         !
         ! --- grab R,Z points ---
         !
         zrho = max(rho(iv),rrmin)
         zphi = phi(iv)

         nregion=0

         call xplasma_ctrans(s,.TRUE.,ierr, &
              rho_in=zrho, theta_in=ztheta, phi_in=zphi, &
              R_out=zR, Z_out=zZ, nregion=nregion)

         if (ierr/=0) then
            write(mline,'(a,f10.5)') &
                 ' ?xplasma_surfgeo: xplasma_ctrans returned error, rho = ',zrho(1)
            call xplasma_errmsg_append(s,trim(mline))
            ierr=9999
            go to 999
         endif
     
         if (maxval(nregion)>1) then
            write(mline,'(a,f10.5)') &
                 ' ?xplasma_surfgeo: surface is beyond magnetic coordinates, rho = ',zrho(1)
            call xplasma_errmsg_append(s,trim(mline))
            ierr=9999
            go to 999
         endif

         zr(itheta) = zr(1)   ! i=itheta is the same point as i=1
         zz(itheta) = zz(1)
     
         ar(:,iv) = zr
         az(:,iv) = zz
      end do

      if(present(auxdata)) then
         call surfgeo(itheta, ivec, ar, az, &
              elong, triang, indent, zmidp, rin_midp, rout_midp, &
              size(auxdata,1), auxdata, &
              jdebug, mline, ierr)
      else
         call surfgeo(itheta, ivec, ar, az, &
              elong, triang, indent, zmidp, rin_midp, rout_midp, &
              size(locdata,1), locdata, &
              jdebug, mline, ierr)
      endif

      if(ierr/=0) then
         call xplasma_errmsg_append(s,mline)
         call xplasma_errmsg_append(s,'error detected in xplasma_rzgeo(xplasma_surfgeo)')
         ierr=9999
         goto 999
      endif

999   continue
      deallocate(ar,az,zrho,ztheta,zphi,zr,zz,nregion,phi)
      if(.NOT.present(auxdata)) then
         deallocate(locdata)
      endif

    end subroutine xplasma_surfgeo

end module xplasma_rzgeo
