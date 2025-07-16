module xplasma_ctran

  !  methods for defining coordinate transformations
  !    [Cartesian] <--> [Cylindric] <--> [Magnetic Flux Coordinates]

  !  related methods:
  !    magnetic axis location
  !    methods for extrema:  min/max R,Z of chosen flux surfaces

  use xplasma_obj
  implicit NONE

  private

  public :: xplasma_ctran1,xplasma_ctrann,xplasma_ctrans
  public :: xplasma_mag_axis
  public :: xplasma_RZminmax,xplasma_Bminmax
  public :: xplasma_RZminmax_plasma,xplasma_RZminmax_extended
  public :: xplasma_bdfind,xplasma_bdfind1,xplasma_bdfindn
  public :: xplasma_rzjac,xplasma_rzjac1,xplasma_rzjacn
  public :: xplasma_vtaub

  integer, parameter, public :: xp_imap_newton = 1
  integer, parameter, public :: xp_imap_polar = 2
  integer, parameter, public :: xp_imap_rzlinear = 3

  interface xplasma_ctrans

     module procedure xplasma_ctran1   ! scalar
     module procedure xplasma_ctrann   ! vector

  end interface

  interface xplasma_bdfind

     module procedure xplasma_bdfind1   ! scalar
     module procedure xplasma_bdfindn   ! vector

  end interface

  interface xplasma_rzjac

     module procedure xplasma_rzjac1   ! scalar
     module procedure xplasma_rzjacn   ! vector

  end interface

  real*8, parameter :: czero = 0.0d0
  real*8, parameter :: cone = 1.0d0
  real*8, parameter :: c2pi = 6.2831853071795862D+00
  real*8, parameter :: cpi = 3.1415926535897931D+00
  real*8, parameter :: ceps4 = 1.0d-4
  real*8, parameter :: ceps13= 1.0d-13

  contains

    subroutine xplasma_mag_axis(s,ier, raxis,zaxis,baxis)

      ! return mag. axis location (scalar, AXISYMMETRY assumed...)
      ! mod(B) can optionally be returned as well

      type (xplasma), pointer :: s
      integer, intent(out) :: ier

      real*8, intent(out), optional :: raxis,zaxis,baxis

      !------------------------------
      integer :: id_axis,id_R,id_Z,id_g,iertmp
      character*32 elemnames(2),elemlabels(2)
      real*8 elemvals(2),zg
      real*8, parameter :: ZERO = 0.0d0
      !------------------------------

      if(present(raxis)) raxis=0
      if(present(zaxis)) zaxis=0
      if(present(baxis)) baxis=0

      call xplasma_common_ids(s,ier, id_axis=id_axis)
      if(ier.ne.0) return

      if(id_axis.eq.0) then

         !  make the mag. axis list object

         call xplasma_ctrans(s,ier, rho_in=czero,theta_in=czero, &
              r_out=elemvals(1),z_out=elemvals(2))

         if(ier.ne.0) return

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         elemnames(1)='Raxis'
         elemnames(2)='Zaxis'

         elemlabels(1)='R of magnetic axis'
         elemlabels(2)='Z of magnetic axis'

         call xplasma_create_list(s,"MAG_AXIS", &
              elemnames,id_axis,ier, &
              label='magnetic axis location', units='m', &
              r8vals=elemvals, chvals=elemlabels)

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)

      else
         
         call xplasma_getList(s,id_axis,elemvals,ier)

      endif

      if(ier.eq.0) then
         if(present(raxis)) raxis=elemvals(1)
         if(present(zaxis)) zaxis=elemvals(2)
         if(present(baxis)) then
            call xplasma_common_ids(s,ier,id_g=id_g)
            if(ier.eq.0) then
               call xplasma_eval_prof(s,id_g,ZERO,zg,ier)
               baxis=zg/elemvals(1)
            endif
         endif
      endif

    end subroutine xplasma_mag_axis

    !======================================================================
    subroutine xplasma_ctran1(s,ier, &
         x_in,y_in,z_in, r_in,phi_in, rho_in,theta_in, &
         tol, maptype, ccw_theta,ccw_phi, i2pi, &
         x_out,y_out,z_out, r_out,phi_out, rho_out,theta_out, &
         nregion,cosphi,sinphi)

      ! SCALAR coordinate mapping interface

      type (xplasma), pointer :: s
      integer, intent(out) :: ier   ! completion code 0=OK

      !------------------------
      ! optional inputs...

      real*8, intent(in), optional :: x_in,y_in,z_in  ! Cartesian Coords in
      real*8, intent(in), optional :: r_in,phi_in   ! Cylindric Coords (z also)
      real*8, intent(in), optional :: rho_in,theta_in ! Flux Coords (phi also)

      !------------------------------------------------------------
      ! the following applies to (R,Z) -> (rho,theta) "inverse" map:

      real*8, intent(in), optional :: tol           ! map tolerance
      
      ! tol applies to Newton map (maptype = xp_imap_newton = default) ONLY

      ! (R_in,Z_in)->(rho_out,theta_out): 
      !         |R_in-R(rho_out,theta_out)| < tol*max(R_in,|Z_in|) 
      !         |Z_in-Z(rho_out,theta_out)| < tol*max(R_in,|Z_in|) 
      ! default: s%bdytol

      integer, intent(in), optional :: maptype

      !  = 1 = xp_imap_newton = default -- iterative Newton root finder
      !        with initial guess generator using polar map (slowest but 
      !        most accurate, accuracy can be controlled with tol.
      !  = 2 = xp_imap_polar -- polar map inversion method
      !        much faster, not as accurate as Newton map
      !  = 3 = xp_imap_rzlinear -- bilinear inverse map -- after 
      !        initialization (using polar map), the fastest, but less
      !        accurate still (depends on R,Z grid resolution).  This 
      !        method requires a user defined (R,Z) rectangle grid; if 
      !        none is available the polar map is used instead.

      !------------------------------------------------------------

      logical, intent(in), optional :: ccw_phi  ! Phi CCW view from above, default=T
      logical, intent(in), optional :: ccw_theta ! Theta CCW in poloidal plane, default=T

      integer, intent(in), optional :: i2pi ! angle coord normalization options
      !  i2pi=1 (default) returned value in range [0,2pi]
      !  i2pi=2 returned value in range [-pi,pi]

      !-------------------------
      ! optional outputs...

      real*8, intent(out), optional :: x_out,y_out,z_out  ! Cartesian Coords in
      real*8, intent(out), optional :: r_out,phi_out   ! Cyl. Coords (z also)
      real*8, intent(out), optional :: rho_out,theta_out ! Flux Coords (phi also)
      real*8, intent(out), optional :: cosphi,sinphi   ! phi sin & cos
      integer, intent(out), optional :: nregion        ! region code
      ! 0=not checked; 1=inside core plasma; 2=outside plasma but inside
      ! mapped region; 3=beyond mapped region.

      !-----------------
      logical :: car_in,cyl_in,flx_in
      logical :: car_out,cyl_out,flx_out
      integer :: intypes_in,intypes_out
      logical :: jccw_phi,jccw_theta

      logical :: car2cyl,cyl2flx
      logical :: flx2cyl,cyl2car

      logical :: sol
      integer :: imap_type

      real*8 :: ztol
      real*8 :: xv(1),yv(1),zv(1),rv(1),phiv(1),rhov(1),thetav(1)
      integer :: iregion(1),j2pi

      integer :: id_R,id_Z
      !-----------------

      !  clear outputs as needed

      ier = 0
      if(present(x_out)) x_out = 0
      if(present(y_out)) y_out = 0
      if(present(z_out)) z_out = 0

      if(present(r_out)) r_out = 0
      if(present(phi_out)) phi_out = 0

      if(present(rho_out)) rho_out = 0
      if(present(theta_out)) theta_out = 0

      !  set CCW flags & tolerance

      jccw_phi=.TRUE.
      jccw_theta=.TRUE.
      if(present(ccw_phi)) jccw_phi=ccw_phi
      if(present(ccw_theta)) jccw_theta=ccw_theta

      imap_type=xp_imap_newton

      call xplasma_global_info(s,ier, scrapeoff=sol)
      if(present(maptype)) then
         if(sol) then
            imap_type=max(xp_imap_newton,min(xp_imap_rzlinear,maptype))
         else
            imap_type=max(xp_imap_newton,min(xp_imap_polar,maptype))
         endif
      endif

      call xplasma_global_info(s,ier, bdytol=ztol)
      if(ier.ne.0) return
      if(present(tol)) ztol=tol
      ztol=max(ceps13,min(ceps4,ztol))

      j2pi=1
      if(present(i2pi)) j2pi=i2pi

      !  check input-- must be from just one of the three possible
      !  coordinate systems.

      !  check output-- there must be some output.

      car_in = .FALSE.
      cyl_in = .FALSE.
      flx_in = .FALSE.

      car_out = .FALSE.
      cyl_out = .FALSE.
      flx_out = .FALSE.

      if(present(x_in).or.present(y_in)) car_in = .TRUE.
      if(present(R_in)) cyl_in = .TRUE.
      if(present(rho_in).or.present(theta_in)) flx_in = .TRUE.

      if((.not.car_in).and.(.not.cyl_in)) then
         if(present(z_in)) cyl_in=.TRUE.
      endif

      if((.not.cyl_in).and.(.not.flx_in)) then
         if(present(phi_in)) cyl_in=.TRUE.
      endif

      if(present(x_out).or.present(y_out)) car_out = .TRUE.
      if(present(R_out)) cyl_out = .TRUE.
      if(present(rho_out).or.present(theta_out)) flx_out = .TRUE.

      if((.not.car_out).and.(.not.cyl_out)) then
         if(present(z_out)) cyl_out=.TRUE.
      endif

      if((.not.cyl_out).and.(.not.flx_out)) then
         if(present(phi_out)) cyl_out=.TRUE.
      endif

      intypes_in = 0
      intypes_out = 0
      if(car_in) then
         intypes_in = intypes_in + 1
         if(cyl_out) intypes_out = intypes_out + 1
         if(flx_out) intypes_out = intypes_out + 1
      endif

      if(cyl_in) then
         intypes_in = intypes_in + 1
         if(car_out) intypes_out = intypes_out + 1
         if(flx_out) intypes_out = intypes_out + 1
      endif

      if(flx_in) then
         intypes_in = intypes_in + 1
         if(cyl_out) intypes_out = intypes_out + 1
         if(car_out) intypes_out = intypes_out + 1
      endif

      if(intypes_in.gt.1) then

         ier=2501

      else if(intypes_out.eq.0) then

         ier=2502

      endif

      if(ier.gt.0) then

         if(car_in) call xplasma_errmsg_append(s, &
              '  cartesian coordinate inputs {x,y,z} detected.')

         if(cyl_in) call xplasma_errmsg_append(s, &
              '  cylindrical coordinate inputs {R,phi,Z} detected.')

         if(flx_in) call xplasma_errmsg_append(s, &
              '  flux coordinate inputs {rho,phi,theta} detected.')

         return

      endif

      !-------------------------------
      ! most error checks done.  See what transforms are needed...

      car2cyl=.FALSE.
      cyl2flx=.FALSE.

      flx2cyl=.FALSE.
      cyl2car=.FALSE.
      
      if(car_in) then
         car2cyl=.TRUE.
         if(flx_out) cyl2flx=.TRUE.
      endif

      if(flx_in) then
         flx2cyl=.TRUE.
         if(car_out) cyl2car=.TRUE.
      endif

      if(cyl_in) then
         if(car_out) cyl2car=.TRUE.
         if(flx_out) cyl2flx=.TRUE.
      endif

      id_R=0; id_Z=0
      if(cyl2flx.or.flx2cyl) then
         call xplasma_common_ids(s,ier,id_R=id_R,id_Z=id_Z)
         if(min(id_R,id_Z).eq.0) then

            ier=2504  ! no flux surfaces...
            return

         endif
      endif

      !-------------------------
      ! load working variables...

      if(car_in) then
         if(present(x_in)) then
            xv(1)=x_in
         else
            xv(1)=0
         endif
         if(present(y_in)) then
            yv(1)=y_in
         else
            yv(1)=0
         endif
         if(present(z_in)) then
            zv(1)=z_in
         else
            zv(1)=0
         endif
      endif

      if(cyl_in) then
         if(present(r_in)) then
            rv(1)=r_in
         else
            rv(1)=0
         endif
         if(present(phi_in)) then
            phiv(1)=phi_in
         else
            phiv(1)=0
         endif
         if(present(z_in)) then
            zv(1)=z_in
         else
            zv(1)=0
         endif
      endif

      if(flx_in) then
         if(present(rho_in)) then
            rhov(1)=rho_in
         else
            rhov(1)=0
         endif
         if(present(phi_in)) then
            phiv(1)=phi_in
         else
            phiv(1)=0
         endif
         if(present(theta_in)) then
            thetav(1)=theta_in
         else
            thetav(1)=0
         endif
      endif


      !-------------------------
      ! compute transformations
      
      call trans(s,id_R,id_Z,car2cyl,cyl2flx,flx2cyl,cyl2car, &
           xv,yv,zv,rv,phiv,rhov,thetav,ztol,imap_type, &
           iregion,jccw_theta,jccw_phi,j2pi, ier)

      !-------------------------
      ! copy results to output

      if(present(nregion)) nregion=iregion(1)

      if(car_out) then
         if(present(x_out)) x_out=xv(1)
         if(present(y_out)) y_out=yv(1)
         if(present(z_out)) z_out=zv(1)
      endif

      if(cyl_out) then
         if(present(r_out)) r_out=rv(1)
         if(present(phi_out)) phi_out=phiv(1)
         if(present(z_out)) z_out=zv(1)
      endif

      if(flx_out) then
         if(present(rho_out)) rho_out=rhov(1)
         if(present(phi_out)) phi_out=phiv(1)
         if(present(theta_out)) theta_out=thetav(1)
      endif

      if(present(cosphi)) cosphi = cos(phiv(1))
      if(present(sinphi)) sinphi = sin(phiv(1))

    end subroutine xplasma_ctran1

    !======================================================================
    subroutine xplasma_ctrann(s,vsize_flag,ier, &
         isize, &
         x_in,y_in,z_in, r_in,phi_in, rho_in,theta_in, &
         tol, maptype, ccw_theta,ccw_phi, i2pi, &
         x_out,y_out,z_out, r_out,phi_out, rho_out,theta_out, &
         nregion,cosphi,sinphi)

      !  VECTOR coordinate mapping interface

      type (xplasma), pointer :: s
      logical, intent(in) :: vsize_flag  ! vector size flag
      !  .TRUE. means size is to be inferred from passed vector sizes,
      !  all of which must match;
      !  .FALSE. means size is to be set by optional argument "isize";
      !  all vectors must size.ge.isize then.

      integer, intent(out) :: ier   ! completion code 0=OK

      !------------------------
      ! optional inputs...

      integer, intent(in), optional :: isize  ! vector size if(vsize_flag)

      real*8, intent(in), dimension(:), optional :: x_in,y_in,z_in  ! Cartesian Coords in
      real*8, intent(in), dimension(:), optional :: r_in,phi_in   ! Cylindric Coords (z also)
      real*8, intent(in), dimension(:), optional :: rho_in,theta_in ! Flux Coords (phi also)

      !------------------------------------------------------------
      ! the following applies to (R,Z) -> (rho,theta) "inverse" map:

      real*8, intent(in), optional :: tol           ! map tolerance
      
      ! tol applies to Newton map (maptype = xp_imap_newton = default) ONLY

      ! (R_in,Z_in)->(rho_out,theta_out): 
      !         |R_in-R(rho_out,theta_out)| < tol*max(R_in,|Z_in|) 
      !         |Z_in-Z(rho_out,theta_out)| < tol*max(R_in,|Z_in|) 
      ! default: s%bdytol; s%bdytol's own default value is 1.0d-10

      integer, intent(in), optional :: maptype

      !  = 1 = xp_imap_newton = default -- iterative Newton root finder
      !        with initial guess generator using polar map (slowest but 
      !        most accurate, accuracy can be controlled with tol.
      !  = 2 = xp_imap_polar -- polar map inversion method
      !        much faster, not as accurate as Newton map
      !  = 3 = xp_imap_rzlinear -- bilinear inverse map -- after 
      !        initialization (using polar map), the fastest, but less
      !        accurate still (depends on R,Z grid resolution).  This 
      !        method requires a user defined (R,Z) rectangle grid; if 
      !        none is available the polar map is used instead.

      !------------------------------------------------------------

      logical, intent(in), optional :: ccw_phi  ! Phi CCW view from above, default=T
      logical, intent(in), optional :: ccw_theta ! Theta CCW in poloidal plane, default=T

      integer, intent(in), optional :: i2pi ! angle coord normalization options
      !  i2pi=1 (default) returned value in range [0,2pi]
      !  i2pi=2 returned value in range [-pi,pi]

      !-------------------------
      ! optional outputs...

      real*8, intent(out), dimension(:), optional :: x_out,y_out,z_out  ! Cartesian Coords in
      real*8, intent(out), dimension(:), optional :: r_out,phi_out   ! Cyl. Coords (z also)
      real*8, intent(out), dimension(:), optional :: rho_out,theta_out ! Flux Coords (phi also)

      real*8, intent(out), dimension(:), optional :: cosphi,sinphi   ! phi sin & cos
      integer, intent(out), dimension(:), optional :: nregion        ! region code
      ! 0=not checked; 1=inside core plasma; 2=outside plasma but inside
      ! mapped region; 3=beyond mapped region.

      !-----------------
      logical :: car_in,cyl_in,flx_in
      logical :: car_out,cyl_out,flx_out
      integer :: intypes_in,intypes_out,isize_min,isize_max,j2pi,jsize
      logical :: jccw_phi,jccw_theta
      character*80 zmsg

      logical :: car2cyl,cyl2flx
      logical :: flx2cyl,cyl2car

      logical :: sol
      integer :: imap_type

      real*8 :: ztol
      real*8, dimension(:), allocatable :: xv,yv,zv,rv,phiv,rhov,thetav
      integer, dimension(:), allocatable :: iregion

      integer :: id_R,id_Z
      !-----------------

      !  clear outputs as needed

      ier = 0
      if(present(x_out)) x_out = 0
      if(present(y_out)) y_out = 0
      if(present(z_out)) z_out = 0

      if(present(r_out)) r_out = 0
      if(present(phi_out)) phi_out = 0

      if(present(rho_out)) rho_out = 0
      if(present(theta_out)) theta_out = 0

      !  set CCW flags & ztol

      jccw_phi=.TRUE.
      jccw_theta=.TRUE.
      if(present(ccw_phi)) jccw_phi=ccw_phi
      if(present(ccw_theta)) jccw_theta=ccw_theta


      imap_type=xp_imap_newton

      call xplasma_global_info(s,ier, scrapeoff=sol)
      if(present(maptype)) then
         if(sol) then
            imap_type=max(xp_imap_newton,min(xp_imap_rzlinear,maptype))
         else
            imap_type=max(xp_imap_newton,min(xp_imap_polar,maptype))
         endif
      endif

      call xplasma_global_info(s,ier, bdytol=ztol)
      if(ier.ne.0) return
      if(present(tol)) ztol=tol
      ztol=max(ceps13,min(ceps4,ztol))

      j2pi=1
      if(present(i2pi)) j2pi=i2pi

      !  check input-- must be from just one of the three possible
      !  coordinate systems.

      !  check output-- there must be some output.

      car_in = .FALSE.
      cyl_in = .FALSE.
      flx_in = .FALSE.

      car_out = .FALSE.
      cyl_out = .FALSE.
      flx_out = .FALSE.

      if(present(x_in).or.present(y_in)) car_in = .TRUE.
      if(present(R_in)) cyl_in = .TRUE.
      if(present(rho_in).or.present(theta_in)) flx_in = .TRUE.

      if((.not.car_in).and.(.not.cyl_in)) then
         if(present(z_in)) cyl_in=.TRUE.
      endif

      if((.not.cyl_in).and.(.not.flx_in)) then
         if(present(phi_in)) cyl_in=.TRUE.
      endif

      if(present(x_out).or.present(y_out)) car_out = .TRUE.
      if(present(R_out)) cyl_out = .TRUE.
      if(present(rho_out).or.present(theta_out)) flx_out = .TRUE.

      if((.not.car_out).and.(.not.cyl_out)) then
         if(present(z_out)) cyl_out=.TRUE.
      endif

      if((.not.cyl_out).and.(.not.flx_out)) then
         if(present(phi_out)) cyl_out=.TRUE.
      endif

      intypes_in = 0
      intypes_out = 0
      if(car_in) then
         intypes_in = intypes_in + 1
         if(cyl_out) intypes_out = intypes_out + 1
         if(flx_out) intypes_out = intypes_out + 1
      endif

      if(cyl_in) then
         intypes_in = intypes_in + 1
         if(car_out) intypes_out = intypes_out + 1
         if(flx_out) intypes_out = intypes_out + 1
      endif

      if(flx_in) then
         intypes_in = intypes_in + 1
         if(cyl_out) intypes_out = intypes_out + 1
         if(car_out) intypes_out = intypes_out + 1
      endif

      if(intypes_in.gt.1) then

         ier=2501

      else if(intypes_out.eq.0) then

         ier=2502

      endif

      if(ier.gt.0) then

         if(car_in) call xplasma_errmsg_append(s, &
              '  cartesian coordinate inputs {x,y,z} detected.')

         if(cyl_in) call xplasma_errmsg_append(s, &
              '  cylindrical coordinate inputs {R,phi,Z} detected.')

         if(flx_in) call xplasma_errmsg_append(s, &
              '  flux coordinate inputs {rho,phi,theta} detected.')

         return

      endif

      !  check sizes of input vectors...

      isize_min=999999999
      isize_max=0

      if(car_in) then
         if(present(x_in)) call ssize(x_in,isize_min,isize_max)
         if(present(y_in)) call ssize(y_in,isize_min,isize_max)
         if(present(z_in)) call ssize(z_in,isize_min,isize_max)
      endif

      if(cyl_in) then
         if(present(R_in)) call ssize(R_in,isize_min,isize_max)
         if(present(phi_in)) call ssize(phi_in,isize_min,isize_max)
         if(present(Z_in)) call ssize(Z_in,isize_min,isize_max)
      endif

      if(flx_in) then
         if(present(rho_in)) call ssize(rho_in,isize_min,isize_max)
         if(present(phi_in)) call ssize(phi_in,isize_min,isize_max)
         if(present(theta_in)) call ssize(theta_in,isize_min,isize_max)
      endif

      if(present(x_out)) call ssize(x_out,isize_min,isize_max)
      if(present(y_out)) call ssize(y_out,isize_min,isize_max)
      if(present(z_out)) call ssize(z_out,isize_min,isize_max)

      if(present(r_out)) call ssize(r_out,isize_min,isize_max)
      if(present(phi_out)) call ssize(phi_out,isize_min,isize_max)
      if(present(cosphi)) call ssize(cosphi,isize_min,isize_max)
      if(present(sinphi)) call ssize(sinphi,isize_min,isize_max)

      if(present(rho_out)) call ssize(rho_out,isize_min,isize_max)
      if(present(theta_out)) call ssize(theta_out,isize_min,isize_max)

      if(present(nregion)) then
         isize_min=min(isize_min,size(nregion))
         isize_max=max(isize_max,size(nregion))
      endif

      jsize=isize_min
      if(.not.vsize_flag) then
         if(present(isize)) then
            jsize=isize
            if(isize_min.lt.jsize) then

               ier=2503
               zmsg=' '
               write(zmsg,*) ' ?xplasma_ctrann: isize argument: ', &
                    jsize,'; min vector size:',isize_min
               call xplasma_errmsg_append(s,zmsg)
            endif
         else
            ier=2503
            zmsg=' '
            write(zmsg,*) ' ?xplasma_ctrann: vsize_flag=.FALSE., but ', &
                 'argument "isize" was not provided.'
            call xplasma_errmsg_append(s,zmsg)
         endif

      else if(vsize_flag) then
         if(isize_min.ne.isize_max) then
            ier=2503
            zmsg=' '
            write(zmsg,*) ' vsize_flag=.TRUE. but input vector sizes ', &
                 'do not match: min, max = ',isize_min, isize_max
            call xplasma_errmsg_append(s,zmsg)
         endif
      endif

      if(ier.ne.0) return

      !-------------------------------
      ! most error checks done.  See what transforms are needed...

      car2cyl=.FALSE.
      cyl2flx=.FALSE.

      flx2cyl=.FALSE.
      cyl2car=.FALSE.
      
      if(car_in) then
         car2cyl=.TRUE.
         if(flx_out) cyl2flx=.TRUE.
      endif

      if(flx_in) then
         flx2cyl=.TRUE.
         if(car_out) cyl2car=.TRUE.
      endif

      if(cyl_in) then
         if(car_out) cyl2car=.TRUE.
         if(flx_out) cyl2flx=.TRUE.
      endif

      id_R=0; id_Z=0
      if(cyl2flx.or.flx2cyl) then
         call xplasma_common_ids(s,ier,id_R=id_R,id_Z=id_Z)
         if(min(id_R,id_Z).eq.0) then

            ier=2504  ! no flux surfaces...
            return

         endif
      endif

      !-------------------------
      ! load working variables...

      allocate(iregion(jsize),xv(jsize),yv(jsize),zv(jsize), &
           rv(jsize),phiv(jsize),rhov(jsize),thetav(jsize))

      if(car_in) then
         if(present(x_in)) then
            xv=x_in(1:jsize)
         else
            xv=0
         endif
         if(present(y_in)) then
            yv=y_in(1:jsize)
         else
            yv=0
         endif
         if(present(z_in)) then
            zv=z_in(1:jsize)
         else
            zv=0
         endif
      endif

      if(cyl_in) then
         if(present(r_in)) then
            rv=r_in(1:jsize)
         else
            rv=0
         endif
         if(present(phi_in)) then
            phiv=phi_in(1:jsize)
         else
            phiv=0
         endif
         if(present(z_in)) then
            zv=z_in(1:jsize)
         else
            zv=0
         endif
      endif

      if(flx_in) then
         if(present(rho_in)) then
            rhov=rho_in(1:jsize)
         else
            rhov=0
         endif
         if(present(phi_in)) then
            phiv=phi_in(1:jsize)
         else
            phiv=0
         endif
         if(present(theta_in)) then
            thetav=theta_in(1:jsize)
         else
            thetav=0
         endif
      endif


      !-------------------------
      ! compute transformations
      
      call trans(s,id_R,id_Z,car2cyl,cyl2flx,flx2cyl,cyl2car, &
           xv,yv,zv,rv,phiv,rhov,thetav,ztol,imap_type, &
           iregion,jccw_theta,jccw_phi,j2pi, ier)

      !-------------------------
      ! copy results to output

      if(present(nregion)) nregion(1:jsize)=iregion

      if(car_out) then
         if(present(x_out)) x_out(1:jsize)=xv
         if(present(y_out)) y_out(1:jsize)=yv
         if(present(z_out)) z_out(1:jsize)=zv
      endif

      if(cyl_out) then
         if(present(r_out)) r_out(1:jsize)=rv
         if(present(phi_out)) phi_out(1:jsize)=phiv
         if(present(z_out)) z_out(1:jsize)=zv
      endif


      if(flx_out) then
         if(present(rho_out)) rho_out(1:jsize)=rhov
         if(present(phi_out)) phi_out(1:jsize)=phiv
         if(present(theta_out)) theta_out(1:jsize)=thetav
      endif

      if(present(cosphi)) cosphi(1:jsize) = cos(phiv)
      if(present(sinphi)) sinphi(1:jsize) = sin(phiv)

      contains
        subroutine ssize(r8vec,isize_min,isize_max)
          real*8, dimension(:), intent(in) :: r8vec
          integer, intent(inout) :: isize_min,isize_max
          integer :: jsize

          jsize=size(r8vec)
          isize_min=min(isize_min,jsize)
          isize_max=max(isize_max,jsize)

        end subroutine ssize

    end subroutine xplasma_ctrann

    !--------------------------------------------------------------------
    !--------------------------------------------------------------------

    subroutine trans(s,id_R,id_Z,car2cyl,cyl2flx,flx2cyl,cyl2car, &
         xv,yv,zv,rv,phiv,rhov,thetav,tol,maptype, &
         iregion,jccw_theta,jccw_phi,j2pi,ier)

      ! **private**
      ! coord transform routine, almost no error checking...

      type(xplasma), pointer :: s
      integer, intent(in) :: id_R,id_Z  ! R(rho,theta) & Z(rho,theta)

      ! flags selecting which transforms to evaluate:
      logical, intent(in) :: car2cyl    
      logical, intent(in) :: cyl2flx
      logical, intent(in) :: flx2cyl
      logical, intent(in) :: cyl2car

      ! coord transform working variables -- mix of input and output depending
      ! on switch settings; all the same size...

      real*8, dimension(:) :: xv,yv,zv,rv,phiv,rhov,thetav

      ! inverse map tolerance

      real*8, intent(in) :: tol
      integer, intent(in) :: maptype  ! type of (R,Z)->(rho,theta) inv. map.

      ! location flag

      integer, intent(out), dimension(:) :: iregion

      ! CCW transforms for angle variables

      logical, intent(in) :: jccw_theta,jccw_phi

      ! range options for angle variables

      integer, intent(in) :: j2pi

      ! completion code -- only set on failure of convergence of inverse
      ! map (generally indicates invalid flux surface geometry)

      integer, intent(out) :: ier

      !--------------------------------------------------
      real*8, dimension(:,:),allocatable :: wkv
      integer :: ii,ids(2)
      !--------------------------------------------------
      
      ier=0
      iregion=0

      if(car2cyl) then

         ! Cartesian -> Cylindrical

         rv = sqrt(xv**2 + yv**2)
         where(rv.eq.czero)
            phiv=czero
         end where
         where(rv.gt.czero)
            phiv = atan2(yv,xv)
         end where

         call angle_trans(phiv,xplasma_phi_coord,jccw_phi,j2pi)

      endif

      if(cyl2flx) then

         if(.not.car2cyl) call angle_trans( &
              phiv,xplasma_phi_coord,jccw_phi,j2pi)

         call cyl2flx_exec(s,rv,phiv,zv,rhov,thetav,tol,maptype,iregion,ier)
         if(ier.ne.0) return

         call angle_trans(thetav,xplasma_theta_coord,jccw_theta,j2pi)
      endif

      if(flx2cyl) then

         call angle_trans(thetav,xplasma_theta_coord,jccw_theta,j2pi)
         call angle_trans(phiv,xplasma_phi_coord,jccw_phi,j2pi)

         ids(1)=id_R
         ids(2)=id_Z

         allocate(wkv(size(rv),2))

         call xplasma_RZeval_2d(s,ids, &
              xplasma_theta_coord,thetav,xplasma_rho_coord,rhov,wkv,ier, &
              force_bounds=.TRUE.)

         rv=wkv(1:size(rv),1)
         zv=wkv(1:size(rv),2)

         deallocate(wkv)

         do ii=1,size(rhov)
            if(rhov(ii).lt.czero) then
               iregion(ii)=5
            else if(rhov(ii).gt.cone) then
               iregion(ii)=4  ! ** refine this **
            else
               iregion(ii)=1
            endif
         enddo

      endif

      if(cyl2car) then
         if(.not.flx2cyl) call angle_trans( &
              phiv,xplasma_phi_coord,jccw_phi,j2pi)

         xv = rv*cos(phiv)
         yv = rv*sin(phiv)

      endif

    end subroutine trans

    subroutine angle_trans(zangl,jcoord,jccw,j2pi)

      type (xplasma), pointer :: s
      real*8, dimension(:) :: zangl  ! angle with range to be normalized...
      integer :: jcoord              ! coordinate id (not currently used)
      logical :: jccw                ! CCW option
      integer :: j2pi                ! range option

      !-------------------
      real*8 :: amin,amax
      integer :: ii,icyc,icheck
      !-------------------

      if((j2pi.le.1).or.(j2pi.gt.2)) then
         amin=czero
         amax=c2pi

      else if(j2pi.eq.2) then
         amin=-cpi
         amax=+cpi

      endif

      do ii=1,size(zangl)

         if(zangl(ii).lt.amin) then
            icyc=(amax-zangl(ii))/c2pi
            zangl(ii)=zangl(ii)+icyc*c2pi
            icheck=1

         else if(zangl(ii).gt.amax) then
            icyc=(zangl(ii)-amin)/c2pi
            zangl(ii)=zangl(ii)-icyc*c2pi
            icheck=1

         else
            icheck=0  ! in range now
         endif

         if(icheck.eq.1) then
            if(zangl(ii).lt.amin) then
               do
                  zangl(ii)=zangl(ii)+c2pi
                  if(zangl(ii).ge.amin) exit
               enddo

            else if(zangl(ii).gt.amax) then
               do 
                  zangl(ii)=zangl(ii)-c2pi
                  if(zangl(ii).le.amax) exit
               enddo
            endif
         endif

         ! now in range...

         if(.not.jccw) then
            zangl(ii)=amax-(zangl(ii)-amin)
         endif

      enddo

    end subroutine angle_trans

    subroutine cyl2flx_exec(s,rv,phiv,zv,rhov,thetav,tol,maptype,iregion,ier)

      use eqi_rzbox_module

      ! *private*
      ! generate initial guess of flux coords (rhov,thetav) corresponding to
      ! real space locations (rv,zv).

      ! then...
      ! refine solution of flux coords (rhov,thetav) corresponding to
      ! real space locations (rv,zv), such that the R & Z of the flux
      ! coordinates returned are within tol*[R of mag. axis] of the
      ! original input values (where possible -- i.e. for all elements
      ! of (rv,zv) in the mapped region.

      ! array dimensions have been set but individual elements of (rv,zv)
      ! could be out of range.

      type(xplasma), pointer :: s

      real*8, dimension(:) :: rv,phiv,zv
      real*8, dimension(:) :: rhov,thetav
      real*8, intent(in) :: tol      ! relative accuracy tolerance
      integer, intent(in) :: maptype ! type of (R,Z)->(rho,theta) map
      integer, dimension(:) :: iregion
      integer, intent(out) :: ier

      !------------------------
      real*8, dimension(:,:), allocatable :: r_orig,z_orig
      real*8 :: rtol,reltol,raxis,zaxis
      real*8, parameter :: czero = 0.0d0
      !------------------------

      !  Method of generating initial guess depends on whether there is
      !  a scrape off layer (SOL) or not...

      sp => s

      if(maptype.eq.xp_imap_rzlinear) then

         call eqi_fastinv(size(rv),rv,zv,phiv,czero,rhov,thetav,iregion,ier)

      else

         iregion=0
         call xplasma_rhoth_inv(s,rv,zv,rhov,thetav,iregion,ier)

      endif
      if(maptype.ne.xp_imap_newton) return

      call xplasma_mag_axis(s,ier, raxis, zaxis)
      if(ier.ne.0) return

      reltol = max(ceps13,min(ceps4,tol))
      rtol = reltol*raxis

      call eqi_inv(size(rv),rv,zv,phiv,reltol,rtol, &
           rhov,thetav,iregion,ier)

    end subroutine cyl2flx_exec

    subroutine xplasma_RZminmax_extended(s, rmin,rmax, zmin,zmax, ier, &
         sol, id_Rgrid, id_Zgrid)

      !  return the R & Z min & max of the extended mapped region; if
      !  there is no such region return the R & Z min & max of the
      !  plasma

      type(xplasma), pointer :: s

      real*8, intent(out) :: Rmin,Rmax
      real*8, intent(out) :: Zmin,Zmax

      integer, intent(out) :: ier

      !----------------------

      logical, intent(out), optional :: sol  ! T if scrape off layer defined
      integer, intent(out), optional :: id_Rgrid,id_Zgrid ! SOL grid ids

      !-------------------------------------
      logical :: isol
      integer :: jd_Rgrid,jd_Zgrid
      real*8, parameter :: zrho_bdy = 1.0d0
      !-------------------------------------

      rmin=0
      rmax=0
      zmin=0
      zmax=0
      if(present(id_Rgrid)) id_Rgrid = 0
      if(present(id_Zgrid)) id_Zgrid = 0

      call xplasma_global_info(s,ier, scrapeoff=isol)
      if(present(sol)) sol = isol
      if(ier.ne.0) return

      if(.not.isol) then

         call xplasma_RZminmax(s,zrho_bdy,ier, &
              rmin=Rmin,rmax=Rmax, zmin=Zmin,zmax=Zmax)

      else

         call xplasma_find_item(s,'__Rgrid',jd_Rgrid,ier, nf_noerr=.TRUE.)
         if(ier.ne.0) return
         if(present(id_Rgrid)) id_Rgrid=jd_Rgrid

         call xplasma_find_item(s,'__Zgrid',jd_Zgrid,ier, nf_noerr=.TRUE.)
         if(ier.ne.0) return
         if(present(id_Zgrid)) id_Zgrid=jd_Rgrid

         call xplasma_grid_info(s,jd_Rgrid,ier, xmin=rmin,xmax=rmax)
         call xplasma_grid_info(s,jd_Zgrid,ier, xmin=zmin,xmax=zmax)
      endif

    end subroutine xplasma_RZminmax_extended

    subroutine xplasma_RZminmax_plasma(s, rmin,rmax, zmin,zmax, ier, &
         ccw_theta,i2pi, thRmin,thRmax, thZmin,thZmax)

      !  return the plasma boundary R & Z min & max

      type(xplasma), pointer :: s

      real*8, intent(out) :: Rmin,Rmax
      real*8, intent(out) :: Zmin,Zmax

      integer, intent(out) :: ier

      !  optionally return the poloidal angle locations

      logical, intent(in), optional :: ccw_theta ! Theta CCW in poloidal plane, default=T

      integer, intent(in), optional :: i2pi ! angle coord normalization options
      !  i2pi=1 (default) returned value in range [0,2pi]
      !  i2pi=2 returned value in range [-pi,pi]

      real*8, intent(out), optional :: thRmin,thRmax
      real*8, intent(out), optional :: thZmin,thZmax

      !-------------------------------------
      real*8, parameter :: zrho_bdy = 1.0d0
      !-------------------------------------

      call xplasma_RZminmax(s,zrho_bdy,ier, &
           rmin=Rmin,rmax=Rmax, zmin=Zmin,zmax=Zmax, &
           ccw_theta=ccw_theta,i2pi=i2pi, &
           thRmin=thRmin,thRmax=thRmax, thZmin=thZmin,thZmax=thZmax)

    end subroutine xplasma_RZminmax_plasma

    subroutine xplasma_RZminmax(s,rho,ier, phi, &
         ccw_theta,i2pi, &
         rmin,rmax, zmin,zmax, &
         thRmin,thRmax, thZmin,thZmax)

      use eqi_rzbox_module

      !  very accurately determine Rmin/max, Zmin/max of a flux surface.

      type (xplasma), pointer :: s
      real*8, intent(in) :: rho       ! required input: flux surface.
      integer, intent(out) :: ier     ! completion code

      real*8, intent(in), optional :: phi   ! toroidal angle location

      real*8, intent(out), optional :: rmin,rmax  ! Rmin & Rmax of surface
      real*8, intent(out), optional :: zmin,zmax  ! Zmin & Zmax of surface

      !  optionally return the poloidal angle locations

      logical, intent(in), optional :: ccw_theta ! Theta CCW in poloidal plane, default=T

      integer, intent(in), optional :: i2pi ! angle coord normalization options
      !  i2pi=1 (default) returned value in range [0,2pi]
      !  i2pi=2 returned value in range [-pi,pi]

      real*8, intent(out), optional :: thRmin,thRmax
      real*8, intent(out), optional :: thZmin,thZmax

      !-------------------
      real*8 :: zrho,zphi,elemvals(8),zrmin,zrmax,zzmin,zzmax,bdytol
      real*8 :: zthRmin,zthRmax,zthZmin,zthZmax
      character*32 elemNames(8),elemLabels(8)
      logical :: ibdy,iwant_R,iwant_Z,jccw
      integer :: id_rzmmx,iertmp,j2pi,inum_elems
      !-------------------

      ier = 0

      if(present(rmin)) rmin=0
      if(present(rmax)) rmax=0

      if(present(thRmin)) thRmin=0
      if(present(thRmax)) thRmax=0

      if(present(zmin)) zmin=0
      if(present(zmax)) zmax=0

      if(present(thZmin)) thZmin=0
      if(present(thZmax)) thZmax=0

      j2pi=1
      if(present(i2pi)) j2pi=i2pi

      jccw=.TRUE.
      if(present(ccw_theta)) jccw=ccw_theta

      !-----------
      !  non-axisymmetry upgrade needed here someday...
      zphi = 0
      if(present(phi)) zphi=phi
      !-----------
      
      call xplasma_global_info(s,ier, bdytol=bdytol)
      if(abs(rho-1.0d0).lt.bdytol) then
         ibdy=.TRUE.
         iwant_R=.TRUE.
         iwant_Z=.TRUE.
         zrho=1.0d0

         call xplasma_common_ids(s,ier, id_RZminmax=id_rzmmx)
         if(ier.ne.0) return

         if(id_rzmmx.gt.0) then
            call xplasma_getList_size(s,id_rzmmx,inum_elems,ier)
            if(inum_elems.ne.8) id_rzmmx=0  ! recompute...
         endif

         if(id_rzmmx.gt.0) then

            ! answer has been saved already: retrieve it.
            call xplasma_getList(s,id_rzmmx,elemvals,ier)

            if(present(rmin)) rmin = elemvals(1)
            if(present(rmax)) rmax = elemvals(2)
            if(present(zmin)) zmin = elemvals(3)
            if(present(zmax)) zmax = elemvals(4)

            if(present(thRmin)) thRmin = set_th(elemvals(5))
            if(present(thRmax)) thRmax = set_th(elemvals(6))
            if(present(thZmin)) thZmin = set_th(elemvals(7))
            if(present(thZmax)) thZmax = set_th(elemvals(8))

            return
         endif
      else
         ibdy=.FALSE.
         iwant_R = present(rmin).or.present(rmax)
         iwant_R = iwant_R.or.present(thRmin).or.present(thRmax)
         iwant_Z = present(zmin).or.present(zmax)
         iwant_Z = iwant_Z.or.present(thZmin).or.present(thZmax)
         zrho=rho
      endif

      if(.not.(iwant_R.or.iwant_Z)) return

      !  compute the extrema

      sp => s   ! module pointer used to pass "s" to root finder argument fcn

      call eqi_rzbox(zrho,zphi, &
           iwant_R,zrmin,zthrmin,iwant_R,zrmax,zthrmax, &
           iwant_Z,zzmin,zthzmin,iwant_Z,zzmax,zthzmax, &
           ier)
      if(ier.ne.0) return

      !  store the result as a list if flag is set

      if(ibdy) then

         elemnames(1)='RMIN'
         elemnames(2)='RMAX'
         elemnames(3)='ZMIN'
         elemnames(4)='ZMAX'

         elemnames(5)='THRMIN'
         elemnames(6)='THRMAX'
         elemnames(7)='THZMIN'
         elemnames(8)='THZMAX'

         elemLabels(1)='minimum R of bdy flux surface'
         elemLabels(2)='maximum R of bdy flux surface'
         elemLabels(3)='minimum Z of bdy flux surface'
         elemLabels(4)='maximum Z of bdy flux surface'

         elemLabels(5)='minimum R of bdy flux surface'
         elemLabels(6)='maximum R of bdy flux surface'
         elemLabels(7)='minimum Z of bdy flux surface'
         elemLabels(8)='maximum Z of bdy flux surface'

         elemvals(1)=zrmin
         elemvals(2)=zrmax
         elemvals(3)=zzmin
         elemvals(4)=zzmax

         elemvals(5)=zthrmin
         elemvals(6)=zthrmax
         elemvals(7)=zthzmin
         elemvals(8)=zthzmax

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         call xplasma_create_list(s,"BDY_RZ_MINMAX", &
              elemnames,id_rzmmx,ier, &
              label='R & Z extrema of bdy flux surface', units='m', &
              r8vals=elemvals, chvals=elemLabels)

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)
      endif

      !  return the requested
      if(present(rmin)) rmin=zrmin
      if(present(rmax)) rmax=zrmax

      if(present(zmin)) zmin=zzmin
      if(present(zmax)) zmax=zzmax

      if(present(thRmin)) thRmin = set_th(zthrmin)
      if(present(thRmax)) thRmax = set_th(zthrmax)

      if(present(thZmin)) thZmin = set_th(zthzmin)
      if(present(thZmax)) thZmax = set_th(zthzmax)

      contains
        real*8 function set_th(zth_in)
          real*8, intent(in) :: zth_in

          real*8 :: zth(1)

          zth=zth_in
          call angle_trans(zth,xplasma_theta_coord,jccw,j2pi)
          set_th=zth(1)

        end function set_th

    end subroutine xplasma_RZminmax

    subroutine xplasma_Bminmax(s,rho,ier, phi, &
         ccw_theta,i2pi, bmin,bmax, thbmin,thbmax)

      use eqi_rzbox_module

      !  very accurately determine Bmin/max of a flux surface

      type (xplasma), pointer :: s
      real*8, intent(in) :: rho       ! required input: flux surface.
      integer, intent(out) :: ier     ! completion code

      real*8, intent(in), optional :: phi   ! toroidal angle location

      logical, intent(in), optional :: ccw_theta ! Theta CCW in poloidal plane, default=T

      integer, intent(in), optional :: i2pi ! angle coord normalization options
      !  i2pi=1 (default) returned value in range [0,2pi]
      !  i2pi=2 returned value in range [-pi,pi]

      real*8, intent(out), optional :: Bmin,Bmax  ! Bmin & Bmax of surface
      real*8, intent(out), optional :: thBmin,thBmax

      !------------------------
      integer :: id_Bmod,j2pi
      real*8 :: zphi,zbmin,zbmax,zthbmin,zthbmax,zth(1)
      logical :: ccwflag
      !------------------------

      call xplasma_common_ids(s,ier, id_Bmod=id_Bmod)
      if(ier.ne.0) return

      if(id_Bmod.eq.0) then
         ier=106
         call xplasma_errmsg_append(sp, &
              '  xplasma_Bminmax: Bmin & Bmax unavailable: no Bmod profile.')
         return
      endif

      if(present(Bmin)) Bmin=0
      if(present(Bmax)) Bmax=0

      if(present(thBmin)) thBmin=0
      if(present(thBmax)) thBmax=0

      zphi=0.0d0
      if(present(phi)) zphi=phi  ! not used now (axisymmetry assumed).

      ccwflag=.TRUE.
      if(present(ccw_theta)) ccwflag=ccw_theta

      j2pi=1
      if(present(i2pi)) j2pi=max(1,min(2,i2pi))

      sp => s

      call eqi_bbox(rho, zbmin,zthbmin, zbmax,zthbmax, ier)

      if(ier.ne.0) return

      if(present(Bmin)) Bmin=zbmin
      if(present(Bmax)) Bmax=zbmax

      if(present(thBmin)) then
         zth=zthbmin
         call angle_trans(zth,xplasma_theta_coord,ccwflag,j2pi)
         thBmin=zth(1)
      endif

      if(present(thBmax)) then
         zth=zthbmax
         call angle_trans(zth,xplasma_theta_coord,ccwflag,j2pi)
         thBmax=zth(1)
      endif

    end subroutine xplasma_Bminmax

    subroutine xplasma_vtaub(s,rho,xb,vtaub,ierr, &
         tol,bcut,indx,bmin,bmax)

      !  evaluate zero banana width approximation to v*tau_bounce

      !  this is the integral over the limits of the orbit of
      !     dtheta*(dlp/dtheta)*(B/Bp)*(1/sqrt(1-B/Brefl))

      !  where Brefl is a trapped orbit's reflection point (vpll=0,v=vperp)
      !  for passing orbits (Brefl > Bmax) the integrand is over the
      !  entire surface contour; for trapped orbits the integrand is
      !  twice the value from theta_lower to theta_upper where these
      !  limits need to be found by a root finder so that
      !     B(theta_lower)=B(theta_upper)=Brefl

      !  the integrand formula follows from the assumptions that v
      !  and mu=vperp**2/B are constants of motion
 
      !  axisymmetry is assumed...

      use eqi_rzbox_module
 
      !-----------
      type (xplasma), pointer :: s

      real*8, intent(in) :: rho  ! surface on which to calculate
      !  note -- must not be magnetic axis; non-singular surface required.

      real*8, dimension(:), intent(in) :: xb  ! evaluation points:
      !   xb(j) = (Brefl(j)-Bmin)/(Bmax-Bmin)
      ! where Bmin and Bmax, the min and mod(B) on the zrho surface, will
      ! be found; xb > 0 is required AND
      !    xb(j+1) > xb(j) required for all j < size(xb).

      real*8, dimension(:), intent(out) :: vtaub ! integral values (m).
      !  (divide by ion velocity (m/sec) to get bounce time, seconds).

      !  size of xb and vtaub must match.

      integer, intent(out) :: ierr   ! completion status, 0=OK

      real*8, intent(in), optional :: tol  ! rel. integration error tolerance
      ! default & recommended value: 1.0d-4; permitted range 1.0d-6 to 1.0d-3

      real*8, intent(in), optional :: bcut ! minimum allowed value for 
      ! approaching singularity at (1-B/Brefl); 
      ! default & recommended value: 1.0d-8; permitted range 1.0d-10 to 1.0d-6.

      integer, intent(in), optional :: indx ! evaluation element index
      ! thread index (for parallelization over surfaces)

      real*8, intent(out), optional :: bmin ! minimum mod(B) on rho surface (T)
      real*8, intent(out), optional :: bmax ! maximum mod(B) on rho surface (T)

      !--------------------------------------
      integer :: i,inum,itrap,indxe
      real*8 :: ztol,zcut,zbmin,zbmax
      real*8, parameter :: ZERO=0.0d0, ONE=1.0d0
      real*8, parameter :: eps3=1.0d-3,eps4=1.0d-4,eps6=1.0d-6
      real*8, parameter :: eps8=1.0d-8,eps10=1.0d-10,eps12=1.0d-12
      character*120 zmsg
      !--------------------------------------

      ztol=eps4
      if(present(tol)) ztol=max(eps6,min(eps3,tol))

      zcut=eps8
      if(present(bcut)) zcut=max(eps12,min(eps6,bcut))

      ierr=0

      if((rho.le.ZERO).or.(rho.gt.ONE)) then
         zmsg=' '
         write(zmsg,*) ' ?xplasma_vtaub: rho out of range 0 < rho <= 1: ',rho
         call xplasma_errmsg_append(s,zmsg)
         ierr=9999
      endif

      inum=size(xb)
      if(inum.ne.size(vtaub)) then
         ierr=9999
         zmsg=' '
         write(zmsg,*) ' ?xplasma_vtaub: size(xb)<>size(vtaub): ', &
              size(xb),' ',size(vtaub)
         call xplasma_errmsg_append(s,zmsg)
      endif
 
      itrap=0
      do i=1,inum
         if(xb(i).le.ZERO) then
            zmsg=' '
            write(zmsg,*) ' ?xplasma_vtaub: "xb" not positive: ', &
                 'xb(',i,')=',xb(i)
            call xplasma_errmsg_append(s,zmsg)
            ierr=9999
         endif
         if(i.lt.inum) then
            if(xb(i).ge.xb(i+1)) then
               zmsg=' '
               write(zmsg,*) ' ?xplasma_vtaub: "xb" not strict ascending: ', &
                    'xb(',i,')=',xb(i),'  xb(',i+1,')=',xb(i+1)
               ierr=9999
            endif
         endif
         if(ierr.ne.0) exit
         if(xb(i).le.ONE) itrap=i
      enddo
      if(ierr.ne.0) return

      ! ----- OK

      indxe=1
      if(present(indx)) indxe=indx

      sp => s
      call eqi_tau_bounce(rho,indxe, &
           itrap,inum,xb,ztol,zcut,zbmin,zbmax,vtaub,ierr)

      if(present(bmin)) bmin=zbmin
      if(present(bmax)) bmax=zbmax

    end subroutine xplasma_vtaub

    subroutine xplasma_bdfind1(s,zR,zZ,ier, &
         maptype, phi_in, rho_in, ccw_theta, ccw_phi, i2pi, outside_only, &
         theta_out, phi_out, dist)

      ! scalar xplasma_bdfind
      ! see xplasma_bdfindn for arg descriptions

      type (xplasma), pointer :: s

      real*8, intent(in) :: zR,zZ
      integer, intent(out) :: ier

      integer, intent(in), optional :: maptype
      real*8, intent(in), optional :: phi_in,rho_in

      logical, intent(in), optional :: ccw_theta,ccw_phi
      integer, intent(in), optional :: i2pi

      logical, intent(in), optional :: outside_only ! flag that caller only
      !  expects to need data from beyond plasma boundary (default .FALSE.)
      !  (this switch only useful for optimization of maptype=3 calls).

      real*8, intent(out), optional :: theta_out,phi_out,dist

      !-----------------
      real*8 :: zR1(1),zZ1(1),phi_in1(1),theta_out1(1),phi_out1(1),dist1(1)
      !-----------------

      zR1 = zR
      zZ1 = zZ

      if(present(phi_in)) then
         phi_in1 = phi_in
      else
         phi_in1 = 0
      endif

      call xplasma_bdfindn(s, 1, zR1, zZ1, ier, &
           maptype=maptype, phi_in=phi_in1, rho_in=rho_in, &
           ccw_theta=ccw_theta, ccw_phi=ccw_phi, i2pi=i2pi, &
           outside_only=outside_only, &
           theta_out=theta_out1, phi_out=phi_out1, dist=dist1)

      if(present(theta_out)) then
         theta_out = theta_out1(1)
      endif

      if(present(phi_out)) then
         phi_out = phi_out1(1)
      endif

      if(present(dist)) then
         dist = dist1(1)
      endif

    end subroutine xplasma_bdfind1

    subroutine xplasma_bdfindn(s,isize,zR,zZ,ier, &
         maptype, phi_in, rho_in, ccw_theta, ccw_phi, i2pi, outside_only, &
         theta_out, phi_out, dist)

      ! accurately find the distance of points (R,Z) from a flux surface
      ! and the closest point on the flux surface

      ! although phi arguments are supplied to support a possible future
      ! upgrade to 3d, these have no effect on axisymmetric calculations.

      use eqi_rzbox_module

      type (xplasma), pointer :: s
      integer, intent(in) :: isize   ! vector size
      !  all array arguments' sizes must be .ge.isize

      !  the target points whose distances from a surface are to be
      !  calculated:

      real*8, dimension(:), intent(in) :: zR   ! R values
      real*8, dimension(:), intent(in) :: zZ   ! Z values

      integer, intent(out) :: ier       ! completion code 0=OK

      !------------
      !  optional arguments

      integer, intent(in), optional :: maptype

      !  = 1 = xp_imap_newton = default -- iterative Newton root finder
      !        with initial guess by searching method, accurate w.r.t.
      !        spline fit of R & Z, to machine precision, but costly...

      !  = 2 = xp_imap_polar -- polar "circle-fit" distance estimate:
      !        The boundary is treated as a piecewise linear interpolation
      !        between circles fit to nearby points on the boundary
      !        much faster, not as accurate as Newton map

      !  = 3 = xp_imap_rzlinear -- bilinear inverse map -- after 
      !        initialization (using polar map), the fastest, but less
      !        accurate still (depends on R,Z grid resolution).  This 
      !        method requires a user defined (R,Z) rectangle grid; if 
      !        none is available the polar map is used instead.

      real*8, dimension(:), intent(in), optional :: phi_in  ! toroidal angle 
      !  (not needed for axisymmetric equilibria)

      real*8, intent(in), optional :: rho_in  ! **SCALAR** -- the surface from
      !  which to calculate distances.  DEFAULT is the boundary surface rho=1.0

      logical, intent(in), optional :: ccw_phi  ! Phi CCW view from above, default=T
      logical, intent(in), optional :: ccw_theta ! Theta CCW in poloidal plane, default=T

      integer, intent(in), optional :: i2pi ! angle coord normalization options
      !  i2pi=1 (default) returned value in range [0,2pi]
      !  i2pi=2 returned value in range [-pi,pi]

      logical, intent(in), optional :: outside_only ! flag that caller only
      !  expects to need data from beyond plasma boundary (default .FALSE.)
      !  (this switch only useful for optimization of maptype=3 calls).

      !-------------
      !  results-- optional but at least one should be present

      real*8, dimension(:), intent(out), optional :: theta_out  ! theta of
      !  nearest point on surface

      real*8, dimension(:), intent(out), optional :: phi_out    ! phi of
      !  nearest point on surface

      real*8, dimension(:), intent(out), optional :: dist       ! distance to
      !  surface.  If dist(j) > 0.0 then (zR(j),zZ(j)) is outside the test
      !  surface; if dist(j) < 0.0 then the point is inside the test surface;
      !  if dist(j) = 0.0 then the point is on the test flux surface.

      !-------------------------
      logical :: axisymm
      integer :: isize_min,j2pi,imap_type,id_Rgrid,iwant
      character*128 zmsg
      real*8, dimension(:), allocatable :: zphi,zth_out,zph_out,zdist
      real*8 :: zrho_in,zrho_bdy,ztol
      logical :: jccw_phi,jccw_theta,iskipi,iskipo,iexact
      !-------------------------

      ier=0
      if(present(theta_out)) theta_out=0
      if(present(phi_out)) phi_out=0
      if(present(dist)) dist=0

      call xplasma_global_info(s,ier, axisymm=axisymm, bdytol=ztol)
      if(ier.ne.0) return

      if(.not.axisymm) then
         ier=2509
         call xplasma_errmsg_append(s, &
              ' xplasma_bdfind requires an axisymmetric equilibrium.')
         return
      endif

      ! check vector sizes

      isize_min = isize

      isize_min = min(isize_min,size(zR))
      isize_min = min(isize_min,size(zZ))
      if(present(phi_in)) isize_min = min(isize_min,size(phi_in))

      if(present(theta_out)) isize_min = min(isize_min,size(theta_out))
      if(present(phi_out)) isize_min = min(isize_min,size(phi_out))
      if(present(dist)) isize_min = min(isize_min,size(dist))

      if(isize_min.lt.isize) then

         ier=2503
         call xplasma_errmsg_append(s,' ...error in xplasma_bdfind:')
         zmsg=' '
         write(zmsg,*) ' isize argument: ',isize,'; min vector size:',isize_min
         call xplasma_errmsg_append(s,zmsg)
         return

      endif

      if(present(rho_in)) then
         zrho_in=rho_in
      else
         zrho_in=1.0d0
      endif

      sp => s
      call eqi_extrap_rho_bdy(zrho_bdy)
      if(zrho_bdy.lt.0.0d0) zrho_bdy=1.0d0

      if((zrho_in.lt.-ztol).or.zrho_in.gt.(zrho_bdy+ztol)) then
         ier=306
         zmsg=' '
         write(zmsg,*) ' xplasma_bdfind: input argument rho_in = ',rho_in, &
              ' not in expected range, typically [0,1].'
         return

      else

         zrho_in=max(ztol,min(zrho_bdy,zrho_in))

      endif

      !  OK...

      j2pi=1
      if(present(i2pi)) j2pi=i2pi

      !  set CCW flags & tolerance

      jccw_phi=.TRUE.
      jccw_theta=.TRUE.
      if(present(ccw_phi)) jccw_phi=ccw_phi
      if(present(ccw_theta)) jccw_theta=ccw_theta

      imap_type = 1
      if(present(maptype)) then
         imap_type=max(1,min(3,maptype))
      endif

      iexact = imap_type.eq.1

      if(imap_type.eq.3) then
         ! fast bilinear map for rho_in=1 only
         if(abs(zrho_in-1.0d0).gt.ztol) then
            imap_type=2
         else
            ! R & Z grids have to have been defined...
            call xplasma_find_item(s,'__ZGRID',id_rgrid,ier,nf_noerr=.TRUE.)
            if(id_rgrid.eq.0) imap_type=2
            call xplasma_find_item(s,'__RGRID',id_rgrid,ier,nf_noerr=.TRUE.)
            if(id_rgrid.eq.0) imap_type=2
         endif
      endif

      !  allocate vectors for eqi_bdfind...

      allocate(zphi(isize),zth_out(isize),zph_out(isize),zdist(isize))

      if(present(phi_in)) then
         zphi=phi_in(1:isize)
      else
         zphi=0
      endif

      if(imap_type.le.2) then

         iskipi=.FALSE.
         iskipo=.FALSE.
         call eqi_bdfind(isize,zR,zZ,zphi,zrho_in,iskipi,iskipo,iexact, &
              zth_out,zph_out,zdist,ier)

      else

         iwant=2
         if(present(outside_only)) then
            if(outside_only) iwant=1
         endif
         zph_out=zphi
         call eqi_bdfast(isize,zR,zZ,iwant,present(theta_out),zdist,zth_out, &
              ier)

      endif

      if(present(theta_out)) then
         call angle_trans(zth_out,xplasma_theta_coord,jccw_theta,j2pi)
         theta_out(1:isize)=zth_out
      endif
      if((.not.axisymm).and.present(phi_out)) then
         call angle_trans(zph_out,xplasma_phi_coord,jccw_phi,j2pi)
         phi_out(1:isize)=zph_out
      endif

      if(present(dist)) dist(1:isize)=zdist

      deallocate(zphi,zth_out,zph_out,zdist)

    end subroutine xplasma_bdfindn

    subroutine xplasma_rzjac1(s,rho_in,theta_in,ier, ccwflag1, &
         r1,z1,drdrho1,dzdrho1,drdtheta1,dzdtheta1, &
         rzdetj1,drhodr1,drhodz1,dthdr1,dthdz1)

      ! evaluate terms associated with 2d Jacobean -- scalar version
      !   [dR/drho dR/dtheta]
      !   [dZ/drho dZ/dtheta]

      type (xplasma), pointer :: s
      real*8 :: rho_in   ! radial coordinate in
      real*8 :: theta_in ! poloidal angle coordinate in
      !  (see ccwflag below)

      integer, intent(out) :: ier      ! completion code 0=OK

      !  Poloidal angle orientation: TRUE (default) if dZ/dtheta is
      !  positive on the large major radius size and negative on the
      !  small major radius side of each flux surface; FALSE for the
      !  reverse.

      logical, intent(in), optional :: ccwflag1

      real*8, intent(out), optional :: r1,z1   ! (R,Z) values if desired
      real*8, intent(out), optional :: drdrho1,dzdrho1 ! rho derivatives
      real*8, intent(out), optional :: drdtheta1,dzdtheta1 ! theta deriv.
      real*8, intent(out), optional :: rzdetj1 ! determinant of 2x2 matrix
      real*8, intent(out), optional :: drhodr1,drhodz1 ! grad(rho)
      real*8, intent(out), optional :: dthdr1,dthdz1   ! grad(theta)
      
      !-----------------
      real*8, dimension(1) :: zwkrho,zwkth,zwkr,zwkz
      real*8, dimension(1) :: zwkrrho,zwkzrho,zwkrth,zwkzth
      real*8, dimension(1) :: zdrhodr,zdrhodz,zdthdr,zdthdz
      logical :: jccwflag,idetj
      real*8 :: zdetj,zdenom
      !-----------------

      ier=0

      jccwflag=.TRUE.
      if(present(ccwflag1)) jccwflag = ccwflag1

      zwkrho=rho_in
      zwkth=theta_in

      if(present(r1)) r1=0
      if(present(r1)) z1=0

      if(present(drdrho1)) drdrho1=0
      if(present(dzdrho1)) dzdrho1=0

      if(present(drdtheta1)) drdtheta1=0
      if(present(dzdtheta1)) dzdtheta1=0

      idetj=.FALSE.

      if(present(rzdetj1)) then
         idetj=.TRUE.
         rzdetj1=0
      endif

      if(present(drhodr1)) then
         idetj=.TRUE.
         drhodr1=0
      endif

      if(present(drhodz1)) then
         idetj=.TRUE.
         drhodz1=0
      endif

      if(present(dthdr1)) then
         idetj=.TRUE.
         dthdr1=0
      endif

      if(present(dthdz1)) then
         idetj=.TRUE.
         dthdz1=0
      endif

      if(present(r1).or.present(z1)) then
         call xplasma_rzjacn(s,zwkrho,zwkth,ier, ccwflag=jccwflag, &
              r=zwkr,z=zwkz)
         if(ier.ne.0) return
         if(present(r1)) r1=zwkr(1)
         if(present(z1)) z1=zwkz(1)
      endif

      if(idetj.or. &
           present(drdrho1).or.present(dzdrho1)) then
         call xplasma_rzjacn(s,zwkrho,zwkth,ier, ccwflag=jccwflag, &
              drdrho=zwkrrho,dzdrho=zwkzrho)
         if(ier.ne.0) return
         if(present(drdrho1)) drdrho1=zwkrrho(1)
         if(present(dzdrho1)) dzdrho1=zwkzrho(1)
      endif

      if(idetj.or. &
           present(drdtheta1).or.present(dzdtheta1)) then
         call xplasma_rzjacn(s,zwkrho,zwkth,ier, ccwflag=jccwflag, &
              drdtheta=zwkrth,dzdtheta=zwkzth)
         if(ier.ne.0) return
         if(present(drdtheta1)) drdtheta1=zwkrth(1)
         if(present(dzdtheta1)) dzdtheta1=zwkzth(1)
      endif

      if(idetj) then

         zdetj=(zwkzrho(1)*zwkrth(1)-zwkrrho(1)*zwkzth(1))
         if(present(rzdetj1)) rzdetj1=zdetj

         if(abs(zdetj).gt.0.0d0) then
            if(present(drhodr1)) drhodr1 = zwkzth(1)/zdetj
            if(present(drhodz1)) drhodz1 =-zwkrth(1)/zdetj
            if(present(dthdr1)) dthdr1 =-zwkrrho(1)/zdetj
            if(present(dthdz1)) dthdz1 = zwkzrho(1)/zdetj
         else
            zdenom=(zwkrrho(1)*zwkrrho(1)+zwkzrho(1)*zwkzrho(1))
            if(present(drhodr1)) drhodr1 = zwkrrho(1)/zdenom
            if(present(drhodz1)) drhodz1 = zwkzrho(1)/zdenom
            if(present(dthdr1)) dthdr1 = 0.0d0
            if(present(dthdz1)) dthdz1 = 0.0d0
         endif
      endif
    end subroutine xplasma_rzjac1

    subroutine xplasma_rzjacn(s,rho_in,theta_in,ier, &
         ccwflag, r,z,drdrho,dzdrho,drdtheta,dzdtheta, &
         rzdetj, drhodr,drhodz, dthdr,dthdz)

      ! evaluate terms associated with 2d Jacobean -- vector version
      !   [dR/drho dR/dtheta]
      !   [dZ/drho dZ/dtheta]

      ! *** all vectors must be of the same length ***

      type (xplasma), pointer :: s
      real*8, dimension(:) :: rho_in   ! vector of radial coordinate in
      real*8, dimension(:) :: theta_in ! vector of poloidal angle coordinate in
      !  (see ccwflag below)

      integer, intent(out) :: ier      ! completion code 0=OK

      !  Poloidal angle orientation: TRUE (default) if dZ/dtheta is
      !  positive on the large major radius size and negative on the
      !  small major radius side of each flux surface; FALSE for the
      !  reverse.

      logical, intent(in), optional :: ccwflag

      real*8, intent(out), dimension(:), optional :: r,z
      ! (R,Z) values if desired

      real*8, intent(out), dimension(:), optional :: drdrho,dzdrho
      ! rho derivatives if desired

      real*8, intent(out), dimension(:), optional :: drdtheta,dzdtheta
      ! theta derivatives if desired

      real*8, intent(out), dimension(:), optional :: rzdetj 
      ! determinant of 2x2 jacobian matrix

      real*8, intent(out), dimension(:), optional :: drhodr,drhodz
      ! grad(rho) components if desired

      real*8, intent(out), dimension(:), optional :: dthdr,dthdz
      ! grad(theta) components if desired

      !--------------------
      real*8, dimension(:,:), allocatable :: zwkbuf
      real*8, dimension(:), allocatable :: zdetj
      real*8 :: zdenom
      integer :: isize_min,isize_max,isize,id_R,id_Z,ids(6),inbuf
      integer :: iddrho(6),iddth(6)
      integer :: ii
      logical :: jccwflag,iaxi,idetj
      integer :: indR,indZ,indRrho,indRth,indZrho,indZth
      !-----------------

      call xplasma_global_info(s,ier, axisymm=iaxi)
      if(ier.ne.0) return
      
      if(.not.iaxi) then
         ier=2509
         call xplasma_errmsg_append(s,'  cannot use xplasma_rzjac...')
         return
      endif

      call xplasma_common_ids(s,ier,id_R=id_R,id_Z=id_Z)
      if(min(id_R,id_Z).eq.0) then

         ier=2504  ! no flux surfaces...
         return

      endif

      inbuf=0
      !-----------------------------

      ! check poloidal angle orientation; set sign term...

      jccwflag=.TRUE.
      if(present(ccwflag)) jccwflag = ccwflag

      !-----------------------------
      ! check vector sizes & clear output vectors that are present...

      isize_min=size(rho_in)
      isize_max=isize_min

      isize_min=min(isize_min,size(theta_in))
      isize_max=max(isize_max,size(theta_in))

      if(present(r)) then
         r=0
         isize_min=min(isize_min,size(r))
         isize_max=max(isize_max,size(r))
         inbuf=inbuf+1
         indR=inbuf
         ids(inbuf)=id_R
         iddrho(inbuf)=0
         iddth(inbuf)=0
      endif

      if(present(z)) then
         z=0
         isize_min=min(isize_min,size(z))
         isize_max=max(isize_max,size(z))
         inbuf=inbuf+1
         indZ=inbuf
         ids(inbuf)=id_Z
         iddrho(inbuf)=0
         iddth(inbuf)=0
      endif

      idetj = present(rzdetj).or. &
           present(drhodr).or.present(drhodz).or. &
           present(dthdr).or.present(dthdz)

      if(present(drdrho)) then
         drdrho=0
         isize_min=min(isize_min,size(drdrho))
         isize_max=max(isize_max,size(drdrho))
      endif

      if(present(drdrho).or.idetj) then
         inbuf=inbuf+1
         indRrho=inbuf
         ids(inbuf)=id_R
         iddrho(inbuf)=1
         iddth(inbuf)=0
      endif

      if(present(dzdrho)) then
         dzdrho=0
         isize_min=min(isize_min,size(dzdrho))
         isize_max=max(isize_max,size(dzdrho))
      endif

      if(present(dzdrho).or.idetj) then
         inbuf=inbuf+1
         indZrho=inbuf
         ids(inbuf)=id_Z
         iddrho(inbuf)=1
         iddth(inbuf)=0
      endif

      if(present(drdtheta)) then
         drdtheta=0
         isize_min=min(isize_min,size(drdtheta))
         isize_max=max(isize_max,size(drdtheta))
      endif

      if(present(drdtheta).or.idetj) then
         inbuf=inbuf+1
         indRth=inbuf
         ids(inbuf)=id_R
         iddrho(inbuf)=0
         iddth(inbuf)=1
      endif

      if(present(dzdtheta)) then
         dzdtheta=0
         isize_min=min(isize_min,size(dzdtheta))
         isize_max=max(isize_max,size(dzdtheta))
      endif

      if(present(dzdtheta).or.idetj) then
         inbuf=inbuf+1
         indZth=inbuf
         ids(inbuf)=id_Z
         iddrho(inbuf)=0
         iddth(inbuf)=1
      endif

      if(present(rzdetj)) then
         rzdetj=0
         isize_min=min(isize_min,size(rzdetj))
         isize_max=max(isize_max,size(rzdetj))
      endif

      if(present(drhodr)) then
         drhodr=0
         isize_min=min(isize_min,size(drhodr))
         isize_max=max(isize_max,size(drhodr))
      endif

      if(present(drhodz)) then
         drhodz=0
         isize_min=min(isize_min,size(drhodz))
         isize_max=max(isize_max,size(drhodz))
      endif

      if(present(dthdr)) then
         dthdr=0
         isize_min=min(isize_min,size(dthdr))
         isize_max=max(isize_max,size(dthdr))
      endif

      if(present(dthdz)) then
         dthdz=0
         isize_min=min(isize_min,size(dthdz))
         isize_max=max(isize_max,size(dthdz))
      endif

      if(isize_min.ne.isize_max) then
         ier=510
         call xplasma_errmsg_append(s, &
              'vector length inconsistency detected in xplasma_rzjac')
         return
      endif

      isize=isize_min

      if(inbuf.eq.0) return

      !-------------------------------------------------------
      !  OK -- allocate workspace

      allocate(zwkbuf(isize,inbuf))

      !  evaluate (R,Z) if requested

      call xplasma_RZeval_2d(s,ids(1:inbuf), &
           xplasma_rho_coord,rho_in,xplasma_theta_coord,theta_in, &
           zwkbuf,ier, ccwflag2=jccwflag, &
           ideriv1s=iddrho(1:inbuf),ideriv2s=iddth(1:inbuf))

      if(present(r)) r = zwkbuf(1:isize,indR)
      if(present(z)) z = zwkbuf(1:isize,indZ)
      if(present(drdrho)) drdrho = zwkbuf(1:isize,indRrho)
      if(present(dzdrho)) dzdrho = zwkbuf(1:isize,indZrho)
      if(present(drdtheta)) drdtheta = zwkbuf(1:isize,indRth)
      if(present(dzdtheta)) dzdtheta = zwkbuf(1:isize,indZth)

      if(idetj) then

         allocate(zdetj(isize))

         zdetj = zwkbuf(1:isize,indZrho)*zwkbuf(1:isize,indRth) - &
              zwkbuf(1:isize,indRrho)*zwkbuf(1:isize,indZth)

         if(present(rzdetj)) rzdetj = zdetj

         do ii=1,isize
            if(abs(zdetj(ii)).gt.0.0d0) then
               if(present(drhodr)) drhodr(ii) =-zwkbuf(ii,indZth)/zdetj(ii)
               if(present(drhodz)) drhodz(ii) = zwkbuf(ii,indRth)/zdetj(ii)
               if(present(dthdr)) dthdr(ii) = zwkbuf(ii,indZrho)/zdetj(ii)
               if(present(dthdz)) dthdz(ii) =-zwkbuf(ii,indRrho)/zdetj(ii)
            else
               zdenom=(zwkbuf(ii,indRrho)**2+zwkbuf(ii,indZrho)**2)
               if(present(drhodr)) drhodr(ii) = zwkbuf(ii,indRrho)/zdenom
               if(present(drhodz)) drhodz(ii) = zwkbuf(ii,indZrho)/zdenom
               if(present(dthdr)) dthdr(ii) = 0.0d0
               if(present(dthdz)) dthdz(ii) = 0.0d0
            endif
         enddo

         deallocate(zdetj)
      endif

      deallocate(zwkbuf)

    end subroutine xplasma_rzjacn

end module xplasma_ctran
