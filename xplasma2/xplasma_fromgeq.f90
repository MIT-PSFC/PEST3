module xplasma_fromgeq

  ! module with routines that help to create various types of xplasma profiles

  use xplasma_obj
  use eqi_rzbox_module

  implicit NONE

  private

  public :: xplasma_fromgeqdsk, xplasma_geqdsk_ppadj
  public :: xplasma_save_psi0, xplasma_get_psi0

  contains

    !-------------------------------------------------------
    subroutine xplasma_fromgeqdsk(s,filename,ier, &
         new_device, nrho, nth, i2pi, lhermite, &
         bdy0, cratio_min, rzspace, &
         kpsi_cur,rho_curbreak, &
         kmom, rhobrk_adj_axis, mom_rtol, &
         equal_arc, &
         gs_errmax)

      ! (re)build equilibrium inverse represenation from EFIT g-eqdsk data
      !    R(rho,theta), Z(rho,theta), g(rho),q(rho),P(rho).

      type (xplasma), pointer :: s

      character*(*), intent(in) :: filename  ! Geqdsk filename or MDS+ path
      !                  ...or...EFIT_INDEX:filename(t=<time>) specification
      !                  see comments in source/xplasma/eqm_fromgeqdsk.f90

      integer, intent(out) :: ier  ! integer status code (returned; 0=OK)

      !---------------
      !  optional...

      logical, intent(in), optional :: new_device ! flag, .TRUE. if reading
      !  equilibrium for a new device.  Default: .TRUE. when reading the
      !  first equilibrium in the history of the xplasma object; .FALSE.
      !  subsequently.  When this is not .TRUE., the limiter data is presumed
      !  already known to xplasma.

      integer, intent(in), optional :: nrho  ! no. of flux surfaces 
      !                                        covering rho = [0:1]
      !  rho = normalized(sqrt(Psi-tor/Psi-tor-at-bdy))
      !        for 20 zones specify nrho=21 -- count includes axis and plasma
      !        boundary-- evenly spaced in rho
      !  DEFAULT: nrho=51 if new_device is .TRUE.; otherwise the value from
      !  the old equilibrium is taken.

      integer, intent(in), optional :: nth   ! no. of theta points -- no. of 
      !        evenly spaced points spanning 0:2pi inclusive.
      !  DEFAULT: nth=101 if new_device is .TRUE.; otherwise the value from
      !  the old equilibrium is taken.

      integer, intent(in), optional :: i2pi  ! set =1 for theta range to
      !  be -pi:pi instead of 0:2pi (the latter is the default).

      logical, intent(in), optional :: lhermite ! flag, .TRUE. for Hermite
      !  interpolation of {R,Z}(rho,theta) -- by default, bicubic spline
      !  interpolation is used.

      real*8, intent(in), optional :: bdy0   ! initial Psi_bdy/Psi_sep
      !  (default 1.0d0) fraction of Psi to enclose in mapped region
      !  The issue here is that a separatrix surface generally cannot be
      !  used as a boundary in an inverse representation, due to q=infinity
      !  and the kink in the surface shape.  So instead one backs off from
      !  the Psi_sep flux, defining a boundary just inside the separatrix.
      !    Recommendation: the default.  "cratio_min" causes the code to
      !  iterate to find a suitable boundary surface.
      !  Accepted range: 0.95d0 to 1.00d0

      real*8, intent(in), optional :: cratio_min  ! minimum curvature ratio
      !  (default 0.080d0) minimum acceptable local curvature radius, 
      !  normalized to midplane half width.  A separatrix kink has a curvature
      !  radius of zero.  This is tested by fitting a circle through successive
      !  sequences of 3 points on the boundary grid [1:nth] and looking at this
      !  circle's radius compared to an estimate of the plasma midplane half-
      !  width.
      !    Recommendation: the default, based on experience.
      !  Accepted range: 0.01d0 to 0.15d0

      logical, intent(in), optional :: rzspace ! .TRUE. to save (R,Z)
      !  grids using the EFIT grids, and compute B(R,Z).  Default = .TRUE.
      !  Setting this .FALSE. may speed some applications.

      integer, intent(in), optional :: kpsi_cur    ! =0 (default) to set
      !  total toroidal current from the EFIT file data (scalar CURRENT_A),
      !  and use EFIT q(psi) data (QPSI) for the q profile in the equilibrium
      !  which affects flux surface spacing.

      !  If instead kpsi_cur = 1 is set, the EFIT total toroidal current 
      !  magnitude and QPSI data are discarded, and the q(psi) and enclosed
      !  current profiles are derived from contour integrals over the 
      !  Psi(R,Z) data.

      !  Advice: if Psi(R,Z) is derived from a high quality free boundary
      !  solution, use kpsi_cur = 1; if unsure, use the defaults.

      real*8, intent(in), optional :: rho_curbreak    ! This number,
      !  constrained to be between 0.75 and 0.95, specifies where to
      !  break the enclosed toroidal current profile and impose a smooth
      !  parabolic interpolation to the edge value (see kpsi_cur for how
      !  this is set).  This is necessary to get smooth dI/drho and <J.B>
      !  profiles out.  Default: 0.90.

      integer, intent(in), optional :: kmom    ! number of Fourier moments
      !  to use in Fourier Spline representation of equilibrium 
      !  (default: not set here; xplasma default is 16).

      real*8, intent(in), optional :: rhobrk_adj_axis ! delta(rho) for
      !  moments trimming operations, and, the maximum rho where flux
      !  surfaces can be touched, in order to assure that d/drho[R0,Z0] --> 0
      !  as rho--> 0, an attributed needed by codes (e.g. RF solvers) that
      !  use the inverted equilibrium.
      !    Default value: 0.15d0.  To disable this feature, set 
      !  rhobrk_adj_axis = -1.0d0 or any other negative number

      real*8, intent(in), optional :: mom_rtol ! moments trimming ratio
      !  (default 1.0d-6).  The EFIT Psi(R,Z) inversion involves a Fourier
      !  representation; this control specifies a method for trimming away
      !  higher order moments that are very small over wide regions of space.
      !  Accepted range: 0.0d0 to 1.0d-4

      logical, intent(in), optional :: equal_arc  ! .TRUE. to use an equal
      !  arc definition for the poloidal angle coordinate; .FALSE. (default)
      !  for VMEC/descur moments power spectrum minimization definition, the
      !  traditional definition for past usage.

      real*8, intent(out), optional :: gs_errmax  ! normalized measure of
      !  Grad Shafranof (GS) error in the EFIT equilibrium.

      !-----------------------------
      real*8 :: zbdy0,zcratio,zmom_rtol,zrhobrk_axis,zerrmax
      integer :: irz,inrho,inth,id,idx,ieq,iertmp,jopt,j2pi
      real*8 :: zdum,xbrki

      logical :: ihermite,inew_device,ienable
      integer :: ipsi_cur

      character*32 :: glist,gelems(7),zauth
      real*8 :: zelems(7)
      character*140 :: ztext(7),zlabel

      real*8 :: rhobrk_min=0.75d0
      real*8 :: rhobrk_max=0.95d0

      !-----------------------------

      ier=0

      call xplasma_global_info(s,ier, initLabel=zlabel, eq_counter=ieq)
      if(ier.ne.0) return

      call xplasma_author_query(s,zauth,ienable,ier)
      if(ier.ne.0) return

      !-----------------------------
      !  check if this an equilibrium for a new device;
      !  if so xplasma re-initialization may be required...

      if(ieq.eq.0) then
         inew_device=.TRUE.
      else 
         inew_device=.FALSE.
         if(present(new_device)) inew_device=new_device
      endif

      if(inew_device.and.(ieq.gt.0)) then
         ! re-initialize xplasma: changing devices!
         call xplasma_init(s,zlabel,iertmp)
         ier=max(ier,iertmp)
         call xplasma_author_set(s,zauth,iertmp)
         ier=max(ier,iertmp)
      endif
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,'error in xplasma_fromgeqdsk.')
         return
      endif

      j2pi=2
      if(present(i2pi)) j2pi=i2pi

      if(present(kmom)) then
         call xplasma_kmom_set(s,kmom,ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s,'error in xplasma_fromgeqdsk.')
            return
         endif
      endif

      ipsi_cur=0
      xbrki=0.9d0

      if(present(kpsi_cur)) then
         ipsi_cur = max(0,min(1,kpsi_cur))
      endif

      if(present(rho_curbreak)) then
         xbrki=min(rhobrk_max,max(rhobrk_min, rho_curbreak))
      endif

      !-----------------------------

      if(inew_device) then
         !  use arguments
         if(.not.present(nrho)) then
            inrho=51
         else
            inrho=max(5,nrho)
         endif

         if(.not.present(nth)) then
            inth=101
         else
            inth=max(21,nth)
         endif

      else
         !  use previous equilibrium grid sizes
         call xplasma_gridId(s,'__RHO',idx)
         call xplasma_grid_size(s,idx,inrho,iertmp)
         ier=max(ier,iertmp)

         call xplasma_gridId(s,'__CHI',idx)
         call xplasma_grid_size(s,idx,inth,iertmp)
         ier=max(ier,iertmp)
      endif
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,'error in xplasma_fromgeqdsk.')
         return
      endif

      ihermite=.FALSE.
      if(present(lhermite)) ihermite=lhermite

      irz=1
      if(inew_device.and.present(rzspace)) then
         if(rzspace) then
            irz=1
         else
            irz=0
         endif
      endif

      glist = '__GFILE_INFO'

      gelems(1) = 'FILENAME'
      zelems(1) = 0
      ztext(1) = filename

      zbdy0=1.0d0
      if(present(bdy0)) then
         zdum=0.95d0
         zbdy0=max(zdum,min(zbdy0,bdy0))  ! enforce range 0.95d0 to 1.00d0
      endif
      gelems(2)='BDY0'
      zelems(2)=zbdy0
      ztext(2)='Relative start point of boundary search'

      zcratio=0.080d0  ! **DEFAULT** curvature ratio limit

      if(present(cratio_min)) then
         zdum=0.01d0
         zcratio=max(zdum,cratio_min)
         zdum=0.15d0
         zcratio=min(zdum,zcratio)        ! enforce range 0.01d0 to 0.15d0
      endif
      gelems(3)='CRATIO'
      zelems(3)=zcratio
      ztext(3)='minimum accepted boundary curvature ratio'

      zrhobrk_axis=0.15d0
      if(present(rhobrk_adj_axis)) then
         zdum=0.30d0
         zrhobrk_axis=min(zdum,rhobrk_adj_axis)
      endif
      gelems(4)='RHOBRK_ADJ_AXIS'
      zelems(4)=zrhobrk_axis
      ztext(4)='width of interior region affected by axis adjustment'

      zmom_rtol=1.0d-6
      if(present(mom_rtol)) then
         zdum=0.0d0
         zmom_rtol=max(zdum,mom_rtol)
         zdum=1.0d-4
         zmom_rtol=min(zdum,zmom_rtol)    ! enforce range 0.0 to 1.0d-4
      endif
      gelems(5)='MOM_RTOL'
      zelems(5)=zmom_rtol
      ztext(5)='trimming parameter for Fourier moments representation'

      gelems(6)='QPSI_OPTION'
      zelems(6)=ipsi_cur
      ztext(6)='value of 1.0 indicates Q(psi) & total current from Psi(R,Z)'

      gelems(7)='RHO_CURBREAK'
      zelems(7)=rho_curbreak
      ztext(7)='break location for matching current profile to edge value.'

      zerrmax=0.0d0
      if(present(gs_errmax)) gs_errmax=zerrmax

      !---------------------------------------------
      ! record geqdsk controls in a list... id not saved...

      call xplasma_create_list(s,glist,gelems,id,ier, &
           label='xplasma_fromgeqdsk inputs list', &
           r8vals=zelems, chvals=ztext)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              'xplasma_fromgeqdsk: '// &
              'error occurred while attempting to define list: '//trim(glist))
         return
      endif

      !---------------------------------------------
      ! OK, build the xplasma equilibrium...

      sp => s

      jopt = 0
      if(present(equal_arc)) then
         if(equal_arc) then
            jopt = 1
            if(.not.present(kmom)) then
               call xplasma_kmom_set(s,64,ier)  ! set FFT moments up to 64
               ! ...for equal arc theta
            endif
         endif
      endif

      call eqi_fromgeqdsk(filename, ihermite, inew_device, &
           inrho, inth, j2pi, zbdy0, zcratio, irz, jopt, &
           ipsi_cur,xbrkI, &
           zrhobrk_axis, zmom_rtol, zerrmax, ier)

      if(ier.ne.0) ier=150

      if(present(gs_errmax)) gs_errmax=zerrmax

    end subroutine xplasma_fromgeqdsk

    subroutine xplasma_geqdsk_ppadj(s,ier,match_auth)

      !  reset G-eqdsk pressure profile "P" using P(bdy) and 
      !  integrating inward using "DPDPSI"; old profile is renamed
      !  "P_ORIG".

      !--------------------------------------------------------
      !  arguments...

      type (xplasma), pointer :: s

      integer, intent(out) :: ier  ! integer status code (returned; 0=OK)

      !---------------
      !  optional...

      logical, intent(in), optional :: match_auth

      !--------------------------------------------------------
      !  local variables:

      logical :: iauth_match
      integer :: id,id_p,id_rho,id_psi,id_dpdpsi,iertmp
      integer :: irho,inrho,irank,icorr
      character*32 :: zauth,zunits
      character*80 :: zlabel

      real*8, dimension(:), allocatable :: zrho,zp,zpsi,zdpdpsi
      real*8 :: zdp,zdpsi,zcorr,zcsum,zpnew,zdelp_norm,zdelp_old,zdelp_new
      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter :: ONE =  1.0d0
      real*8, parameter :: HALF = 0.5d0

      !--------------------------------------------------------
      !  first,check if necessary data is available:

      call xplasma_profId(s,'P',id_p)
      call xplasma_profId(s,'dPdPsi',id_dpdpsi)
      call xplasma_profId(s,'Psi',id_psi)

      ier=0
      if(id_p.eq.0) then
         call xplasma_errmsg_append(s, &
              '?xplasma_geqdsk_ppadj: "P" profile not found.')
         ier=ier+1
      endif
      if(id_psi.eq.0) then
         call xplasma_errmsg_append(s, &
              '?xplasma_geqdsk_ppadj: "Psi" profile not found.')
         ier=ier+1
      endif
      if(id_dpdpsi.eq.0) then
         call xplasma_errmsg_append(s, &
              '?xplasma_geqdsk_ppadj: "dPdPsi" profile not found.')
         ier=ier+1
      endif

      if(ier.ne.0) then
         ier=9999
         return
      endif

      !--------------------------------------------------------
      !  reset authorship if requested

      iauth_match = .FALSE.
      if(present(match_auth)) iauth_match = match_auth

      if(iauth_match) then
         call xplasma_get_item_info(s,id_P,iertmp, author=zauth)
         call xplasma_author_set(s,zauth,iertmp)
      endif

      do

         !--------------------------------------------------------
         !  copy old P profile

         call xplasma_copy_item(s,id_P,'P_ORIG',id,ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 '?xplasma_geqdsk_ppadj: xplasma_copy_item failed.')
            exit
         endif

         !--------------------------------------------------------
         !  get grid of P profile; allocate working arrays

         call xplasma_prof_info(s,id_P,ier, rank=irank, gridId1=id_rho, &
              label=zlabel, units=zunits)
         if(ier.ne.0) exit

         if(irank.ne.1) then
            call xplasma_errmsg_append(s, &
                 '?xplasma_geqdsk_ppadj: "P" profile rank is not 1.')
            ier=9999
            exit
         endif

         call xplasma_grid_info(s, id_rho, ier, size=inrho)
         if(ier.ne.0) exit

         allocate(zrho(inrho),zp(inrho),zpsi(inrho),zdpdpsi(inrho))

         call xplasma_grid(s, id_rho, zrho, ier)
         if(ier.ne.0) exit

         call xplasma_eval_prof(s,id_p,zrho,zp,ier)
         if(ier.ne.0) exit

         call xplasma_eval_prof(s,id_psi,zrho,zpsi,ier)
         if(ier.ne.0) exit

         call xplasma_eval_prof(s,id_dpdpsi,zrho,zdpdpsi,ier)
         if(ier.ne.0) exit

         !--------------------------------------------------------
         !  compute new P profile in two passes; 1st pass, compute
         !  an averaged correction factor for the integration to reach
         !  the original P(0)...

         zcorr=ONE

         icorr=0
         zcsum=ZERO
         zdelp_norm=(maxval(zp)-zp(inrho))*HALF

         do irho=inrho-1,1,-1
            zdpsi = zpsi(irho+1)-zpsi(irho)
            zdp = zdpsi * (zdpdpsi(irho) + zdpdpsi(irho+1))*HALF

            zdelp_old = zp(irho)-zp(inrho)

            !  minus sign due to inward direction of integration
            zp(irho) = zp(irho+1) - zdp

            if(zdelp_old.gt.zdelp_norm) then
               zdelp_new = zp(irho)-zp(inrho)
               icorr = icorr + 1
               zcsum = zcsum + (zdelp_new/zdelp_old)
            endif
         enddo

         if(icorr.gt.0) then
            zcorr = ONE/(zcsum/icorr)

            do irho=inrho-1,1,-1
               zdpsi = zpsi(irho+1)-zpsi(irho)
               zdp = zdpsi * (zdpdpsi(irho) + zdpdpsi(irho+1))*HALF

               !  minus sign due to inward direction of integration
               zp(irho) = zp(irho+1) - zcorr*zdp
            enddo
         endif
            
         !--------------------------------------------------------
         !  store new P profile

         call xplasma_create_1dprof(s,'P',id_rho,zp,id,ier, &
              ispline=2,ibca=1,zbca=ZERO,ibcb=1, &
              label=trim(zlabel), units=trim(zunits))

         exit
      enddo
      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
                 '?xplasma_geqdsk_ppadj: an error occurred.')
      endif

      !--------------------------------------------------------
      !  cleanup

      if(allocated(zp)) deallocate(zrho,zp,zpsi,zdpdpsi)

      if(iauth_match) then
         call xplasma_author_clear(s,zauth,iertmp)
      endif

    end subroutine xplasma_geqdsk_ppadj

    subroutine xplasma_save_psi0(s,psi0_val,ierr)

      ! save Psi0 offset to machine axis
      !   Psi(R,Z) is saved in xplasma with Psi=0 @ the magnetic axis
      !   Psi0+psi(R,Z) restores the EFIT conventional Psi = R*A_phi,
      !     which has psi->0 as R->0

      type (xplasma), pointer :: s

      real*8, intent(in) :: psi0_val  ! Psi value to save (pol. flux, Wb/rad).
      integer, intent(out) :: ierr    ! status code returned, 0=OK

      !--------------

      integer :: ivals(1),id_list
      real*8 :: r8vals(1)
      character*64 :: chvals(1)
      character*32 :: enames(1)

      !--------------

      enames(1) = 'Psi0'
      ivals(1)=1
      r8vals(1)=Psi0_val
      chvals(1)='Psi offset to machine axis'

      call xplasma_create_list(s,'Psi0_container',enames,id_list,ierr, &
           label='Container: Psi offset to machine axis', units='Wb/rad', &
           ivals=ivals, r8vals=r8vals, chvals=chvals)

    end subroutine xplasma_save_psi0

    subroutine xplasma_get_psi0(s,psi0_val,ifind,ierr)

      ! retrieve Psi0 offset to machine axis, if available
      ! This is the value saved by xplasma_save_Psi0

      type (xplasma), pointer :: s

      real*8, intent(out) :: psi0_val ! Psi value to save (pol. flux, Wb/rad).
      integer, intent(out) :: ifind   ! =1 if Psi0 value found, 0 otherwise.
      integer, intent(out) :: ierr    ! status code returned, 0=OK

      ! if Psi0 value is not found, ierr=0 and ifind=0 are returned.

      !--------------------------

      integer :: indx_psi0,indx,nlist
      real*8, dimension(:), allocatable :: psi0a
      character*32, dimension(:), allocatable :: enames
      integer id_list                   ! xplasma id:  Psi0 container list

      !--------------------------

      ifind=0
      ierr=0
      psi0_val = 0.0d0

      call xplasma_listID(s,'Psi0_container',id_list)
      if(id_list.eq.0) return

      call xplasma_getList_size(s,id_list,nlist,ierr)
      if(ierr.ne.0) return
  
      allocate(enames(nlist),psi0a(nlist))

      do
         call xplasma_getList_names(s,id_list,enames,ierr)
         if(ierr.ne.0) exit

         call xplasma_getList_r8vals(s,id_list,psi0a,ierr)
         if(ierr.ne.0) exit

         indx_psi0=0
         do indx=1,nlist
            call uupper(enames(indx))
            if(enames(indx).eq.'PSI0') then
               indx_psi0=indx
               exit
            endif
         enddo

         if(indx_psi0.gt.0) then
            psi0_val = psi0a(indx_psi0)
            ifind=1
         endif

         exit
      enddo

      deallocate(enames,psi0a)

    end subroutine xplasma_get_psi0

end module xplasma_fromgeq
