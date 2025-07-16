subroutine eqm_irhofun(id_axis,zlbl,inprof,zprof,iflag,id,ierr)

  use xplasma_obj_instance
  use eq_module

  !  make XPLASMA profile -- integrated quantity; smooth by 1/2 zone width
  !  to insure smooth derivative from spline

  implicit NONE
  integer, intent(in) :: id_axis    ! axis id -- must be rho or akin to rho
  character*(*), intent(in) :: zlbl ! name for profile function to be created
  integer, intent(in) :: inprof     ! size of the integrated data profile
  real*8, intent(in) :: zprof(inprof) ! the integrated data provided
                         ! if inprof = size(x_axis) zprof(1)=0 is expected
                         ! if inprof = size(x_axis)-1 the axial data point
                         !    is presumed to be omitted.
  integer, intent(in) :: iflag      ! =1 -- volume normalization;
                         ! derivative evaluations -> W/m^3, #/sec/m^3, etc.
                                    ! =2 -- area normalization;
                         ! derivative evaluations -> A/m^2 (current density).

  integer, intent(out) :: id        ! id of stored profile (if successful)
  integer, intent(out) :: ierr      ! completion code, 0=OK.

  !  if an error occurs and ierr is set, id=0 will be returned.

  !-------------------------
  integer :: iertmp
  !-------------------------

  call xoi_author_set(iertmp)

  call xplasma_irhofun(s,id_axis,zlbl,inprof,zprof,iflag,id,ierr)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
  endif

  call xoi_author_clear(iertmp)

end subroutine eqm_irhofun

subroutine eqm_irhofun_sm(id_axis,zlbl,inprof,zprof,iflag,zsm,id,ierr)

  use xplasma_obj_instance
  use eq_module

  !  make XPLASMA profile -- integrated quantity; smooth by at least 1/2 
  !  zone width or increased amount to insure smooth derivative from spline.

  !  the multiplier on the amount of smoothing is zsm/(x(2)-x(1)), where
  !  (x(2)-x(1)) is the spacing of the id_axis identified rho grid by rho=0.

  implicit NONE
  integer, intent(in) :: id_axis    ! axis id -- must be rho or akin to rho
  character*(*), intent(in) :: zlbl ! name for profile function to be created
  integer, intent(in) :: inprof     ! size of the integrated data profile
  real*8, intent(in) :: zprof(inprof) ! the integrated data provided
                         ! if inprof = size(x_axis) zprof(1)=0 is expected
                         ! if inprof = size(x_axis)-1 the axial data point
                         !    is presumed to be omitted.
  integer, intent(in) :: iflag      ! =1 -- volume normalization;
                         ! derivative evaluations -> W/m^3, #/sec/m^3, etc.
                                    ! =2 -- area normalization;
                         ! derivative evaluations -> A/m^2 (current density).

  real*8, intent(in) :: zsm    ! smoothing width, in units of id_axis (rho)
  !  this gets converted into a multiplier based on the zone spacing at rho=0
  !  the minimal smoothing is basedon the zone spacing; the maximum smoothing
  !  corresponds to a multiplier that yields (1/4) the profile width based
  !  on the zone spacing at rho=0.

  integer, intent(out) :: id        ! id of stored profile (if successful)
  integer, intent(out) :: ierr      ! completion code, 0=OK.

  !  if an error occurs and ierr is set, id=0 will be returned.

  !-------------------------
  integer :: iertmp
  !-------------------------

  call xoi_author_set(iertmp)

  call xplasma_irhofun(s,id_axis,zlbl,inprof,zprof,iflag,id,ierr, smdelx=zsm)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
  endif

  call xoi_author_clear(iertmp)

end subroutine eqm_irhofun_sm
