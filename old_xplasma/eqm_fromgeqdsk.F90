subroutine eqm_fromgeqdsk(label,zfile,ns,nt1,fbdy,cbdy,irz, &
     ier)

  !  make an xplasma equilibrium representation from contents of
  !  a G-EQDSK file or MDSplus EFIT data.  This is done by finding an
  !  appropriately spaced set of fieldline contours psi=constant from
  !  EFIT's psi(R,Z) data.

  !  dmc Feb. 2001

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  character*(*), intent(in) :: label     ! xplasma eq. label
  character*(*), intent(in) :: zfile     ! path to G-EQDSK file / MDSplus
 
  integer, intent(in) :: ns              ! no. of surfaces incl. axis
         ! abs(ns)= no. of surfaces; if ns.lt.0 use Hermite fit for surfaces
 
  integer, intent(in) :: nt1             ! abs(nt1) = no. of theta pts /surface
  !  sign(nt1): +1 for traditional moments; -1 for moments based on equal arc
  !                                            theta definition
 
  integer, intent(in) :: irz             ! =1 for B(R,Z) calculation
 
  real*8, intent(in) :: fbdy             ! bdy back-off parameter (to start)
  real*8, intent(in) :: cbdy             ! min. curvature parameter
 
  integer, intent(out) :: ier            ! completion code, 0=OK

  !----------------------------------------------

  integer :: icur_psi = 0
  real*8, parameter :: xbrki = 0.9d0

  !----------------------------------------------

  call eqm_fromgeqdsk_iadj(label,zfile,ns,nt1,fbdy,cbdy,icur_psi,xbrki,irz, &
     ier)

end subroutine eqm_fromgeqdsk

subroutine eqm_fromgeqdsk_iadj(label,zfile,ns,nt1,fbdy,cbdy, &
     icur_psi,xbrki,irz, &
     ier)

  !  make an xplasma equilibrium representation from contents of
  !  a G-EQDSK file or MDSplus EFIT data.  This is done by finding an
  !  appropriately spaced set of fieldline contours psi=constant from
  !  EFIT's psi(R,Z) data.

  !  dmc Feb. 2001

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  character*(*), intent(in) :: label     ! xplasma eq. label
  character*(*), intent(in) :: zfile     ! path to G-EQDSK file / MDSplus
 
  integer, intent(in) :: ns              ! no. of surfaces incl. axis
         ! abs(ns)= no. of surfaces; if ns.lt.0 use Hermite fit for surfaces
 
  integer, intent(in) :: nt1             ! abs(nt1) = no. of theta pts /surface
  !  sign(nt1): +1 for traditional moments; -1 for moments based on equal arc
  !                                            theta definition
 
  integer, intent(in) :: icur_psi        ! =0: use EFIT file data QPSI for
  ! q(psi) and EFIT file data CURRENT_A for total enclosed toroidal current;
  ! =1:  discard file values and re-derive q(psi) and enclosed current from
  ! analysis of the Psi(R,Z) data.

  ! Advice: if the Psi(R,Z) data in the EFIT input is from a high quality
  ! free boundary GS solution, icur_psi=1 is better; if not, icur_psi=0 is
  ! better.  EFIT datasets extracted from TRANSP time slices are not based
  ! on a free boundary solution, and, icur_psi=0 should be used.

  !  Most experimental EFIT data can be safely processed either way.

  real*8, intent(in) :: xbrki            ! break location for current profile;
  !  enclosed toroidal current is a smooth parabolic interpolation to the 
  !  edge value, matching point and slope with the current profile computed
  !  from Psi(R,Z), at normalized toroidal flux (xbrki).  Recommended value
  !  is: xbrki=0.9; constraint: 0.75 <= xbrki <= 0.95.

  integer, intent(in) :: irz             ! =1 for B(R,Z) calculation
 
  real*8, intent(in) :: fbdy             ! bdy back-off parameter (to start)
  real*8, intent(in) :: cbdy             ! min. curvature parameter
 
  integer, intent(out) :: ier            ! completion code, 0=OK

  !----------------------------------------------

  !  input arguments
  !  ---------------

  !  label ~ 30 characters -- a label for xplasma

  !  zfile - a file path, e.g. /u/efit/nstx/g123456_6800.eqdsk
  !         i.e. a path to a single EFIT G-EQDSK file...

  !               --or--

  !     EFIT_INDEX:<file-path>(<index>)

  !       where <file-path> is [<directory-path>/]<filename>
  !       gives an "EFIT INDEX ascii file" that associates a set of time
  !       points with a set of separate G-EQDSK files; the file format
  !       is illustrated in this example; "|" denotes start of each text
  !       lines:

  !         |ntimes=4
  !         |path='/pub/outgoing/pshare'
  !         |time=0.15   filename=g008500.00150
  !         |time=0.20   filename=g008500.00200
  !         |time=0.25   filename=g008500.00250
  !         |time=0.294  filename=g008500.00294

  !      and the listed files are G-EQDSK files associated with each time.
  !      The times *must* be in ascending order.

  !      The full path to each G-EQDSK file will be <directory-path>/<filename>
  !      where <directory-path> is what leads to the index file itself (i.e.
  !      the index file and the G-EQDSK files are expected all to be in the
  !      same directory; the path record in the index file is *ignored*).

  !      The field <index> should decode as integer and refers to the
  !      time index; 1 => 1st time, 2 => 2nd time, etc.

  !               --or--

  !          an MDSplus path:

  !     MDS+[/reduce|/noreduce]:<server[:port]|LOCAL>:<tree>(<shot>;t=<time>)

  !    where the initial MDS+ specifies MDSplus not a file;

  !    /reduce means just read the specified timeslice (this is the default)
  !    /noreduce means read and cache data from all available EFIT times

  !    <server...> specifies the MDSplus server, can be "LOCAL"
  !                or "CONNECTED", which means the server connection has
  !                already been made

  !    <tree> specifies the MDSplus tree, e.g. "EFIT01" or "EFIT02"
  !           value can be "OPENED", which means, the tree is already open.
  !           if value "OPENED" is provided, the correct shot number should
  !           still be provided, for information/identification purposes.

  !    <shot> specifies the MDSplus (or experiment) shot number
  !    <time> is an ascii encoding of the currently desired time of interest

  !       examples of valid MDSplus paths:

  !       1.  MDS+:LOCAL:EFIT01(104370;t=0.15)
  !       2.  MDS+:CONNECTED:OPENED(104370;t=0.15)
  !       3.  MDS+/reduce:europa.pppl.gov:8501:EFIT01(104370;t=1.0e-1)
  !       4.  MDS+/noreduce:cmoda.psfc.mit.edu:EFIT09(123456789;t=2.3)

  !     notes:  => in cases 3 & 4, the code connects to the server, opens
  !             the tree, reads the data, closes the tree, and disconnects
  !             from the server.
  !             => in case 1, the code opens the tree, reads the data, and
  !             closes the tree.
  !             => in case 2, the code just reads the data; the caller is
  !             expected to have taken care of any MDS+ connect or open 
  !             operations.

  !  ns -- abs(ns)=number of surfaces, spaced in equal steps of sqrt(tor.flux),
  !       to be generated in the xplasma representation, from the EFIT data
  !       if ns.lt.0, use Hermite instead of Spline fits for R(rho,chi),
  !       Z(rho,chi), and B(rho,chi) components computed in xplasma
  !       NOTE: spline is still calculated; it is replaced with a Hermite
  !       representation.

  !  nt1 -- abs(nt1) = number of theta points per surface in the xplasma 
  !         representation; sign(nt1) = -1 for moments based on equal arc
  !         definition of theta; +1 for traditional VMEC/descur power
  !         spectrum minimization definition
  !         

  !  fbdy -- delta_psi(use)/delta_psi(EFIT-bdy), btw 0.95 and 1.00
  !         use this to specify a plasma boundary "slightly in from" the
  !         plasma boundary given by EFIT.

  !    note that eqi_fromgeqdsk will also back off from the nominal EFIT
  !    boundary, if when calculating its fieldline contour it is found not
  !    to close.  So:  fbdy=1.0d0 as input is safe and usually desirable.

  !  cbdy -- minimum allowed boundary contour local curvature, normalized
  !    to boundary midplane half radius.  Setting this .gt.0 is another way
  !    to bind the boundary surface away from the separatrix-- the software
  !    will automatically reduce delta_psi(use) as needed to satisfy this
  !    condition.

  !  irz -- =1: extend the xplasma representation and coordinate scheme to
  !    cover EFIT's R,Z box; define psi(R,Z), B(R,Z) as xplasma objects
  !         =0: just form R(theta,rho),Z(theta,rho), psi(rho), g(rho), q(rho)
  !    P(rho) but do not save the (R,Z) information.
  !    ...irz=1 requires more work

  !----------------

  !  output argument:  ier=0 denotes success

  !    if successful, xplasma interpolation routines are available for
  !    evaluation of metrics, field, g, q, etc.,

  !----------------
  integer :: iertmp,iwarn,ikmom,ibrki
  real*8, parameter :: mom_rtol = 1.0d-6
  real*8, parameter :: xbrk_ax0 = 0.15d0
  !----------------
  !  access the EFIT data; build xplasma

  call xmoments_kmom_get(ikmom)

  call eqm_select(label,1)   !! this deletes old xplasma contents !!

  call xoi_author_set(iertmp)

  call xplasma_fromgeqdsk(s,zfile,ier, &
       new_device=.TRUE., nrho=abs(ns), nth=abs(nt1), &
       bdy0=fbdy, cratio_min=cbdy, rzspace = (irz.eq.1), &
       kpsi_cur=icur_psi, rho_curbreak=xbrki, &
       kmom = ikmom, mom_rtol=mom_rtol, rhobrk_adj_axis = xbrk_ax0, &
       equal_arc = (nt1.lt.0), lhermite = (ns.lt.0) )

  call xoi_author_clear(iertmp)

  call eqi_fromg_warn(iwarn)
  if((iwarn.gt.0).or.(ier.ne.0)) call xplasma_error(s,ier,lunerr)

#ifdef __DEBUG
  call eqdbg_plot
#endif

end subroutine eqm_fromgeqdsk_iadj

subroutine eqm_geqdsk_ppadj(ier)

  !  Replace presure profile (based on P(psi) from G-eqdsk file originally)
  !  Use P(bdy) and integrate inward using dP/dPsi with a finite difference
  !  approximation.  The profile "DPDPSI" must be defined.  The original
  !  "P" profile is renamed "P_ORIG".

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  integer, intent(out) :: ier
  !----------------------------------
  integer :: iertmp
  !----------------------------------

  call xoi_author_set(iertmp)

  call xplasma_geqdsk_ppadj(s,ier)
  if(ier.ne.0) then
     write(lunerr,*) ' ?eqm_geqdsk_ppadj: xplasma routine reports error:'
     call xplasma_error(s,ier,lunerr)
  endif

  call xoi_author_clear(iertmp)

#ifdef __DEBUG
  call eqdbg_plot
#endif

end subroutine eqm_geqdsk_ppadj
