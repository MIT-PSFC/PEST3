module xplasma_obj
  !
  ! *** FORTRAN-95 REQUIRED ***
  !   self-initializing compound data types defined with allocatable array
  !   elements
  !
  !-----------------------------------------------------------------
  ! maintenance note:  any change to definition of XPLASMA object
  !   must result in changes to XPLASMA pack/unpack and read/write
  !   methods, to assure completeness of these operations!
  ! In some cases e.g. moments data, the approach is to use XPLASMA's
  !   own black box method to create storage for data in a form 
  !   already known to the I/O routines.
  !
  !-----------------------------------------------------------------
  ! Mod DMC Jan 2010: xplasma_trace_name & global counter added, for
  !    debugging.
  !
  ! Mod DMC (following Pest-III developer Dylan Brennan) ajac_maxvar
  !    default value reduced from 1.0e-2 to 1.0e-3 -- allows more poloidal
  !    variation in det(J) to be accepted, by default
  !
  !-----------------------------------------------------------------
  !
  ! define the XPLASMA data type-- a container for a sorted list of items
  ! (i.e. a dictionary) each element of which consists of:
  !   a name
  !   a label
  !   a data type
  !   an index pointing to an element in an array of the indicated data
  !     types
  !
  ! the index points to the particular instance of a data type to which
  ! the dictionary entry refers.
  !
  ! the container maintains a 1d list of each supported data type
  !
  ! submodules define the supported data types.
  !
  ! supported data types:
  !
  !   dtype=1 -- list datatype
  !   dtype=2 -- coordinate datatype (container for grids).
  !   dtype=3 -- grid datatype (each one a discretization of a coordinate).
  !   dtype=4 -- profile datatype (splines, etc., defined over some grids).
  !   dtype=5 -- "black box" datatype (each black box contains an integer
  !              buffer of some length and a real*8 buffer of some length).
  !
  ! to add a data type:
  !   (1) define the (private) data type.
  !   (2) add private methods (read, write, copy, init).
  !   (3) add dictionary support (expandable pointer array).
  !   (4) add public support (create, remove, query, use, ...)
  ! in general: use existing types as a model...
  !
  ! modification DMC Aug 2006-- reduced use of pointers; memory management
  ! modified.  Pointers are used only for expandable lists; there is no
  ! garbage collection; "free lists" are used to track availability of slots
  ! for new data items; when free lists are empty, item lists are expanded;
  ! entries for deleted items are returned to the free lists.
  !-----------------------------------------

  use ezcdf
  use logmod   ! jec 08Jan2010 - see xplist_free
  use execsystem
  !
  implicit NONE

  !----------------
  ! The container type and the "xpeval" type both are public...

  private
  public :: xplasma,xpeval
  
  ! container contains additional derived types but these are accessible
  ! only via the public interface
  !----------------
  ! generic methods -- xplasma container; dictionary

  public :: xplasma_init, xplasma_free, xplasma_relabel
  public :: xplasma_write, xplasma_read, xplasma_copy, xplasma_avg_elements
  public :: xplasma_pack, xplasma_unpack
  public :: xpdebug

  public ::  xplasma_label_item, &
       xplasma_copy_item, xplasma_remove_item, &
       xplasma_find_item, xplasma_get_item_info, &
       xplasma_num_items, xplasma_contents, &
       xplasma_find_grids, xplasma_find_profs, xplasma_find_blkbxs, &
       xplasma_get_counter, xplasma_reset_author

  public :: &
       xplasma_listId,xplasma_coordId,xplasma_gridId,xplasma_profId, &
       xplasma_blkbxId

  public :: &
       xplasma_author_set, xplasma_author_pause, &
       xplasma_author_resume, xplasma_author_clear, &
       xplasma_author_query, xplasma_remove_author
  
  public :: &
       xplasma_time_set, xplasma_time_get, &
       xplasma_bdytol_set, xplasma_rzOrder_set, &
       xplasma_ajac_maxvar_set, xplasma_field_ccw_set, &
       xplasma_kmom_set, xplasma_misc_set, &
       xplasma_eqCount_incr, &
       xplasma_global_info, xplasma_common_ids, xplasma_set_pressure_id

  public :: xplasma_error,xplasma_errmsg_append,xplasma_clear_msgs

  !----------------
  ! methods for Fourier Spline

  public :: xpmom_init,xpmom_setup,xpmom_info,xpmom_get1,xpmom_free
  public :: xpmom_protect,xpmom_unprotect

  !----------------
  ! methods for *list* items in the container
  public :: xplasma_create_list,xplasma_merge_lists

  public :: xplasma_list_info
  public :: xplasma_getList_size
  public :: xplasma_getList_names

  public :: xplasma_setList_ivals,xplasma_setList_r8vals
  public :: xplasma_setList_chvals,xplasma_setList_r4vals
  public :: xplasma_getList_ivals,xplasma_getList_r8vals
  public :: xplasma_getList_chvals,xplasma_getList_r4vals
  public :: xplasma_getList,xplasma_setList

  interface xplasma_setList
     module procedure &
          xplasma_setList_ivals, xplasma_setList_r8vals, &
          xplasma_setList_chvals, xplasma_setList_r4vals
  end interface

  interface xplasma_getList
     module procedure xplasma_getList_size, &
          xplasma_getList_ivals, xplasma_getList_r8vals, &
          xplasma_getList_chvals, xplasma_getList_r4vals
  end interface

  !----------------
  ! public methods for *coordinate* data type

  public :: xplasma_create_coord
  public :: xplasma_coord_info
  public :: xplasma_coord_gridId
  public :: xplasma_coord_isPeriodic

  !----------------
  ! public methods for *grid* data type

  public :: xplasma_create_grid
  public :: xplasma_grid_size
  public :: xplasma_grid_info
  public :: xplasma_grid

  ! careful, this is powerful:
  public :: xplasma_remove_grid  ! remove grid *AND* profiles that use grid

  !----------------
  ! public methods for *profile* data type
  !  create methods...

  public :: xplasma_create_1dprof
  public :: xplasma_create_2dprof
  public :: xplasma_create_3dprof

  public :: xplasma_prof_lintrans

  public :: xplasma_create_prof
  interface xplasma_create_prof
     module procedure xplasma_create_1dprof, &
          xplasma_create_2dprof, xplasma_create_3dprof
  end interface

  !  query/info method
  public :: xplasma_prof_info
  public :: xplasma_prof_gridInfo

  !  evaluation methods
  public :: xplasma_getProf_1data
  public :: xplasma_getProf_2data
  public :: xplasma_getProf_3data

  public :: xplasma_eval_1dprofx,xplasma_eval_1dprofxs
  public :: xplasma_eval_1dprofsx,xplasma_eval_1dprofsxs

  public :: xplasma_eval_2dprofx,xplasma_eval_2dprofxs
  public :: xplasma_eval_2dprofsx,xplasma_eval_2dprofsxs

  public :: xplasma_eval_3dprofx,xplasma_eval_3dprofxs
  public :: xplasma_eval_3dprofsx,xplasma_eval_3dprofsxs

  public :: xplasma_xpeval_1dprof, xplasma_xpeval_2dprof, xplasma_xpeval_3dprof

  public :: xplasma_psi_range
  public :: xplasma_eval_prof

  interface xplasma_eval_prof
     module procedure xplasma_getProf_1data, &
          xplasma_getProf_2data, xplasma_getProf_3data, &
          xplasma_eval_1dprofx,xplasma_eval_1dprofxs, &
          xplasma_eval_1dprofsx,xplasma_eval_1dprofsxs, &
          xplasma_eval_2dprofx,xplasma_eval_2dprofxs, &
          xplasma_eval_2dprofsx,xplasma_eval_2dprofsxs, &
          xplasma_eval_3dprofx,xplasma_eval_3dprofxs, &
          xplasma_eval_3dprofsx,xplasma_eval_3dprofsxs, &
          xplasma_xpeval_1dprof, xplasma_xpeval_2dprof, xplasma_xpeval_3dprof
  end interface

  public :: xpeval_info

  public :: xplasma_diff_prof
  public xplasma_diff_1dprof, xplasma_diff_2dprof, xplasma_diff_3dprof

  interface xplasma_diff_prof
     module procedure xplasma_diff_1dprof, xplasma_diff_2dprof, &
          xplasma_diff_3dprof
  end interface

  !-----------------
  !  this routine allows conversion of Splines to Hermite representation
  !  the C2 property of the splines is preserved

  public :: xplasma_hermitize

  !-----------------
  !  this routine gives evaluation of R(rho,theta),Z(rho,theta) and
  !  derivatives, with automatic failover to R & Z defined over
  !  extrapolated flux coordinates as needed:

  public :: xplasma_RZeval_2d
  !-----------------

  public :: xplasma_x_lookup1, xplasma_x_lookupn, xplasma_x_lookup
  public :: xpeval_free
  interface xplasma_x_lookup
     module procedure xplasma_x_lookup1,xplasma_x_lookupn
  end interface

  !----------------
  !  specialized...  an (R,Z) -> (rho,theta) inverse map estimator

  public :: xplasma_rhoth_inv

  !----------------
  ! methods for creating numerical integration capability
  !   (xplasma_integ_dva is public for technical reasons but
  !    need not be advertised to users...)

  public :: xplasma_create_integ,xplasma_augment_integ
  public :: xplasma_integ_dva,xplasma_integ_info,xplasma_integ_scache

  !  integrand cache #fcns parameter -- should match specification in
  !  eqi_intg_module: eqi_intg_nfi

  integer, parameter :: eqi_intg_neval = 10 ! dmc Mar 2006 -- this is verified.

  !----------------
  ! generic "black box" methods

  public :: xplasma_create_blackBox
  public :: xplasma_blackBox_info
  public :: xplasma_blackBox_retrieve

  !----------------------------------------------
  !  xplasma known data type id codes...

  integer, parameter, public :: xplasma_listType = 1 ! LIST data type code
  integer, parameter, public :: xplasma_coordType= 2 ! COORDinate data type
  integer, parameter, public :: xplasma_gridType = 3 ! GRID data type code
  integer, parameter, public :: xplasma_profType = 4 ! PROFILE data type code
  integer, parameter, public :: xplasma_blkbxType= 5 ! BLACK BOX data type code

  !----------------------------------------------
  logical, parameter, public :: xp_use_coord_vecsize = .TRUE.
  logical, parameter, public :: xp_set_coord_vecsize = .FALSE.

  !----------------------------------------------
  !  standard grid coordinate id codes
  !  NOTE the ordering must be consistent with the order of definition
  !  of standard coordinates in xplasma_init

  integer, parameter, public :: xplasma_rho_coord = 1   ! radial flux coord.
  integer, parameter, public :: xplasma_theta_coord = 2 ! poloidal angle coord.
  integer, parameter, public :: xplasma_phi_coord = 3   ! toroidal angle coord.
  integer, parameter, public :: xplasma_R_coord =   4   ! R coord
  integer, parameter, public :: xplasma_Z_coord =   5   ! Z coord
  integer, parameter, public :: xplasma_rhox_coord =6   ! "rho" extrapolation
  integer, parameter, public :: xplasma_thx_coord = 7   ! "theta" extrapolation
  integer, parameter, public :: xplasma_B_coord =   8   ! mod(B) grid
  integer, parameter, public :: xplasma_vpv_coord = 9   ! vpll/v coord.

  integer, parameter, public :: xplasma_flux_coords = 1 ! generic flux coord.
  integer, parameter, public :: xplasma_cyl_coords = 2  ! generic cyl. coord.

  character*32, parameter, public :: xplasma_root = '__XPLASMA_SYSTEM'
  character*32, parameter, public :: xplasma_xmhd = '__XPLASMA_MHDEQ_DERIVED'

  !----------------------------------------------
  !  Constants (private):
  !  Range of known type codes: xplasma_minType:xplasma_maxType

  integer, parameter :: xplasma_minType = 1  ! lowest known type code
  integer, parameter :: xplasma_maxType = 5  ! highest known type code

  integer, parameter :: maxmsgs = 50

  real*8, parameter :: c2pi = 6.2831853071795862D+00
  real*8, parameter :: cpi = 3.1415926535897931D+00
  real*8, parameter :: ceps4 = 1.0d-4
  real*8, parameter :: ceps7 = 1.0d-7
  real*8, parameter :: ceps10 = 1.0d-10
  real*8, parameter :: czero = 0.0d0
  real*8, parameter :: cone = 1.0d0

  integer, parameter :: ilun_dbg=6  ! private, for debugging only

  !----------------------------------------------
  !  SAVEd global variables

  integer :: internal_counter = 0     ! private internal instance counter

  character*32, public :: xplasma_trace_name = ' '   ! (debug) data trace name

  !
  !  xplasma_trace_name can be set for debugging: to trace the creation
  !  or deletion of named data items within any xplasma object
  !
  !----------------------------------------------
  !  type definitions -- private to this module
  !
  !  (1) supported data types:

  type :: xplist  ! LIST

     ! pointers used: lists are expandable...

     integer :: size = 0       ! length of list

     character*32, dimension(:), pointer :: enames => NULL()
        ! names of list elements

     integer, dimension(:,:), pointer :: chindx => NULL()
        ! character data indices

     real*8, dimension(:), pointer :: r8vals => NULL()
        ! element values (real*8)

     integer, dimension(:), pointer :: ivals => NULL()
        ! element values (integer)

     character, dimension(:), pointer :: chbuf => NULL()
        ! character data buffer

  end type xplist

  type :: xpcoord ! COORDinate

     logical :: periodic = .FALSE.  ! T if it is a periodic "angle" coordinate
     integer :: ngridc = 0     ! no. of grids discretizing this coordinate
     integer :: ngrid_max = 0  ! current max. no. of grids for this coordinate

     integer, pointer, dimension(:) :: grid_ids => NULL()
        ! xplasma ids of grid objects that use this coordinate (expandable)

  end type xpcoord

  type :: xpgrid  ! GRID

     integer :: size = 0     ! no. of grid points

     real*8, dimension(:,:), allocatable :: xpkg
                             ! the grid (x(1)<x(2)<...<x(size))
                             !  (as a "pspline" lookup package) (1:size,1:4)

     integer :: coord = 0                 ! grid coord (R Z rho theta phi ...)
     integer :: nrefs = 0                 ! no. of references to this grid

     !  each grid is a particular discretization of a physical coordinate
     !  (system defined: R Z rho theta phi... or [user defined]...)

     !  coord  -- identifies which coordinate is discretized
     !  nrefs  -- identifies number of references from profile objects

  end type xpgrid

  integer, parameter, private :: maxRank=3  ! MAX rank of profiles
  ! if this ever changes care will be needed to maintain backwards
  ! compatibility with existing XPLASMA NetCDF files.

  integer, parameter, private :: maxApro=1  ! MAX #associated profiles
  ! Ids of associated profiles can be used interchangeably-- E.g. if we
  ! have B(R,Z) and B(rho,theta) the id of either of these can be used
  ! for interpolation of either of these profiles; if this ever changes
  ! care will be needed to maintain backwards compatibility with XPLASMA
  ! NetCDF files.

  integer, parameter, public :: xplasma_prof_maxRank = maxRank

  type :: xpprof   ! PROFILE (spline/hermite/linear interp); max rank = 3

     integer :: rank = 0       ! rank of profile (i.e. #of indep. coordinates)
     integer :: gridIds(maxRank) ! associated grids
     integer :: profIds(maxApro) ! associated profiles
     integer :: kspline = 0      ! interpolation type of profile
     !  kspline=0, piecewise linear; =1, Hermite; =2, cubic spline
     integer :: prof_counter = 0 ! profile creation/update counter (~version no.)
     integer :: size = 0       ! currently allocated buffer size
     real*8, dimension(:), allocatable :: buf      ! data buffer

  end type xpprof

  type :: xpblkbx  ! BLACK BOX -- data for user defined purposes
     ! (xplasma itself uses this for numerical integration data structures)
     !   elements:  a TYPE CODE, user defined types >= 100
     !              an integer vector of some length
     !              a real*8 vector of some length

     integer :: type = 0
     integer, dimension(:), allocatable :: iwork
     real*8, dimension(:), allocatable :: r8work

  end type xpblkbx
  
  !----------------------------------------------
  ! (2)   *** dictionary entry for a single data item ***

  type :: xplasma_item
 
     character*32 :: name = " "
     character*32 :: author = " "
     character, dimension(:), allocatable :: label
     character, dimension(:), allocatable :: units

     integer :: dtype = 0   ! data type referred to by this entry
     integer :: dindex = 0  ! index to data type instance

  end type xplasma_item

  !----------------------------------------------
  !  the following entire type is private to xplasma_obj -- aid in 
  !  evaluating Fourier representation of equilibrium

  type :: xmom2d

     integer :: kmom = 0 ! # of moments 0:kmom
     integer :: ksurf= 0 ! # of rho surfaces
     integer :: kmaxe= 1 ! # max exponent N for /x**N scaling of moments

     real*8, dimension(:,:), allocatable :: xrho   ! rho surfaces (ksurf)
     real*8, dimension(:,:,:), allocatable :: &
          sxrmmc,sxrmms,sxzmmc,sxzmms  ! Fourier coef splines (4,ksurf,0:kmom)

  end type xmom2d

  !----------------------------------------------
  ! (3)   *** container -- public object, private contents ***

  type :: xplasma

     !  ** when this definition changes please update COPY READ and WRITE
     !     methods!!

     private  ! this type definition is public but the contents are private...

     character, dimension(:), allocatable :: label ! for entire xplasma dataset

     !---------------------------
     ! author stack -- only "authors" make changes to xplasma state
     integer :: max_author = 0  ! max no. of authors
     integer :: cur_author = 0 ! current author

     character*32, dimension(:), pointer :: authors => NULL() ! author names
     integer, dimension(:), pointer :: write_enable => NULL()
         ! 1 for TRUE, 0 for FALSE

     !---------------------------
     ! DICTIONARY -- catalog of all data items in xplasma

     ! the number of known items in the dictionary:
     integer :: nitems = 0
     ! current allocated dictionary capacity (expandable):
     integer :: nmax = 0

     ! sort order
     integer, dimension(:), pointer :: order => NULL()
     ! the dictionary entries
     type (xplasma_item), dimension(:), pointer :: dict => NULL()

     !---------------------------
     ! top level data items  (not in the dictionary)

     integer :: nguard = 0           ! initialization indicator (0=NO).

     ! axisymmetry flag:
     logical :: axisymm = .TRUE.     ! T means axisymmetric (tokamak) geometry

     ! scrapeoff flag:
     logical :: scrapeoff = .FALSE.  ! T means scrapeoff region is defined

     real*8 :: time = 0.0d0   ! a time (s) may be associated with the object.

     integer :: bphi_ccw = 0  ! B_phi counter-clockwise (CCW) info
     ! 0: unspecified;  1: CCW;  -1: CW, clockwise

     integer :: jphi_ccw = 0  ! J_phi counter-clockwise (CCW) info
     ! 0: unspecified;  1: CCW;  -1: CW, clockwise

     !---------------------------
     ! xplasma constituent data types
     ! (each instance also has a dictionary entry)
     !
     !  (1) -- lists

     integer :: nlists = 0
     integer :: nlist_free = 0
     type (xplist), dimension(:), allocatable :: lists
     integer, dimension(:), allocatable :: list_freelist

     !  (2) -- coordinates

     integer :: ncoords = 0
     integer :: ncoord_free = 0
     type (xpcoord), dimension(:), allocatable :: coords
     integer, dimension(:), allocatable :: coord_freelist

     !  (3) -- grids

     integer :: ngrids = 0
     integer :: ngrid_free = 0
     type (xpgrid), dimension(:), allocatable :: grids
     integer, dimension(:), allocatable :: grid_freelist
     integer, dimension(:), allocatable :: grid_equiv

     !  (4) -- profiles

     integer :: prof_counter = 0 ! incremented for each profile create/update
     integer :: eq_counter = 0   ! incremented for each equilibrium update
     real*8 :: bdytol = 1.0d-10  ! out-of-bounds tol. for prof. evaluations
     real*8 :: ajac_maxvar = 1.0d-3
               ! max rel. variation in det(J) on a flux surface
     integer :: rzOrder = 2      ! spline order (1=C1 Hermite or 2=C2 spline)
     !  preset for R(rho,theta),Z(rho,theta) equilibrium representation
     integer :: kmom = 16        ! #Fourier moments for spectral representation

     integer :: nprofs = 0
     integer :: nprof_free = 0
     type (xpprof), dimension(:), allocatable :: profs
     integer, dimension(:), allocatable :: prof_freelist

     !  (5) -- Black Boxes -- Opaque objects, store/retrieve/read/write/copy 
     !         ONLY (a type code is recorded; xplasma uses some types for
     !         numerical integration control).

     integer :: nblkbxs = 0
     integer :: nblkbx_free = 0
     type (xpblkbx), dimension(:), allocatable :: blkbxs
     integer, dimension(:), allocatable :: blkbx_freelist

     !-------------------------------------------
     !  private, derived data -- not in NetCDF file
     !  this is for speed & convenience-- pointers to commonly used dictionary
     !  entries...

     integer :: id_g = 0     ! G profile
     integer :: id_psi = 0   ! Psi profile
     integer :: id_P = 0     ! MHD Pressure profile (optional)
     integer :: id_R = 0     ! R(x,theta)
     integer :: id_Rx = 0    ! R(x,theta)-- extrapolated flux coordinate space
     integer :: id_Z = 0     ! Z(x,theta)
     integer :: id_Zx = 0    ! Z(x,theta)-- extrapolated flux coordinate space
     integer :: id_Bmod = 0  ! |B|(x,theta)
     integer :: id_BR = 0    ! BR(x,theta)
     integer :: id_BZ = 0    ! BZ(x,theta)

     integer :: id_axis = 0  ! list: {Raxis,Zaxis}
     integer :: id_RZminmax = 0   ! list: {Rmin,Rmax,Zmin,Zmax of bdy surface}

     !-------------------------------------------

     !  Fourier representation of equilibrium
     !  this is not explicitly saved to / restored from NetCDF file.  Instead
     !  a blackbox storage block is used

     type (xmom2d) :: xm

     !-------------------------------------------
     !  error message buffers
     !   ... if set these are printed and cleared by a call to "xplasma_error".

     integer :: nmsgs = 0 
     character*128 msgs(maxmsgs)

     !  NOTE: error message buffers not saved to / restored from NetCDF file.

  end type xplasma

  !----------------------------------------------
  !  (4)  utility data types (private for this module)
  !

  type :: xpeval

     private

     !  gather together results of xplookup & normalization of x values
     !  for profile function evaluation

     integer :: idx     = 0 ! grid ID
     integer :: inum    = 0 ! number of interpolations sought
     integer :: inx     = 0 ! number of x values in grid (idx)
     integer :: icoord  = 0 ! coordinate discretized by grid (idx)
     logical :: ccwflag = .TRUE. ! CCW option for angle variables (default T)

     integer, dimension(:), pointer :: iderivs => NULL() ! deriv. controls:
     !  one specified for each evaluation: 0->f, 1->df/dx, etc.

     integer, dimension(:), pointer :: ix => NULL()  ! zone lookup results
     real*8, dimension(:), pointer  :: dxn => NULL() ! normalized displacement
     real*8, dimension(:), pointer  :: hx => NULL()  ! zone widths
     real*8, dimension(:), pointer  :: hxi => NULL() ! inverse zone widths

  end type xpeval

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  contains

    !========================================
    !  generic methods for container object
    !========================================

    subroutine xplasma_init(s,zlabel,ier,total)

      !  called once at the beginning of the life of an xplasma
      !  instance...

      type (xplasma), pointer :: s
      character*(*), intent(in) :: zlabel
      integer, intent(out) :: ier  ! completion code (0=OK)

      ! the following argument has an effect if an xplasma_init call is
      ! applied to an object that has already been used...

      logical, intent(in), optional :: total ! .TRUE. to delete list structures
      ! (by default, list elements are emptied but the structures left intact
      ! for performance reasons).

      !--------------------/local
      integer :: i,ilen ! non-blank length of label
      integer :: id
      logical, parameter :: not_periodic = .FALSE.
      logical, parameter :: periodic = .TRUE.
      integer :: luntty = 6
      !--------------------

      ier = 0

      if(.not.associated(s)) then
         ier=101
         return
      endif

      if(s%nguard.ne.0) call xplasma_freesub(s,ier,total)

      internal_counter = internal_counter + 1
      s%nguard=internal_counter

      if(xplasma_trace_name.ne.' ') then
         write(luntty,*) ' %xplasma_obj(xplasma_INIT): object initialization: '
         write(luntty,*) '  label: "'//trim(zlabel)//'"'
         write(luntty,*) '  internal ID: ',s%nguard
      endif

      if(ier.ne.0) return

      s%axisymm=.TRUE.
      s%scrapeoff=.FALSE.
      s%bphi_ccw=0
      s%jphi_ccw=0

      s%time = 0

      s%prof_counter = 0
      s%eq_counter = 0

      s%bdytol = ceps10

!!! resets on this allowed to persist across re-initialization:  
!!!       s%ajac_maxvar = 1.0d-3
!!! (dmc -- support legacy usage)

      s%rzOrder = 2
      s%kmom = 16

      !  set label...

      call xplasma_relabel(s,zlabel)

      ! OK -- start author stack

      call xplasma_author_set(s,xplasma_root,ier)
      if(ier.ne.0) return

      ! OK -- add the standard coordinates...

      call xplasma_create_coord(s,'__rho_coord',not_periodic,id,ier, &
           label='rho',units='-')
      if(ier.ne.0) return
      if(id.ne.xplasma_rho_coord) then
         ier=88
         call xplasma_errmsg_append(s, &
              'xplasma initialization (xplasma_init) code error detected!')
         return
      endif

      call xplasma_create_coord(s,'__theta_coord',periodic,id,ier, &
           label='theta',units='radians')
      if(ier.ne.0) return
      if(id.ne.xplasma_theta_coord) then
         ier=88
         call xplasma_errmsg_append(s, &
              'xplasma initialization (xplasma_init) code error detected!')
         return
      endif

      call xplasma_create_coord(s,'__phi_coord',periodic,id,ier, &
           label='phi',units='radians')
      if(ier.ne.0) return
      if(id.ne.xplasma_phi_coord) then
         ier=88
         call xplasma_errmsg_append(s, &
              'xplasma initialization (xplasma_init) code error detected!')
         return
      endif

      call xplasma_create_coord(s,'__R_coord',not_periodic,id,ier, &
           label='R',units='m')
      if(ier.ne.0) return
      if(id.ne.xplasma_R_coord) then
         ier=88
         call xplasma_errmsg_append(s, &
              'xplasma initialization (xplasma_init) code error detected!')
         return
      endif

      call xplasma_create_coord(s,'__Z_coord',not_periodic,id,ier, &
           label='Z',units='m')
      if(ier.ne.0) return
      if(id.ne.xplasma_Z_coord) then
         ier=88
         call xplasma_errmsg_append(s, &
              'xplasma initialization (xplasma_init) code error detected!')
         return
      endif

      call xplasma_create_coord(s,'__rhox_coord',not_periodic,id,ier, &
           label='rho(extrapolated)',units='-')
      if(ier.ne.0) return
      if(id.ne.xplasma_rhox_coord) then
         ier=88
         call xplasma_errmsg_append(s, &
              'xplasma initialization (xplasma_init) code error detected!')
         return
      endif

      call xplasma_create_coord(s,'__thx_coord',periodic,id,ier, &
           label='theta(extrapolated space)',units='-')
      if(ier.ne.0) return
      if(id.ne.xplasma_thx_coord) then
         ier=88
         call xplasma_errmsg_append(s, &
              'xplasma initialization (xplasma_init) code error detected!')
         return
      endif

      call xplasma_create_coord(s,'__B_coord',not_periodic,id,ier, &
           label='mod(B)',units='T')
      if(ier.ne.0) return
      if(id.ne.xplasma_B_coord) then
         ier=88
         call xplasma_errmsg_append(s, &
              'xplasma initialization (xplasma_init) code error detected!')
         return
      endif

      call xplasma_create_coord(s,'__vpv_coord',not_periodic,id,ier, &
           label='vpll/v',units='-')
      if(ier.ne.0) return
      if(id.ne.xplasma_vpv_coord) then
         ier=88
         call xplasma_errmsg_append(s, &
              'xplasma initialization (xplasma_init) code error detected!')
         return
      endif

      call xplasma_author_pause(s,xplasma_root,ier)

      ! clear private IDs

      s%id_g = 0
      s%id_psi = 0
      s%id_P = 0
      s%id_R = 0
      s%id_Z = 0
      s%id_Rx = 0
      s%id_Zx = 0
      s%id_BMOD= 0
      s%id_BR = 0
      s%id_BZ = 0
      s%id_axis = 0
      s%id_RZminmax = 0

    end subroutine xplasma_init

    !========================================
    subroutine xplasma_relabel(s,zlabel)

      !  set or update the global label of the xplasma object

      type (xplasma), pointer :: s
      character*(*), intent(in) :: zlabel

      !-------------
      integer :: i,ilen
      !-------------

      ilen = len(trim(zlabel))

      !  dmc -- min length now 2 to avoid future NetCDF "scalar" error

      if(allocated(s%label)) deallocate(s%label)

      allocate(s%label(max(2,ilen))); s%label(1:2)='  '
      if(ilen.gt.0) then
         do i=1,ilen
            s%label(i) = zlabel(i:i)
         enddo
      endif
    end subroutine xplasma_relabel

    !========================================
    subroutine xplasma_clear_msgs(s)

      ! private method: clear all error messages

      type (xplasma), pointer :: s

      !----------------
      integer :: iertmp
      !----------------

      if(.not.associated(s)) then
         return
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init via xplasma_clear_msgs call)',iertmp)
      endif

      s%msgs = ' '
      s%nmsgs = 0

    end subroutine xplasma_clear_msgs

    !========================================
    subroutine xplasma_free(s,ier,total)

      !  free memory associated with xplasma object

      type (xplasma), pointer :: s
      integer, intent(out) :: ier

      logical, intent(in), optional :: total ! .TRUE. to delete list structures
      ! (by default, list elements are emptied but the structures left intact
      ! for performance reasons).

      if(.not.associated(s)) then
         ier=101
         return
      endif

      call xplasma_init(s,'(init due to xplasma_free call)',ier,total)

    end subroutine xplasma_free
      
    !========================================
    subroutine xplasma_freesub(s,ier,total)

      !  free memory associated with xplasma object (private implementation)

      type (xplasma), pointer :: s
      integer, intent(out) :: ier

      logical, intent(in), optional :: total ! .TRUE. to delete list structures
      ! (by default, list elements are emptied but the structures left intact
      ! for performance reasons).

      type (xmom2d), pointer :: m

      !--------------------/local
      integer :: i,itchk
      logical :: itotal
      integer :: luntty = 6
      character*32 :: znam1
      !--------------------

      itotal = .FALSE.
      if(present(total)) itotal = total

      ier=0

      if(xplasma_trace_name.ne.' ') then
         znam1=xplasma_trace_name
         call uupper(znam1)
         itchk=1
      else
         znam1=' '
         itchk=0
      endif

      if(allocated(s%label)) deallocate(s%label)

      !  remove author list

      if(associated(s%authors)) deallocate(s%authors)
      if(associated(s%write_enable)) deallocate(s%write_enable)
      s%max_author=0
      s%cur_author=0

      !  remove dictionary

      do i=1,s%nmax
         if(allocated(s%dict(i)%label)) deallocate(s%dict(i)%label)
         if(allocated(s%dict(i)%units)) deallocate(s%dict(i)%units)
         if((itchk.eq.1).and.(s%dict(i)%dindex.gt.0)) then
            if(znam1.eq.s%dict(i)%name) then
               call tchk(s,'deleted in xplasma_freesub',znam1)
            endif
         endif
      enddo
      if(allocated(s%label)) deallocate(s%label)
      if(associated(s%order)) deallocate(s%order)
      if(associated(s%dict)) deallocate(s%dict)
      s%nitems=0
      s%nmax=0

      !  clean up constituent data item lists...
      !  all slots are re-initialized empty but not deleted unless TOTAL is set

      if(allocated(s%lists)) then
         do i=1,s%nlists
            call xplist_free(s%lists(i))
         enddo
         if(itotal) then
            deallocate(s%list_freelist,s%lists)
            s%nlists=0
            s%nlist_free=0
         else
            s%nlist_free = s%nlists
            do i=1,s%nlists
               s%list_freelist(i)=s%nlists + 1 - i
            enddo
         endif
      endif

      !  clean up constituent data item coords...

      if(allocated(s%coords)) then
         do i=1,s%ncoords
            call xpcoord_free(s%coords(i))
         enddo
         if(itotal) then
            deallocate(s%coord_freelist,s%coords)
            s%ncoords=0
            s%ncoord_free=0
         else
            s%ncoord_free = s%ncoords
            do i=1,s%ncoords
               s%coord_freelist(i)=s%ncoords + 1 - i
            enddo
         endif
      endif

      !  clean up constituent data grids...

      if(allocated(s%grids)) then
         do i=1,s%ngrids
            call xpgrid_free(s%grids(i))
         enddo
         if(itotal) then
            deallocate(s%grid_freelist,s%grids,s%grid_equiv)
            s%ngrids=0
            s%ngrid_free=0
         else
            s%ngrid_free = s%ngrids
            do i=1,s%ngrids
               s%grid_freelist(i)=s%ngrids + 1 - i
               s%grid_equiv(i)=0
            enddo
         endif
      endif

      !  clean up constituent profiles...

      if(allocated(s%profs)) then
         do i=1,s%nprofs
            call xpprof_free(s%profs(i))
         enddo
         if(itotal) then
            deallocate(s%prof_freelist,s%profs)
            s%nprofs=0
            s%nprof_free=0
         else
            s%nprof_free = s%nprofs
            do i=1,s%nprofs
               s%prof_freelist(i)=s%nprofs + 1 - i
            enddo
         endif
      endif

      !  clean up constituent BLACK BOX datasets...

      if(allocated(s%blkbxs)) then
         do i=1,s%nblkbxs
            call xpblkbx_free(s%blkbxs(i))
         enddo
         if(itotal) then
            deallocate(s%blkbx_freelist,s%blkbxs)
            s%nblkbxs=0
            s%nblkbx_free=0
         else
            s%nblkbx_free = s%nblkbxs
            do i=1,s%nblkbxs
               s%blkbx_freelist(i)=s%nblkbxs + 1 - i
            enddo
         endif
      endif

      !  miscellaneous...
      
      m => s%xm
      if(allocated(m%xrho)) call xpmom_free(s)

    end subroutine xplasma_freesub

    !========================================
    subroutine xplasma_author_set(s,author,ier)
      
      !  push new author onto stack -- e.g. before writing some xplasma
      !  output...

      type (xplasma), pointer :: s
      character*(*), intent(in) :: author
      integer, intent(out) :: ier

      !------------------
      integer :: innew
      character*40 zauthor
      !------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_author_set call)',ier)
      endif
      if(ier.ne.0) return

      call xplasma_checkName(author,ier)
      if(ier.ne.0) then
         ier=54
         call xplasma_errmsg_append(s, &
              '  passed author name string was: "'//trim(author)//'".')
         return
      endif

      zauthor = author
      call uupper(zauthor)

      ! check for dynamic stack growth...

      innew = s%cur_author + 1
      call xplasma_authorStack_expand(s,innew)

      s%authors(innew)=zauthor
      s%write_enable(innew)=1
      s%cur_author = innew

    end subroutine xplasma_author_set

    !========================================
    subroutine xplasma_author_pause(s,author,ier)
      
      !  disable write access for current author
      !  e.g. as precaution prior to passing control to 
      !  different author in called subprocess...

      type (xplasma), pointer :: s
      character*(*), intent(in) :: author
      integer, intent(out) :: ier

      !------------------
      character*32 zauthor
      !------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_author_pause call)',ier)
      endif
      if(ier.ne.0) return

      zauthor = author
      call uupper(zauthor)

      if(s%authors(s%cur_author).ne.zauthor) then
         ier=52
         call xplasma_errmsg_append(s, &
              ' xplasma_author_pause call failed.')
         call xplasma_errmsg_append(s, &
              ' passed author name was: '//trim(zauthor))
         call xplasma_errmsg_append(s, &
              ' current author for write access: '// &
              trim(s%authors(s%cur_author)))
         return
      else if (s%cur_author.lt.1) then
         ier=53
         return
      endif

      s%write_enable(s%cur_author)=0

    end subroutine xplasma_author_pause

    !========================================
    subroutine xplasma_author_resume(s,author,ier)
      
      !  re-enable write access for current author
      !  e.g. after after control is returned from subprocess...

      type (xplasma), pointer :: s
      character*(*), intent(in) :: author
      integer, intent(out) :: ier

      !------------------
      character*32 zauthor
      !------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_author_resume call)',ier)
      endif
      if(ier.ne.0) return

      zauthor = author
      call uupper(zauthor)

      if(s%authors(s%cur_author).ne.zauthor) then
         ier=52
         call xplasma_errmsg_append(s, &
              ' xplasma_author_resume call failed.')
         call xplasma_errmsg_append(s, &
              ' passed author name was: '//trim(zauthor))
         call xplasma_errmsg_append(s, &
              ' current author for write access: '// &
              trim(s%authors(s%cur_author)))
         return
      else if (s%cur_author.le.1) then
         ier=53
         return
      endif

      s%write_enable(s%cur_author)=1

    end subroutine xplasma_author_resume

    !========================================
    subroutine xplasma_author_clear(s,author,ier)
      
      !  pop author off stack -- return to prior author
      !  passed name must match 

      type (xplasma), pointer :: s
      character*(*), intent(in) :: author
      integer, intent(out) :: ier

      !------------------
      character*32 zauthor
      !------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_author_clear call)',ier)
      endif
      if(ier.ne.0) return

      zauthor = author
      call uupper(zauthor)

      if(s%authors(s%cur_author).ne.zauthor) then
         ier=52
         call xplasma_errmsg_append(s, &
              ' xplasma_author_clear call failed.')
         call xplasma_errmsg_append(s, &
              ' passed author name was: '//trim(zauthor))
         call xplasma_errmsg_append(s, &
              ' current author for write access: '// &
              trim(s%authors(s%cur_author)))
         return
      else if (s%cur_author.le.1) then
         ier=53
         return
      endif

      s%authors(s%cur_author)=' '
      s%write_enable(s%cur_author)=0

      s%cur_author = s%cur_author - 1

    end subroutine xplasma_author_clear

    !========================================
    subroutine xplasma_global_info(s,ier, init_check, &
         initLabel,time,axisymm,scrapeoff,bphi_ccw,jphi_ccw,nitems, &
         prof_counter,eq_counter,bdytol,ajac_maxVaR,rzOrder,kmom)

      ! this routine returns global xplasma object information, according
      ! to which optional arguments are provided...

      type (xplasma), pointer :: s
      integer, intent(out) :: ier   ! error set iff s is uninitialized.

      logical, intent(in), optional :: init_check   ! .TRUE. to report error
      ! if xplasma object is not marked as initialized (default .FALSE.)

      ! this is the label passed in an xplasma_(re)Init call...
      character*(*), intent(out), optional :: initLabel

      ! can also use xplasma_time_get for this:
      real*8, intent(out), optional :: time

      logical, intent(out), optional :: axisymm   ! axisymmetry flag

      logical, intent(out), optional :: scrapeoff ! SOL flag

      integer, intent(out), optional :: bphi_ccw  ! Bphi counter-clockwise flag
      ! 0 means: never specified; 1: counter-clockwise (CCW); -1: CW

      integer, intent(out), optional :: jphi_ccw  ! Jphi counter-clockwise flag
      ! 0 means: never specified; 1: counter-clockwise (CCW); -1: CW

      integer, intent(out), optional :: nitems    ! number of items in xplasma
      ! dictionary at time of call.

      integer, intent(out), optional :: prof_counter  ! profile counter
      ! global counter incremented for each creation or update of each profile

      integer, intent(out), optional :: eq_counter  ! profile counter
      ! global counter incremented for each update of MHD equilibrium

      real*8, intent(out), optional :: bdytol   ! out-of-range error tolerance
      ! relative tolerance applied for interpolations over non-periodic 
      ! coordinates; for periodic coordinates shifts of 2*n*pi are applied
      ! as necessary to bring interpolations in range

      real*8, intent(out), optional :: ajac_maxVar   ! det[J] variation tol.
      ! relative tolerance for variations of det[J] within a flux surface;
      ! used to check for singular or near-singular metric tensor in 
      ! R(rho,theta),Z(rho,theta) flux surface definitions.

      integer, intent(out), optional :: rzOrder ! spline fit order for MHD
      ! equilibrium geometry R(rho,theta) z(rho,theta) data
      ! 1 for C1 Hermite; 2 for C2 Spline.

      integer, intent(out), optional :: kmom    ! #moments for Fourier
      ! representation of equilibrium

      !--------------------/local
      integer :: isize,i
      logical :: icheck
      !--------------------

      icheck=.FALSE.
      if(present(init_check)) icheck = init_check

      if(present(initLabel)) initLabel=' '
      if(present(time)) time=0
      if(present(axisymm)) axisymm=.FALSE.
      if(present(scrapeoff)) scrapeoff=.FALSE.
      if(present(bphi_ccw)) bphi_ccw=0
      if(present(jphi_ccw)) jphi_ccw=0
      if(present(prof_counter)) prof_counter=0
      if(present(eq_counter)) eq_counter=0
      if(present(bdytol)) bdytol=0
      if(present(ajac_maxVar)) ajac_maxVar=0
      if(present(rzOrder)) rzOrder=0
      if(present(kmom)) kmom=0
      if(present(nitems)) nitems = 0

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         if(icheck) then
            ier=101
         else
            call xplasma_init(s,'(init due to xplasma_global_info call)',ier)
         endif
      endif
      if(ier.ne.0) return

      if(present(initLabel)) then
         isize=min(size(s%label),len(initLabel))
         do i=1,isize
            initLabel(i:i)=s%label(i)
         enddo
      endif

      if(present(time)) time = s%time

      if(present(axisymm)) axisymm = s%axisymm

      if(present(scrapeoff)) scrapeoff = s%scrapeoff

      if(present(bphi_ccw)) bphi_ccw = s%bphi_ccw
      if(present(jphi_ccw)) jphi_ccw = s%jphi_ccw

      if(present(nitems)) nitems = s%nitems

      if(present(prof_counter)) prof_counter = s%prof_counter
      if(present(eq_counter)) eq_counter = s%eq_counter

      if(present(bdytol)) bdytol = s%bdytol

      if(present(ajac_maxVar)) ajac_maxVar = s%ajac_maxVar

      if(present(rzOrder)) rzOrder = s%rzOrder

      if(present(kmom)) kmom = s%kmom

    end subroutine xplasma_global_info

    !========================================
    subroutine xplasma_common_ids(s,ier, &
         id_g,id_psi,id_P,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
         id_axis,id_RZminmax)

      !  return commonly requested profile IDs
      !  g(rho), psi(rho), R(rho,theta), Z(rho,theta)

      type (xplasma), pointer :: s
      integer, intent(out) :: ier

      integer, intent(out), optional :: id_g
      integer, intent(out), optional :: id_psi
      integer, intent(out), optional :: id_P
      integer, intent(out), optional :: id_R
      integer, intent(out), optional :: id_Z
      integer, intent(out), optional :: id_BR
      integer, intent(out), optional :: id_BZ
      integer, intent(out), optional :: id_BMOD


      integer, intent(out), optional :: id_axis
      integer, intent(out), optional :: id_RZminmax

      ! return common IDs for standard profiles or 0 if the profile
      ! is not currently available

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_common_ids call)',ier)
      endif
      if(ier.ne.0) return

      if(present(id_g)) id_g = s%id_g
      if(present(id_psi)) id_psi = s%id_psi
      if(present(id_P)) id_P = s%id_P

      if(present(id_R)) id_R = s%id_R
      if(present(id_Z)) id_Z = s%id_Z

      if(present(id_BR)) id_BR = s%id_BR
      if(present(id_BZ)) id_BZ = s%id_BZ
      if(present(id_Bmod)) id_Bmod = s%id_Bmod

      if(present(id_axis)) id_axis=s%id_axis
      if(present(id_RZminmax)) id_RZminmax=s%id_RZminmax

    end subroutine xplasma_common_ids

    !========================================
    subroutine xplasma_set_pressure_id(s,id_pressure,ier)

      !  set the id for pressure

      type (xplasma), pointer :: s
      integer, intent(in) :: id_pressure
      integer, intent(out) :: ier  ! completion code,  0=OK

      !--------------------------------
      integer :: irank,igrid,icoord
      !--------------------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_set_pressure_id call)',ier)
      endif
      if(ier.ne.0) return

      call xplasma_internal_ids(s)
      if(s%id_g.eq.0) then
         ier=9999
         call xplasma_errmsg_append(s, &
              ' ?xplasma_set_pressure_id: no "G" profile.')
         call xplasma_errmsg_append(s, &
              '  Need to acquire MHD equilibrium before calling routine.')
         return
      endif

      call xplasma_prof_info(s,id_pressure,ier, rank=irank,gridId1=igrid)

      if(ier.eq.0) then
         if(irank.ne.1) then
            ier=9999
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_set_pressure_id: only 1d f(rho) profile can be used; rank > 1.')
         else
            call xplasma_grid_info(s,igrid,ier, coord=icoord)
            if(icoord.ne.xplasma_rho_coord) then
               if(ier.ne.0) ier=9999
               call xplasma_errmsg_append(s, &
                    ' ?xplasma_set_pressure_id: only 1d f(rho) profile can be used; x coordinate not "rho".')
            endif
         endif
      endif

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_set_pressure_id: invalid pressure profile id.')
         return
      endif

      !  OK

      s%id_P = id_pressure

    end subroutine xplasma_set_pressure_id

    !========================================
    subroutine xplasma_time_set(s,ztime,ier)

      !  set the time (seconds) affiliated with xplasma object "s".

      type (xplasma), pointer :: s
      real*8, intent(in) :: ztime
      integer, intent(out) :: ier

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_time_set call)',ier)
      endif
      if(ier.ne.0) return
      
      s%time = ztime

    end subroutine xplasma_time_set

    !========================================
    subroutine xplasma_time_get(s,ztime,ier)

      !  fetch the time (seconds) affiliated with xplasma object "s".

      type (xplasma), pointer :: s
      real*8, intent(out) :: ztime
      integer, intent(out) :: ier

      ztime = 0

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_time_get call)',ier)
      endif
      if(ier.ne.0) return

      ztime = s%time

    end subroutine xplasma_time_get

    !========================================
    subroutine xplasma_bdytol_set(s,zbtol,ier)

      !  set the out-of-bounds tolerance for profile interpolations
      !  a maximum value of 1.0d-4 is allowed...

      type (xplasma), pointer :: s
      real*8, intent(in) :: zbtol
      integer, intent(out) :: ier

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_bdytol_set call)',ier)
      endif
      if(ier.ne.0) return

      s%bdytol = min(1.0d-4,zbtol)

    end subroutine xplasma_bdytol_set

    !========================================
    subroutine xplasma_kmom_set(s,kmom,ier)

      !  set the # of Fourier moments for the optional Fourier representation
      !  of the current equilibrium.  Default value is 16 not counting the
      !  0'th moment.  Minimum value is 2, max is 64.

      type (xplasma), pointer :: s
      integer, intent(in) :: kmom
      integer, intent(out) :: ier

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_kmom_set call)',ier)
      endif
      if(ier.ne.0) return

      s%kmom = max(2,min(64,kmom))

    end subroutine xplasma_kmom_set

    !========================================
    subroutine xplasma_ajac_maxvar_set(s,zvar,ier)

      !  set the amount of variation of det[J] to be allowed on a flux
      !  surface-- numerical control for detection of singular equilibria

      !  det J is det  | dR/drho   dR/dtheta |
      !                | dZ/drho   dZ/dtehta |

      type (xplasma), pointer :: s
      real*8, intent(in) :: zvar
      integer, intent(out) :: ier

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_ajac_maxvar_set call)',ier)
      endif
      if(ier.ne.0) return

      if(zvar.gt.czero) then
         s%ajac_maxVar = min(0.1d0,max(1.0d-7,zvar))
      else
         s%ajac_maxVar = 1.0d-3
      endif

    end subroutine xplasma_ajac_maxvar_set

    subroutine xplasma_misc_set(s,ier, scrapeoff)

      !  SET miscellaneous flags inside xplasma object.
      !  this routine should not be called from user applications code
      
      !  out the moment only the SOL flag is accessible.

      type (xplasma), pointer :: s
      integer, intent(out) :: ier

      logical, intent(in), optional :: scrapeoff

      !-----------------
      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         ier=101
      endif
      if(ier.ne.0) return

      if(present(scrapeoff)) then
         s%scrapeoff = scrapeoff
      endif

    end subroutine xplasma_misc_set

    !========================================
    subroutine xplasma_eqCount_incr(s,ier)

      !  increment the equilibrium counter
      !  (public within the xplasma library, but probably should not
      !  be directly accessed by user).

      type (xplasma), pointer :: s
      integer, intent(out) :: ier
 
      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         ier=101
      endif
      if(ier.ne.0) return

      s%eq_counter = s%eq_counter + 1

    end subroutine xplasma_eqCount_incr

    !========================================
    subroutine xplasma_field_ccw_set(s,ier, bphi_ccw,jphi_ccw)

      !  set the Bphi and/or Jphi orientations (for tokamak fields)
      !  CCW -- counter-clockwise
      !  1 means (B or J is pointed) CCW looking down on machine from above
      ! -1 means CW, clockwise orientation.

      type (xplasma), pointer :: s
      integer, intent(out) :: ier   ! status code, 0=OK

      !  for the toroidal magnetic field:
      integer, intent(in), optional :: bphi_ccw  ! CCW value for Bphi

      !  for the toroidal plasma current:
      integer, intent(in), optional :: jphi_ccw  ! CCW value for Jphi

      !------------------------
      integer :: icount
      !------------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_field_ccw_set call)',ier)
      endif
      if(ier.ne.0) return

      if(present(bphi_ccw)) then
         if(abs(bphi_ccw).ne.1) then
            call xplasma_errmsg_append(s, &
                 ' xplasma_field_ccw_set:  +1 or -1 required for bphi_ccw.')
            ier=60
         else
            s%bphi_ccw = bphi_ccw
         endif
      endif

      if(present(jphi_ccw)) then
         if(abs(jphi_ccw).ne.1) then
            call xplasma_errmsg_append(s, &
                 ' xplasma_field_ccw_set:  +1 or -1 required for jphi_ccw.')
            ier=60
         else
            s%jphi_ccw = jphi_ccw
         endif
      endif

    end subroutine xplasma_field_ccw_set

    !========================================
    subroutine xplasma_rzOrder_set(s,irzOrder,ier)

      !  set the default fit order for R(rho,theta),Z(rho,theta)
      !  equilibrium specification
      !    min value:  1 -- C1 Hermite fit will be used
      !    max value:  2 -- C2 Spline fit will be used (default)

      type (xplasma), pointer :: s
      integer, intent(in) :: irzOrder
      integer, intent(out) :: ier

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_rzOrder_set call)',ier)
      endif
      if(ier.ne.0) return

      s%rzOrder = max(1,min(2,irzOrder))

    end subroutine xplasma_rzOrder_set

    subroutine xplasma_get_counter(s,icount,ier)

      !  get profile counter value-- current global value for this
      !  xplasma object.
      !  ... the next profile added or replaced will receive this
      !  value + 1 as its profile counter value.

      type (xplasma), pointer :: s
      integer, intent(out) :: icount ! profile counter for this profile
      integer, intent(out) :: ier   ! completion code, 0=OK

      !------------------------------------------------

      icount = 0

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_get_counter call)',ier)
      endif
      if(ier.ne.0) return

      icount = s%prof_counter

    end subroutine xplasma_get_counter

    !========================================
    subroutine xplasma_copy_item(s,id_old,zname,id_new,ier)

      !  make a copy of an item-- this is a full memory copy, not
      !  just a pointer copy

      type (xplasma), pointer :: s
      integer, intent(in) :: id_old       ! id of old item to copy from
      character*(*), intent(in) :: zname  ! name of new item to copy to
      integer, intent(out) :: id_new      ! id of new item created by copy
      integer, intent(out) :: ier         ! status code (0=OK)

      !---------------------------
      integer :: itype,i1,i2
      !---------------------------

      id_new=0

      call xplasma_get_item_info(s,id_old,ier,itype=itype)
      if(ier.ne.0) return

      call xplasma_add_item(s,zname,itype,id_new,ier)
      if(ier.ne.0) return

      i1=s%dict(id_old)%dindex
      i2=s%dict(id_new)%dindex

      if(itype.eq.xplasma_listType) then
         call xplist_fullcopy(s%lists(i1),s%lists(i2))

      else if(itype.eq.xplasma_coordType) then
         call xpcoord_fullcopy(s%coords(i1),s%coords(i2))

      else if(itype.eq.xplasma_gridType) then
         call xpgrid_fullcopy(s%grids(i1),s%grids(i2))

      else if(itype.eq.xplasma_profType) then
         call xpprof_fullcopy(s%profs(i1),s%profs(i2))

      else if(itype.eq.xplasma_blkbxType) then
         call xpblkbx_fullcopy(s%blkbxs(i1),s%blkbxs(i2))

      endif

    end subroutine xplasma_copy_item

    !========================================
    !  private now...
    subroutine xplasma_add_item(s,zname,itype,id,ier)

      ! add an (initialized, but empty) item to the xplasma dictionary
      ! item has name and type assigned here; use xplasma_add_label to
      ! add a longer descriptive label.

      type (xplasma), pointer :: s
      character*(*), intent(in) :: zname  ! name of new item
      integer, intent(in) :: itype       ! type of new item
      integer, intent(out) :: id         ! id of new item (0 if error)
      integer, intent(out) :: ier        ! status code returned

      ! on exit, ier>0 means there was an error
      ! ier=0 means a normal exit

      !--------------------/local
      integer :: i,inp,iloc,iexact,iertmp,jtype,idum,indx
      logical :: have_slot
      character*128 msgbuf
      !--------------------

      id=0

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         ier=101
      endif
      if(ier.ne.0) return

      if(s%write_enable(s%cur_author).ne.1) then
         ier=50
         return
      endif

      if(itype.ne.xplasma_profType) then
         call xplasma_reserve_check(s,itype,0,zname,idum,ier)
         if(ier.ne.0) return
      endif

      call xplasma_checkName(zname,ier)
      if(ier.ne.0) return

      call xplasma_checkType(itype,ier)
      if(ier.ne.0) return

      ! see if item name is already in use and also where the item belongs
      ! in alphabetic order...

      call xplasma_locate(s,zname,iloc,iexact)

      if(iexact.eq.1) then
         id=s%order(iloc)

         jtype=s%dict(id)%dtype
         if(jtype.gt.0) then
            if((s%dict(id)%author.ne.s%authors(s%cur_author)).and. &
                 (s%authors(s%cur_author).ne.xplasma_root)) then
               ier=51
               write(msgbuf,*) ' item name: "',trim(s%dict(id)%name), &
                    '"; item author: ',trim(s%dict(id)%author)
               call xplasma_errmsg_append(s,msgbuf)
               write(msgbuf,*) ' name of author attempting modification: ', &
                    trim(s%authors(s%cur_author))
               call xplasma_errmsg_append(s,msgbuf)
               return
            endif
         endif

         ! delete contents of old item, if any...
         if(jtype.gt.0) then

            if(jtype.eq.itype) then

               !  item being replaced by another item of the same type;
               !  label is retained; free memory associated with old item

               indx = s%dict(id)%dindex
               call xplasma_remove_item(s, id, ier)

               !  s%dict(id)%dindex pointed to an initialized but empty item...
               have_slot = .TRUE.

               !  pop off free list
               if(itype.eq.xplasma_listType) then
                  s%list_freelist(s%nlist_free)=0
                  s%nlist_free = s%nlist_free-1

               else if(itype.eq.xplasma_coordType) then
                  s%coord_freelist(s%ncoord_free)=0
                  s%ncoord_free = s%ncoord_free-1

               else if(itype.eq.xplasma_gridType) then
                  s%grid_freelist(s%ngrid_free)=0
                  s%ngrid_free = s%ngrid_free-1

               else if(itype.eq.xplasma_profType) then
                  s%prof_freelist(s%nprof_free)=0
                  s%nprof_free = s%nprof_free-1

               else if(itype.eq.xplasma_blkbxType) then
                  s%blkbx_freelist(s%nblkbx_free)=0
                  s%nblkbx_free = s%nblkbx_free-1

               endif
            else

               !  replacement item not of the same type...

               call xplasma_remove_item(s, id, ier)

               have_slot = .FALSE.

            endif

         else

            have_slot = .FALSE.  ! referenced a deleted item...

         endif

      else

         have_slot = .FALSE.

         if(s%nitems.eq.s%nmax) then
            !  make room for new item...

            call xplasma_expand(s,s%nitems+1,ier)
            if(ier.ne.0) return
         endif

         inp=s%nitems
         s%nitems = s%nitems + 1
         id = s%nitems

         !  maintain alphabetic ordering...

         do i=inp,iloc,-1
            s%order(i+1)=s%order(i)
         enddo
         s%order(iloc)=id

      endif

      s%dict(id)%dtype = itype
      s%dict(id)%author = s%authors(s%cur_author)

      if(have_slot) then
         s%dict(id)%dindex = indx
         ! name is OK -- was matched.

      else

         s%dict(id)%name = zname
         call uupper(s%dict(id)%name)

         if(xplasma_trace_name.ne.' ') then
            call tchk(s,'created in xplasma_add_item',s%dict(id)%name)
         endif

         ! initialize an empty list item; save its address

         if(itype.eq.xplasma_listType) then
            ! new list item...
            if(s%nlist_free.eq.0) then
               call xplasma_lists_expand(s,s%nlists+1,iertmp)
            endif
            s%dict(id)%dindex = s%list_freelist(s%nlist_free)
            s%list_freelist(s%nlist_free)=0
            s%nlist_free = s%nlist_free - 1 
         endif

         ! initialize an empty coordinate item; save its address

         if(itype.eq.xplasma_coordType) then
            ! new coord item...
            if(s%ncoord_free.eq.0) then
               call xplasma_coords_expand(s,s%ncoords+1,iertmp)
            endif
            s%dict(id)%dindex = s%coord_freelist(s%ncoord_free)
            s%coord_freelist(s%ncoord_free)=0
            s%ncoord_free = s%ncoord_free - 1 
         endif

         ! initialize an empty grid item; save its address

         if(itype.eq.xplasma_gridType) then
            ! new grid item...
            if(s%ngrid_free.eq.0) then
               call xplasma_grids_expand(s,s%ngrids+1,iertmp)
            endif
            s%dict(id)%dindex = s%grid_freelist(s%ngrid_free)
            s%grid_freelist(s%ngrid_free)=0
            s%ngrid_free = s%ngrid_free - 1 
         endif

         ! initialize an empty prof item; save its address

         if(itype.eq.xplasma_profType) then
            ! new prof item...
            if(s%nprof_free.eq.0) then
               call xplasma_profs_expand(s,s%nprofs+1,iertmp)
            endif
            s%dict(id)%dindex = s%prof_freelist(s%nprof_free)
            s%prof_freelist(s%nprof_free)=0
            s%nprof_free = s%nprof_free - 1 
         endif

         ! initialize an empty integrator dataset; save its address

         if(itype.eq.xplasma_blkbxType) then
            ! new blkbx item...
            if(s%nblkbx_free.eq.0) then
               call xplasma_blkbxs_expand(s,s%nblkbxs+1,iertmp)
            endif
            s%dict(id)%dindex = s%blkbx_freelist(s%nblkbx_free)
            s%blkbx_freelist(s%nblkbx_free)=0
            s%nblkbx_free = s%nblkbx_free - 1 
         endif

      endif

      if(id.le.0) then
         ier=9999
         call xplasma_errmsg_append(s, &
              ' ?? unexpected xplasma error -- code error detected!')
         return  ! this should not happen!
      endif

    end subroutine xplasma_add_item

    subroutine xplasma_reserve_check(s,itype,irank,zname,iflag,ier)
      ! private routine, check for "reserved" names:
      !
      !   f(rho):  G, Psi
      !   f(rho,theta):  R, Z

      type (xplasma), pointer :: s
      integer, intent(in) :: itype  ! type of item
      integer, intent(in) :: irank  ! 0 for non-profile, 1 for 1d profile
      !   2 for 2d profile, 3 for 3d profile

      character*(*), intent(in) :: zname  ! the profile or object name
      integer, intent(out) :: iflag       ! =1 or 2 if a reserved name
      !  =0 if not; =1 if reserved for MHD equilibrium i.e. {G,PSI,R,Z}.
      integer, intent(out) :: ier         ! =0 if OK; =89 if reserved name
      !                                               AND of wrong type...

      !-------------------
      character*32 ztest
      !-------------------

      ztest = zname
      call uupper(ztest)

      iflag=0
      ier=0

      if((ztest.eq.'R').or.(ztest.eq.'Z').or. &
           (ztest.eq.'BR').or.(ztest.eq.'BZ').or.(ztest.eq.'BMOD').or. &
           (ztest.eq.'__R_EXTRAP').or.(ztest.eq.'__Z_EXTRAP')) then
         if((itype.eq.xplasma_profType).and.(irank.eq.2)) then
            if((ztest.eq.'R').or.(ztest.eq.'Z')) then
               iflag=1
            else
               iflag=2  ! not R or Z
            endif
         else
            call xplasma_error_reserve(s,xplasma_profType,2,ztest,ier)
            ! this sets ier=89
         endif

      else if((ztest.eq.'G').or.(ztest.eq.'PSI')) then
         if((itype.eq.xplasma_profType).and.(irank.eq.1)) then
            if((ztest.eq.'G').or.(ztest.eq.'PSI')) then
               iflag=1
            else
               iflag=2  ! not G or PSI
            endif
         else
            call xplasma_error_reserve(s,xplasma_profType,1,ztest,ier)  
            ! this sets ier=89
         endif

      else if((ztest.eq.'MAG_AXIS').or.(ztest.eq.'BDY_RZ_MINMAX')) then
         if(itype.eq.xplasma_listType) then
            iflag=2
         else
            call xplasma_error_reserve(s,xplasma_listType,0,ztest,ier)  
            ! this sets ier=89
         endif
      endif
    end subroutine xplasma_reserve_check

    subroutine xplasma_error_reserve(s,itype,irank,zname,ier)

      !  set error code and append message: reserved name error
      !    R(rho,theta) Z(rho,theta)
      !    G(rho), Psi(rho)  --------------- all are reserved profiles.
      !    MAG_AXIS  BDY_RZ_MINMAX

      type (xplasma), pointer :: s
      integer, intent(in) :: itype,irank
      character*(*), intent(in) :: zname
      integer, intent(out) :: ier


      ier=89

      call xplasma_errmsg_append(s, &
           '  The name "'//trim(zname)//'" is reserved and must be a')

      if((itype.eq.xplasma_profType).and.(irank.eq.1)) then

         call xplasma_errmsg_append(s, &
              '  profile of the form f(rho).')

      else if((itype.eq.xplasma_profType).and.(irank.eq.2)) then

         call xplasma_errmsg_append(s, &
              '  profile of the form f(rho,theta).')

      else if(itype.eq.xplasma_listType) then

         call xplasma_errmsg_append(s, &
              '  a list of data items.')

      endif

      call xplasma_errmsg_append(s, &
           '  ...reserved profile names are: "R" "Z" "G" "PSI" "BR" "BZ" "BMOD".')
      call xplasma_errmsg_append(s, &
           '  ...reserved list names are:  "MAG_AXIS" "BDY_RZ_MINMAX"')

    end subroutine xplasma_error_reserve

    !========================================
    subroutine xplasma_label_item(s,id,ier, &
         label,units)

      !  label an existing item as identified by "id"...

      type (xplasma), pointer :: s
      integer, intent(in) :: id          ! id of existing item
      integer, intent(out) :: ier        ! status code returned

      character*(*), intent(in), optional :: label ! label for item
      character*(*), intent(in), optional :: units ! units label for item

      !--------------------/local
      integer :: i,ilen ! non-blank length of label
      character*32 zunits
      integer :: iertmp
      !--------------------

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      call xplasma_author_check(s,id,ier)
      if(ier.ne.0) return

      if(present(label)) then
         if(allocated(s%dict(id)%label)) deallocate(s%dict(id)%label)

         ilen = len(trim(label))
         allocate(s%dict(id)%label(max(1,ilen)))
         if(ilen.gt.0) then
            do i=1,ilen
               s%dict(id)%label(i) = label(i:i)
            enddo
         else
            s%dict(id)%label = ' '
         endif
      endif

      if(present(units)) then
         if(s%dict(id)%dtype.eq.xplasma_gridType) then
            ier= 309
            i=s%dict(id)%dindex
            i=s%grids(i)%coord
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_label_item: cannot change units label for grid')
            call xplasma_errmsg_append(s, &
                 '  object-- units determined by coordinate associated')
            call xplasma_errmsg_append(s, &
                 '  with the grid: '//trim(s%dict(id)%name)//' -> '// &
                 trim(s%dict(i)%name))
            call xplasma_get_item_info(s,i,iertmp, units=zunits)
            call xplasma_errmsg_append(s, &
                 '  units string: "'//trim(zunits)//'".')
         else

            if(allocated(s%dict(id)%units)) deallocate(s%dict(id)%units)

            ilen = len(trim(units))
            allocate(s%dict(id)%units(max(1,ilen)))
            if(ilen.gt.0) then
               do i=1,ilen
                  s%dict(id)%units(i) = units(i:i)
               enddo
            else
               s%dict(id)%units = ' '
            endif

         endif
      endif

    end subroutine xplasma_label_item

    !========================================
    subroutine xplasma_reset_author(s,id,old_author,new_author,ier)

      !  change the author ID of an xplasma item.

      type (xplasma), pointer :: s
      integer, intent(in) :: id          ! id of existing item

      character*(*), intent(in) :: old_author   ! old author id: must match!
      character*(*), intent(in) :: new_author   ! new author id: syntax checked

      integer, intent(out) :: ier        ! status code returned

      !--------------------
      character*40 zauth
      !--------------------

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      zauth = old_author
      call uupper(zauth)

      if(zauth.ne.s%dict(id)%author) then
         ier=52
         call xplasma_errmsg_append(s, &
              ' item '//trim(s%dict(id)%name)//' author id: '//trim(s%dict(id)%author))
         call xplasma_errmsg_append(s, &
              ' xplasma_reset_author: old_author name received: '//trim(zauth))
         return
      endif

      zauth = new_author
      call uupper(zauth)
      call xplasma_checkName(zauth,ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' xplasma_reset_author: error in new_author name: '//trim(zauth))
         return
      endif

      ! OK

      s%dict(id)%author = zauth

    end subroutine xplasma_reset_author

    !========================================
    subroutine xplasma_remove_item(s,id,ier)

      !  remove an existing item as identified by "id"...
      !  => the item remains in the dictionary, but its data type
      !     is set to zero -- indicating no associated data.
      !  
      !  the item is deallocated, but the name remains in the dictionary.
      !
      !  some items cannot be removed.  For example, a grid cannot be
      !  deleted, if a profile object exists which is dependent upon it.
      !
      !  see also xplasma_add_item -- this can be used to *replace* an
      !  an item, without the requirement that it be removed first.

      type (xplasma), pointer :: s
      integer, intent(in) :: id          ! id of existing item
      integer, intent(out) :: ier        ! status code returned

      !--------------------/local
      integer :: itype,iindex,i,ii,inp,icoord,igrid,irank
      integer :: inum,imatch,iprof,itchk
      integer :: luntty = 6
      !--------------------

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      call xplasma_author_check(s,id,ier)
      if(ier.ne.0) return

      if(xplasma_trace_name.ne.' ') then
         itchk=1
      else
         itchk=0
      endif

      !  save type information

      itype = s%dict(id)%dtype
      iindex= s%dict(id)%dindex

      !  error checks

      if(itype.eq.xplasma_coordType) then
         if(s%coords(iindex)%ngridc .gt. 0) then
            ier=401
            call xplasma_errmsg_append(s, &
                 ' coordinate name: '//trim(s%dict(id)%name))
         endif
      endif

      if(itype.eq.xplasma_gridType) then
         if(s%grids(iindex)%nrefs .gt. 0) then
            ier=301
            call xplasma_errmsg_append(s, &
                 ' grid name: '//trim(s%dict(id)%name))
         endif
      endif

      if(ier.ne.0) return

      !  clean up master list
      !  named slot & author name & labeling remains in dictionary...

      if(itchk.eq.1) then
         call tchk(s,'deleted in xplasma_remove_item',s%dict(id)%name)
      endif

      s%dict(id)%dtype=0
      s%dict(id)%dindex=0

      ! clean up constituent data item list

      if(itype.eq.xplasma_listType) then
         call xplist_free(s%lists(iindex))
         s%nlist_free = s%nlist_free + 1
         s%list_freelist(s%nlist_free) = iindex
      endif

      ! clean up coordinate object
      ! list is compressed now... no need for future garbage collection.

      if(itype.eq.xplasma_coordType) then
         call xpcoord_free(s%coords(iindex))
         s%ncoord_free = s%ncoord_free + 1
         s%coord_freelist(s%ncoord_free) = iindex
      endif

      ! clean up constituent grid
      ! remove reference in coordinate record

      if(itype.eq.xplasma_gridType) then
         icoord = s%grids(iindex)%coord  ! xplasma id
         icoord = s%dict(icoord)%dindex  ! actual index

         inum = s%coords(icoord)%ngridc
         imatch = 0
         do i=1,inum
            if(s%coords(icoord)%grid_ids(i).eq.id) then
               imatch=i
               exit
            endif
         enddo
         if(imatch.eq.0) then
            ier=9999  ! this indicates an error in the code
            call xplasma_errmsg_append(s, &
                 ' ?? unexpected xplasma error -- code error detected!')
            return    ! a grid must have its id recorded in its coordinate list
         endif

         !  remove grid reference
         
         do i=imatch+1,inum
            s%coords(icoord)%grid_ids(i-1) = s%coords(icoord)%grid_ids(i)
         enddo
         s%coords(icoord)%grid_ids(inum) = 0
         s%coords(icoord)%ngridc = inum - 1

         !  OK... now delete the grid description itself; item entry
         !  remains for later garbage collection...

         call xpgrid_free(s%grids(iindex))
         s%ngrid_free = s%ngrid_free + 1
         s%grid_freelist(s%ngrid_free) = iindex
         s%grid_equiv(iindex) = 0
      endif

      ! cleanup up constituent profile

      if(itype.eq.xplasma_profType) then

         ! decrement grid reference counters

         irank = s%profs(iindex)%rank
         do ii=1,irank
            igrid = s%profs(iindex)%gridIds(ii)  ! xplasma grid id
            igrid = s%dict(igrid)%dindex         ! actual index

            s%grids(igrid)%nrefs = s%grids(igrid)%nrefs - 1

         enddo

         ! at present a profile can have up to maxApro = 1 assocated profile.
         ! associations point both ways, so, fixup on delete...
         do ii=1,maxApro
            iprof = s%profs(iindex)%profIds(ii)
            if(iprof.gt.0) then
               iprof=s%dict(igrid)%dindex
               s%profs(iprof)%profIds(maxApro)=0 ! maxApro = 1 assumed here!
            endif
         enddo

         !  OK... now delete the profile itself, and compress the profile list.

         call xpprof_free(s%profs(iindex))
         s%nprof_free = s%nprof_free + 1
         s%prof_freelist(s%nprof_free) = iindex
      endif

      ! clean up numerical integrator control dataset

      if(itype.eq.xplasma_blkbxType) then
         call xpblkbx_free(s%blkbxs(iindex))
         s%nblkbx_free = s%nblkbx_free + 1
         s%blkbx_freelist(s%nblkbx_free) = iindex
      endif

      ! update (clear) internal pointers

      call xplasma_internal_ids(s)

    end subroutine xplasma_remove_item

    subroutine xplasma_internal_ids(s)
      
      type (xplasma), pointer :: s

      integer :: iertmp ! items not found-- not an error here

      call xplasma_find_item(s,"G",s%id_g,iertmp,nf_noerr=.TRUE.)
      call xplasma_find_item(s,"PSI",s%id_psi,iertmp,nf_noerr=.TRUE.)
      if(s%id_g.eq.0) s%id_P=0

      call xplasma_find_item(s,"R",s%id_R,iertmp,nf_noerr=.TRUE.)
      call xplasma_find_item(s,"Z",s%id_Z,iertmp,nf_noerr=.TRUE.)
      call xplasma_find_item(s,"__R_extrap",s%id_Rx,iertmp,nf_noerr=.TRUE.)
      call xplasma_find_item(s,"__Z_extrap",s%id_Zx,iertmp,nf_noerr=.TRUE.)

      call xplasma_find_item(s,"BR",s%id_BR,iertmp,nf_noerr=.TRUE.)
      call xplasma_find_item(s,"BZ",s%id_BZ,iertmp,nf_noerr=.TRUE.)
      call xplasma_find_item(s,"BMOD",s%id_BMOD,iertmp,nf_noerr=.TRUE.)

      call xplasma_find_item(s,"MAG_AXIS",s%id_axis,iertmp,nf_noerr=.TRUE.)
      call xplasma_find_item(s,"BDY_RZ_MINMAX",s%id_RZminmax,iertmp,nf_noerr=.TRUE.)

    end subroutine xplasma_internal_ids

    !========================================
    subroutine xplasma_copy(s1,s2,ier)
      !  full memory copy (not just pointer copy) s1=>s2
      !  s1 & s2 both need to have been initialized...

      type (xplasma), pointer :: s1,s2
      integer, intent(out) :: ier

      !--------------------------------------------
      integer :: isize,i,jsize,itchk
      integer :: luntty = 6
      character*32 :: tname
      !--------------------------------------------

      ier=0
      if((.not.associated(s1)).or.(.not.associated(s2))) then
         ier=101
      else if(s1%nguard.eq.0) then
         ier=101
      endif
      if(ier.ne.0) return

      if(s2%nguard.ne.0) call xplasma_freesub(s2,ier,total=.TRUE.)
      if(ier.ne.0) return

      internal_counter = internal_counter + 1
      s2%nguard=internal_counter

      itchk=0
      if(xplasma_trace_name.ne.' ') then
         write(luntty,*) ' %xplasma_obj(xplasma_COPY): object init by copy: '
         write(luntty,*) '  internal ID: ',s2%nguard
         tname=xplasma_trace_name
         call uupper(tname)
         itchk=1
      endif

      if(allocated(s2%label)) deallocate(s2%label)

      allocate(s2%label(size(s1%label)))
      s2%label = s1%label

!  copy AUTHOR LIST

      s2%max_author = s1%max_author
      s2%cur_author = s1%cur_author
      if(s2%max_author.gt.0) then
         allocate(s2%authors(s2%max_author)); s2%authors = ' '
         s2%authors(1:s2%cur_author) = s1%authors(1:s1%cur_author)
         allocate(s2%write_enable(s2%max_author)); s2%write_enable = 0
         s2%write_enable(1:s2%cur_author) = s1%write_enable(1:s1%cur_author)
      endif

!  copy DICTIONARY with sort order and array sizes preserved

      s2%nitems = s1%nitems
      s2%nmax   = s1%nmax

      isize = s2%nmax

      if(isize.gt.0) then
         allocate(s2%order(isize))
         allocate(s2%dict(isize))
         s2%order = s1%order
      endif

      do i=1,isize
         call xpltest_alloc(s2%dict(i))
         if(i.gt.s2%nitems) then
            s2%order(i)=0
            s2%dict(i)%name=' '
            s2%dict(i)%author=' '
            s2%dict(i)%dtype=0
            s2%dict(i)%dindex=0
         else
            s2%dict(i)%name = s1%dict(i)%name
            s2%dict(i)%author = s1%dict(i)%author
            s2%dict(i)%dtype = s1%dict(i)%dtype
            s2%dict(i)%dindex = s1%dict(i)%dindex
            if(itchk.gt.0) then
               if(s2%dict(i)%dindex.gt.0) then
                  if(s2%dict(i)%name.eq.tname) then
                     call tchk(s2,'created by xplasma_copy',tname)
                  endif
               endif
            endif
            if(allocated(s1%dict(i)%label)) then
               jsize = size(s1%dict(i)%label)
               allocate(s2%dict(i)%label(jsize))
               s2%dict(i)%label = s1%dict(i)%label
            endif
            if(allocated(s1%dict(i)%units)) then
               jsize = size(s1%dict(i)%units)
               allocate(s2%dict(i)%units(jsize))
               s2%dict(i)%units = s1%dict(i)%units
            endif
         endif
      enddo

      ! s2%nguard set by initialization already (otherwise xplasma_free returns
      ! an error code).

      !  flags...

      s2%axisymm = s1%axisymm
      s2%scrapeoff = s1%scrapeoff
      s2%bphi_ccw = s1%bphi_ccw
      s2%jphi_ccw = s1%jphi_ccw

      !  time:

      s2%time = s1%time

      !  LISTS...

      isize=s1%nlists
      s2%nlists = isize
      s2%nlist_free = s1%nlist_free

      if(isize.gt.0) then
         allocate(s2%lists(isize))
         allocate(s2%list_freelist(isize))
      endif

      do i=1,isize
         call xplist_fullcopy(s1%lists(i),s2%lists(i))
         s2%list_freelist(i)=s1%list_freelist(i)
      enddo

      !  COORDS...

      isize=s1%ncoords
      s2%ncoords = isize
      s2%ncoord_free = s1%ncoord_free

      if(isize.gt.0) then
         allocate(s2%coords(isize))
         allocate(s2%coord_freelist(isize))
      endif

      do i=1,isize
         call xpcoord_fullcopy(s1%coords(i),s2%coords(i))
         s2%coord_freelist(i)=s1%coord_freelist(i)
      enddo

      !  GRIDS...

      isize=s1%ngrids
      s2%ngrids = isize
      s2%ngrid_free = s1%ngrid_free

      if(isize.gt.0) then
         allocate(s2%grids(isize))
         allocate(s2%grid_freelist(isize))
         allocate(s2%grid_equiv(isize))
      endif

      do i=1,isize
         call xpgrid_fullcopy(s1%grids(i),s2%grids(i))
         s2%grid_freelist(i)=s1%grid_freelist(i)
         s2%grid_equiv(i)=s1%grid_equiv(i)
      enddo

      !  PROFILES...

      s2%eq_counter = s1%eq_counter
      s2%prof_counter = s1%prof_counter

      s2%bdytol = s1%bdytol
      s2%ajac_maxvar = s1%ajac_maxvar
      s2%rzOrder = s1%rzOrder
      s2%kmom = s1%kmom

      isize=s1%nprofs
      s2%nprofs = isize
      s2%nprof_free = s1%nprof_free

      if(isize.gt.0) then
         allocate(s2%profs(isize))
         allocate(s2%prof_freelist(isize))
      endif

      do i=1,isize
         call xpprof_fullcopy(s1%profs(i),s2%profs(i))
         s2%prof_freelist(i)=s1%prof_freelist(i)
      enddo
            
      !  BLACK BOXES... --also, control data for numerical integrations...

      isize=s1%nblkbxs
      s2%nblkbxs = isize
      s2%nblkbx_free = s1%nblkbx_free

      if(isize.gt.0) then
         allocate(s2%blkbxs(isize))
         allocate(s2%blkbx_freelist(isize))
      endif

      do i=1,isize
         call xpblkbx_fullcopy(s1%blkbxs(i),s2%blkbxs(i))
         s2%blkbx_freelist(i)=s1%blkbx_freelist(i)
      enddo

      call xplasma_internal_ids(s2)
 
    end subroutine xplasma_copy

    subroutine xplasma_avg_elements(s,zwght,savg,zauth,ier)

      !  for elements owned by author (zauth), with floating point data d_s
      !  from s and d_savg from savg,
      !     d_savg <-- (1-zwght)*d_savg + wght*d_s
      !  and replace the floating point data in savg.

      !  for profiles to be replaced, the profile name must match,
      !  the author must match AND the grids must match;

      !  for list element data to be replaced, the listname and the
      !  element name must match

      !  for blackbox data, the black box name and the R8 sizes must match.

      type (xplasma), pointer :: s
      real*8, intent(in) :: zwght     ! BTW 0 and 1 !!
      type (xplasma), pointer :: savg
      character*(*), intent(in) :: zauth
      integer, intent(out) :: ier

      !----------------------------------
      integer :: i,ii,jj,isize,isizg,isiza,isizl,isizp,id,ida,ier1,ier2,idiff
      integer :: id1,id2,id3,id1a,id2a,id3a,inum1,inum2,irank1,irank2
      integer :: iap1,iap2,isiz1,isiz2,irc
      integer :: ityp1,ityp2,htyp1(10),htyp2(10)
      character*32 :: zauthor,zname
      real*8 :: zwghti,ztol,zval
      integer, dimension(:), allocatable :: idmatch,ids,idd
      real*8, dimension(:), allocatable :: g1,g2,vals,vala
      real*8, dimension(:), pointer :: r8n,r8a
      integer, dimension(:), pointer :: i4n,i4a
      character*32, dimension(:), allocatable :: elems,elema
      logical :: allmatch
      !----------------------------------

      ier=0

      !-------------------------------
      !  authorname in uppercase

      zauthor=zauth
      call uupper(zauthor)

      !-------------------------------
      !  basic error check: 0 <= zwght <= 1

      zwghti=zwght
      if(zwghti.lt.czero) then
         if(zwghti.ge.-ceps10) then
            zwghti=czero
         else
            ier=9999
            call xplasma_errmsg_append(s, &
                 "xplasma_avg_elements: averaging weight not btw 0.0 and 1.0: less than zero.")
         endif

      else if(zwghti.gt.cone) then
         if(zwghti.le.cone+ceps10) then
            zwghti=cone
         else
            ier=9999
            call xplasma_errmsg_append(s, &
                 "xplasma_avg_elements: averaging weight not btw 0.0 and 1.0: greater than one.")
         endif

      endif
      if(ier.ne.0) return

      !-------------------------------
      ! check for grid matches
      ! for each grid ids(i) in s, find a matching grid idmatch(i)

      call xplasma_find_grids(s,ier, num_grids=isize)
      if(ier.ne.0) return

      allocate(idmatch(isize),ids(isize)); idmatch = 0

      ! do loop to catch deallocation at end
      do
         call xplasma_find_grids(s,ier, id_grids=ids)
         if(ier.ne.0) exit

         do i=1,isize
            id=ids(i)
            call xplasma_grid_info(s, id, ier, name=zname, size=isizg)
            if(ier.ne.0) exit

            ! find corresponding grid in savg
            call xplasma_gridId(savg, zname, ida)

            if(ida.eq.0) cycle  ! no match: grid not found in savg
            call xplasma_grid_info(savg, ida, ier, size=isiza)
            if(ier.ne.0) exit

            if(isizg.ne.isiza) cycle  ! no match: sizes differ

            allocate(g1(isiza),g2(isiza))

            call xplasma_grid(s,id,g1,ier1)
            call xplasma_grid(s,id,g2,ier2)
            ier=max(ier1,ier2)
            if(ier.eq.0) then
               idiff=0
               ztol=ceps10*max(abs(g1(1)),abs(g1(isiza)))
               do ii=1,isiza
                  if(abs(g1(ii)-g2(ii)).gt.ztol) then
                     idiff=ii
                     exit
                  endif
               enddo
            endif
            deallocate(g1,g2)

            if(ier.ne.0) exit
            if(idiff.ne.0) cycle

            idmatch(i)=ida  ! no error; grid s[ids(i)] matchs savg[idmatch(i)]
         enddo
         if(ier.ne.0) exit

         !-------------------------------
         ! find profiles in s owned by requested author

         call xplasma_find_profs(s, ier, author_only=zauthor, num_profs=isizp)
         if(ier.ne.0) exit

         allocate(idd(isizp))
         call xplasma_find_profs(s, ier, author_only=zauthor, id_profs=idd)
         if(ier.ne.0) then
            deallocate(idd)
            exit
         endif

         ! for each found profile in s, find corresponding profile in savg

         do i=1,isizp
            id=idd(i)
            call xplasma_prof_info(s,id,ier, &
                 rank=irank1, splineType=ityp1, hybridType=htyp1, &
                 gridId1=id1, gridId2=id2, gridId3=id3,  name=zname)
            if(ier.ne.0) exit

            call xplasma_profId(savg,zname,ida)
            if(ida.eq.0) cycle  ! no profile found

            call xplasma_prof_info(savg,ida,ier, &
                 rank=irank2, splineType=ityp2, hybridType=htyp2, &
                 gridId1=id1a, gridId2=id2a, gridId3=id3a)
            if(ier.ne.0) exit

            if(irank1.ne.irank2) cycle  ! check rank match
            if(ityp1.ne.ityp2) cycle    ! check spline type match

            idiff=0                     ! and check hybrid spline type match
            do ii=1,irank1
               if(htyp1(ii).ne.htyp2(ii)) idiff=idiff+1
            enddo
            if(idiff.gt.0) cycle

            !  check the grids

            irc=0
            do ii=1,isize
               if(id1.eq.ids(ii)) then
                  if(id1a.ne.idmatch(ii)) idiff=idiff+1
                  irc=irc+1
                  if(irc.eq.irank1) exit
               endif
               if(irank1.gt.1) then
                  if(id2.eq.ids(ii)) then
                     if(id2a.ne.idmatch(ii)) idiff=idiff+1
                     irc=irc+1
                     if(irc.eq.irank1) exit
                  endif
               endif
               if(irank1.gt.2) then
                  if(id3.eq.ids(ii)) then
                     if(id3a.ne.idmatch(ii)) idiff=idiff+1
                     irc=irc+1
                     if(irc.eq.irank1) exit
                  endif
               endif
            enddo

            if(idiff.gt.0) cycle

            ! OK: name, rank, interp. type, and grids ALL match, so...

            iap1 = s%dict(id)%dindex
            iap2 = savg%dict(ida)%dindex

            isiz1 = s%profs(iap1)%size
            isiz2 = savg%profs(iap2)%size

            isiz1=min(isiz1,isiz2) ! these should match, though

            do ii=1,isiz1
               zval = s%profs(iap1)%buf(ii)*zwght
               zval = zval + savg%profs(iap2)%buf(ii)*(cone-zwght)
               savg%profs(iap2)%buf(ii)=zval
            enddo
         enddo
         deallocate(idd)
         if(ier.ne.0) exit

         exit
      enddo

      deallocate(idmatch,ids)

      !------------------------------------
      ! average matching elements on scalar lists

      call xplasma_num_items(s,ier, author_only=zauthor, num_lists=isize)
      if(ier.ne.0) return

      allocate(ids(isize))
      call xplasma_contents(s,ier, author_only=zauthor, id_lists=ids)
      if(ier.ne.0) then
         deallocate(ids)
         return
      endif

      do i=1,isize
         id=ids(i)
         call xplasma_list_info(s,id,ier, name=zname,nelems=isizl)
         if(ier.ne.0) exit

         call xplasma_listId(savg,zname,ida)
         if(ida.eq.0) cycle

         call xplasma_list_info(s,ida,ier, nelems=isiza)
         if(ier.ne.0) exit

         allocate(elems(isizl),vals(isizl),elema(isiza),vala(isiza))

         call xplasma_getList_names(s,id,elems,ier1, r8vals=vals)
         call xplasma_getList_names(savg,ida,elema,ier2, r8vals=vala)
         ier=max(ier1,ier2)
         if(ier.eq.0) then
            do ii=1,isizl
               do jj=1,isiza
                  if(elems(ii).eq.elema(jj)) then
                     zval=zwght*vals(ii) + (1-zwght)*vala(jj)
                     vala(jj)=zval
                     exit
                  endif
               enddo
            enddo

            call xplasma_author_set(savg,zauthor,ier1)
            call xplasma_setList_r8vals(savg,ida,vala,ier)
            call xplasma_author_clear(savg,zauthor,ier1)

         endif

         deallocate(elems,vals,elema,vala)
         if(ier.ne.0) exit
      enddo

      deallocate(ids)
      if(ier.ne.0) return

      !------------------------------------
      ! average blackbox R8 data

      call xplasma_num_items(s,ier, author_only=zauthor, num_blkbxs=isize)
      if(ier.ne.0) return

      allocate(ids(isize))
      call xplasma_contents(s,ier, author_only=zauthor, id_blkbxs=ids)
      if(ier.ne.0) then
         deallocate(ids)
         return
      endif

      do i=1,isize
         id=ids(i)
         call xplasma_blackbox_info(s,id,ier, name=zname, r8size=isizl)
         if(ier.ne.0) exit

         call xplasma_blkbxId(savg,zname,ida)
         if(ida.eq.0) cycle

         call xplasma_blackbox_info(savg,ida,ier, r8size=isiza)
         if(ier.ne.0) exit

         if(isiza.ne.isizl) cycle  ! size mismatch

         call xplasma_blackBox_retrieve(s,id,ier, ia_ptr=i4n, r8a_ptr=r8n)
         if(ier.ne.0) exit

         call xplasma_blackBox_retrieve(savg,ida,ier, ia_ptr=i4a, r8a_ptr=r8a)
         if(ier.ne.0) exit

         do ii=1,isiza
            zval = zwght*r8n(ii) + (1-zwght)*r8a(ii)
            r8a(ii) = zval
         enddo
      enddo

      deallocate(ids)
      if(ier.ne.0) return
      
    end subroutine xplasma_avg_elements
    
    !========================================
    subroutine xplasma_write(s,filename,ier)

      ! write NetCDF file of xplasma contents

      type (xplasma), pointer :: s
      character*(*), intent(in) :: filename  ! file to write...
      integer, intent(out) :: ier        ! status code returned

      !--------------------------------------
      integer :: icdf   ! NetCDF file handle
      !--------------------------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         ier=101
      endif
      if(ier.ne.0) return

      call ezcdf_open(icdf,filename,'w',ier)
      if(ier.ne.0) then
         ier=10
         call xplasma_errmsg_append(s, &
              ' ezcdf_open (write) failure on file: '//trim(filename))
         return
      endif

      call xplasma_ncdefine(s,icdf,ier)
      call xplasma_ncwrite(s,icdf,ier)

      call ezcdf_close(icdf)

    end subroutine xplasma_write

    !========================================
    subroutine xplasma_pack(s,iptr,r8ptr,ier)

      ! encode entire xplasma contents into two dynamically allocated
      ! arrays, iptr(...) -- integer; r8ptr(...) -- real*8

      ! CAUTION -- it is up to caller to assure that no previously allocated
      ! memory is on these pointers -- they are reallocated here...

      type (xplasma), pointer :: s
      integer, dimension(:), pointer :: iptr
      real*8, dimension(:), pointer :: r8ptr
      integer, intent(out) :: ier        ! status code returned

      !--------------------------------------
      integer :: iictr,r8ctr  ! size counters
      !--------------------------------------

      nullify(iptr)
      nullify(r8ptr)

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         ier=101
      endif
      if(ier.ne.0) return

      allocate(iptr(1000),r8ptr(1000))  ! initial allocation

      iptr(1)=123456789   ! guard values
      iptr(2)=-iptr(1)    ! if xplasma_unpack sees these at the start of
      iptr(3:4)=iptr(1:2) ! an integer array, it assumes input is valid.

      iptr(5)=6           ! #integers in use
      iptr(6)=0           ! #floats in use

      iptr(7:size(iptr)) = 0

      r8ptr = 0.0
      r8ctr = 0

      call xplasma_packmaster(s,iptr,r8ptr)

    end subroutine xplasma_pack

    subroutine xplasma_unpack(s,iarray,r8array,ier)

      !  unpack previously packed integer and float arrays to
      !  reconstruct an xplasma object

      type (xplasma), pointer :: s       ! object contents replaced
      integer, dimension(:), intent(in) :: iarray  ! integer packed data
      real*8, dimension(:), intent(in) :: r8array  ! real*8 packed data
      integer, intent(out) :: ier        ! status code returned

      !-----------------------------
      integer, parameter :: itest = 123456789
      integer :: isize,r8size,id_xmom2d
      !-----------------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'xplasma_read:initialization',ier)
      endif
      if(ier.ne.0) return

      if(size(iarray).lt.1000) then
         ier=9999
         call xplasma_errmsg_append(s,'?xplasma_unpack: iarray too small.')
      endif

      if(size(r8array).lt.1000) then
         ier=9999
         call xplasma_errmsg_append(s,'?xplasma_unpack: r8array too small.')
      endif

      if(ier.ne.0) return

      if((iarray(1).ne.itest).or.(iarray(2).ne.-itest).or. &
           (iarray(3).ne.itest).or.(iarray(4).ne.-itest)) then
         ier=9999
         call xplasma_errmsg_append(s, &
              '?xplasma_unpack: invalid unpack data (guard words invalid).')
         return
      endif

      isize = size(iarray)
      r8size = size(r8array)

      if(iarray(5).gt.isize) then
         call xplasma_errmsg_append(s, &
              '?xplasma_unpack: iarray smaller than size data found.')
         ier=9999
      endif

      if(iarray(6).gt.r8size) then
         call xplasma_errmsg_append(s, &
              '?xplasma_unpack: r8array smaller than size data found.')
         ier=9999
      endif

      if(ier.ne.0) return

      !  OK -- pack arguments pass all tests;
      !  delete old object contents & replace...

      call xplasma_freesub(s,ier,total=.TRUE.)

      call xplasma_unpackmaster(s,iarray,r8array,ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s, '?xplasma_unpackmaster error.')
         ier=9999

      else
         call xplasma_internal_ids(s)

         call xplasma_find_item(s,'__XMOM2D',id_xmom2d,ier,nf_noerr=.TRUE.)
         if(id_xmom2d.ne.0) call xpmom_fetch(s)
      endif

    end subroutine xplasma_unpack

    !========================================
    subroutine xplasma_read(s,filename,ier)

      ! read NetCDF file of xplasma contents
      ! mod DMC Sept 2006 -- if filename is URL (starts with "http://") try
      ! to download...

      type (xplasma), pointer :: s
      character*(*), intent(in) :: filename  ! file to read...
      integer, intent(out) :: ier        ! status code returned

      !--------------------------------------
      integer :: icdf   ! NetCDF file handle

!      integer :: jsystem,ii,ilen,istat,iflag
      integer :: ii,ilen,istat,iflag
      integer :: luntty = 6
      character*8 ztest
      character*200 zcmd
      !--------------------------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'xplasma_read:initialization',ier)
      endif
      if(ier.ne.0) return

      iflag=0
      ilen=len(trim(filename))
      if(ilen.gt.7) then
         ztest=filename(1:8)
         call uupper(ztest)
         if((ztest(1:7).eq.'HTTP://').and.(ztest(8:8).ne.'/')) then
            write(luntty,*) ' '
            write(luntty,*) ' %xplasma_read: attempting URL download: ', &
                 filename(1:ilen)
            iflag=1
         endif
      endif

      if(iflag.eq.1) then
         do 
            if(filename(ilen:ilen).eq.'/') then
               ilen=ilen-1
            else
               exit
            endif
         enddo
         do ii=ilen-1,1,-1
            if(filename(ii:ii).eq.'/') exit
         enddo
         zcmd = '/usr/bin/curl '//trim(filename)//' > ./'//filename(ii+1:ilen)
         write(luntty,*) trim(zcmd)
         istat = jsystem(zcmd)
         if(istat.ne.0) then
            ier=11
            write(luntty,*)  ' download failed! '
         else

            write(luntty,*)  ' download to : '//filename(ii+1:ilen)//' complete.'
            call ezcdf_open(icdf,filename(ii+1:ilen),'r',ier)

         endif
         write(luntty,*) ' '

      else
         !  normal file read...
         call ezcdf_open(icdf,filename,'r',ier)
      endif

      if(ier.ne.0) then
         ier=11
         call xplasma_errmsg_append(s, &
              ' ezcdf_open (read) failure on file: '//trim(filename))
         return
      endif

      call xplasma_ncread(s,icdf,ier)

      call ezcdf_close(icdf)

    end subroutine xplasma_read

    !========================================
    subroutine xplasma_find_item(s,zname,id,ier, nf_noerr)
      
      !  find named item -- return id

      type (xplasma), pointer :: s
      character*(*), intent(in) :: zname  ! name of item to find
      integer, intent(out) :: id         ! id of item found
      integer, intent(out) :: ier        ! status code returned

      logical, intent(in), optional :: nf_noerr  ! .TRUE. to suppress error

      !  if there is no such item, id=0 and ier=106 or ier=108 are returned.
      !  if (nf_noerr) is .TRUE. id=0 and ier=0 are returned if not found.

      !--------------------/local
      integer :: iloc,iexact
      logical :: inf_noerr
      !--------------------

      id=0

      inf_noerr=.FALSE.
      if(present(nf_noerr)) inf_noerr = nf_noerr

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_find_item call)',ier)
      endif
      if(ier.ne.0) return

      call xplasma_locate(s,zname,iloc,iexact)
      if(iexact.eq.1) then
         id = s%order(iloc)
         if(s%dict(id)%dtype.eq.0) then
            if(.not.inf_noerr) then
               ier = 108  ! the item once existed but was deleted...
               call xplasma_errmsg_append(s, &
                    'item "'//trim(s%dict(id)%name)//'" was deleted.')
            endif
            id = 0
         endif
      else
         if(.not.inf_noerr) then
            ier = 106
         endif
      endif

    end subroutine xplasma_find_item

    !========================================
    subroutine xplasma_get_item_info(s,id,ier, &
         nf_noerr, &
         name,label,units,itype,author)

      ! return info on item-- what info determined by presence of
      ! optional arguments

      type (xplasma), pointer :: s
      integer, intent(in) :: id    ! id of item for which info is requested
      integer, intent(out) :: ier  ! return status code  (0=OK)

      !  optional control:
      logical, intent(in), optional :: nf_noerr  ! .TRUE. to allow ivalid id

      !  here is what can be returned...

      character*(*), intent(out), optional :: name     ! name of item
      character*(*), intent(out), optional :: label    ! item label
      character*(*), intent(out), optional :: units    ! item units label
      integer, intent(out), optional :: itype          ! item type
      character*(*), intent(out), optional :: author   ! name of author

      !--------------
      integer :: ic,iskip
      !--------------

      iskip=0
      ier=0
      if(present(nf_noerr)) then
         if(nf_noerr) then
            if((id.lt.1).or.(id.gt.s%nitems)) then
               iskip=1
            endif
         endif
      endif

      if(iskip.eq.0) then
         call xplasma_errck0(s,id,ier)
      endif

      if((iskip.eq.1).or.(ier.ne.0)) then
         if(present(name)) name=' '
         if(present(label)) label=' '
         if(present(units)) units=' '
         if(present(itype)) itype=0
         if(present(author)) author=' '
         return
      endif

      if(present(name)) then
         name=s%dict(id)%name
      endif

      if(present(label)) then
         if(.not.allocated(s%dict(id)%label)) then
            label=' '
         else
            do ic=1,len(label)
               if(ic.le.size(s%dict(id)%label)) then
                  label(ic:ic)=s%dict(id)%label(ic)
               else
                  label(ic:ic)=' '
               endif
            enddo
         endif
      endif

      if(present(units)) then
         if(.not.allocated(s%dict(id)%units)) then
            units=' '
         else
            do ic=1,len(units)
               if(ic.le.size(s%dict(id)%units)) then
                  units(ic:ic)=s%dict(id)%units(ic)
               else
                  units(ic:ic)=' '
               endif
            enddo
         endif
      endif

      if(present(itype)) then
         itype=s%dict(id)%dtype
      endif

      if(present(author)) then
         author=s%dict(id)%author
      endif

    end subroutine xplasma_get_item_info

    !========================================
    subroutine xplasma_listId(s,name,id)

      ! return list ID given name of list -- no error code; id=0 if not
      ! found or if there is an error...

      type(xplasma), pointer :: s
      character*(*), intent(in) :: name
      integer, intent(out) :: id

      integer :: ier,itype

      call xplasma_find_item(s,name,id,ier,nf_noerr=.TRUE.)
      if(id.eq.0) return

      if(s%dict(id)%dtype.ne.xplasma_listType) then
         id=0  ! not a list...
      endif

    end subroutine xplasma_listId

    !========================================
    subroutine xplasma_coordId(s,name,id)

      ! return coord ID given name of coord -- no error code; id=0 if not
      ! found or if there is an error...

      type(xplasma), pointer :: s
      character*(*), intent(in) :: name
      integer, intent(out) :: id

      integer :: ier,itype

      call xplasma_find_item(s,name,id,ier,nf_noerr=.TRUE.)
      if(id.eq.0) return

      if(s%dict(id)%dtype.ne.xplasma_coordType) then
         id=0  ! not a coord...
      endif

    end subroutine xplasma_coordId

    !========================================
    subroutine xplasma_gridId(s,name,id)

      ! return grid ID given name of grid -- no error code; id=0 if not
      ! found or if there is an error...

      type(xplasma), pointer :: s
      character*(*), intent(in) :: name
      integer, intent(out) :: id

      integer :: ier,itype

      call xplasma_find_item(s,name,id,ier,nf_noerr=.TRUE.)
      if(id.eq.0) return

      if(s%dict(id)%dtype.ne.xplasma_gridType) then
         id=0  ! not a grid...
      endif

    end subroutine xplasma_gridId

    !========================================
    subroutine xplasma_profId(s,name,id)

      ! return prof ID given name of prof -- no error code; id=0 if not
      ! found or if there is an error...

      type(xplasma), pointer :: s
      character*(*), intent(in) :: name
      integer, intent(out) :: id

      integer :: ier,itype

      call xplasma_find_item(s,name,id,ier,nf_noerr=.TRUE.)
      if(id.eq.0) return

      if(s%dict(id)%dtype.ne.xplasma_profType) then
         id=0  ! not a prof...
      endif

    end subroutine xplasma_profId

    !========================================
    subroutine xplasma_blkbxId(s,name,id)

      ! return blkbx ID given name of blkbx -- no error code; id=0 if not
      ! found or if there is an error...

      type(xplasma), pointer :: s
      character*(*), intent(in) :: name
      integer, intent(out) :: id

      integer :: ier,itype

      call xplasma_find_item(s,name,id,ier,nf_noerr=.TRUE.)
      if(id.eq.0) return

      if(s%dict(id)%dtype.ne.xplasma_blkbxType) then
         id=0  ! not a blkbx...
      endif

    end subroutine xplasma_blkbxId

    !========================================
    subroutine xplasma_author_query(s,author,ienable,ier)

      ! return information on current author with xplasma write access

      type (xplasma), pointer :: s
      character*(*), intent(out) :: author  ! current author with write access
      logical, intent(out) :: ienable       ! access enable flag T/F
      integer, intent(out) :: ier           ! status code, 0=OK

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_author_query call)',ier)
      endif
      if(ier.ne.0) return
      
      author = s%authors(s%cur_author)
      ienable = (s%write_enable(s%cur_author) .eq. 1)

    end subroutine xplasma_author_query
    
    !========================================
    subroutine xplasma_remove_author(s,authorName,ier)

      ! remove ALL ITEMS associated with author = authorName

      type (xplasma), pointer :: s
      character*(*), intent(in) :: authorName
      integer, intent(out) :: ier

      !----------------------
      character*32 :: zauthor
      integer :: i,iertmp
      !----------------------

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_author_query call)',ier)
      endif
      if(ier.ne.0) return

      zauthor = authorName
      call uupper(zauthor)

      call xplasma_author_set(s,zauthor,iertmp)

      !  remove profiles first
      do i=1,s%nitems
         if(zauthor.eq.s%dict(i)%author) then
            if(s%dict(i)%dindex.gt.0) then
               if(s%dict(i)%dtype.eq.xplasma_profType) then
                  call xplasma_remove_item(s,i,iertmp)
               endif
            endif
         endif
      enddo

      !  remove grids next
      do i=1,s%nitems
         if(zauthor.eq.s%dict(i)%author) then
            if(s%dict(i)%dindex.gt.0) then
               if(s%dict(i)%dtype.eq.xplasma_gridType) then
                  call xplasma_remove_item(s,i,iertmp)
               endif
            endif
         endif
      enddo

      !  remove the rest...
      do i=1,s%nitems
         if(zauthor.eq.s%dict(i)%author) then
            if(s%dict(i)%dindex.gt.0) then
               call xplasma_remove_item(s,i,iertmp)
            endif
         endif
      enddo

      call xplasma_author_clear(s,zauthor,iertmp)

    end subroutine xplasma_remove_author

    !========================================
    subroutine xplasma_error(s,ier,ilun)

      type (xplasma), pointer :: s          ! msg details may be in s msg buffer(s)
      integer, intent(in) :: ier  ! error code (returned by xplasma routine)
      integer, intent(in) :: ilun ! Fortran I/O unit number where to write msg

      !  error reporting routine
      !  ier < 0 -- warning
      !  ier = 0 -- no error, but write trailer messages...
      !  ier > 0 -- error

      integer :: i

      !------------------------
      !  only enter main body of routine, if the xplasma object appears to
      !  be properly initialized.

      if(.not.associated(s)) then
         write(ilun,*) ' ?xplasma_error: object pointer (s) not associated.'
         return

      else if(s%nguard.eq.0) then
         write(ilun,*) ' ?xplasma_error: xplasma object was not initialized.'
         return

      endif

      if(ier.eq.11) then
         write(ilun,*) ' ?xplasma_write: NetCDF file open failure.'
      else if(ier.eq.12) then
         write(ilun,*) ' ?xplasma_read: NetCDF file open failure.'
      else if(ier.eq.13) then
         write(ilun,*) ' ?xplasma_read: NetCDF file read error.'

      else if(ier.eq.50) then
         write(ilun,*) &
              ' ?xplasma:  xplasma object locked, no modifications permitted.'
      else if(ier.eq.51) then
         write(ilun,*) &
              ' ?xplasma:  attempt to modify item owned by another author.'
      else if(ier.eq.52) then
         write(ilun,*) ' ?xplasma: author name verification failed (mismatch).'
      else if(ier.eq.53) then
         write(ilun,*) ' ?xplasma: stack underflow -- no known authors.'
      else if(ier.eq.54) then
         write(ilun,*) ' ?xplasma: author name contains illegal character or is too long.'

      else if(ier.eq.60) then
         write(ilun,*) ' %xplasma_contents: array too small, contents truncated.'
      else if(ier.eq.61) then
         write(ilun,*) ' %xplasma_find_*: array too small, contents truncated.'

      else if(ier.eq.88) then
         write(ilun,*) ' ?xplasma_init: code error, coordinate id mismatch.'

      else if(ier.eq.89) then
         write(ilun,*) ' ?xplasma_create_prof: reserved name error.'

      else if(ier.eq.90) then
         write(ilun,*) &
              ' ?xplasma_create_prof: mixed (inconsistent) spline selection.'
         write(ilun,*) &
              '  --both "ispline" and "hspline" optional arguments are present.'

      else if(ier.eq.101) then
         write(ilun,*) ' ?received NULL xplasma object pointer.'

      else if(ier.eq.102) then
         write(ilun,*) ' ?xplasma item name argument contains invalid character.'
         write(ilun,*) '  -- all characters must be alphanumeric or "$" or "_".'
         write(ilun,*) '  -- first character must not be a numeric digit.'
         write(ilun,*) '  -- no imbedded whitespace; name must not be blank.'
      else if(ier.eq.103) then
         write(ilun,*) ' ?xplasma item name too long: >32 characters.'
      else if(ier.eq.104) then
         write(ilun,*) ' ?xplasma item type argument value is invalid.'
      else if(ier.eq.105) then
         write(ilun,*) ' ?xplasma item id is invalid.'
      else if(ier.eq.106) then
         write(ilun,*) ' ?xplasma name not valid-- no such item.'
      else if(ier.eq.107) then
         write(ilun,*) ' ?xplasma -- non-axisymmetry not implemented!'
      else if(ier.eq.108) then
         write(ilun,*) ' ?xplasma name not valid-- item was deleted.'
      else if(ier.eq.109) then
         write(ilun,*) ' ?xplasma -- scrapeoff region (limiter & grids) undefined.'

      else if(ier.eq.150) then
         write(ilun,*) &
              ' ?xplasma -- error in G-eqdsk MHD eq. mapping (eqi_fromgeqdsk)'

      else if(ier.eq.200) then
         write(ilun,*) ' ?xplist: specified item is not a list.'
      else if(ier.eq.201) then
         write(ilun,*) ' ?xplist: attempt to access unitialized list object.'
      else if(ier.eq.202) then
         write(ilun,*) ' ?xplist: illegal character in list element name.'
      else if(ier.eq.203) then
         write(ilun,*) ' ?xplist: list element name too long (max 32 chars).'
      else if(ier.eq.204) then
         write(ilun,*) ' ?xplist: size passed does not match actual list size.'
      else if(ier.eq.205) then
         write(ilun,*) ' ?xplist: list of element names contains duplicates.'

      else if(ier.eq.300) then
         write(ilun,*) ' ?xpgrid: specified item is not a grid.'
      else if(ier.eq.301) then
         write(ilun,*) ' ?xpgrid: grid item in use, cannot be deleted.'
      else if(ier.eq.302) then
         write(ilun,*) ' ?xpgrid: passed grid has nx < 2.'
      else if(ier.eq.303) then
         write(ilun,*) ' ?xpgrid: passed grid not strict ascending.'
      else if(ier.eq.304) then
         write(ilun,*) ' ?xpgrid: pspline r8genxpkg failure (unexpected).'
      else if(ier.eq.305) then
         write(ilun,*) ' ?xpgrid: passed grid size does not match stored grid size.'
      else if(ier.eq.306) then
         write(ilun,*) ' ?xpgrid: range of rho (or x) grid must be [0,1]'
      else if(ier.eq.307) then
         write(ilun,*) ' ?xpgrid: angle coordinate max-min = 2pi required.'
      else if(ier.eq.308) then
         write(ilun,*) ' ?xpgrid: [xmin,xmax] not consistent with existing discretizations.'
         write(ilun,*) '  use xplasma_coord_info to get expected range.'

      else if(ier.eq.309) then
         write(ilun,*) '  grid UNITS label cannot be directly modified.'
      else if(ier.eq.400) then
         write(ilun,*) ' ?xpcoord: specified item is not a coordinate.'
      else if(ier.eq.401) then
         write(ilun,*) ' ?xpcoord: grid discretizations exist, cannot delete coordinate.'
      else if(ier.eq.402) then
         write(ilun,*) ' ?xpcoord: coordinate grid index out of range.'
      else if(ier.eq.403) then
         write(ilun,*) ' ?xpcoord: coordinate contains no defined grids.'

      else if(ier.eq.500) then
         write(ilun,*) ' ?xpprof: specified item is not an Xplasma profile.'
      else if(ier.eq.501) then
         write(ilun,*) ' ?xpprof: profile x grid id invalid.'
      else if(ier.eq.502) then
         write(ilun,*) ' ?xpprof: explicit boundary condition data size is incorrect.'
      else if(ier.eq.503) then
         write(ilun,*) ' ?xpprof: profile contents have not been defined.'
      else if(ier.eq.504) then
         write(ilun,*) ' ?xpprof: invalid boundary condition control argument.'
      else if(ier.eq.505) then
         write(ilun,*) ' ?xpprof: xplasma_define_prof spline type error.'
      else if(ier.eq.506) then
         write(ilun,*) ' ?xpprof: size of profile data is incorrect.'
      else if(ier.eq.507) then
         write(ilun,*) ' ?xpprof: not a 1d profile.'
      else if(ier.eq.508) then
         write(ilun,*) ' ?xpprof: profile has more x axes than space provided.'
      else if(ier.eq.510) then
         write(ilun,*) ' ?xpprof: input vector / output vector sizes inconsistent.'
      else if(ier.eq.511) then
         write(ilun,*) ' ?xpprof: profile independent coordinates inconsistent.'
      else if(ier.eq.512) then
         write(ilun,*) ' ?xpprof: profile rank inconsistent with evaluation routine.'
      else if(ier.eq.513) then
         write(ilun,*) ' ?xpprof: interpolating function cannot evaluate derivative.'
      else if(ier.eq.514) then
         write(ilun,*) ' ?xpprof: x arguments for interpolation are out of bounds.'
      else if(ier.eq.515) then
         write(ilun,*) ' ?xpprof: in a call involving a profile over multiple'
         write(ilun,*) '  coordinates, the same coordinate is used more than once.'
      else if(ier.eq.516) then
         write(ilun,*) ' ?xpprof: scalar and vector derivative controls both present.'
      else if(ier.eq.517) then
         write(ilun,*) ' ?xpprof: x grid lookup information missing or from wrong grid.'
      else if(ier.eq.518) then
         write(ilun,*) ' ?xpprof: size mismatch x grid lookup information; return vector size.'
      else if(ier.eq.519) then
         write(ilun,*) ' ?xpprof: rank of profile not correct.'
      else if(ier.eq.520) then
         write(ilun,*) ' ?xpprof: associated profile ID is not valid.'
      else if(ier.eq.521) then
         write(ilun,*) ' ?xplasma_getProf: non-existent coefficient data requested.'
      else if(ier.eq.522) then
         write(ilun,*) ' ?xplasma_getProf: dimensioning of coefficient data array not correct.'
      else if(ier.eq.523) then
         write(ilun,*) ' ?xplasma_create_prof: coeff optional input only valid if spline or hermite is present.'
      else if(ier.eq.524) then
         write(ilun,*) ' ?xplasma_create_prof: dimensioning of coeff array is not correct.'
      else if(ier.eq.525) then
         write(ilun,*) ' ?xplasma_create_prof: coeff array data and BC data cannot BOTH be present.'

      else if(ier.eq.600) then
         write(ilun,*) ' ?xpblkbx: passed id does not refer to valid x grid or coordinate.'
      else if(ier.eq.601) then
         write(ilun,*) ' ?xpblkbx: passed id is not "rho".'
      else if(ier.eq.602) then
         write(ilun,*) ' ?xpblkbx: passed id is not "theta".'
      else if(ier.eq.603) then
         write(ilun,*) ' ?xpblkbx: passed integration grid not strict ascending.'
      else if(ier.eq.604) then
         write(ilun,*) ' ?xpblkbx: optional subroutine arguments error.'
      else if(ier.eq.605) then
         write(ilun,*) ' ?xpblkbx: passed id not a 1d integrator dataset.'
      else if(ier.eq.606) then
         write(ilun,*) ' ?xpblkbx: specified item is not a "black box".'
      else if(ier.eq.607) then
         write(ilun,*) ' ?xpblkbx: "black box" contents were never defined.'
      else if(ier.eq.608) then
         write(ilun,*) ' ?xpblkbx: {R,Z} flux surfaces must be defined first.'
      else if(ier.eq.609) then
         write(ilun,*) ' ?xpblkbx: integration range exceeds coordinate range.'
      else if(ier.eq.610) then
         write(ilun,*) ' ?xpblkbx: integration result array size incorrect.'
      else if(ier.eq.611) then
         write(ilun,*) ' ?xpblkbx: invalid name for numerical integration.'
      else if(ier.eq.612) then
         write(ilun,*) ' ?xplasma: passed array is of incorrect length.'
      else if(ier.eq.613) then
         write(ilun,*) ' ?xpblkbx: dVol or dVol/drho was already stored.'
      else if(ier.eq.614) then
         write(ilun,*) ' ?xpblkbx: no cache available in this integrator.'
      else if(ier.eq.615) then
         write(ilun,*) ' ?xpblkbx: no [R,Z] grid rectangle defined.'
      else if(ier.eq.616) then
         write(ilun,*) ' ?xplasma_create_blackbox: error with optional arguments.'

      else if(ier.eq.701) then
         write(ilun,*) ' ?xplasma_prof: grid not on "R" coordinate.'
      else if(ier.eq.702) then
         write(ilun,*) ' ?xplasma_prof: grid not on "Z" coordinate.'

      else if(ier.eq.2001) then
         write(ilun,*) ' ?xplasma_rzgeo module error.'
      else if(ier.eq.2002) then
         write(ilun,*) ' ?xplasma_rzgeo: Jacobian singularity check failure.'
      else if(ier.eq.2003) then
         write(ilun,*) ' ?xplasma_rzgeo: Boundary shape check failure.'
         write(ilun,*) '  atan((Zbdy(theta)-Zaxis)/(Rbdy(theta)-Raxis)) '// &
              'must be monotonic with margin to spare.'

      else if(ier.eq.2501) then
         write(ilun,*) ' ?xplasma_ctrans: ambiguous input.'

      else if(ier.eq.2502) then
         write(ilun,*) ' ?xplasma_ctrans: no output specified.'

      else if(ier.eq.2503) then
         write(ilun,*) ' ?xplasma_ctrans: an argument vector length is too small.'

      else if(ier.eq.2504) then
         write(ilun,*) ' ?xplasma_ctrans: {R,Z} flux surfaces have not been defined.'

      else if(ier.eq.2509) then
         write(ilun,*) ' ?xplasma_ctrans: an axisymmetric geometry was expected.'
      else if(ier.eq.2510) then
         write(ilun,*) ' ?xplasma_rzjac: actual or near singularity in metric Jacobian'
      else if(ier.eq.2511) then
         write(ilun,*) ' ?xplasma_ctrans: 2x2 Newton search failure.'

      else if(ier.eq.3000) then
         write(ilun,*) ' ?xplasma_sol: error setting up contour limiter.'
      else if(ier.eq.3001) then
         write(ilun,*) ' ?xplasma_sol: an equilibrium must be defined before.'
         write(ilun,*) '  a limiter based on plasma boundary location or on'
         write(ilun,*) '  "circles and lines" can be specified.'
      else if(ier.eq.3002) then
         write(ilun,*) ' ?xplasma_mklim_special: error with optional arguments.'
      else if(ier.eq.3003) then
         write(ilun,*) ' ?xplasma_limcheck: plasma boundary intersects limiter.'
      else if(ier.eq.3004) then
         write(ilun,*) ' ?xplasma_limcheck: no limiter has been defined.'
      else if(ier.eq.3005) then
         write(ilun,*) ' ?xplasma_limcheck: no plasma equilibrium has been defined.'
         write(ilun,*) '  plasma {R,Z}(rho,theta) flux surfaces not defined.'

      else if(ier.eq.3006) then
         write(ilun,*) '  R&Z rectangle grid range error.'

      else if(ier.ne.0) then
         write(ilun,*) ' ?xplasma object error (unspecified).'
      endif

      !  except for "NULL pointer" error, write out additional
      !  data on error if available; then clear the msg* buffers
      !  (so that stale messages won't be seen on reuse of object)

      if(ier.eq.101) return

      do i=1,s%nmsgs
         write(ilun,*) trim(s%msgs(i))
      enddo

      call xplasma_clear_msgs(s)

    end subroutine xplasma_error

    subroutine xplasma_errmsg_append(s,zmsg)
      !  append message to xplasma object error message buffers

      type (xplasma), pointer :: s
      character*(*), intent(in) :: zmsg

      if(s%nmsgs.eq.maxmsgs) then
         s%msgs(s%nmsgs) = ' ** ERROR REPORT BUFFER OVERFLOW **'

      else
         s%nmsgs=s%nmsgs+1
         s%msgs(s%nmsgs) = zmsg
      endif

    end subroutine xplasma_errmsg_append

    !------------------------------------------------------
    subroutine xplasma_find_grids(s,ier, &
         author_only, icoord, &
         num_grids, id_grids)

      ! find number of grids & grid ids

      type (xplasma), pointer :: s
      integer, intent(out) :: ier ! error code -- only set if s not initialized

      character*(*), intent(in), optional :: author_only
      !  set this to restrict the list to items owned by the specified author

      integer, intent(in), optional :: icoord
      !  set this to restrict the list to grids that discretize the coordinate
      !  indicated by icoord; if icoord=0 (the default) all grids are returned.

      integer, intent(out), optional :: num_grids
      !  number of grids found

      integer, dimension(:), intent(out), optional :: id_grids
      !  the list of ids of the grids found...

      !-----------------------------
      integer :: inum,isize
      integer, dimension(:), allocatable :: ids
      integer :: jcoord,ii,jj,itmp,iertmp
      !-----------------------------

      ier = 0

      jcoord = 0
      if(present(icoord)) jcoord=icoord

      !------------------

      if(present(num_grids)) num_grids=0
      if(present(id_grids)) id_grids=0

      !------------------

      call xplasma_num_items(s,ier, author_only, num_grids=inum)
      if(ier.ne.0) return
      if(inum.eq.0) return

      allocate(ids(inum))
      call xplasma_contents(s,ier, author_only, id_grids=ids)
      if(ier.ne.0) then
         deallocate(ids)
         return
      endif

      if(present(id_grids)) isize=size(id_grids)

      if(jcoord.eq.0) then
         if(present(num_grids)) num_grids=inum
         if(present(id_grids)) then
            if(isize.lt.inum) then
               ier=61
               id_grids = ids(1:isize)
            else
               id_grids(1:inum)=ids(1:inum)
            endif
         endif

      else

         ii=0
         do jj=1,inum
            call xplasma_grid_info(s,ids(jj),iertmp, coord=itmp)
            if(itmp.eq.jcoord) then

               ii=ii+1
               if(present(id_grids)) then
                  if(ii.gt.isize) then
                     ier=61
                  else
                     id_grids(ii)=ids(jj)
                  endif
               endif
            endif
         enddo

         if(present(num_grids)) num_grids=ii

      endif

      deallocate(ids)

    end subroutine xplasma_find_grids

    !------------------------------------------------------
    subroutine xplasma_find_profs(s, ier, &
         author_only, icoord1, icoord2, icoord3, &
         num_profs, id_profs)

      ! return the total number and/or the ids of profiles which
      ! are defined over the specified coordinate(s).

      type (xplasma), pointer :: s
      integer, intent(out) :: ier ! error code -- only set if s not initialized

      character*(*), intent(in), optional :: author_only
      !  set this to restrict the list to items owned by the specified author

      integer, intent(in), optional :: icoord1
      integer, intent(in), optional :: icoord2
      integer, intent(in), optional :: icoord3
      !  set these to restrict the list to profiles over the specified
      !  coordinates:  if NONE are set, a list of all profiles are returned;
      !  if only icoord1 is set, profiles vs. (icoord1) are returned
      !  if icoord1,icoord2 are set, profiles vs. (icoord1,icoord2) 
      !  (or icoord2,icoord1) are returned; if all 3 are set, only profiles
      !  vs. all three coordinates in any order are returned.

      integer, intent(out), optional :: num_profs
      !  number of profs found

      integer, dimension(:), intent(out), optional :: id_profs
      !  the list of ids of the profs found...

      !-----------------------------
      integer :: inum,isize
      integer, dimension(:), allocatable :: ids
      integer :: jcoord(3),irank,ii,jj,k,kk,itmp,iertmp
      integer :: jrank,gridId1,gridId2,gridId3,imatch
      integer :: gridCoord(3)
      !-----------------------------

      ier = 0

      jcoord = 0
      if(present(icoord1)) jcoord(1)=icoord1
      if(present(icoord2)) jcoord(2)=icoord2
      if(present(icoord3)) jcoord(3)=icoord3

      irank=0
      do kk=1,size(jcoord)
         if(jcoord(kk).eq.0) exit
         irank=irank+1
      enddo

      !------------------

      if(present(num_profs)) num_profs=0
      if(present(id_profs)) id_profs=0

      !------------------

      call xplasma_num_items(s,ier, author_only, num_profs=inum)
      if(ier.ne.0) return
      if(inum.eq.0) return

      allocate(ids(inum))
      call xplasma_contents(s,ier, author_only, id_profs=ids)
      if(ier.ne.0) then
         deallocate(ids)
         return
      endif

      if(present(id_profs)) isize=size(id_profs)

      if(irank.eq.0) then

         ! return all...

         if(present(num_profs)) num_profs=inum
         if(present(id_profs)) then
            if(isize.lt.inum) then
               ier=61
               id_profs=ids(1:isize)
            else
               id_profs(1:inum)=ids
            endif
         endif
               
      else

         ii=0
         do jj=1,inum

            call xplasma_prof_info(s,ids(jj),iertmp, rank=jrank, &
                 gridId1=gridId1,gridId2=gridId2,gridId3=gridId3)

            if(irank.ne.jrank) cycle

            call xplasma_grid_info(s,gridId1,iertmp, coord=gridCoord(1))

            if(irank.ge.2) then
               call xplasma_grid_info(s,gridId2,iertmp, coord=gridCoord(2))
            endif

            if(irank.ge.3) then
               call xplasma_grid_info(s,gridId3,iertmp, coord=gridCoord(3))
            endif

            do k=1,irank
               imatch=0
               do kk=1,irank
                  if(gridCoord(kk).eq.jcoord(k)) then
                     imatch=kk
                     exit
                  endif
               enddo
               if(imatch.eq.0) exit
            enddo

            if(imatch.eq.0) cycle

            ! profile matches all coordinates & has correct rank...

            ii=ii+1
            if(present(id_profs)) then
               if(ii.le.isize) then
                  id_profs(ii)=ids(jj)
               else
                  ier=61
               endif
            endif
         enddo

         if(present(num_profs)) num_profs=ii

      endif
      deallocate(ids)

    end subroutine xplasma_find_profs

    !------------------------------------------------------
    subroutine xplasma_find_blkbxs(s,ier, &
         author_only, itype, &
         num_blkbxs, id_blkbxs)

      ! find number of black boxes and/or black box ids
      ! only count those whose type matches "itype"-- if specified.

      type (xplasma), pointer :: s
      integer, intent(out) :: ier ! error code -- only set if s not initialized

      character*(*), intent(in), optional :: author_only
      !  set this to restrict the list to items owned by the specified author

      integer, intent(in), optional :: itype
      !  set this to restrict the list to black boxes of the specified type.
      !  if itype=0 (the default) all available black boxes are returned.

      integer, intent(out), optional :: num_blkbxs
      !  number of black boxes found

      integer, dimension(:), intent(out), optional :: id_blkbxs
      !  the list of ids of the black boxes found...

      !-----------------------------
      integer :: inum,isize
      integer, dimension(:), allocatable :: ids
      integer :: jtype,ii,jj,itmp,iertmp
      !-----------------------------

      ier = 0

      jtype = 0
      if(present(itype)) jtype=itype

      !------------------

      if(present(num_blkbxs)) num_blkbxs=0
      if(present(id_blkbxs)) id_blkbxs=0

      !------------------

      call xplasma_num_items(s,ier, author_only, num_blkbxs=inum)
      if(ier.ne.0) return
      if(inum.eq.0) return

      allocate(ids(inum))
      call xplasma_contents(s,ier, author_only, id_blkbxs=ids)
      if(ier.ne.0) then
         deallocate(ids)
         return
      endif

      if(present(id_blkbxs)) isize=size(id_blkbxs)

      if(jtype.eq.0) then
         if(present(num_blkbxs)) num_blkbxs=inum
         if(present(id_blkbxs)) then
            if(isize.lt.inum) then
               ier=61
               id_blkbxs = ids(1:isize)
            else
               id_blkbxs(1:inum)=ids(1:inum)
            endif
         endif

      else

         ii=0
         do jj=1,inum
            call xplasma_blackbox_info(s,ids(jj),iertmp, type=itmp)
            if(itmp.eq.jtype) then

               ii=ii+1
               if(present(id_blkbxs)) then
                  if(ii.gt.isize) then
                     ier=61
                  else
                     id_blkbxs(ii)=ids(jj)
                  endif
               endif
            endif
         enddo

         if(present(num_blkbxs)) num_blkbxs=ii

      endif

      deallocate(ids)

    end subroutine xplasma_find_blkbxs

    !------------------------------------------------------
    subroutine xplasma_num_items(s,ier, &
         author_only, &
         num_lists, num_plists, num_coords, num_grids, num_profs, num_blkbxs)

      !  return the total numbers of items of various types in xplasma 
      !  container object.

      type (xplasma), pointer :: s
      integer, intent(out) :: ier ! error code -- only set if s not initialized

      character*(*), intent(in), optional :: author_only
      !  set this to restrict the number to items owned by the specified author

      integer, intent(out), optional :: num_lists   ! no. of lists
      integer, intent(out), optional :: num_plists  ! no. lists refering
                                                    !     to profiles
      integer, intent(out), optional :: num_coords  ! no. of coordinates
      integer, intent(out), optional :: num_grids   ! no. of grids
      integer, intent(out), optional :: num_profs   ! no. of profiles
      integer, intent(out), optional :: num_blkbxs  ! no. of black boxes

      !---------------------------
      integer :: i,itype
      character*32 :: zauth
      logical :: iauth
      !---------------------------

      if(present(num_lists)) num_lists=0
      if(present(num_plists)) num_plists=0
      if(present(num_coords)) num_coords=0
      if(present(num_grids)) num_grids=0
      if(present(num_profs)) num_profs=0
      if(present(num_blkbxs)) num_blkbxs=0

      if(present(author_only)) then
         iauth=.TRUE.
         zauth=author_only
         call uupper(zauth)
      else
         iauth=.FALSE.
      endif

      ier=0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_num_items call)',ier)
      endif
      if(ier.ne.0) return

      do i=1,s%nitems

         if(iauth) then
            if(s%dict(i)%author.ne.zauth) cycle
         endif

         itype = s%dict(i)%dtype
         if(present(num_lists)) then
            if(itype.eq.xplasma_listType) num_lists=num_lists+1
         endif
         if(present(num_plists)) then
            if(itype.eq.xplasma_listType) then
               if(xp_is_plist(s,i)) num_plists=num_plists+1
            endif
         endif
         if(present(num_coords)) then
            if(itype.eq.xplasma_coordType) num_coords=num_coords+1
         endif
         if(present(num_grids)) then
            if(itype.eq.xplasma_gridType) num_grids=num_grids+1
         endif
         if(present(num_profs)) then
            if(itype.eq.xplasma_profType) num_profs=num_profs+1
         endif
         if(present(num_blkbxs)) then
            if(itype.eq.xplasma_blkbxType) num_blkbxs=num_blkbxs+1
         endif
      enddo

    end subroutine xplasma_num_items

    subroutine xplasma_contents(s,ier, &
         author_only, &
         id_lists, id_plists, id_coords, id_grids, id_profs, id_blkbxs)

      !  return sorted list(s) of ids of objects of indicated type(s).
      !    lists of lists, lists of profiles, etc.

      type (xplasma), pointer :: s
      integer, intent(out) :: ier ! error code, 0=OK

      character*(*), intent(in), optional :: author_only
      !  set this to restrict the list to items owned by the specified author

      integer, dimension(:), intent(out), optional :: id_lists  ! lists
      integer, dimension(:), intent(out), optional :: id_plists ! lists that
                                                      ! refer to profiles
      integer, dimension(:), intent(out), optional :: id_coords ! coordinates
      integer, dimension(:), intent(out), optional :: id_grids  ! grids
      integer, dimension(:), intent(out), optional :: id_profs  ! profiles
      integer, dimension(:), intent(out), optional :: id_blkbxs ! black boxes

      !------------------------------
      integer :: i,ii,itype,ilists,iplists,icoords,igrids,iprofs,iblkbxs
      character*128 msgbuf
      character*32 :: zauth
      logical :: iauth
      !------------------------------

      if(present(id_lists)) id_lists = 0
      if(present(id_plists)) id_plists = 0
      if(present(id_coords)) id_coords = 0
      if(present(id_grids)) id_grids = 0
      if(present(id_profs)) id_profs = 0
      if(present(id_blkbxs)) id_blkbxs = 0

      if(present(author_only)) then
         iauth=.TRUE.
         zauth=author_only
         call uupper(zauth)
      else
         iauth=.FALSE.
      endif

      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_contents call)',ier)
      endif
      if(ier.ne.0) return

      ilists=0
      iplists=0
      icoords=0
      igrids=0
      iprofs=0
      iblkbxs=0

      do i=1,s%nitems
         call xplasma_order(s,i,ii,ier)

         if(iauth) then
            if(s%dict(ii)%author.ne.zauth) cycle
         endif

         if(present(id_lists)) then
            if(s%dict(ii)%dtype.eq.xplasma_listType) then
               ilists=ilists+1
               if(ilists.le.size(id_lists)) then
                  id_lists(ilists)=ii
               else
                  ier=60
                  if(ilists.eq.(size(id_lists)+1)) then
                     call xplasma_errmsg_append(s, &
                          '  ...list of lists truncated.')
                  endif
               endif
            endif
         endif

         if(present(id_plists)) then
            if(s%dict(ii)%dtype.eq.xplasma_listType) then
               if(xp_is_plist(s,ii)) then
                  iplists=iplists+1
                  if(iplists.le.size(id_plists)) then
                     id_plists(iplists)=ii
                  else
                     ier=60
                     if(iplists.eq.(size(id_plists)+1)) then
                        call xplasma_errmsg_append(s, &
                             '  ...list of profile lists truncated.')
                     endif
                  endif
               endif
            endif
         endif

         if(present(id_coords)) then
            if(s%dict(ii)%dtype.eq.xplasma_coordType) then
               icoords=icoords+1
               if(icoords.le.size(id_coords)) then
                  id_coords(icoords)=ii
               else
                  ier=60
                  if(icoords.eq.(size(id_coords)+1)) then
                     call xplasma_errmsg_append(s, &
                          '  ...list of coordinates truncated.')
                  endif
               endif
            endif
         endif

         if(present(id_grids)) then
            if(s%dict(ii)%dtype.eq.xplasma_gridType) then
               igrids=igrids+1
               if(igrids.le.size(id_grids)) then
                  id_grids(igrids)=ii
               else
                  ier=60
                  if(igrids.eq.(size(id_grids)+1)) then
                     call xplasma_errmsg_append(s, &
                          '  ...list of grids truncated.')
                  endif
               endif
            endif
         endif

         if(present(id_profs)) then
            if(s%dict(ii)%dtype.eq.xplasma_profType) then
               iprofs=iprofs+1
               if(iprofs.le.size(id_profs)) then
                  id_profs(iprofs)=ii
               else
                  ier=60
                  if(iprofs.eq.(size(id_profs)+1)) then
                     call xplasma_errmsg_append(s, &
                          '  ...list of profiles truncated.')
                  endif
               endif
            endif
         endif

         if(present(id_blkbxs)) then
            if(s%dict(ii)%dtype.eq.xplasma_blkbxType) then
               iblkbxs=iblkbxs+1
               if(iblkbxs.le.size(id_blkbxs)) then
                  id_blkbxs(iblkbxs)=ii
               else
                  ier=60
                  if(iblkbxs.eq.(size(id_blkbxs)+1)) then
                     call xplasma_errmsg_append(s, &
                          '  ...list of black  boxes truncated.')
                  endif
               endif
            endif
         endif

      enddo

    end subroutine xplasma_contents

    !---------------------------
    ! private logical function

    logical function xp_is_plist(s,id)
      !  return TRUE if list contains profiles ONLY
      !  otherwise FALSE

      type (xplasma), pointer :: s
      integer, intent(in) :: id  ! list ID (not checked)

      !-------------------------------
      integer :: iertmp,isize,ielem,itype,idpro
      character*32, dimension(:), allocatable :: pnames
      character*32 ztest1,ztest2
      integer, dimension(:), allocatable :: ids
      !-------------------------------

      call xplasma_getList(s,id,isize,iertmp)
      if((isize.eq.0).or.(iertmp.ne.0)) then
         xp_is_plist=.FALSE.
         return
      endif

      xp_is_plist = .TRUE.

      allocate(pnames(isize),ids(isize))

      do
         call xplasma_getlist_names(s,id,pnames,iertmp)
         if(iertmp.ne.0) exit
         call xplasma_getlist_ivals(s,id,ids,iertmp)
         if(iertmp.ne.0) exit

         do ielem=1,isize
            ztest1=pnames(ielem)

            idpro=abs(ids(ielem))
            call xplasma_get_item_info(s,idpro,iertmp, &
                 itype=itype, name=ztest2)

            if(itype.ne.xplasma_profType) then
               xp_is_plist=.FALSE.
               exit
            else if(iertmp.ne.0) then
               exit
            else
               if(ztest1.ne.ztest2) then
                  xp_is_plist=.FALSE.
                  exit
               endif
            endif
         enddo

         exit
      enddo

      if(iertmp.ne.0) xp_is_plist = .FALSE.

      deallocate(pnames,ids)
    end function xp_is_plist

    !---------------------------
    ! private sorting routine

    subroutine xplasma_order(s,i,isort,ier)

      !  return the index to the item in order position (i)
      !  return 0 if (i) is out of range

      type (xplasma), pointer :: s
      integer, intent(in) :: i
      integer, intent(out) :: isort
      integer, intent(out) :: ier ! error code -- only set if s not initialized

      isort=0
      ier = 0

      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then

         ier = 101

      else
         if((i.gt.0).and.(i.le.s%nitems)) then

            isort = s%order(i)

         endif
      endif

    end subroutine xplasma_order

    !========================================
    !  Public methods dealing with lists...
    !========================================

    subroutine xplasma_create_list(s,listname,enames,id,ier, &
         label, units, ivals, r8vals, chvals, change_status)
      !
      !  create list & define list elements
      !
      type (xplasma), pointer :: s
      character*(*), intent(in) :: listname        ! name of list being created
      character*(*), intent(in), dimension(:) :: enames ! names of elements
      integer, intent(out) :: id          ! id of list item
      integer, intent(out) :: ier         ! completion code

      character*(*), intent(in), optional :: label  ! label for entire list
      character*(*), intent(in), optional :: units  ! units label for list

      integer,intent(in), dimension(:), optional :: ivals   ! integer values
      real*8, intent(in), dimension(:), optional :: r8vals  ! real*8 values
      character*(*), intent(in), dimension(:), optional :: chvals  ! str values
      integer, intent(out), optional :: change_status

      !  change_status optional output:
      !    0 -- list contents unmodified;
      !    1 -- elements preserved but data values changed;
      !    2 -- elements were added or deleted;
      !    3 -- list did not previously exist; was created.

      !---------------
      integer :: iertmp,nelems,ielems,i,idiff,jstat,inum_changed
      character*32 zname,tmp_units
      character*128 tmp_label
      character*32 :: chknames(size(enames))
      logical :: val_only
      !---------------

      if(present(change_status)) change_status = -1 ! error value

      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_create_list call)',ier)
      endif
      if(ier.ne.0) return

      !  see if list is already known and only values are being changed...

      val_only = .FALSE.

      call xplasma_listId(s,listname,id)
      if(id.eq.0) then
         jstat = 3   ! a new list
      else
         jstat = 2   ! for now assume changes to elements

         call xplasma_list_info(s,id,iertmp, nelems=ielems)
         if(ielems.eq.size(enames)) then
            call xplasma_getList_names(s,id,chknames,iertmp)

            idiff=0
            do i=1,ielems
               zname=enames(i)
               call uupper(zname)
               if(chknames(i).ne.zname) then
                  idiff=1
                  exit
               endif
            enddo

            if(idiff.eq.0) then
               jstat = 0 ! for now assume NO changes
               val_only = .TRUE.
            endif
         endif
      endif
               
      nelems=size(enames)

      if(.not.val_only) then

         id=0
         call xplasma_add_item(s,listname,xplasma_listType,id,ier)
         if(ier.ne.0) return

         call xplasma_define_list(s,id,nelems,enames,ier)
         if(ier.ne.0) return

      endif

      iertmp=0

      if(present(label)) then
         call xplasma_label_item(s,id,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_label_item(s,id,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      if(present(ivals)) then
         call xplasma_setList(s,id,ivals,ier, nchange=inum_changed)
         if(inum_changed.gt.0) jstat=max(1,jstat)
         iertmp=max(iertmp,ier)
      endif

      if(present(r8vals)) then
         call xplasma_setList(s,id,r8vals,ier, nchange=inum_changed)
         if(inum_changed.gt.0) jstat=max(1,jstat)
         iertmp=max(iertmp,ier)
      endif

      if(present(chvals)) then
         call xplasma_setList(s,id,chvals,ier, nchange=inum_changed)
         if(inum_changed.gt.0) jstat=max(1,jstat)
         iertmp=max(iertmp,ier)
      endif

      ier=iertmp

      if(ier.eq.0) then
         zname=listname
         call uupper(zname)
         if(zname.eq.'MAG_AXIS') then
            s%id_axis = id
         endif
         if(zname.eq.'BDY_RZ_MINMAX') then
            s%id_RZminmax = id
         endif
      endif

      if(present(change_status)) change_status = jstat

      return
    end subroutine xplasma_create_list

    subroutine xplasma_define_list(s,id,nelems,enames,ier)
      !
      !  **private** now...
      !  define list elements -- item created by xplasma_add_item call
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      integer, intent(in) :: nelems       ! no. of list elements
      character*(*), intent(in) :: enames(nelems)  ! names of list elements
      integer, intent(out) :: ier         ! completion code

      !------------
      character*32 :: enamu,enamu2
      integer :: i,j,iadr
      !------------

      call xplasma_errck0_list(s,id,ier)
      if(ier.ne.0) return

      call xplasma_author_check(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex

      call xplist_free(s%lists(iadr))

      ier=0
      if(max(0,nelems).eq.0) return

      do i=1,nelems
         call xplasma_checkname_sub(enames(i),ier)
         if(ier.eq.1) then
            ier=202
            call xplasma_errmsg_append(s, &
                 'element name string was: "'//trim(enames(i))//'".')
         else if(ier.eq.2) then
            ier=203
            call xplasma_errmsg_append(s, &
                 'element name string was: "'//trim(enames(i))//'".')
         else if(ier.ne.0) then
            ier=9999
            call xplasma_errmsg_append(s, &
                 'element name string was: "'//trim(enames(i))//'".')
         endif
         if(ier.ne.0) return
      enddo

      ! check against name duplication within list

      do i=1,nelems-1
         enamu=enames(i)
         call uupper(enamu)
         do j=i+1,nelems
            enamu2=enames(j)
            call uupper(enamu2)

            if(enamu.eq.enamu2) then
               ier = 205
               call xplasma_errmsg_append(s, &
                    'list element name "'//trim(enamu)// &
                    '" occurs more than once.')
            endif
         enddo
      enddo

      if(ier.ne.0) return

      ! OK...

      s%lists(iadr)%size=nelems
      allocate(s%lists(iadr)%enames(nelems))
      allocate(s%lists(iadr)%chindx(2,nelems))
      allocate(s%lists(iadr)%r8vals(nelems))
      allocate(s%lists(iadr)%ivals(nelems))
      do i=1,nelems
         enamu=enames(i)
         call uupper(enamu)
         s%lists(iadr)%enames(i)=enamu
         s%lists(iadr)%chindx(1:2,i)=0
         s%lists(iadr)%r8vals(i)=0
         s%lists(iadr)%ivals(i)=0
      enddo

    end subroutine xplasma_define_list

    !========================================
    subroutine xplasma_merge_lists(s,id1,id2,ier, idelete)
      !
      !  modify list (id1) by merging in elements from 2nd list (id2)
      !  optionally delete list (id2) when done.
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id1  ! list merged into
      integer, intent(in) :: id2  ! list merged from
      integer, intent(out) :: ier ! status completion code

      logical, intent(in), optional :: idelete ! .TRUE. to delete list (id2)
      !  when merge is done -- default is .FALSE.

      !------------------------
      integer :: i,j,inold,innew,intot,ifound,iertmp
      integer :: iadr1,iadr2,isize1,isize2,icmax,ic,ic1,ic2,ie1,ie2,ik

      character*32 :: enamu
      
      !  from the original list...

      character*32, dimension(:), pointer :: ep      ! names of list elements
      real*8, dimension(:), pointer :: r8vals        ! element values (real*8)
      integer, dimension(:), pointer :: ivals        ! element values (integer)
      integer, dimension(:,:), pointer :: chindx     ! ch indices for merger
      character, dimension(:), pointer :: chbuf      ! ch buffer for merger

      logical, allocatable, dimension(:) :: ldup
      integer, allocatable, dimension(:) :: kdup

      logical :: jdelete
      !------------------------

      jdelete = .FALSE.
      if(present(idelete)) jdelete=idelete

      call xplasma_errck0_list(s,id1,ier)
      if(ier.ne.0) return

      call xplasma_errck0_list(s,id2,ier)
      if(ier.ne.0) return

      call xplasma_author_check(s,id1,ier)
      if(ier.ne.0) return

      if(jdelete) then
         call xplasma_author_check(s,id2,ier)
         if(ier.ne.0) return
      endif

      ! first look for replacement elements (names matching names already
      ! in list); zero out any corresponding values; count the number that
      ! are *not* replacements.

      iadr1 = s%dict(id1)%dindex
      iadr2 = s%dict(id2)%dindex

      isize1 = s%lists(iadr1)%size
      isize2 = s%lists(iadr2)%size

      allocate(kdup(isize1)); kdup = 0
      allocate(ldup(isize2))

      inold = isize1

      innew=0
      do i=1,isize2
         ifound=0
         enamu=s%lists(iadr2)%enames(i)
         do j=1,inold
            if(enamu.eq.s%lists(iadr1)%enames(j)) then
               ifound=j
               exit
            endif
         enddo
         if(ifound.eq.0) then
            innew=innew+1
            ldup(i)=.FALSE.
         else
            ldup(i)=.TRUE.
            kdup(ifound)=i
            s%lists(iadr1)%r8vals(ifound)=s%lists(iadr2)%r8vals(i)
            s%lists(iadr1)%ivals(ifound)=s%lists(iadr2)%ivals(i)
            !  deferred: character data transfer...
         endif
      enddo

      intot=inold+innew

      ! update character data...

      allocate(chindx(2,intot))
      icmax = s%lists(iadr1)%chindx(2,isize1) + s%lists(iadr2)%chindx(2,isize2)
      allocate(chbuf(icmax))

      ie2=0

      do i=1,isize1
         if(kdup(i).gt.0) then
            ik=kdup(i)
            ic1=s%lists(iadr2)%chindx(1,ik)
            ic2=s%lists(iadr2)%chindx(2,ik)

            ie1=ie2+1
            ie2=ie2+(ic2-ic1)+1

            do ic=ic1,ic2
               chbuf(ie1+ic-ic1) = s%lists(iadr2)%chbuf(ic)
            enddo

         else
            ic1=s%lists(iadr1)%chindx(1,i)
            ic2=s%lists(iadr1)%chindx(2,i)

            ie1=ie2+1
            ie2=ie2+(ic2-ic1)+1

            do ic=ic1,ic2
               chbuf(ie1+ic-ic1) = s%lists(iadr1)%chbuf(ic)
            enddo
         endif
         chindx(1,i) = ie1
         chindx(2,i) = ie2
      enddo

      if(innew.gt.0) then

         ! append list to existing, non-empty list.

         s%lists(iadr1)%size = intot

         ep => s%lists(iadr1)%enames
         r8vals => s%lists(iadr1)%r8vals
         ivals => s%lists(iadr1)%ivals

         allocate(s%lists(iadr1)%enames(intot))
         s%lists(iadr1)%enames(1:inold) = ep(1:inold)

         allocate(s%lists(iadr1)%r8vals(intot))
         s%lists(iadr1)%r8vals(1:inold) = r8vals(1:inold)

         allocate(s%lists(iadr1)%ivals(intot))
         s%lists(iadr1)%ivals(1:inold)  = ivals(1:inold)

         j=inold
         do i=1,isize2
            if(.not.ldup(i)) then
               j=j+1
               s%lists(iadr1)%enames(j)=s%lists(iadr2)%enames(i)
               s%lists(iadr1)%r8vals(j)=s%lists(iadr2)%r8vals(i)
               s%lists(iadr1)%ivals(j)=s%lists(iadr2)%ivals(i)

               ic1=s%lists(iadr2)%chindx(1,i)
               ic2=s%lists(iadr2)%chindx(2,i)

               ie1=ie2+1
               ie2=ie2+(ic2-ic1)+1

               do ic=ic1,ic2
                  chbuf(ie1+ic-ic1) = s%lists(iadr2)%chbuf(ic)
               enddo
               chindx(1,j) = ie1
               chindx(2,j) = ie2

            endif
         enddo

      endif

      if(associated(s%lists(iadr1)%chindx)) deallocate(s%lists(iadr1)%chindx)
      allocate(s%lists(iadr1)%chindx(2,intot))
      s%lists(iadr1)%chindx = chindx

      if(associated(s%lists(iadr1)%chbuf)) deallocate(s%lists(iadr1)%chbuf)
      if(ie2.gt.0) then
         allocate(s%lists(iadr1)%chbuf(ie2))
         s%lists(iadr1)%chbuf(1:ie2) = chbuf(1:ie2)
      endif

      deallocate(ldup,kdup)

      if(jdelete) then
         call xplasma_remove_item(s,id2,iertmp)
      endif

    end subroutine xplasma_merge_lists

    !========================================
    subroutine xplasma_list_info(s,id,ier, &
         nelems,name,label,units,author)
      !
      !  get the size of an existing list
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      integer, intent(out) :: ier         ! completion code

      integer, intent(out), optional :: nelems  ! no. of list elements in list

      character*(*), intent(out), optional :: name
      character*(*), intent(out), optional :: label
      character*(*), intent(out), optional :: units
      character*(*), intent(out), optional :: author

      !---------------------
      integer :: iadr,iertmp
      !---------------------

      if(present(nelems)) nelems=0

      if(present(name)) name=' '
      if(present(label)) label=' '
      if(present(units)) units=' '
      if(present(author)) author=' '

      call xplasma_errck0_list(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex

      if(present(nelems)) nelems = s%lists(iadr)%size

      iertmp = 0

      if(present(name)) then
         call xplasma_get_item_info(s,id,ier, name=name)
         iertmp=max(iertmp,ier)
      endif

      if(present(label)) then
         call xplasma_get_item_info(s,id,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_get_item_info(s,id,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      if(present(author)) then
         call xplasma_get_item_info(s,id,ier, author=author)
         iertmp=max(iertmp,ier)
      endif

      ier=iertmp

    end subroutine xplasma_list_info

    !========================================
    subroutine xplasma_getList_size(s,id,nelems,ier)
      !
      !  get the size of an existing list
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      integer, intent(out) :: nelems      ! no. of list elements in list
      integer, intent(out) :: ier         ! completion code


      !---------------------
      integer :: iadr
      !---------------------

      call xplasma_errck0_list(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      nelems = s%lists(iadr)%size

    end subroutine xplasma_getList_size

    !========================================
    subroutine xplasma_getList_names(s,id,enames,ier, &
         name,label,units,author,ivals,r8vals,chvals)
      !
      !  get the element names from an existing list
      !   optional arguments allow additional data to be fetched as well.
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      character*(*), dimension(:), intent(out) :: enames ! the names...
      integer, intent(out) :: ier         ! completion code

      character*(*), intent(out), optional :: name   ! name (of whole list)
      character*(*), intent(out), optional :: label  ! label (for whole list)
      character*(*), intent(out), optional :: units  ! units label (whole list)
      character*(*), intent(out), optional :: author ! author

      integer, intent(out), dimension(:), optional :: ivals   ! integer values
      real*8, intent(out), dimension(:), optional :: r8vals   ! real*8 values
      character*(*), intent(out), dimension(:), optional :: chvals  ! str vals

      !---------------------
      integer :: iadr,iertmp,nelems
      !---------------------

      nelems = size(enames)

      call xplasma_errck0_list_adr(s,id,nelems,iadr,ier)
      if(ier.ne.0) return

      enames(1:nelems)=s%lists(iadr)%enames(1:nelems)

      iertmp=0

      if(present(name)) then
         call xplasma_get_item_info(s,id,ier, name=name)
         iertmp=max(iertmp,ier)
      endif

      if(present(label)) then
         call xplasma_get_item_info(s,id,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_get_item_info(s,id,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      if(present(author)) then
         call xplasma_get_item_info(s,id,ier, author=author)
         iertmp=max(iertmp,ier)
      endif

      if(present(ivals)) then
         call xplasma_getList(s,id,ivals,ier)
         iertmp=max(iertmp,ier)
      endif
      
      if(present(r8vals)) then
         call xplasma_getList(s,id,r8vals,ier)
         iertmp=max(iertmp,ier)
      endif
      
      if(present(chvals)) then
         call xplasma_getList(s,id,chvals,ier)
         iertmp=max(iertmp,ier)
      endif

      ier=iertmp

    end subroutine xplasma_getList_names

    !========================================
    subroutine xplasma_setList_r4vals(s,id,r4vals,ier, nchange)
      !
      !  set floating point values in list (r4 precision)
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      real, dimension(:), intent(in) :: r4vals  ! the floating point values
      integer, intent(out) :: ier         ! completion code

      integer, intent(out), optional :: nchange  ! number of changed values

      !----------------
      real*8 :: r8vals(size(r4vals))
      !----------------

      r8vals = r4vals
      call xplasma_setList_r8vals(s,id,r8vals,ier, nchange)

    end subroutine xplasma_setList_r4vals

    !========================================
    subroutine xplasma_getList_r4vals(s,id,r4vals,ier)
      !
      !  get floating point values in list (r4 precision)
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      real, dimension(:), intent(out) :: r4vals ! the floating point values
      integer, intent(out) :: ier         ! completion code

      !----------------
      real*8 :: r8vals(size(r4vals))
      integer :: nelems
      !----------------

      nelems = size(r4vals)

      call xplasma_getList_r8vals(s,id,r8vals,ier)
      if(ier.eq.0) then
         r4vals(1:nelems) = r8vals(1:nelems)
      else
         r4vals(1:nelems) = 0
      endif

    end subroutine xplasma_getList_r4vals

    !========================================
    subroutine xplasma_setList_ivals(s,id,ivals,ier, nchange)
      !
      !  set integer values in list
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      integer, dimension(:), intent(in) :: ivals   ! the integer values
      integer, intent(out) :: ier         ! completion code

      integer, intent(out), optional :: nchange  ! number of changed values

      !----------------
      integer :: iadr,nelems,ielem,iiv
      !----------------

      nelems = size(ivals)

      call xplasma_errck0_list_adr(s,id,nelems,iadr,ier)
      if(ier.ne.0) return

      call xplasma_author_check(s,id,ier)

      if(ier.eq.0) then
         if(present(nchange)) then
            nchange = 0
            do ielem=1,nelems
               iiv=s%lists(iadr)%ivals(ielem)
               if(iiv.ne.ivals(ielem)) nchange = nchange + 1
            enddo
         endif
         s%lists(iadr)%ivals(1:nelems) = ivals(1:nelems)
      endif

    end subroutine xplasma_setList_ivals

    !========================================
    subroutine xplasma_getList_ivals(s,id,ivals,ier)
      !
      !  get integer values in list (r4 precision)
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      integer, dimension(:), intent(out) :: ivals  ! the integer values
      integer, intent(out) :: ier         ! completion code

      !----------------
      integer :: iadr,nelems
      !----------------

      nelems=size(ivals)

      call xplasma_errck0_list_adr(s,id,nelems,iadr,ier)
      if(ier.eq.0) then
         ivals(1:nelems) = s%lists(iadr)%ivals(1:nelems)
      else
         ivals(1:nelems) = 0
      endif

    end subroutine xplasma_getList_ivals

    !========================================
    subroutine xplasma_setList_r8vals(s,id,r8vals,ier, nchange)
      !
      !  set real*8 values in list
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      real*8, dimension(:), intent(in) :: r8vals  ! the real*8 values
      integer, intent(out) :: ier         ! completion code

      integer, intent(out), optional :: nchange  ! number of changed values

      !----------------
      integer :: iadr,nelems,ielem
      real*8 :: zv
      !----------------

      nelems = size(r8vals)

      call xplasma_errck0_list_adr(s,id,nelems,iadr,ier)
      if(ier.ne.0) return

      call xplasma_author_check(s,id,ier)

      if(ier.eq.0) then
         if(present(nchange)) then
            nchange = 0
            do ielem=1,nelems
               zv=s%lists(iadr)%r8vals(ielem)
               if(zv.ne.r8vals(ielem)) nchange = nchange + 1
            enddo
         endif
         s%lists(iadr)%r8vals(1:nelems) = r8vals(1:nelems)
      endif

    end subroutine xplasma_setList_r8vals

    !========================================
    subroutine xplasma_getList_r8vals(s,id,r8vals,ier)
      !
      !  get real*8 values in list
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      real*8, dimension(:), intent(out) :: r8vals ! the real*8 values
      integer, intent(out) :: ier         ! completion code

      !----------------
      integer :: iadr,nelems
      !----------------

      nelems=size(r8vals)

      call xplasma_errck0_list_adr(s,id,nelems,iadr,ier)
      if(ier.eq.0) then
         r8vals(1:nelems) = s%lists(iadr)%r8vals(1:nelems)
      else
         r8vals(1:nelems) = 0
      endif

    end subroutine xplasma_getList_r8vals

    !========================================
    subroutine xplasma_setList_chvals(s,id,chvals,ier, nchange)
      !
      !  set character string values in list
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      character*(*), dimension(:), intent(in) :: chvals 
                                          ! the character string values

      integer, intent(out) :: ier         ! completion code

      integer, intent(out), optional :: nchange  ! number of changed values

      !----------------------------
      integer :: iadr,i,ic,ic1,ic2,ici,ilen,ilena,nelems,ielem,ichange
      character*1 cch
      !----------------------------

      if(present(nchange)) nchange = -1   ! error value

      nelems = size(chvals)

      call xplasma_errck0_list_adr(s,id,nelems,iadr,ier)
      if(ier.ne.0) return

      call xplasma_author_check(s,id,ier)

      if(ier.eq.0) then

         if(.not.associated(s%lists(iadr)%chbuf)) then
            ichange=nelems  ! no prior character values so all are changes
         else
            ichange = 0
            do ielem = 1,nelems
               ilen=len(trim(chvals(ielem)))
               ilena=max(1,ilen)
               ic1=s%lists(iadr)%chindx(1,ielem)
               ic2=s%lists(iadr)%chindx(2,ielem)
               if((ic2-ic1+1).ne.ilena) then
                  ichange = ichange + 1  ! length mismatch
               else if(ilen.eq.0) then
                  cch = s%lists(iadr)%chbuf(ic1)
                  if(cch.ne.' ') then
                     ichange = ichange + 1  ! empty value mismatch
                  endif
               else
                  i=0
                  do ic=ic1,ic2
                     i=i+1
                     cch = s%lists(iadr)%chbuf(ic)
                     if(cch.ne.chvals(ielem)(i:i)) then
                        ichange = ichange + 1  ! value mismatch
                        exit
                     endif
                  enddo
               endif
            enddo
         endif

         if(present(nchange)) nchange = ichange

         if(ichange.gt.0) then
            if(associated(s%lists(iadr)%chbuf)) deallocate(s%lists(iadr)%chbuf)
            ic=0
            do i=1,nelems
               ilen=max(1,len(trim(chvals(i))))
               s%lists(iadr)%chindx(1,i)=ic+1
               s%lists(iadr)%chindx(2,i)=ic+ilen
               ic=ic+ilen
            enddo

            allocate(s%lists(iadr)%chbuf(ic))

            do i=1,nelems
               ic1=s%lists(iadr)%chindx(1,i)
               ic2=s%lists(iadr)%chindx(2,i)
               ilen=len(trim(chvals(i)))
               if(ilen.gt.0) then
                  do ic=ic1,ic2
                     ici=ic-ic1+1
                     s%lists(iadr)%chbuf(ic) = chvals(i)(ici:ici)
                  enddo
               else
                  s%lists(iadr)%chbuf(ic1:ic2)=' '
               endif
            enddo
         endif
      endif

    end subroutine xplasma_setList_chvals

    !========================================
    subroutine xplasma_getList_chvals(s,id,chvals,ier)
      !
      !  get character string values in list
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of list item
      character*(*), dimension(:), intent(out) :: chvals ! the character string values
      integer, intent(out) :: ier         ! completion code

      !--------------------------------------
      integer :: iadr,i,ic,ic1,ic2,ici,ilen,nelems
      !--------------------------------------

      nelems = size(chvals)

      chvals = ' '

      call xplasma_errck0_list_adr(s,id,nelems,iadr,ier)
      if(ier.ne.0) return

      if(.not.associated(s%lists(iadr)%chindx)) return  ! OK -- all blank

      do i=1,nelems
         ic1=s%lists(iadr)%chindx(1,i)
         ic2=s%lists(iadr)%chindx(2,i)
         do ic=ic1,ic2
            ici=ic-ic1+1
            chvals(i)(ici:ici) = s%lists(iadr)%chbuf(ic)
         enddo
      enddo

    end subroutine xplasma_getList_chvals

    !========================================
    !  Public methods dealing with coordinates
    !========================================

    subroutine xplasma_create_coord(s,coordname,periodic,id,ier, &
      label,units)
    
      !  create a coordinate with the specified name; return new coord id

      type (xplasma), pointer :: s
      character*(*), intent(in) :: coordname
      logical, intent(in) :: periodic  ! T if named coordinate is periodic
      integer, intent(out) :: id  ! id of coordinate object
      integer, intent(out) :: ier ! completion code

      character*(*), intent(in), optional :: label
      character*(*), intent(in), optional :: units

      !------------------------------------------------
      integer :: iertmp
      !------------------------------------------------

      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_create_coord call)',ier)
      endif
      if(ier.ne.0) return

      id=0
      call xplasma_add_item(s,coordname,xplasma_coordType,id,ier)
      if(ier.ne.0) return

      call xplasma_define_coord(s,id,periodic,ier)
      if(ier.ne.0) return

      iertmp=0

      if(present(label)) then
         call xplasma_label_item(s,id,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_label_item(s,id,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      ier = iertmp

    end subroutine xplasma_create_coord

    !========================================
    subroutine xplasma_define_coord(s,id,periodic,ier)

      !  **private** 
      !  define a new coordinate, at the specified id

      type (xplasma), pointer :: s
      integer, intent(in) :: id   ! id of coordinate object
      logical, intent(in) :: periodic  ! T if named coordinate is periodic
      integer, intent(out) :: ier ! completion code

      !------------------------------------------------
      integer :: iadr
      !------------------------------------------------

      call xplasma_errck0_coord(s,id,ier)
      if(ier.ne.0) return

      call xplasma_author_check(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex

      call xpcoord_free(s%coords(iadr))
      if(ier.ne.0) return

      s%coords(iadr)%periodic = periodic

    end subroutine xplasma_define_coord

    !========================================
    subroutine xplasma_coord_info(s,id,ier, &
         ngrids,periodic,xmin,xmax,name,label,units,author)

      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! coordinate object id
      integer, intent(out) :: ier   ! completion code, 0=OK

      integer, intent(out), optional :: ngrids    ! # of grid discretizations
      logical, intent(out), optional :: periodic  ! periodicity flag
      real*8, intent(out), optional :: xmin       ! min value
      real*8, intent(out), optional :: xmax       ! max value

      character*(*), intent(out), optional :: name
      character*(*), intent(out), optional :: label
      character*(*), intent(out), optional :: units
      character*(*), intent(out), optional :: author

      !--------------------
      real*8 :: zmin,zmax
      integer :: iadr,iertmp
      !--------------------

      if(present(ngrids)) ngrids=0
      if(present(periodic)) periodic=.FALSE.
      if(present(xmin)) xmin=0
      if(present(xmax)) xmax=0

      if(present(name)) name=' '
      if(present(label)) label=' '
      if(present(units)) units=' '
      if(present(author)) author=' '

      call xplasma_errck0_coord(s,id,ier)
      if(ier.ne.0) return

      iertmp=0

      if(present(ngrids)) then
         iadr = s%dict(id)%dindex
         ngrids = s%coords(iadr)%ngridc
      endif

      if(present(periodic)) then
         call xplasma_coord_isPeriodic(s,id,periodic,ier)
         iertmp=max(iertmp,ier)
      endif

      if(present(xmin).or.present(xmax)) then
         call xplasma_coord_minmax(s,id,zmin,zmax,ier)
         iertmp=max(iertmp,ier)

         if(present(xmin)) then
            if(ier.eq.0) then
               xmin=zmin
            else
               xmin=0
            endif
         endif
            
         if(present(xmax)) then
            if(ier.eq.0) then
               xmax=zmax
            else
               xmax=0
            endif
         endif

      endif

      if(present(name)) then
         call xplasma_get_item_info(s,id,ier, name=name)
         iertmp=max(iertmp,ier)
      endif

      if(present(label)) then
         call xplasma_get_item_info(s,id,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_get_item_info(s,id,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      if(present(author)) then
         call xplasma_get_item_info(s,id,ier, author=author)
         iertmp=max(iertmp,ier)
      endif

      ier=iertmp

    end subroutine xplasma_coord_info

    !========================================
    subroutine xplasma_coord_gridId(s,id,indx,idgrid,ier)

      !  get grid id from coordinate grid list index

      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! coordinate object id
      integer, intent(in) :: indx   ! which grid is desired
      integer, intent(out) :: idgrid  ! grid id of desired grid
      integer, intent(out) :: ier   ! completion code, 0=OK

      !------------------------------------------------
      integer :: iadr
      character*128 msgbuf
      !------------------------------------------------

      call xplasma_errck0_coord(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex

      if((indx.le.0).or.(indx.gt.s%coords(iadr)%ngridc)) then
         write(msgbuf,*) ' grid index argument value: ',indx
         call xplasma_errmsg_append(s,msgbuf)
         ier=402
         return
      endif

      idgrid = s%coords(iadr)%grid_ids(indx)

    end subroutine xplasma_coord_gridId

    !========================================
    subroutine xplasma_coord_isPeriodic(s,id,periodic,ier)

      !  get information whether coordinate object is a periodic angle coord.

      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! coordinate object id
      logical, intent(out) :: periodic ! T if this is a periodic coordinate
      integer, intent(out) :: ier   ! completion code, 0=OK

      !------------------------------------------------
      integer :: iadr,itype,icoord
      !------------------------------------------------

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      itype = s%dict(id)%dtype
      
      if(itype.eq.xplasma_coordType) then
         periodic = s%coords(iadr)%periodic
      else if(itype.eq.xplasma_gridType) then
         icoord = s%grids(iadr)%coord
         iadr = s%dict(icoord)%dindex
         periodic = s%coords(iadr)%periodic
      else 
         call xplasma_errmsg_append(s, &
              'item name: "'//trim(s%dict(id)%name)//'".')
         ier=400
      endif

    end subroutine xplasma_coord_isPeriodic

    !========================================
    subroutine xplasma_coord_minmax(s,id,xmin,xmax,ier)

      ! **private**
      ! return coordinate min & max
      ! always [0,2pi] for angle coordinates
      ! always [0,1] for rho

      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! coordinate object id
      real*8, intent(out) :: xmin,xmax  ! range returned
      integer, intent(out) :: ier   ! completion code, 0=OK

      !------------------------------------------------
      integer :: iadr,iadg,inum,inumx,idgrid
      logical :: periodic
      !------------------------------------------------

      call xplasma_errck0_coord(s,id,ier)
      if(ier.ne.0) return

      if(id.eq.xplasma_rho_coord) then
         xmin=0
         xmax=1
         return
      endif

      iadr = s%dict(id)%dindex
      periodic = s%coords(iadr)%periodic
      
      if(periodic) then
         xmin=0
         xmax=c2pi
         return
      endif

      inum = s%coords(iadr)%ngridc
      if(inum.eq.0) then
         ier=403
      else
         !  all grids must cover the same range; so take the 1st one...
         idgrid = s%coords(iadr)%grid_ids(1)
         iadg = s%dict(idgrid)%dindex
         inumx = s%grids(iadg)%size
         xmin = s%grids(iadg)%xpkg(1,1)
         xmax = s%grids(iadg)%xpkg(inumx,1)
      endif

    end subroutine xplasma_coord_minmax

    !========================================
    !  Public methods dealing with grids...
    !========================================

    subroutine xplasma_create_grid(s,gridname,icoord,x,id,ier, &
         ccwflag,label)
    
      !  create a new grid with the specified name; return new grid id

      type (xplasma), pointer :: s
      character*(*), intent(in) :: gridname
      integer, intent(in) :: icoord  ! id of coordinate object to which the
                                     ! grid object belongs
      !  icoord can also be a grid id, in which case the coordinate associated
      !  with this new grid is the one associated with grid (icoord)

      real*8, intent(in) :: x(:)  ! the grid, x(j).lt.x(j+1) required.
      integer, intent(out) :: id  ! id of grid object
      integer, intent(out) :: ier ! completion code

      !  specify if toroidal or poloidal angle grid is drawn counterclockwise:
      logical, intent(in), OPTIONAL :: ccwflag
      !    theta: when viewing plasma cross section to right of machine
      !           centerline
      !    phi:   when viewing machine from above
      !    Default: .TRUE.

      character*(*), intent(in), optional :: label
      !  units label inherited from associated coordinate object

      !  NOTE: xplasma saves all angle grids with CCW orientation; a CW
      !  oriented grid will be redefined with its direction reversed.

      !------------------------------------------------
      integer :: iersave,jcoord,isize
      logical :: iccw
      integer :: nx  ! = size(x)
      !------------------------------------------------

      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_create_grid call)',ier)
      endif
      if(ier.ne.0) return

      id=0
      call xplasma_add_item(s,gridname,xplasma_gridType,id,ier)
      if(ier.ne.0) return

      nx=size(x)
      if(nx.lt.2) then
         ier=302
         call xplasma_errmsg_append(s, &
              ' ?xplasma_create_grid: grid "'//trim(gridname)// &
              ' has fewer than 2 points; 2 is the minimum.')
         return
      endif

      iccw=.TRUE.
      if(present(ccwflag)) iccw = ccwflag

      call xplasma_define_grid(s,id,icoord,nx,x,ier, iccw)
      if(ier.ne.0) return

      call xplasma_grid_info(s,id,ier, coord=jcoord)
      if(ier.ne.0) return

      if(allocated(s%dict(jcoord)%units)) then
         isize=size(s%dict(jcoord)%units)
         if(allocated(s%dict(id)%units)) deallocate(s%dict(id)%units)
         allocate(s%dict(id)%units(isize))
         s%dict(id)%units = s%dict(jcoord)%units
      endif

      if(present(label)) then
         call xplasma_label_item(s,id,ier, label=label)
      else
         if(allocated(s%dict(jcoord)%label)) then
            isize=size(s%dict(jcoord)%label)
            if(allocated(s%dict(id)%label)) deallocate(s%dict(id)%label)
            allocate(s%dict(id)%label(isize))
            s%dict(id)%label = s%dict(jcoord)%label
         endif
      endif

      return
    end subroutine xplasma_create_grid

    !========================================
    subroutine xplasma_define_grid(s,id,icoord,nx,x,ier, &
         ccwflag)

      !  define a new grid, at the specified grid id

      type (xplasma), pointer :: s
      integer, intent(in) :: id   ! id of grid object
      integer, intent(in) :: icoord  ! id of coordinate object to which the
                                     ! grid object belongs
      !  icoord can also be a grid id, in which case the coordinate associated
      !  with this new grid is the one associated with grid (icoord)

      integer, intent(in) :: nx   ! no. of grid points
      real*8, intent(in) :: x(nx) ! the grid, x(j).lt.x(j+1) required.
      integer, intent(out) :: ier ! completion code

      !  specify if toroidal or poloidal angle grid is drawn counterclockwise:
      logical, intent(in), OPTIONAL :: ccwflag
      !    theta: when viewing plasma cross section to right of machine
      !           centerline
      !    phi:   when viewing machine from above
      !    Default: .TRUE.

      !  NOTE: xplasma saves all angle grids with CCW orientation; a CW
      !  oriented grid will be redefined with its direction reversed.

      !------------------------------------------------
      integer :: jcoord,iadr,iadc,iadg,itype,i,ii,iper
      integer :: imsg,itol,ialg,igrids,igridp,inx,imax
      integer, dimension(:), allocatable :: inums
      logical :: iccw,imatch
      real*8 :: xtmp(nx),xmin,xmax,ztol,ztol0,ztol2,zx1,zx2,znorm

      logical :: l_tmp

      !------------------------------------------------
      !  check if id is valid
      call xplasma_errck0_grid(s,id,ier)
      if(ier.ne.0) return

      call xplasma_author_check(s,id,ier)
      if(ier.ne.0) return

      !--------------------------
      !  check for optional argument

      iccw=.TRUE.
      if(present(ccwflag)) iccw=ccwflag

      !--------------------------
      !  check if icoord is a grid id-- in which case, associate this new
      !  grid with the same coordinate as the indicated old grid.

      call xplasma_errck0(s,icoord,ier)
      if(ier.eq.0) then
         itype=s%dict(icoord)%dtype
         if((itype.ne.xplasma_gridType).and.(itype.ne.xplasma_coordType)) then
            call xplasma_errck0_coord(s,icoord,ier)  ! gen error
         endif
      endif

      !  if neither, it is an error...
      if(ier.ne.0) return

      !--------------------------
      !  get coordinate data item index
      if(itype.eq.xplasma_coordType) then
         jcoord=icoord
      else
         jcoord=s%dict(icoord)%dindex  ! grid object index
         jcoord=s%grids(jcoord)%coord  ! corresponding coordinate id
      endif

      !--------------------------
      !  "jcoord" points to coordinate of which current grid is a 
      !  particular discretization...

      iadc = s%dict(jcoord)%dindex     ! adr of coordinate data object

      iadr = s%dict(id)%dindex         ! adr of grid data object

      call xpgrid_free(s%grids(iadr))

      !  check size and ordering of x

      if(nx.lt.2) then
         ier=302
         return
      endif

      do i=1,nx-1
         if(x(i+1).le.x(i)) then
            ier=303
            exit
         endif
      enddo

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              'name of grid: "'//trim(s%dict(id)%name)//'".')
         return
      endif

      !--------------------------
      !  range check...
      
      l_tmp = s%coords(iadc)%periodic 
      if(jcoord.eq.xplasma_rho_coord) then

         !  rho: [0,1]

         if((abs(x(1)).gt.s%bdytol).or.(abs(x(nx)-1).gt.s%bdytol)) then
            ier=306
         endif

      else if( l_tmp ) then

         !  angle coord: [x,x+2pi], x arbitrary but of order unity...

         if(abs(x(nx)-x(1)-c2pi).gt.s%bdytol) then
            ier=307
         endif

      else

         !  other coordinates: all discretizations must cover same range:
         !  check this...

         if(s%coords(iadc)%ngridc.gt.0) then
            iadg=s%coords(iadc)%grid_ids(1)   ! grid id
            iadg=s%dict(iadg)%dindex         ! => grid data item address
            inx =s%grids(iadg)%size
            xmin=s%grids(iadg)%xpkg(1,1)
            xmax=s%grids(iadg)%xpkg(inx,1)

            znorm = max(abs(xmin),abs(xmax)) ! at least one is guaranteed > 0

            if((abs(xmin-x(1)).gt.ceps10*znorm).or. &
                 (abs(xmax-x(nx)).gt.ceps10*znorm)) then
               ier=308
            endif

         else
            continue  ! this is the first grid -- it sets the range.
         endif

      endif

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              'name of grid: "'//trim(s%dict(id)%name)//'".')
         return
      endif

      !--------------------------
      !  looks OK...

      s%grids(iadr)%size = nx
      s%grids(iadr)%coord = jcoord
      s%grids(iadr)%nrefs = 0

      !--------------------------
      !  deal with possible non-CCW angle coordinate

      if(.not.s%coords(iadc)%periodic) then
         iper=0
         iccw=.TRUE.
      else
         iper=1
      endif

      !  non-ccw periodic coordinate
      if(iccw) then
         xtmp = x
      else
         xmin=x(1)
         xmax=x(nx)
         do i=2,nx-1
            xtmp(i)=xmin+(xmax-x(nx+1-i))
         enddo
      endif

      !--------------------------
      !  create XPKG for interpolation (pspline library call)

      allocate(s%grids(iadr)%xpkg(nx,4))

      s%grids(iadr)%xpkg=czero !MG, I added this line to please Intel compiler
      imsg=0
      itol=0
      ztol=itol
      ialg=3

      call r8genxpkg(nx,xtmp,s%grids(iadr)%xpkg, &
           iper,imsg,itol,ztol,ialg,ier)

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              'name of grid: "'//trim(s%dict(id)%name)//'".')
         ier=304            ! this shouldn't happen -- probable xplasma bug...
      endif

      if(ier.ne.0) return

      !--------------------------
      !  record grid as belonging to coordinate (iadc)

      igrids=s%coords(iadc)%ngridc

      call xplasma_coord_expand(s%coords(iadc),igrids+1)

      igrids=igrids+1
      s%coords(iadc)%ngridc = igrids
      s%coords(iadc)%grid_ids(igrids) = id

      !--------------------------
      !  check if grid is equivalent to any other

      !  first compress equiv list numbers

      !  find the equiv numbers in use
      allocate(inums(s%ngrids)); inums=0
      do i=1,s%ngrids
         ii=s%grid_equiv(i)
         if(ii.gt.0) then
            if(ii.gt.s%ngrids) then
               call xplasma_errmsg_append(s, &
                    'xplasma internal error: grid_equiv overflow.')
               ier=999
               exit
            else
               inums(ii)=ii
            endif
         endif
      enddo
      if(ier.ne.0) return

      !  "compression" map for equiv numbers
      ii=0
      do i=1,s%ngrids
         if(inums(i).gt.0) then
            ii=ii+1
            inums(i)=ii
         endif
      enddo
      imax=ii  ! save for possible later use

      !  apply map
      do i=1,s%ngrids
         ii=s%grid_equiv(i)
         if(ii.gt.0) then
            s%grid_equiv(i)=inums(ii)
         endif
      enddo
      deallocate(inums)

      !  now see if current grid (at index iadr) matches any other

      ztol0=1.0d-12*max(abs(x(1)),abs(x(nx)))

      imatch=.FALSE.
      do i=1,s%ngrids
         if(i.eq.iadr) cycle
         ii=s%grids(i)%size
         if(ii.eq.0) cycle
         if(ii.ne.nx) cycle

         ! size match found

         zx1=s%grids(i)%xpkg(1,1)
         zx2=s%grids(i)%xpkg(nx,1)
         ztol2=1.0d-12*max(abs(zx1),abs(zx2))
         ztol=min(ztol0,ztol2)

         imatch=.TRUE.
         do ii=1,nx
            zx1=x(ii)
            zx2=s%grids(i)%xpkg(ii,1)
            if(abs(zx1-zx2).gt.ztol) then
               imatch=.FALSE.
               exit
            endif
         enddo

         if(.not.imatch) cycle

         ! full match found

         s%grid_equiv(iadr) = s%grid_equiv(i)
         exit
      enddo

      if(.not.imatch) then
         s%grid_equiv(iadr) = imax+1
      endif

    end subroutine xplasma_define_grid

    subroutine xplasma_coord_expand(d,ineed)
      !  private: expand grid_ids list in 

      type (xpcoord) :: d
      integer, intent(in) :: ineed
      integer, dimension(:), pointer :: itmp

      integer :: igridn,igridp

      igridn=d%ngrid_max
      if(ineed.gt.igridn) then
         igridp=igridn
         if(igridp.gt.0) then
            itmp => d%grid_ids
         endif
         do
            igridn=max(4,2*igridn)
            if(igridn.ge.ineed) exit
         enddo
         allocate(d%grid_ids(igridn))
         d%grid_ids(igridp+1:igridn)=0
         if(igridp.gt.0) then
            d%grid_ids(1:igridp)=itmp(1:igridp)
            deallocate(itmp)
         endif
         d%ngrid_max = igridn
      endif
    end subroutine xplasma_coord_expand

    !========================================
    subroutine xplasma_grid_size(s,id,nx,ier)

      !  get grid size

      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! grid object id
      integer, intent(out) :: nx    ! grid size
      integer, intent(out) :: ier   ! completion code, 0=OK

      !------------------------------------------------
      integer :: iadr
      !------------------------------------------------

      call xplasma_errck0_grid(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex

      nx = s%grids(iadr)%size

    end subroutine xplasma_grid_size

    !========================================
    subroutine xplasma_grid(s,id,x,ier, ccwflag)

      !  get grid

      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! grid id object
      real*8, intent(out) :: x(:)   ! the grid returned
      integer, intent(out) :: ier   ! completion code, 0=OK

      logical, intent(in), optional :: ccwflag  ! .FALSE. to reverse grid
      ! (default is .TRUE.) (ignored unless grid is periodic)

      !------------------------------------------------
      integer :: iadr,inx,ix,nx,iertmp
      logical :: jccwflag,periodic
      real*8 :: xmin,xmax
      !------------------------------------------------

      call xplasma_errck0_grid(s,id,ier)
      if(ier.ne.0) return

      nx = size(x)

      iadr = s%dict(id)%dindex
      inx = s%grids(iadr)%size
      if(inx.ne.nx) then
         call xplasma_errmsg_append(s, &
              'name of grid: "'//trim(s%dict(id)%name)//'".')
         ier=305
         return
      endif

      jccwflag=.TRUE.
      if(present(ccwflag)) then
         call xplasma_grid_info(s,id,iertmp, perio=periodic)
         if(periodic) jccwflag = ccwflag
      endif

      if(jccwflag) then
         x(1:nx)=s%grids(iadr)%xpkg(1:nx,1)
      else
         xmin=s%grids(iadr)%xpkg(1,1)
         xmax=s%grids(iadr)%xpkg(nx,1)
         do ix=1,nx
            x(ix) = xmin + (xmax - s%grids(iadr)%xpkg(nx+1-ix,1))
         enddo
      endif

    end subroutine xplasma_grid

    !========================================
    subroutine xplasma_remove_grid(s,id,ier)

      ! remove grid AND all profiles which use the grid

      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! grid id object
      integer, intent(out) :: ier   ! completion code, 0=OK

      !------------------------------------------------
      integer :: iertmp,ii,inum,idg1,idg2,idg3
      integer, dimension(:), allocatable :: ids
      !------------------------------------------------

      call xplasma_errck0_grid(s,id,ier)  ! check have valid grid ID
      if(ier.ne.0) return

      call xplasma_author_check(s,id,ier) ! check that author code matches
      if(ier.ne.0) return

      !  OK; set author to root, so all xplasma_remove_item actions work.

      call xplasma_author_set(s, xplasma_root, iertmp)
      ier=max(ier,iertmp)

      !  find all profiles...

      call xplasma_num_items(s,iertmp, num_profs=inum)
      ier=max(ier,iertmp)
      allocate(ids(inum))

      call xplasma_find_profs(s,iertmp, id_profs=ids)
      ier=max(ier,iertmp)

      !  remove profiles that use the specified grid...

      do ii=1,inum
         call xplasma_prof_info(s,ids(ii),iertmp, &
              gridId1=idg1,gridId2=idg2,gridId3=idg3)
         ier=max(ier,iertmp)
         
         if((idg1.eq.id).or.(idg2.eq.id).or.(idg3.eq.id)) then
            call xplasma_remove_item(s,ids(ii),iertmp)
            ier=max(ier,iertmp)
         endif
      enddo

      !  remove the grid itself

      call xplasma_remove_item(s,id,iertmp)
      ier=max(ier,iertmp)

      call xplasma_author_clear(s, xplasma_root, iertmp)
      ier=max(ier,iertmp)

    end subroutine xplasma_remove_grid

    !========================================
    subroutine xplasma_grid_info(s,id,ier, &
         size,xmin,xmax,perio,coord,nrefs,name,label,units,author)

      !  get grid reference count (# of profile items using this grid)

      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! grid id object
      integer, intent(out) :: ier   ! completion code, 0=OK

      integer, intent(out), optional :: size  ! size of grid
      real*8, intent(out), optional :: xmin   ! minimum grid value
      real*8, intent(out), optional :: xmax   ! maximum grid value
      logical, intent(out), optional :: perio ! T for periodic grids
      integer, intent(out), optional :: coord ! coordinate associated w/grid
      integer, intent(out), optional :: nrefs ! no. of references to grid

      character*(*), intent(out), optional :: name
      character*(*), intent(out), optional :: label
      character*(*), intent(out), optional :: units
      character*(*), intent(out), optional :: author

      !------------------------------------------------
      integer :: iadr,iertmp,isize,icoord,iadc
      !------------------------------------------------

      if(present(size)) size=0
      if(present(coord)) coord=0
      if(present(nrefs)) nrefs=0
      if(present(xmin)) xmin=0
      if(present(xmax)) xmax=0
      if(present(perio)) perio=.FALSE.

      if(present(name)) name=' '
      if(present(label)) label=' '
      if(present(units)) units=' '
      if(present(author)) author=' '

      call xplasma_errck0_grid(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex

      isize = s%grids(iadr)%size
      if(present(size)) then
         size = isize
      endif

      icoord = s%grids(iadr)%coord
      if(present(coord)) then
         coord = icoord
      endif

      if(present(xmin)) then
         xmin = s%grids(iadr)%xpkg(1,1)
      endif

      if(present(xmax)) then
         xmax = s%grids(iadr)%xpkg(isize,1)
      endif

      if(present(perio)) then
         iadc = s%dict(icoord)%dindex
         perio = s%coords(iadc)%periodic
      endif

      if(present(nrefs)) then
         nrefs = s%grids(iadr)%nrefs
      endif

      iertmp=0

      if(present(name)) then
         call xplasma_get_item_info(s,id,ier, name=name)
         iertmp=max(iertmp,ier)
      endif

      if(present(label)) then
         call xplasma_get_item_info(s,id,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_get_item_info(s,id,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      if(present(author)) then
         call xplasma_get_item_info(s,id,ier, author=author)
         iertmp=max(iertmp,ier)
      endif

      ier=iertmp

    end subroutine xplasma_grid_info

    !========================================
    !  Public methods dealing with profiles...
    !========================================

    subroutine xplasma_prof_info(s,id,ier, &
         rank,splineType,hybridType,gridId1,gridId2,gridId3, &
         profId1,gridType,counter, &
         name,label,units,author)

      !  base argauments in/out
      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! profile object id
      integer, intent(out) :: ier   ! completion code, 0=OK

      !---------------------------------------------------
      !  select desired information with optional arguments...

      integer, intent(out), optional :: rank  ! rank

      integer, intent(out), optional :: splineType  ! 0 for piecewise linear,
      ! 1 for C1 Hermite, 2 for C2 cubic spline

      integer, dimension(:), intent(out), optional :: hybridType ! type code
      ! for interpolation along each dimension: 0 for piecewise linear,
      ! 1 for C1 Hermite, 2 for C2 cubic spline.  For non-hybrid interpolants
      ! the same value is returned for each element up to the rank of the
      ! profile object.  If this argument is present, its dimensioning must
      ! match or exceed the profile rank.

      integer, intent(out), optional :: gridId1  ! id of 1st grid
      integer, intent(out), optional :: gridId2  ! id of 2nd grid (or zero)
      integer, intent(out), optional :: gridId3  ! id of 3rd grid (or zero)

      integer, intent(out), optional :: profId1  ! id of associated profile
      !  (sometimes f(R,Z) has an associated f(rho,theta) profile available,
      !  and vice versa).

      integer, intent(out), optional :: gridType ! generic grid type code
      ! 0=unknown or mixed
      ! 1=flux coords (rho) or (rho,theta) [phi someday]
      ! 2=cylindric coordinates (R,Z) [phi someday]

      integer, intent(out), optional :: counter  ! profile counter
      ! each new version of this profile gets a new counter value.

      character*(*), intent(out), optional :: name
      character*(*), intent(out), optional :: label
      character*(*), intent(out), optional :: units
      character*(*), intent(out), optional :: author

      !------------------------------------------------
      integer :: iadr,iertmp,i,igid,iadg,icoord,itype,irank,ispline
      character*128 msgbuf
      !------------------------------------------------

      if(present(rank)) rank = 0
      if(present(splineType)) splineType = 0
      if(present(hybridType)) hybridType = 0
      if(present(gridId1)) gridId1 = 0
      if(present(gridId2)) gridId2 = 0
      if(present(gridId3)) gridId3 = 0
      if(present(profId1)) profId1 = 0
      if(present(gridType)) gridType = 0
      if(present(counter)) counter = 0

      if(present(name)) name=' '
      if(present(label)) label=' '
      if(present(units)) units=' '
      if(present(author)) author=' '

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex

      irank = s%profs(iadr)%rank

      if(present(rank)) then
         rank = irank
      endif

      if(present(splineType)) then
         splineType = s%profs(iadr)%kspline
      endif

      if(present(hybridType)) then
         if(size(hybridType).lt.irank) then
            ier=512
            msgbuf=' '
            write(msgbuf,*) &
                 ' ?xplasma_prof_info: dimension of optional argument '// &
                 ' "hybridType": ',size(hybridType)
            call xplasma_errmsg_append(s,msgbuf)
            msgbuf=' '
            write(msgbuf,*) &
                 '  rank of referenced profile (',id,'): ',irank
            call xplasma_errmsg_append(s,msgbuf)
            return
         endif

         ispline = s%profs(iadr)%kspline
         call spline_type_decode(ispline,hybridType(1:irank),irank)
      endif

      if(present(counter)) then
         counter = s%profs(iadr)%prof_counter
      endif

      if(present(gridId1)) then
         gridId1 = s%profs(iadr)%gridIds(1)
      endif

      if(present(gridId2)) then
         gridId2 = s%profs(iadr)%gridIds(2)
      endif

      if(present(gridId3)) then
         gridId3 = s%profs(iadr)%gridIds(3)
      endif

      if(present(profId1)) then
         profId1 = s%profs(iadr)%profIds(1)
      endif

      if(present(gridType)) then
         do i=1,3
            igid=s%profs(iadr)%gridIds(i)
            if(igid.eq.0) cycle

            iadg = s%dict(igid)%dindex
            icoord = s%grids(iadg)%coord

            itype = 0

            if(icoord.eq.xplasma_rho_coord) itype=xplasma_flux_coords
            if(icoord.eq.xplasma_theta_coord) itype=xplasma_flux_coords
            if(icoord.eq.xplasma_rhox_coord) itype=xplasma_flux_coords
            if(icoord.eq.xplasma_thx_coord) itype=xplasma_flux_coords

            if(icoord.eq.xplasma_R_coord) itype=xplasma_cyl_coords
            if(icoord.eq.xplasma_Z_coord) itype=xplasma_cyl_coords

            if(icoord.eq.xplasma_B_coord) then
               gridType = 0   ! answer is: neither!
               exit
            endif

            if(gridType.eq.0) then
               gridType = itype

            else if(gridType.ne.itype) then
               gridType = 0   ! a mix?
               exit
            endif

         enddo
      endif


      iertmp = 0

      if(present(name)) then
         call xplasma_get_item_info(s,id,ier, name=name)
         iertmp=max(iertmp,ier)
      endif

      if(present(label)) then
         call xplasma_get_item_info(s,id,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_get_item_info(s,id,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      if(present(author)) then
         call xplasma_get_item_info(s,id,ier, author=author)
         iertmp=max(iertmp,ier)
      endif

      ier=iertmp

    end subroutine xplasma_prof_info

    !==========================================================
    subroutine xplasma_prof_gridInfo(s,id,icoord,idgrid,ierr, &
         id_actual)

      type (xplasma), pointer :: s
      integer, intent(in) :: id     ! id of profile about which info is sought
      integer, intent(in) :: icoord ! coordinate sought

      integer, intent(out) :: idgrid   ! grid over this coordinate, or zero
      !  zero means profile identified by "id" is not a function of the
      !  specified coordinate; this is not considered an error.

      integer, intent(out) :: ierr  ! error status code returned, 0=OK
      !  errors: profile id is invalid; coordinate id is invalid...

      integer, intent(out), optional :: id_actual  ! "actual" profile id
      !  explanation: for some quantities there are two representations,
      !  e.g. Bphi(R,Z) and Bphi(rho,theta).  The passed id can refer to
      !  either of these; the coordinate sought is found for the appropriate
      !  representation.  Thus, if id points to Bphi(R,Z) and info on a rho
      !  grid is sought, id_actual (if present) will be set to the id of 
      !  Bphi(rho,theta).

      !-----------------------

      integer :: id_funs(2),idf,iertmp
      integer :: idi,ingrids,idg1,idg2,idg3,jcoord

      !-----------------------

      id_funs(1)=id
      call xplasma_prof_info(s,id_funs(1),ierr, profId1=id_funs(2))
      if(ierr.ne.0) return

      call xplasma_coord_info(s,icoord,ierr, ngrids=ingrids)
      if(ierr.ne.0) return

      idgrid=0
      if(present(id_actual)) id_actual=0

      do idi=1,2
         idf=id_funs(idi)
         if(idf.eq.0) cycle

         call xplasma_prof_info(s,idf,iertmp, &
              gridId1=idg1, gridId2=idg2, gridId3=idg3)

         if(idg1.gt.0) then
            call xplasma_grid_info(s,idg1,ierr, coord=jcoord)
            if(jcoord.eq.icoord) then
               idgrid=idg1
               if(present(id_actual)) id_actual=idf
               exit
            endif
         endif

         if(idg2.gt.0) then
            call xplasma_grid_info(s,idg2,ierr, coord=jcoord)
            if(jcoord.eq.icoord) then
               idgrid=idg2
               if(present(id_actual)) id_actual=idf
               exit
            endif
         endif

         if(idg3.gt.0) then
            call xplasma_grid_info(s,idg3,ierr, coord=jcoord)
            if(jcoord.eq.icoord) then
               idgrid=idg3
               if(present(id_actual)) id_actual=idf
               exit
            endif
         endif

      enddo

    end subroutine xplasma_prof_gridInfo

    !==========================================================
    subroutine xplasma_create_1dprof(s,profname,idx,data,idf,ier, &
         ispline, ibca,zbca,ibcb,zbcb, coeffs, ccwflag,assoc_id,label,units)
    
      !  create a new profile with the specified name

      !--------
      !  standard INPUT arguments...

      type (xplasma), pointer :: s
      character*(*), intent(in) :: profname  ! name of profile
      integer, intent(in) :: idx         ! id of grid over which it is defined
      real*8, dimension(:), intent(in) :: data  ! the profile data
      !  (NOTE: size of data must match size of grid;
      !         for step function:  must match grid size - 1)

      !--------
      !  standard OUTPUT arguments...

      integer, intent(out) :: idf        ! id of profile now being (re)defined.
      integer, intent(out) :: ier        ! completion code (0=normal)

      !--------
      !  optional input arguments...
      !  defaults: pc linear; for splines: "not a knot" BC
      !            angle grid: normal CCW orientation.

      integer, intent(in), optional :: ispline  ! order of fit
      ! -1 -- step fcn 0 -- pc linear; 1 -- Hermite; 2 -- cubic spline;
      ! default depends on array dimensioning-- if data array dimensions
      ! match grid size, 0 is default; if dimensions match ([grid size]-1),
      ! -1 is the default and the only valid option.

      integer, intent(in), optional :: ibca     ! BC control @ LHS of grid
      real*8, intent(in), optional :: zbca      ! BC data @ LHS of grid
      integer, intent(in), optional :: ibcb     ! BC control @ RHS of grid
      real*8, intent(in), optional :: zbcb      ! BC data @ RHS of grid

      ! ibca=ibcb=0 = default: "not a knot" for splines; Akima for Hermite
      ! for periodic functions BC args are ignored; periodic BC is used.

      ! ibc[*]=1 -- df/dx specified in zbc[*] (0 if absent)
      ! ibc[*]=2 -- d^2f/dx^2 specified in zbc[*] (0 if absent)

      ! INSTEAD of boundary conditions the entire set of nodal 1st derivatives
      ! (if Hermite) or spline coefficients (if Spline) can be provided; if
      ! this is done it is up to the caller to make sure these are correct!

      real*8, intent(in), dimension(:), optional :: coeffs

      ! it is an error to input coeffs for piecewise linear or zonal data...

      logical, intent(in), optional :: ccwflag  ! T means data vs. CCW angle
      ! grid; F means CW.  Default=T.  ccwflag is ignored if the x grid of
      ! the profile is not an angle coordinate.

      integer, intent(in), optional :: assoc_id

      character*(*), intent(in), optional :: label
      character*(*), intent(in), optional :: units

      !------------------------------------------
      integer :: igrids(maxRank),imatch,idimx1,idimx1a,iertmp,imul,isize
      integer :: iadr,jadr,iadx,i,ipx,idum,iflag,iada,istat,ii
      logical :: perio1,mccwflag1
      integer :: jbcx1a,jbcx1b,jbcmin,jbcmax,maxok
      integer :: jspline,icount,jsign
      real*8 :: zzbcx1a,zzbcx1b,zdx,zcoeff
      real*8, dimension(:), allocatable :: dtmp
      character*128 msgbuf
      character*32 zname

      integer :: k_tmp,idba
      !------------------------------------------

      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_create_1dprof call)',ier)
      endif
      if(ier.ne.0) return

      idf = 0

      call xplasma_reserve_check(s,xplasma_profType,1,profname,iflag,ier)
      if(ier.ne.0) return

      call xplasma_add_item(s,profname,xplasma_profType,idf,ier)
      if(ier.ne.0) return

      iadr = s%dict(idf)%dindex
      call xpprof_free(s%profs(iadr))  ! delete any old data if necessary...

      !-----------------
      !  check x grid

      call xplasma_grid_size(s,idx,idimx1,ier)
      if(ier.ne.0) then
         msgbuf=' '
         write(msgbuf,*) 'invalid x1 grid id (xplasma_create_1dprof): ',idx
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)
         ier=501
      endif
      if(ier.ne.0) return

      if(iflag.eq.1) then
         !  require f(rho)
         call xplasma_grid_info(s,idx,ier, coord=idum)
         if(idum.ne.xplasma_rho_coord) then
            call xplasma_error_reserve(s,xplasma_profType,1,profname,ier)
            ier=89
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
            return
         endif
      endif

      !-----------------
      if(present(assoc_id)) then
         if(assoc_id.ne.0) then
            k_tmp = 0
            if((assoc_id.ge.1).and.(assoc_id.le.s%nitems)) then
               k_tmp = s%dict(assoc_id)%dtype
            endif
            if((assoc_id.lt.0).or.(assoc_id.gt.s%nitems)) then
               ier=520
               msgbuf=' '
               write(msgbuf,*) '  assoc_id = ',assoc_id,' value out of range.'
               call xplasma_errmsg_append(s,msgbuf)
            else if( k_tmp .ne.xplasma_profType) then
               ier=520
               msgbuf=' '
               write(msgbuf,*) '  assoc_id = ',assoc_id,' does not point to a profile, it points to: ',trim(s%dict(assoc_id)%name)
               call xplasma_errmsg_append(s,msgbuf)
            else
               iada=s%dict(assoc_id)%dindex
               k_tmp = s%profs(iada)%rank
               idba = s%profs(iada)%profIds(1)  ! back ref must be 0 or "self"
               if((idba.gt.0).and.(idba.ne.idf)) then
                  ier=520
                  msgbuf=' '
                  write(msgbuf,*) '  assoc_id = ',assoc_id,' points to ', &
                       trim(s%dict(assoc_id)%name)
                  call xplasma_errmsg_append(s,msgbuf)
                  call xplasma_errmsg_append(s,' which is already associated.')

               else if( k_tmp .eq.1) then
                  ier=520
                  msgbuf=' '
                  write(msgbuf,*) '  assoc_id = ',assoc_id,' points to ', &
                       trim(s%dict(assoc_id)%name)
                  call xplasma_errmsg_append(s,msgbuf)
                  call xplasma_errmsg_append(s,'  two 1d profiles cannot be associated.')
               else
                  call ck_coords(s,idx,0,0,assoc_id,istat)
                  if(istat.ne.4) then
                     ier=520
                     msgbuf=' '
                     write(msgbuf,*) '  assoc_id = ',assoc_id,' points to ', &
                          trim(s%dict(assoc_id)%name)
                     call xplasma_errmsg_append(s,msgbuf)
                     call xplasma_errmsg_append(s, &
                          '  association not possible-- coordinates not disjoint.')
                  endif
               endif
            endif
            if(ier.eq.520) call xplasma_errmsg_append(s, &
                 ' error in xplasma_create_1dprof, profname='//trim(profname))
         endif
      endif
      if(ier.ne.0) return

      !-----------------
      !  check for CW/-coordinate reversal

      call xplasma_coord_isPeriodic(s,idx,perio1,iertmp)

      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag)) mccwflag1 = ccwflag

      !------------------
      if(size(data).eq.idimx1-1) then
         jspline=-1
      else
         jspline=0
      endif

      if(iflag.eq.1) jspline=s%rzOrder  ! for equilibrium piece

      if(present(ispline)) jspline=ispline

      if((jspline.lt.-1).or.(jspline.gt.2)) then
         ier=505
         write(msgbuf,*) &
              'interpolation fit method argument value invalid: ',jspline
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) &
              '  use 0 for piecewise linear, 1 for Hermite, 2 for spline.'
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)
      endif
      if(ier.ne.0) return

      jsign=1
      if(jspline.le.0) then
         maxok=0
         if(present(coeffs)) then
            ier=523
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
            return
         endif

      else if(jspline.eq.1) then
         maxok=1
         if(.not.mccwflag1) jsign = -1  ! for Hermite coeffs only

      else if(jspline.eq.2) then
         maxok=2
      endif

      !------------------
      if(present(coeffs)) then
         if(present(ibca).or.present(ibcb).or. &
              present(zbca).or.present(zbcb)) then
            ier=525
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
            return
         endif
      endif

      !------------------
      jbcx1a=0
      jbcx1b=0
      if(present(ibca)) jbcx1a=ibca
      if(present(ibcb)) jbcx1b=ibcb

      jbcmin=min(jbcx1a,jbcx1b)
      jbcmax=max(jbcx1a,jbcx1b)

      if((jbcmin.lt.0).or.(jbcmax.gt.maxok)) then
         ier=504
         write(msgbuf,*) &
              'boundary condition type argument value invalid:',jbcmin,jbcmax
         call xplasma_errmsg_append(s,msgbuf)
         if(maxok.eq.0) then
            write(msgbuf,*) ' for piecewise linear, no BC options allowed.'
         else if(maxok.eq.1) then
            write(msgbuf,*) &
                 ' for Hermite, only default (0) or 1st derivative (1) BCs are supported.'
         else if(maxok.eq.2) then
            write(msgbuf,*) &
                 ' BC option argument values between 0 and 2 expected.'
         endif
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)
      endif
      if(ier.ne.0) return

      !------------------
      zzbcx1a=0
      zzbcx1b=0
      if(perio1) then
         jbcx1a=-1
         jbcx1b=-1
      else
         if(present(zbca)) zzbcx1a=zbca
         if(present(zbcb)) zzbcx1b=zbcb
      endif

      !------------------
      !  check data size

      if(jspline.ge.0) then
         idimx1a=idimx1
      else
         idimx1a=idimx1-1
      endif
      if(size(data).ne.idimx1a) then
         ier=506
         write(msgbuf,*) ' ?xpprof: size of profile data = ',size(data)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   does not match x grid size  = ',idimx1, &
              ' x grid name: ',s%dict(idx)%name
         call xplasma_errmsg_append(s,msgbuf)
         if(jspline.eq.-1) call xplasma_errmsg_append(s, &
              '    NOTE: step function data: need size = [grid size] - 1')
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)

      endif

      if(present(coeffs)) then
         if(size(coeffs).ne.idimx1a) then
            ier=524
            write(msgbuf,*) ' ?xpprof: size of coeffs data = ',size(coeffs)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '   does not match x grid size  = ',idimx1, &
                 ' x grid name: ',s%dict(idx)%name
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
         endif
      endif
      if(ier.ne.0) return

      !------------------
      ! create buffer space for spline coefficients

      if(jspline.le.0) then
         imul=1
      else if(jspline.eq.1) then
         imul=2
      else if(jspline.eq.2) then
         imul=2
      endif

      isize=idimx1a*imul

      s%profs(iadr)%rank = 1
      s%profs(iadr)%gridIds = 0
      s%profs(iadr)%gridIds(1) = idx
      s%profs(iadr)%profIds = 0
      if(present(assoc_id)) then
         s%profs(iadr)%profIds(1) = assoc_id
         if(assoc_id.gt.0) then
            iada = s%dict(assoc_id)%dindex
            s%profs(iada)%profIds(1) = idf
         endif
      endif
      s%profs(iadr)%kspline = jspline
      s%profs(iadr)%size = isize
      allocate(s%profs(iadr)%buf(isize))

      !------------------
      ! allocate tmp space for data; reverse if flag is set

      allocate(dtmp(idimx1a))
      if(mccwflag1) then
         dtmp = data
      else
         do i=1,idimx1a
            dtmp(i)=data(idimx1a+1-i)
         enddo
      endif

      !------------------
      ! compute spline coefficients as needed

      iadx=s%dict(idx)%dindex

      if(present(coeffs)) then
         ii=0
         do i=1,idimx1a
            ii=ii+1
            s%profs(iadr)%buf(ii)=dtmp(i)
            if(mccwflag1) then
               zcoeff = coeffs(i)*jsign
            else
               zcoeff = coeffs(idimx1a+1-i)*jsign
            endif
            ii=ii+1
            s%profs(iadr)%buf(ii)=zcoeff
         enddo
            
      else if(jspline.le.0) then
         s%profs(iadr)%buf = dtmp

      else if(jspline.eq.1) then

         s%profs(iadr)%buf(1:isize:imul) = dtmp

         if(jbcx1a.eq.-1) then
            ipx=1
         else if((jbcx1a.eq.1).or.(jbcx1b.eq.1)) then
            ipx=2
            if(jbcx1a.eq.1) then
               s%profs(iadr)%buf(2)=zzbcx1a
            else
               zdx = s%grids(iadx)%xpkg(2,1) - s%grids(iadx)%xpkg(1,1)
               s%profs(iadr)%buf(2)=(dtmp(2)-dtmp(1))/zdx
            endif
            if(jbcx1b.eq.1) then
               s%profs(iadr)%buf(isize)=zzbcx1b
            else
               zdx=s%grids(iadx)%xpkg(idimx1,1)-s%grids(iadx)%xpkg(idimx1-1,1)
               s%profs(iadr)%buf(isize)=(dtmp(idimx1)-dtmp(idimx1-1))/zdx
            endif
         else
            ipx=0
         endif
         call r8akherm1p(s%grids(iadx)%xpkg(1:idimx1,1),idimx1, &
              s%profs(iadr)%buf,idum,ipx,ier)
         if(ier.ne.0) then
            ier=9999
            msgbuf = &
                 ' unexpected failure in r8akherm1p, xplasma_create_1dprof'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
         endif

      else if(jspline.eq.2) then
         s%profs(iadr)%buf(1:isize:imul) = dtmp
         call r8mkspline(s%grids(iadx)%xpkg(1:idimx1,1),idimx1, &
              s%profs(iadr)%buf, &
              jbcx1a,zzbcx1a,jbcx1b,zzbcx1b, &
              idum,ier)
         if(ier.ne.0) then
            ier=9999
            msgbuf = &
                 ' unexpected failure in r8cspline, xplasma_create_1dprof'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
         endif

      endif

      deallocate(dtmp)

      if(ier.le.0) then
         s%prof_counter = s%prof_counter + 1
         s%profs(iadr)%prof_counter = s%prof_counter
         s%grids(iadx)%nrefs = s%grids(iadx)%nrefs + 1
      endif

      iertmp=0

      if(present(label)) then
         call xplasma_label_item(s,idf,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_label_item(s,idf,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      ier = iertmp

      if((iflag.gt.0).and.(ier.eq.0)) then
         call xplasma_internal_ids(s)
      endif

      if(ier.ne.0) then
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)
      endif

    end subroutine xplasma_create_1dprof

    !==========================================================
    subroutine xplasma_create_2dprof(s,profname, &
         idx1,idx2,data,idf,ier, &
         ispline,hspline, &
         ibcx1a,zbcx1a,ibcx1b,zbcx1b, &
         ibcx2a,zbcx2a,ibcx2b,zbcx2b, &
         coeffs, ccwflag1,ccwflag2, &
         assoc_id, label,units)
    
      !  create a new profile with the specified name

      !--------
      !  standard INPUT arguments...

      type (xplasma), pointer :: s
      character*(*), intent(in) :: profname  ! name of profile
      integer, intent(in) :: idx1    ! id of 1st grid over which it is defined
      integer, intent(in) :: idx2    ! id of 2nd grid over which it is defined
      real*8, dimension(:,:), intent(in) :: data  ! the profile data
      !  (NOTE: sizes of dimensions MUST match sizes of identified grids)
      !  (For step function: dimensions must be ([grid size] - 1) EXACTLY).

      !--------
      !  standard OUTPUT arguments...

      integer, intent(out) :: idf       ! id of container for this profile.
      integer, intent(out) :: ier       ! completion code (0=normal)

      !--------
      !  optional input arguments...
      !  defaults: pc linear; for splines: "not a knot" BC
      !            angle grid: normal CCW orientation.

      !-----------------------------------
      !  ispline or hspline can be present but not both...

      integer, intent(in), optional :: ispline  ! order of fit
      ! -1 -- step fcn 0 -- pc linear; 1 -- Hermite; 2 -- cubic spline;
      ! default depends on array dimensioning-- if data array dimensions
      ! match grid size, 0 is default; if dimensions match ([grid size]-1),
      ! -1 is the default and the only valid option.
      !    (ispline control applies to all grid dimensions)

      integer, intent(in), dimension(:), optional :: hspline  ! hyrbrid spline
      ! control: hspline(1) = 1st dimension spline order; hspline(2) = 2nd
      !   dimension spline order; {-1,0,1,2} for {zonal, piecewise linear,
      !   Hermite, C2 Spline}; any combination legal, except that Hermite and
      !   C2 splines cannot currently be mixed.

      !  use hspline to make a hybrid spline -- e.g. pc-linear in 1 dimension,
      !   Hermite or cubic spline in the other.
      !-----------------------------------

      integer, intent(in), optional :: ibcx1a     ! BC control @ LHS of grid x1
      real*8, intent(in), dimension(:), optional :: zbcx1a ! BC data @ LHS of grid x1
      integer, intent(in), optional :: ibcx1b     ! BC control @ RHS of grid x1
      real*8, intent(in), dimension(:), optional :: zbcx1b ! BC data @ RHS of grid x1

      integer, intent(in), optional :: ibcx2a     ! BC control @ LHS of grid x2
      real*8, intent(in), dimension(:), optional :: zbcx2a ! BC data @ LHS of grid x2
      integer, intent(in), optional :: ibcx2b     ! BC control @ RHS of grid x2
      real*8, intent(in), dimension(:), optional :: zbcx2b ! BC data @ RHS of grid x2

      ! ibc*a=ibc*b=0 = default: "not a knot" for splines; Akima for Hermite
      ! for periodic functions BC args are ignored; periodic BC is used.

      ! ibc[*]=1 -- df/dx specified in zbc[*] (0 if absent)
      ! ibc[*]=2 -- d^2f/dx^2 specified in zbc[*] (0 if absent)

      ! INSTEAD of boundary conditions the entire set of nodal 1st derivatives
      ! (if Hermite) or spline coefficients (if Spline) can be provided; if
      ! this is done it is up to the caller to make sure these are correct!

      real*8, intent(in), dimension(:,:,:), optional :: coeffs

      !  1st dimension of coeffs: #coeffs/nodal point
      !     e.g. for BiHermite: 3: d/dx1, d/dx2, d2/dx1dx2
      !  2nd dimension of coeffs: matches 1st dimension of data
      !  3nd dimension of coeffs: matches 2nd dimension of data
      ! it is an error to input coeffs for piecewise linear or zonal data...

      logical, intent(in), optional :: ccwflag1  ! T means data vs. CCW angle
      ! grid; F means CW.  Default=T.

      logical, intent(in), optional :: ccwflag2  ! T means data vs. CCW angle
      ! grid; F means CW.  Default=T.

      ! ccwflag1 applies to x1 *if* it is periodic; otherwise it is ignored.
      ! ccwflag2 applies to x2 *if* it is periodic; otherwise it is ignored.

      integer, intent(in), optional :: assoc_id  ! id of associated profile

      character*(*), intent(in), optional :: label
      character*(*), intent(in), optional :: units

      !------------------------------------------
      integer :: jspline(2),jbcx1a,jbcx1b,jbcx2a,jbcx2b,icount,inrho,idrho,iada
      integer :: kspline
      logical :: mccwflag1,mccwflag2
      real*8, dimension(:), allocatable :: zzbcx1a,zzbcx1b,zzbcx2a,zzbcx2b
      real*8, dimension(:,:), allocatable :: dtmp
      real*8, dimension(:,:,:), allocatable :: ctmp
      integer :: idimx1,idimx2,iertmp,i,ii,j,jj,kk,iadr,ipx1,ipx2,iflag,istat
      integer :: idimx1a,idimx2a
      integer :: jbcmin,jbcmax,maxok,imul,isize,iadx1,iadx2,idum1,idum2
      logical :: perio1,perio2
      integer :: jcoord1,jcoord2
      integer :: ksign(3),isize_exp
      character*128 msgbuf
      character*32 zname

      integer :: k_tmp,idba
      !------------------------------------------

      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_create_2dprof call)',ier)
      endif
      if(ier.ne.0) return

      idf = 0

      call xplasma_reserve_check(s,xplasma_profType,2,profname,iflag,ier)
      if(ier.ne.0) return

      call xplasma_add_item(s,profname,xplasma_profType,idf,ier)
      if(ier.ne.0) return

      iadr = s%dict(idf)%dindex
      call xpprof_free(s%profs(iadr))  ! delete any old data if necessary...

      !-----------------
      !  check x axes

      iertmp=0
      call xplasma_grid_size(s,idx1,idimx1,ier)
      if(ier.ne.0) then
         ier=501
         write(msgbuf,*) 'invalid x1 grid id (xplasma_create_2dprof): ',idx1
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         iertmp=max(iertmp,ier)
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)
      endif

      call xplasma_grid_size(s,idx2,idimx2,ier)
      if(ier.ne.0) then
         ier=501
         write(msgbuf,*) 'invalid x2 grid id (xplasma_create_2dprof): ',idx2
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         iertmp=max(iertmp,ier)
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)
      endif

      ier=iertmp
      if(ier.ne.0) return

      !-----------------
      if(present(assoc_id)) then
         if(assoc_id.ne.0) then
            k_tmp = 0
            if((assoc_id.ge.1).and.(assoc_id.le.s%nitems)) then
               k_tmp = s%dict(assoc_id)%dtype
            endif
            if((assoc_id.lt.0).or.(assoc_id.gt.s%nitems)) then
               ier=520
               msgbuf=' '
               write(msgbuf,*) '  assoc_id = ',assoc_id,' value out of range.'
               call xplasma_errmsg_append(s,msgbuf)
            else if( k_tmp .ne.xplasma_profType) then
               ier=520
               msgbuf=' '
               write(msgbuf,*) '  assoc_id = ',assoc_id,' does not point to a profile, it points to: ',trim(s%dict(assoc_id)%name)
               call xplasma_errmsg_append(s,msgbuf)
            else
               iada=s%dict(assoc_id)%dindex
               idba = s%profs(iada)%profIds(1)  ! back ref must be 0 or "self"
               if((idba.gt.0).and.(idba.ne.idf)) then
                  ier=520
                  msgbuf=' '
                  write(msgbuf,*) '  assoc_id = ',assoc_id,' points to ', &
                       trim(s%dict(assoc_id)%name)
                  call xplasma_errmsg_append(s,msgbuf)
                  call xplasma_errmsg_append(s,' which is already associated.')

               else
                  call ck_coords(s,idx1,idx2,0,assoc_id,istat)
                  if(istat.ne.4) then
                     ier=520
                     msgbuf=' '
                     write(msgbuf,*) '  assoc_id = ',assoc_id,' points to ', &
                          trim(s%dict(assoc_id)%name)
                     call xplasma_errmsg_append(s,msgbuf)
                     call xplasma_errmsg_append(s, &
                          '  association not possible-- coordinates not disjoint.')
                  endif
               endif
            endif
            if(ier.eq.520) call xplasma_errmsg_append(s, &
                 ' error in xplasma_create_2dprof, profname='//trim(profname))
         endif
      endif
      if(ier.ne.0) return

      !-----------------
      !  check for CW/-coordinate reversal

      call xplasma_coord_isPeriodic(s,idx1,perio1,iertmp)
      call xplasma_coord_isPeriodic(s,idx2,perio2,iertmp)

      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1 = ccwflag1

      mccwflag2=.TRUE.
      if(perio2.and.present(ccwflag2)) mccwflag2 = ccwflag2

      iadx1=s%dict(idx1)%dindex
      iadx2=s%dict(idx2)%dindex

      jcoord1=s%grids(iadx1)%coord
      jcoord2=s%grids(iadx2)%coord

      if(jcoord1.eq.jcoord2) then
         ier=515
         write(msgbuf,*) '  ',trim(s%dict(idx1)%name),' coordinate is ', &
              trim(s%dict(jcoord1)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  ',trim(s%dict(idx2)%name),' coordinate is ', &
              trim(s%dict(jcoord2)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  profile name is: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      else
         if(iflag.eq.1) then
            idum2=0
            if((jcoord1.eq.xplasma_rho_coord).or. &
                 (jcoord2.eq.xplasma_rho_coord)) idum2 = idum2 + 1
            if((jcoord1.eq.xplasma_theta_coord).or. &
                 (jcoord2.eq.xplasma_theta_coord)) idum2 = idum2 + 1
            if(idum2.ne.2) then
               call xplasma_error_reserve(s,xplasma_profType,2,profname,ier)
               ier=89
               return
            endif
         endif
      endif

      !------------------
      if(present(ispline).and.present(hspline)) then
         ier=90
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      if(size(data,1).eq.(idimx1-1)) then
         jspline=-1
      else
         jspline=0
      endif

      if(iflag.eq.1) then
         jspline=s%rzOrder
      endif

      if(present(ispline)) jspline=ispline

      if(present(hspline)) then
         if(size(hspline).ne.2) then
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_create_2dprof: size(hspline).ne.2')
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
            ier=91
            return
         else
            jspline=hspline
         endif
      endif

      call spline_type_ck_encode(s,'xplasma_create_2dprof',s%dict(idf)%name, &
           jspline,kspline,ier)
      if(ier.ne.0) return

      if(jspline(1).lt.0) then
         idimx1a=idimx1-1
      else
         idimx1a=idimx1
      endif

      if(jspline(2).lt.0) then
         idimx2a=idimx2-1
      else
         idimx2a=idimx2
      endif

      if(maxval(jspline).le.0) then
         if(present(coeffs)) then
            ier=523
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
            return
         endif
      endif

      if(present(coeffs)) then
         if(present(ibcx1a).or.present(ibcx1b).or. &
              present(zbcx1a).or.present(zbcx1b).or. &
              present(ibcx2a).or.present(ibcx2b).or. &
              present(zbcx2a).or.present(zbcx2b)) then
            ier=525
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
            return
         endif
      endif

      !------------------
      jbcx1a=0
      jbcx1b=0
      if(present(ibcx1a)) jbcx1a=ibcx1a
      if(present(ibcx1b)) jbcx1b=ibcx1b

      jbcx2a=0
      jbcx2b=0
      if(present(ibcx2a)) jbcx2a=ibcx2a
      if(present(ibcx2b)) jbcx2b=ibcx2b

      isize_exp=1
      do ii=1,2
         if(jspline(ii).le.0) then
            maxok=0
         else if(jspline(ii).eq.1) then
            maxok=1
            isize_exp=isize_exp*2
         else if(jspline(ii).eq.2) then
            maxok=2
            isize_exp=isize_exp*2
         endif

         if(ii.eq.1) then
            jbcmin=min(jbcx1a,jbcx1b)
            jbcmax=max(jbcx1a,jbcx1b)
         else
            jbcmin=min(jbcx2a,jbcx2b)
            jbcmax=max(jbcx2a,jbcx2b)
         endif

         if((jbcmin.lt.0).or.(jbcmax.gt.maxok)) then
            if(ier.eq.0) then
               write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
               call xplasma_errmsg_append(s,msgbuf)
            endif
            ier=504
            write(msgbuf,*) &
                 'boundary condition type argument value invalid:',jbcmin,jbcmax
            call xplasma_errmsg_append(s,msgbuf)
            if(maxok.eq.0) then
               write(msgbuf,*) ' for zone or piecewise linear, no BC options allowed.'
            else if(maxok.eq.1) then
               write(msgbuf,*) &
                    ' for Hermite, only default (0) or 1st derivative (1) BCs are supported.'
            else if(maxok.eq.2) then
               write(msgbuf,*) &
                    ' BC option argument values between 0 and 2 expected.'
            endif
            call xplasma_errmsg_append(s,msgbuf)
         endif
      enddo
      if(ier.ne.0) return

      isize_exp=isize_exp-1
      if(present(coeffs)) then
         call ksign_2d(mccwflag1,mccwflag2,jspline,isize_exp,ksign)
      endif

      !------------------
      !  if explicit boundary conditions are indicated, make sure that
      !  the supplied data is correctly sized...

      allocate(zzbcx1a(idimx2a),zzbcx1b(idimx2a))
      zzbcx1a=0
      zzbcx1b=0

      allocate(zzbcx2a(idimx1a),zzbcx2b(idimx1a))
      zzbcx2a=0
      zzbcx2b=0

      if(perio1) then
         jbcx1a=-1
         jbcx1b=-1
      else
         if(jbcx1a.gt.0) then
            if(present(zbcx1a)) then
               if(size(zbcx1a).ne.idimx2a) then
                  ier=502
                  write(msgbuf,*) 'size of x1(1) bdy condition data does',&
                       ' not match size of dimension (x2): ',idimx2a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do ii=1,idimx2a
                     if(mccwflag2) then
                        jj=ii
                     else
                        jj=idimx2a+1-ii
                     endif
                     zzbcx1a(jj)=zbcx1a(ii)
                  enddo
               endif
            endif
         endif
         if(jbcx1b.gt.0) then
            if(present(zbcx1b)) then
               if(size(zbcx1b).ne.idimx2a) then
                  ier=502
                  write(msgbuf,*) 'size of x1(Nx1) bdy condition data does',&
                       ' not match size of dimension (x2): ',idimx2a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do ii=1,idimx2a
                     if(mccwflag2) then
                        jj=ii
                     else
                        jj=idimx2a+1-ii
                     endif
                     zzbcx1b(jj)=zbcx1b(ii)
                  enddo
               endif
            endif
         endif
      endif

      if(perio2) then
         jbcx2a=-1
         jbcx2b=-1
      else
         if(jbcx2a.gt.0) then
            if(present(zbcx2a)) then
               if(size(zbcx2a).ne.idimx1a) then
                  ier=502
                  write(msgbuf,*) 'size of x2(1) bdy condition data does',&
                       ' not match size of dimension (x1): ',idimx1a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do ii=1,idimx1a
                     if(mccwflag1) then
                        jj=ii
                     else
                        jj=idimx1a+1-ii
                     endif
                     zzbcx2a(jj)=zbcx2a(ii)
                  enddo
               endif
            endif
         endif
         if(jbcx2b.gt.0) then
            if(present(zbcx2b)) then
               if(size(zbcx2b).ne.idimx1a) then
                  ier=502
                  write(msgbuf,*) 'size of x2(Nx2) bdy condition data does',&
                       ' not match size of dimension (x1): ',idimx1a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do ii=1,idimx1a
                     if(mccwflag1) then
                        jj=ii
                     else
                        jj=idimx1a+1-ii
                     endif
                     zzbcx2b(jj)=zbcx2b(ii)
                  enddo
               endif
            endif
         endif
      endif
      if(ier.ne.0) then
         deallocate(zzbcx1a,zzbcx1b,zzbcx2a,zzbcx2b)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      !------------------
      !  check data size

      if(size(data,1).ne.idimx1a) then
         ier=506
         write(msgbuf,*) ' ?xpprof: size of profile 1st dim = ',size(data,1)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   does not match x1 grid size  = ',idimx1, &
              ' x1 grid name: ',trim(s%dict(idx1)%name)
         call xplasma_errmsg_append(s,msgbuf)
         if(jspline(1).eq.-1) call xplasma_errmsg_append(s, &
              '    NOTE: step function data: need size = [grid size] - 1')
      endif

      if(size(data,2).ne.idimx2a) then
         ier=506
         write(msgbuf,*) ' ?xpprof: size of profile 2nd dim = ',size(data,2)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   does not match x2 grid size  = ',idimx2, &
              ' x2 grid name: ',trim(s%dict(idx2)%name)
         call xplasma_errmsg_append(s,msgbuf)
         if(jspline(2).eq.-1) call xplasma_errmsg_append(s, &
              '    NOTE: step function data: need size = [grid size] - 1')
      endif

      if(present(coeffs)) then
         if(size(coeffs,1).ne.isize_exp) then
            ier=524
            write(msgbuf,*) ' ?xpprof: 1st dim size of coeffs: ',size(coeffs,1)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  does not match expected coeffs/node: ',isize_exp
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  expected: 2**N-1, where N = #of spline or hermite dimensions.'
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(size(coeffs,2).ne.idimx1a) then
            ier=524
            write(msgbuf,*) ' ?xpprof: size of 2nd dim of coeffs data = ',size(coeffs,2)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '   does not match x1 grid size  = ',idimx1, &
                 ' x1 grid name: ',trim(s%dict(idx1)%name)
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(size(coeffs,3).ne.idimx2a) then
            ier=524
            write(msgbuf,*) ' ?xpprof: size of 3rd dim of coeffs data = ',size(coeffs,3)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '   does not match x2 grid size  = ',idimx2, &
                 ' x2 grid name: ',trim(s%dict(idx2)%name)
            call xplasma_errmsg_append(s,msgbuf)
         endif
      endif

      if(ier.ne.0) then
         deallocate(zzbcx1a,zzbcx1b,zzbcx2a,zzbcx2b)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      !------------------
      ! create buffer space for spline coefficients

      imul=1
      do ii=1,2
         if(jspline(ii).gt.0) then
            imul=imul*2
         endif
      enddo

      isize=idimx1a*idimx2a*imul

      s%profs(iadr)%rank = 2
      s%profs(iadr)%gridIds = 0
      s%profs(iadr)%gridIds(1) = idx1
      s%profs(iadr)%gridIds(2) = idx2
      s%profs(iadr)%profIds = 0
      if(present(assoc_id)) then
         s%profs(iadr)%profIds(1) = assoc_id
         if(assoc_id.gt.0) then
            iada = s%dict(assoc_id)%dindex
            s%profs(iada)%profIds(1) = idf
         endif
      endif
      s%profs(iadr)%kspline = kspline
      s%profs(iadr)%size = isize
      allocate(s%profs(iadr)%buf(isize))

      !------------------
      ! allocate tmp space for data; reverse if flag is set

      allocate(dtmp(idimx1a,idimx2a))
      if(mccwflag1.and.mccwflag2) then
         dtmp = data
      else
         if(present(coeffs)) allocate(ctmp(isize_exp,idimx1a,idimx2a))
         do j=1,idimx2a
            if(mccwflag2) then
               jj=j
            else
               jj=idimx2a+1-j
            endif
            do i=1,idimx1a
               if(mccwflag1) then
                  ii=i
               else
                  ii=idimx1a+1-i
               endif
               dtmp(ii,jj)=data(i,j)
               if(present(coeffs)) then
                  do kk=1,isize_exp
                     ctmp(kk,ii,jj)=ksign(kk)*coeffs(kk,i,j)
                  enddo
               endif
            enddo
         enddo
      endif

      !------------------
      ! compute spline coefficients as needed

      ii=-imul+1
      do j=1,idimx2a
         do i=1,idimx1a
            ii=ii+imul
            s%profs(iadr)%buf(ii)=dtmp(i,j)
         enddo
      enddo

      if(imul.eq.1) then
         continue   ! pclin or zone data already in place
         ! no coefficients to be evaluated

      else if(present(coeffs)) then

         ii=-imul+1
         do j=1,idimx2a
            do i=1,idimx1a
               ii=ii+imul
               if(mccwflag1.and.mccwflag2) then
                  s%profs(iadr)%buf(ii+1:ii+isize_exp)=coeffs(1:isize_exp,i,j)
               else
                  s%profs(iadr)%buf(ii+1:ii+isize_exp)=ctmp(1:isize_exp,i,j)
               endif
            enddo
         enddo

      else

         ! create spline or Hermite coefficients

         call r8mkintrp2d(s%grids(iadx1)%xpkg(1:idimx1,1),idimx1, &
              s%grids(iadx2)%xpkg(1:idimx2,1),idimx2, &
              jspline, &
              s%profs(iadr)%buf,imul,idimx1a,idimx2a, &
              jbcx1a,zzbcx1a,jbcx1b,zzbcx1b, &
              jbcx2a,zzbcx2a,jbcx2b,zzbcx2b, &
              ier)

         if(ier.ne.0) then
            ier=9999
            msgbuf = &
                 ' unexpected failure in r8mkintrp, xplasma_create_2dprof'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         
      endif

      deallocate(dtmp,zzbcx1a,zzbcx1b,zzbcx2a,zzbcx2b)
      if(allocated(ctmp)) deallocate(ctmp)

      if(ier.le.0) then
         s%prof_counter = s%prof_counter + 1
         s%profs(iadr)%prof_counter = s%prof_counter
         s%grids(iadx1)%nrefs = s%grids(iadx1)%nrefs + 1
         s%grids(iadx2)%nrefs = s%grids(iadx2)%nrefs + 1
      endif

      iertmp=0

      if(present(label)) then
         call xplasma_label_item(s,idf,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_label_item(s,idf,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      ier = iertmp

      if((iflag.gt.0).and.(ier.eq.0)) then
         call xplasma_internal_ids(s)
      endif

      if(ier.ne.0) then
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)
      endif

    end subroutine xplasma_create_2dprof

    !==========================================================
    subroutine xplasma_create_3dprof(s,profname, &
         idx1,idx2,idx3,data,idf,ier, &
         ispline,hspline,ibcx1a,zbcx1a,ibcx1b,zbcx1b, &
         ibcx2a,zbcx2a,ibcx2b,zbcx2b, &
         ibcx3a,zbcx3a,ibcx3b,zbcx3b, &
         coeffs, ccwflag1,ccwflag2,ccwflag3, &
         assoc_id,label,units)
    
      !  create a new profile with the specified name

      !--------
      !  standard INPUT arguments...

      type (xplasma), pointer :: s
      character*(*), intent(in) :: profname  ! name of profile
      integer, intent(in) :: idx1    ! id of 1st grid over which it is defined
      integer, intent(in) :: idx2    ! id of 2nd grid over which it is defined
      integer, intent(in) :: idx3    ! id of 3rd grid over which it is defined
      real*8, dimension(:,:,:), intent(in) :: data  ! the profile data
      !  (NOTE: sizes of dimensions MUST match sizes of identified grids)
      !  (For step function: dimensions must be ([grid size] - 1) EXACTLY).

      !--------
      !  standard OUTPUT arguments...

      integer, intent(out) :: idf     ! id of newly created profile.
      integer, intent(out) :: ier     ! completion code (0=normal)

      !--------
      !  optional input arguments...
      !  defaults: pc linear; for splines: "not a knot" BC
      !            angle grid: normal CCW orientation.

      !-----------------------------------
      !  ispline or hspline can be present but not both...

      integer, intent(in), optional :: ispline  ! order of fit
      ! -1 -- step fcn 0 -- pc linear; 1 -- Hermite; 2 -- cubic spline;
      ! default depends on array dimensioning-- if data array dimensions
      ! match grid size, 0 is default; if dimensions match ([grid size]-1),
      ! -1 is the default and the only valid option.
      !    (ispline control applies to all grid dimensions)

      integer, intent(in), dimension(:), optional :: hspline  ! hyrbrid spline
      ! control: hspline(1) = 1st dimension spline order; hspline(2) = 2nd
      !   dimension spline order; hspline(3) = 3rd dimension spline order;
      !   each dimension set to one of {-1,0,1,2} for {zonal, piecewise linear,
      !   Hermite, C2 Spline}; any combination legal, except that Hermite and
      !   C2 splines cannot currently be mixed.

      !  use hspline to make a hybrid spline -- e.g. pc-linear in 1 dimension,
      !   Hermite or cubic spline in the other.
      !-----------------------------------

      integer, intent(in), optional :: ibcx1a     ! BC control @ LHS of grid x1
      real*8, intent(in), dimension(:,:), optional :: zbcx1a ! BC data @ LHS of grid x1
      integer, intent(in), optional :: ibcx1b     ! BC control @ RHS of grid x1
      real*8, intent(in), dimension(:,:), optional :: zbcx1b ! BC data @ RHS of grid x1

      integer, intent(in), optional :: ibcx2a     ! BC control @ LHS of grid x2
      real*8, intent(in), dimension(:,:), optional :: zbcx2a ! BC data @ LHS of grid x2
      integer, intent(in), optional :: ibcx2b     ! BC control @ RHS of grid x2
      real*8, intent(in), dimension(:,:), optional :: zbcx2b ! BC data @ RHS of grid x2

      integer, intent(in), optional :: ibcx3a     ! BC control @ LHS of grid x3
      real*8, intent(in), dimension(:,:), optional :: zbcx3a ! BC data @ LHS of grid x3
      integer, intent(in), optional :: ibcx3b     ! BC control @ RHS of grid x3
      real*8, intent(in), dimension(:,:), optional :: zbcx3b ! BC data @ RHS of grid x3

      ! ibc*a=ibc*b=0 = default: "not a knot" for splines; Akima for Hermite
      ! for periodic functions BC args are ignored; periodic BC is used.

      ! ibc[*]=1 -- df/dx specified in zbc[*] (0 if absent)
      ! ibc[*]=2 -- d^2f/dx^2 specified in zbc[*] (0 if absent)

      ! INSTEAD of boundary conditions the entire set of nodal 1st derivatives
      ! (if Hermite) or spline coefficients (if Spline) can be provided; if
      ! this is done it is up to the caller to make sure these are correct!

      real*8, intent(in), dimension(:,:,:,:), optional :: coeffs

      !  1st dimension of coeffs: #coeffs/nodal point
      !     e.g. for TriHermite: 7: d/dx1, d/dx2, d/dx3,
      !                             d2/dx1dx2, d2/dx1dx3, d2/dx2dx3,
      !                             d3/dx1dx2dx3
      !  2nd dimension of coeffs: matches 1st dimension of data
      !  3nd dimension of coeffs: matches 2nd dimension of data
      !  4th dimension of coeffs: matches 3rd dimension of data
      ! it is an error to input coeffs for piecewise linear or zonal data...

      logical, intent(in), optional :: ccwflag1  ! T means data vs. CCW angle
      ! grid; F means CW.  Default=T.

      logical, intent(in), optional :: ccwflag2  ! T means data vs. CCW angle
      ! grid; F means CW.  Default=T.

      logical, intent(in), optional :: ccwflag3  ! T means data vs. CCW angle
      ! grid; F means CW.  Default=T.

      ! ccwflag1 applies to x1 if it is periodic; otherwise it is ignored.
      ! ccwflag2 applies to x2 if it is periodic; otherwise it is ignored.
      ! ccwflag3 applies to x3 if it is periodic; otherwise it is ignored.

      integer, intent(in), optional :: assoc_id  ! id of associated profile

      character*(*), intent(in), optional :: label
      character*(*), intent(in), optional :: units

      !------------------------------------------
      integer :: jspline(3),jbcx1a,jbcx1b,jbcx2a,jbcx2b,jbcx3a,jbcx3b,icount
      integer :: jbcmin,jbcmax,maxok,iada,kspline
      logical :: mccwflag1,mccwflag2,mccwflag3
      real*8, dimension(:,:), allocatable :: zzbcx1a,zzbcx1b,zzbcx2a,zzbcx2b
      real*8, dimension(:,:), allocatable :: zzbcx3a,zzbcx3b
      real*8, dimension(:,:,:), allocatable :: dtmp
      real*8, dimension(:,:,:,:), allocatable :: ctmp
      integer :: iadr,i,j,k,ii,jj,kk,iadx1,iadx2,iadx3,istat
      integer :: idimx1,idimx2,idimx3,iertmp,iertmp2,iflag
      integer :: idimx1a,idimx2a,idimx3a
      integer :: imul,isize,ipx1,ipx2,ipx3,idum1,idum2,idum3
      logical :: perio1,perio2,perio3
      integer :: jcoord1,jcoord2,jcoord3
      integer :: ksign(7),isize_exp
      character*128 msgbuf

      integer :: k_tmp,idba
      !------------------------------------------

      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_create_3dprof call)',ier)
      endif
      if(ier.ne.0) return

      idf = 0

      call xplasma_reserve_check(s,xplasma_profType,3,profname,iflag,ier)
      if(ier.ne.0) return

      call xplasma_add_item(s,profname,xplasma_profType,idf,ier)
      if(ier.ne.0) return

      iadr = s%dict(idf)%dindex
      call xpprof_free(s%profs(iadr))  ! delete any old data if necessary...

      !-----------------
      !  check x axes

      iertmp=0
      call xplasma_grid_size(s,idx1,idimx1,ier)
      if(ier.ne.0) then
         ier=501
         write(msgbuf,*) 'invalid x1 grid id (xplasma_create_3dprof): ',idx1
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         iertmp=max(iertmp,ier)
      endif

      call xplasma_grid_size(s,idx2,idimx2,ier)
      if(ier.ne.0) then
         ier=501
         write(msgbuf,*) 'invalid x2 grid id (xplasma_create_3dprof): ',idx2
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         iertmp=max(iertmp,ier)
      endif

      call xplasma_grid_size(s,idx3,idimx3,ier)
      if(ier.ne.0) then
         ier=501
         write(msgbuf,*) 'invalid x3 grid id (xplasma_create_3dprof): ',idx3
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         iertmp=max(iertmp,ier)
      endif

      ier=iertmp
      if(ier.ne.0) return

      !-----------------
      if(present(assoc_id)) then
         if(assoc_id.ne.0) then
            k_tmp = 0
            if((assoc_id.ge.1).and.(assoc_id.le.s%nitems)) then
               k_tmp = s%dict(assoc_id)%dtype
            endif
            if((assoc_id.lt.0).or.(assoc_id.gt.s%nitems)) then
               ier=520
               msgbuf=' '
               write(msgbuf,*) '  assoc_id = ',assoc_id,' value out of range.'
               call xplasma_errmsg_append(s,msgbuf)
            else if( k_tmp .ne.xplasma_profType) then
               ier=520
               msgbuf=' '
               write(msgbuf,*) '  assoc_id = ',assoc_id,' does not point to a profile, it points to: ',trim(s%dict(assoc_id)%name)
               call xplasma_errmsg_append(s,msgbuf)
            else
               iada=s%dict(assoc_id)%dindex
               idba = s%profs(iada)%profIds(1)  ! back ref must be 0 or "self"
               if((idba.gt.0).and.(idba.ne.idf)) then
                  ier=520
                  msgbuf=' '
                  write(msgbuf,*) '  assoc_id = ',assoc_id,' points to ', &
                       trim(s%dict(assoc_id)%name)
                  call xplasma_errmsg_append(s,msgbuf)
                  call xplasma_errmsg_append(s,' which is already associated.')

               else
                  call ck_coords(s,idx1,idx2,idx3,assoc_id,istat)
                  if(istat.ne.4) then
                     ier=520
                     msgbuf=' '
                     write(msgbuf,*) '  assoc_id = ',assoc_id,' points to ', &
                          trim(s%dict(assoc_id)%name)
                     call xplasma_errmsg_append(s,msgbuf)
                     call xplasma_errmsg_append(s, &
                          '  association not possible-- coordinates not disjoint.')
                  endif
               endif
            endif
            if(ier.eq.520) call xplasma_errmsg_append(s, &
                 ' error in xplasma_create_3dprof, profname='//trim(profname))
         endif
      endif
      if(ier.ne.0) return

      !-----------------
      !  check for CW/-coordinate reversal

      call xplasma_coord_isPeriodic(s,idx1,perio1,iertmp)
      call xplasma_coord_isPeriodic(s,idx2,perio2,iertmp)
      call xplasma_coord_isPeriodic(s,idx3,perio3,iertmp)

      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1 = ccwflag1

      mccwflag2=.TRUE.
      if(perio2.and.present(ccwflag2)) mccwflag2 = ccwflag2

      mccwflag3=.TRUE.
      if(perio3.and.present(ccwflag3)) mccwflag3 = ccwflag3

      iadx1=s%dict(idx1)%dindex
      iadx2=s%dict(idx2)%dindex
      iadx3=s%dict(idx3)%dindex

      jcoord1=s%grids(iadx1)%coord
      jcoord2=s%grids(iadx2)%coord
      jcoord3=s%grids(iadx3)%coord

      if((jcoord1.eq.jcoord2).or.(jcoord1.eq.jcoord3)) then
         ier=515
         write(msgbuf,*) '  ',trim(s%dict(idx1)%name),' coordinate is ', &
              trim(s%dict(jcoord1)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  ',trim(s%dict(idx2)%name),' coordinate is ', &
              trim(s%dict(jcoord2)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  ',trim(s%dict(idx3)%name),' coordinate is ', &
              trim(s%dict(jcoord3)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  profile name is: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      !------------------
      if(present(ispline).and.present(hspline)) then
         ier=90
         return
      endif

      if(size(data,1).eq.(idimx1-1)) then
         jspline=-1
      else
         jspline=0
      endif
      if(present(ispline)) jspline=ispline

      if(present(hspline)) then
         if(size(hspline).ne.3) then
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_create_3dprof: size(hspline).ne.3')
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
            ier=91
            return
         else
            jspline=hspline
         endif
      endif

      call spline_type_ck_encode(s,'xplasma_create_3dprof',s%dict(idf)%name, &
           jspline,kspline,ier)
      if(ier.ne.0) return

      if(jspline(1).lt.0) then
         idimx1a=idimx1-1
      else
         idimx1a=idimx1
      endif

      if(jspline(2).lt.0) then
         idimx2a=idimx2-1
      else
         idimx2a=idimx2
      endif

      if(jspline(3).lt.0) then
         idimx3a=idimx3-1
      else
         idimx3a=idimx3
      endif

      if(maxval(jspline).le.0) then
         if(present(coeffs)) then
            ier=523
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
            return
         endif
      endif

      if(present(coeffs)) then
         if(present(ibcx1a).or.present(ibcx1b).or. &
              present(zbcx1a).or.present(zbcx1b).or. &
              present(ibcx2a).or.present(ibcx2b).or. &
              present(zbcx2a).or.present(zbcx2b).or. &
              present(ibcx3a).or.present(ibcx3b).or. &
              present(zbcx3a).or.present(zbcx3b)) then
            ier=525
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
            return
         endif
      endif

      !------------------
      jbcx1a=0
      jbcx1b=0
      if(present(ibcx1a)) jbcx1a=ibcx1a
      if(present(ibcx1b)) jbcx1b=ibcx1b

      jbcx2a=0
      jbcx2b=0
      if(present(ibcx2a)) jbcx2a=ibcx2a
      if(present(ibcx2b)) jbcx2b=ibcx2b

      jbcx3a=0
      jbcx3b=0
      if(present(ibcx3a)) jbcx3a=ibcx3a
      if(present(ibcx3b)) jbcx3b=ibcx3b

      isize_exp=1
      do ii=1,3
         if(jspline(ii).le.0) then
            maxok=0
         else if(jspline(ii).eq.1) then
            maxok=1
            isize_exp=isize_exp*2
         else if(jspline(ii).eq.2) then
            maxok=2
            isize_exp=isize_exp*2
         endif

         if(ii.eq.1) then
            jbcmin=min(jbcx1a,jbcx1b)
            jbcmax=max(jbcx1a,jbcx1b)
         else if(ii.eq.2) then
            jbcmin=min(jbcx2a,jbcx2b)
            jbcmax=max(jbcx2a,jbcx2b)
         else
            jbcmin=min(jbcx3a,jbcx3b)
            jbcmax=max(jbcx3a,jbcx3b)
         endif

         if((jbcmin.lt.0).or.(jbcmax.gt.maxok)) then
            if(ier.eq.0) then
               write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
               call xplasma_errmsg_append(s,msgbuf)
            endif
            ier=504
            write(msgbuf,*) &
                 'boundary condition type argument value invalid:',jbcmin,jbcmax
            call xplasma_errmsg_append(s,msgbuf)
            if(maxok.eq.0) then
               write(msgbuf,*) ' for piecewise linear, no BC options allowed.'
            else if(maxok.eq.1) then
               write(msgbuf,*) &
                    ' for Hermite, only default (0) or 1st derivative (1) BCs are supported.'
            else if(maxok.eq.2) then
               write(msgbuf,*) &
                    ' BC option argument values between 0 and 2 expected.'
            endif
            call xplasma_errmsg_append(s,msgbuf)
         endif
      enddo
      if(ier.ne.0) return

      isize_exp=isize_exp-1
      if(present(coeffs)) then
         call ksign_3d(mccwflag1,mccwflag2,mccwflag3,jspline,isize_exp,ksign)
      endif

      !------------------
      !  if explicit boundary conditions are indicated, make sure that
      !  the supplied data is correctly sized...

      allocate(zzbcx1a(idimx2a,idimx3a),zzbcx1b(idimx2a,idimx3a))
      zzbcx1a=0
      zzbcx1b=0

      allocate(zzbcx2a(idimx1a,idimx3a),zzbcx2b(idimx1a,idimx3a))
      zzbcx2a=0
      zzbcx2b=0

      allocate(zzbcx3a(idimx1a,idimx2a),zzbcx3b(idimx1a,idimx2a))
      zzbcx3a=0
      zzbcx3b=0

      if(perio1) then
         jbcx1a=-1
         jbcx1b=-1
      else
         if(jbcx1a.gt.0) then
            if(present(zbcx1a)) then
               if((size(zbcx1a,1).ne.idimx2a).or. &
                    (size(zbcx1a,2).ne.idimx3a)) then
                  ier=502
                  write(msgbuf,*) 'size of x1(1) bdy condition data does',&
                       ' not match sizes of dimensions 2 & 3: ',idimx2a,idimx3a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do j=1,idimx3a
                     if(mccwflag3) then
                        jj=j
                     else
                        jj=idimx3a+1-j
                     endif
                     do i=1,idimx2a
                        if(mccwflag2) then
                           ii=i
                        else
                           ii=idimx2a+1-i
                        endif

                        zzbcx1a(ii,jj)=zbcx1a(i,j)
                     enddo
                  enddo
               endif
            endif
         endif
         if(jbcx1b.gt.0) then
            if(present(zbcx1b)) then
               if((size(zbcx1b,1).ne.idimx2a).or. &
                    (size(zbcx1b,2).ne.idimx3a)) then
                  ier=502
                  write(msgbuf,*) 'size of x1(Nx1) bdy condition data does',&
                       ' not match sizes of dimensions 2 & 3: ',idimx2a,idimx3a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do j=1,idimx3a
                     if(mccwflag3) then
                        jj=j
                     else
                        jj=idimx3a+1-j
                     endif
                     do i=1,idimx2a
                        if(mccwflag2) then
                           ii=i
                        else
                           ii=idimx2a+1-i
                        endif

                        zzbcx1b(ii,jj)=zbcx1b(i,j)
                     enddo
                  enddo
               endif
            endif
         endif
      endif

      if(perio2) then
         jbcx2a=-1
         jbcx2b=-1
      else
         if(jbcx2a.gt.0) then
            if(present(zbcx2a)) then
               if((size(zbcx2a,1).ne.idimx1a).or. &
                    (size(zbcx2a,2).ne.idimx3a)) then
                  ier=502
                  write(msgbuf,*) 'size of x2(1) bdy condition data does',&
                       ' not match sizes of dimensions 1 & 3: ',idimx1a,idimx3a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do j=1,idimx3a
                     if(mccwflag3) then
                        jj=j
                     else
                        jj=idimx3a+1-j
                     endif
                     do i=1,idimx1a
                        if(mccwflag1) then
                           ii=i
                        else
                           ii=idimx1a+1-i
                        endif

                        zzbcx2a(ii,jj)=zbcx2a(i,j)
                     enddo
                  enddo
               endif
            endif
         endif
         if(jbcx2b.gt.0) then
            if(present(zbcx2b)) then
               if((size(zbcx2b,1).ne.idimx1a).or. &
                    (size(zbcx2b,2).ne.idimx3a)) then
                  ier=502
                  write(msgbuf,*) 'size of x2(Nx2) bdy condition data does',&
                       ' not match sizes of dimensions 1 & 3: ',idimx1a,idimx3a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do j=1,idimx3a
                     if(mccwflag3) then
                        jj=j
                     else
                        jj=idimx3a+1-j
                     endif
                     do i=1,idimx1a
                        if(mccwflag1) then
                           ii=i
                        else
                           ii=idimx1a+1-i
                        endif

                        zzbcx2b(ii,jj)=zbcx2b(i,j)
                     enddo
                  enddo
               endif
            endif
         endif
      endif

      if(perio3) then
         jbcx3a=-1
         jbcx3b=-1
      else
         if(jbcx3a.gt.0) then
            if(present(zbcx3a)) then
               if((size(zbcx3a,1).ne.idimx1a).or. &
                    (size(zbcx3a,2).ne.idimx2a)) then
                  ier=502
                  write(msgbuf,*) 'size of x3(1) bdy condition data does',&
                       ' not match sizes of dimensions 1 & 2: ',idimx1a,idimx2a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do j=1,idimx2a
                     if(mccwflag2) then
                        jj=j
                     else
                        jj=idimx2a+1-j
                     endif
                     do i=1,idimx1a
                        if(mccwflag1) then
                           ii=i
                        else
                           ii=idimx1a+1-i
                        endif

                        zzbcx3a(ii,jj)=zbcx3a(i,j)
                     enddo
                  enddo
               endif
            endif
         endif
         if(jbcx3b.gt.0) then
            if(present(zbcx3b)) then
               if((size(zbcx3b,1).ne.idimx1a).or. &
                    (size(zbcx3b,2).ne.idimx2a)) then
                  ier=502
                  write(msgbuf,*) 'size of x3(Nx3) bdy condition data does',&
                       ' not match sizes of dimensions 1 & 2: ',idimx1a,idimx2a
                  call xplasma_errmsg_append(s,msgbuf)
               else
                  do j=1,idimx2a
                     if(mccwflag2) then
                        jj=j
                     else
                        jj=idimx2a+1-j
                     endif
                     do i=1,idimx1a
                        if(mccwflag1) then
                           ii=i
                        else
                           ii=idimx1a+1-i
                        endif

                        zzbcx3b(ii,jj)=zbcx3b(i,j)
                     enddo
                  enddo
               endif
            endif
         endif
      endif

      if(ier.ne.0) then
         deallocate(zzbcx1a,zzbcx1b,zzbcx2a,zzbcx2b,zzbcx3a,zzbcx3b)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      !------------------
      !  check data size

      if(size(data,1).ne.idimx1a) then
         ier=506
         write(msgbuf,*) ' ?xpprof: size of profile 1st dim = ',size(data,1)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   does not match x1 grid size  = ',idimx1, &
              ' x grid name: ',trim(s%dict(idx1)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         if(jspline(1).eq.-1) call xplasma_errmsg_append(s, &
              '    NOTE: step function data: need size = [grid size] - 1')
      endif

      if(size(data,2).ne.idimx2a) then
         ier=506
         write(msgbuf,*) ' ?xpprof: size of profile 2nd dim = ',size(data,2)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   does not match x2 grid size  = ',idimx2, &
              ' x grid name: ',trim(s%dict(idx2)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         if(jspline(2).eq.-1) call xplasma_errmsg_append(s, &
              '    NOTE: step function data: need size = [grid size] - 1')
      endif

      if(size(data,3).ne.idimx3a) then
         ier=506
         write(msgbuf,*) ' ?xpprof: size of profile 3rd dim = ',size(data,3)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   does not match x3 grid size  = ',idimx3, &
              ' x grid name: ',trim(s%dict(idx3)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         if(jspline(3).eq.-1) call xplasma_errmsg_append(s, &
              '    NOTE: step function data: need size = [grid size] - 1')
      endif

      if(present(coeffs)) then
         if(size(coeffs,1).ne.isize_exp) then
            ier=524
            write(msgbuf,*) ' ?xpprof: 1st dim size of coeffs: ',size(coeffs,1)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  does not match expected coeffs/node: ',isize_exp
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  expected: 2**N-1, where N = #of spline or hermite dimensions.'
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(size(coeffs,2).ne.idimx1a) then
            ier=524
            write(msgbuf,*) ' ?xpprof: size of 2nd dim of coeffs data = ',size(coeffs,2)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '   does not match x1 grid size  = ',idimx1, &
                 ' x1 grid name: ',trim(s%dict(idx1)%name)
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(size(coeffs,3).ne.idimx2a) then
            ier=524
            write(msgbuf,*) ' ?xpprof: size of 3rd dim of coeffs data = ',size(coeffs,3)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '   does not match x2 grid size  = ',idimx2, &
                 ' x2 grid name: ',trim(s%dict(idx2)%name)
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(size(coeffs,4).ne.idimx3a) then
            ier=524
            write(msgbuf,*) ' ?xpprof: size of 4th dim of coeffs data = ',size(coeffs,4)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '   does not match x3 grid size  = ',idimx3, &
                 ' x3 grid name: ',trim(s%dict(idx3)%name)
            call xplasma_errmsg_append(s,msgbuf)
         endif
      endif

      if(ier.ne.0) then
         deallocate(zzbcx1a,zzbcx1b,zzbcx2a,zzbcx2b,zzbcx3a,zzbcx3b)
         write(msgbuf,*) '   profile name: ',trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      !------------------
      ! create buffer space for spline coefficients

      imul=1
      do ii=1,3
         if(jspline(ii).gt.0) then
            imul=imul*2
         endif
      enddo

      isize=idimx1a*idimx2a*idimx3a*imul

      s%profs(iadr)%rank = 3
      s%profs(iadr)%gridIds = 0
      s%profs(iadr)%gridIds(1) = idx1
      s%profs(iadr)%gridIds(2) = idx2
      s%profs(iadr)%gridIds(3) = idx3
      s%profs(iadr)%profIds = 0
      if(present(assoc_id)) then
         s%profs(iadr)%profIds(1) = assoc_id
         if(assoc_id.gt.0) then
            iada = s%dict(assoc_id)%dindex
            s%profs(iada)%profIds(1) = idf
         endif
      endif
      s%profs(iadr)%kspline = kspline
      s%profs(iadr)%size = isize
      allocate(s%profs(iadr)%buf(isize))

      !------------------
      ! allocate tmp space for data; reverse if flag is set

      allocate(dtmp(idimx1a,idimx2a,idimx3a))
      if(mccwflag1.and.mccwflag2.and.mccwflag3) then
         dtmp = data
      else
         if(present(coeffs)) allocate(ctmp(isize_exp,idimx1a,idimx2a,idimx3a))
         do k=1,idimx3a
            if(mccwflag3) then
               kk=k
            else
               kk=idimx3a+1-k
            endif
            do j=1,idimx2a
               if(mccwflag2) then
                  jj=j
               else
                  jj=idimx2a+1-j
               endif
               do i=1,idimx1a
                  if(mccwflag1) then
                     ii=i
                  else
                     ii=idimx1a+1-i
                  endif
                  dtmp(ii,jj,kk)=data(i,j,k)
                  if(present(coeffs)) then
                     ctmp(1:isize_exp,ii,jj,kk) = &
                          ksign(1:isize_exp)*coeffs(1:isize_exp,i,j,k)
                  endif
               enddo
            enddo
         enddo
      endif

      !------------------
      ! compute spline coefficients as needed

      ii=-imul+1
      do k=1,idimx3a
         do j=1,idimx2a
            do i=1,idimx1a
               ii=ii+imul
               s%profs(iadr)%buf(ii)=dtmp(i,j,k)
            enddo
         enddo
      enddo

      if(imul.eq.1) then
         continue   ! pclin or zone data already in place
         ! no coefficients to be evaluated

      else if(present(coeffs)) then

         ii=-imul+1
         do k=1,idimx3a
            do j=1,idimx2a
               do i=1,idimx1a
                  ii=ii+imul
                  if(mccwflag1.and.mccwflag2.and.mccwflag3) then
                     s%profs(iadr)%buf(ii+1:ii+isize_exp) = &
                          coeffs(1:isize_exp,i,j,k)
                  else
                     s%profs(iadr)%buf(ii+1:ii+isize_exp) = &
                          ctmp(1:isize_exp,i,j,k)
                  endif
               enddo
            enddo
         enddo

      else

         call r8mkintrp3d(s%grids(iadx1)%xpkg(1:idimx1,1),idimx1, &
              s%grids(iadx2)%xpkg(1:idimx2,1),idimx2, &
              s%grids(iadx3)%xpkg(1:idimx3,1),idimx3, &
              jspline, &
              s%profs(iadr)%buf,imul,idimx1a,idimx2a,idimx3a, &
              jbcx1a,zzbcx1a,jbcx1b,zzbcx1b, &
              jbcx2a,zzbcx2a,jbcx2b,zzbcx2b, &
              jbcx3a,zzbcx3a,jbcx3b,zzbcx3b, &
              ier)
         if(ier.ne.0) then
            ier=9999
            msgbuf = &
                 ' unexpected failure in r8mkbicub, xplasma_create_3dprof'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
            call xplasma_errmsg_append(s,msgbuf)
         endif

      endif

      deallocate(dtmp,zzbcx1a,zzbcx1b,zzbcx2a,zzbcx2b,zzbcx3a,zzbcx3b)
      if(allocated(ctmp)) deallocate(ctmp)

      if(ier.le.0) then
         s%prof_counter = s%prof_counter + 1
         s%profs(iadr)%prof_counter = s%prof_counter
         s%grids(iadx1)%nrefs = s%grids(iadx1)%nrefs + 1
         s%grids(iadx2)%nrefs = s%grids(iadx2)%nrefs + 1
         s%grids(iadx3)%nrefs = s%grids(iadx3)%nrefs + 1
      endif

      iertmp=0

      if(present(label)) then
         call xplasma_label_item(s,idf,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_label_item(s,idf,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      ier = iertmp

      if((iflag.gt.0).and.(ier.eq.0)) then
         call xplasma_internal_ids(s)
      endif

      if(ier.ne.0) then
         write(msgbuf,*) ' ?xpprof: profname: '//trim(profname)
         call xplasma_errmsg_append(s,msgbuf)
      endif

    end subroutine xplasma_create_3dprof

    subroutine spline_type_ck_encode(s,slbl,name,jspline,kspline,ier)
      
      ! PRIVATE
      ! check & encode hybrid spline codes
      !   jspline(1) = spline code for 1st dimension
      !   jspline(2) = spline code for 2nd dimension
      !     ...

      !     -1: zonal step function
      !      0: piecewise linear
      !      1: Hermite (C1) spline
      !      2: Cubic (C2) spline

      ! if not all jspline elements are the same, then, 2 values can be present
      ! one drawn from the set {-1,0} and the other from the set {1,2}.

      type (xplasma), pointer :: s
      character*(*), intent(in) :: slbl               ! subroutine label
      character*(*), intent(in) :: name               ! profile name
      integer, dimension(:), intent(in) :: jspline    ! hybrid spline codes
      integer, intent(out) :: kspline                 ! combined spline code
      integer, intent(out) :: ier                     ! status code (0=OK)

      !---------------------------------
      logical :: present(-1:2)
      integer :: ii,imin,imax,inum
      character*128 msgbuf
      !---------------------------------

      present = .FALSE.
      ier = 0
      kspline = -99

      do ii=1,size(jspline)
         if((jspline(ii).lt.-1).or.(jspline(ii).gt.2)) then
            ier=505
            write(msgbuf,*) ' '//trim(slbl)//': '// &
                 'interpolation fit method argument value(s) invalid: ',jspline
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  use -1 for step fcn; '// &
                 '0 for piecewise linear, 1 for Hermite, 2 for spline.'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '   profile name: ',trim(name)
            call xplasma_errmsg_append(s,msgbuf)
            exit
         endif
         present(jspline(ii))=.TRUE.
      enddo
      if(ier.ne.0) return

      inum=0
      imin=3
      imax=-2

      do ii=-1,2
         if(present(ii)) then
            inum=inum+1
            imin=min(imin,ii)
            imax=max(imax,ii)
         endif
      enddo

      if(inum.eq.1) then
         kspline=imin  ! =imax; only 1 code found.
         return

      else if(present(1).and.present(2)) then
         ier=505
         msgbuf=' '
         write(msgbuf,*) ' '//trim(slbl)//': '// &
              'interpolation fit method argument values invalid: ',jspline
         call xplasma_errmsg_append(s,msgbuf)
         call xplasma_errmsg_append(s,'  Hermite/Spline hybrid not supported.')
      endif

      if(ier.ne.0) then
         msgbuf=' '
         write(msgbuf,*) ' '//trim(slbl)//': '// &
              'interpolation fit combination invalid: ',jspline
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '   profile name: ',trim(name)
      else
         ! encode
         kspline=0
         do ii=1,size(jspline)
            kspline = 100*kspline + (4+jspline(ii))
         enddo
      endif

    end subroutine spline_type_ck_encode

    subroutine spline_type_decode(kspline,jspline,irank)
      ! PRIVATE
      !  decode kspline to get fit order along each dimension

      integer, intent(in) :: kspline
      integer, intent(in) :: irank
      integer, intent(out) :: jspline(irank)

      !-----------------------
      ! local
      integer :: idenom(irank),ii,ispline
      !-----------------------

      if(kspline.lt.5) then
         jspline=kspline
         return
      endif

      !  Hybrid

      ispline=kspline

      idenom(irank)=1
      do ii=irank-1,1,-1
         idenom(ii)=idenom(ii+1)*100
      enddo

      do ii=1,irank
         jspline(ii) = ispline/idenom(ii)
         ispline = ispline - jspline(ii)*idenom(ii)
         jspline(ii) = jspline(ii) - 4
      enddo
    end subroutine spline_type_decode

    !---------------------------------------------------------
    subroutine xplasma_prof_lintrans(s,id,ierr, factor, offset)

      ! apply an in-place linear transformation to a profile:
      !   f -> factor*f + offset

      ! data values have offset applied; derivatives only have factor applied.

      type (xplasma), pointer :: s
      integer, intent(in) :: id    ! profile ID
      integer, intent(out) :: ierr ! status code returned 0 = normal

      real*8, intent(in), optional :: factor   ! multiplicative factor
      real*8, intent(in), optional :: offset   ! additive offset

      !----------------------------------
      integer :: iadr,isize,imul,ii,ia,irank,ispline
      integer, dimension(:), allocatable :: hybridType
      real*8 :: zfactor,zoffset
      !----------------------------------

      call xplasma_errck1_prof(s,id,ierr)
      if(ierr.ne.0) return

      call xplasma_author_check(s,id,ierr)
      if(ierr.ne.0) return

      if(.not.(present(factor).or.present(offset))) then
         return  ! no transformation to be done
      endif

      zoffset = 0.0d0
      zfactor = 1.0d0

      if(present(offset)) zoffset = offset
      if(present(factor)) zfactor = factor

      ! OK: execute: f -> zfactor*f + zoffset 

      iadr = s%dict(id)%dindex

      irank = s%profs(iadr)%rank
      isize = s%profs(iadr)%size
      ispline = s%profs(iadr)%kspline

      ! OK have info needed to apply transformation

      if(zfactor.ne.1.0d0) then
         ! apply multiplicative factor
         do ii=1,isize
            s%profs(iadr)%buf(ii) = zfactor*s%profs(iadr)%buf(ii)
         enddo
      endif

      if(zoffset.ne.0.0d0) then
         ! apply offset
         ! first find spacing of data points

         allocate(hybridType(irank))
         call spline_type_decode(ispline,hybridType,irank)

         imul=1
         do ii=1,irank
            if(hybridType(ii).gt.0) imul=imul*2
         enddo

         deallocate(hybridType)

         do ii=1,isize,imul
            s%profs(iadr)%buf(ii) = s%profs(iadr)%buf(ii) + zoffset
         enddo
      endif

    end subroutine xplasma_prof_lintrans

    !---------------------------------------------------------
    subroutine xplasma_psi_range(s,psimin,psimax)

      !  return poloidal flux psimin (psi on grid) and psimax (psi at edge)
      !  Wb/rad.  Psimin < Psimax -- profile is really |Psi|.

      !  If an error occurs, Psimin=Psimax=0 is returned.

      !  Sign information:
      !  The direction of the plasma current (jphi_ccw) can be
      !  had via an xplasma_global_info call.

      type (xplasma), pointer :: s
      real*8, intent(out) :: psimin,psimax  ! Psi range Wb/rad -- poloidal flux

      !-----------------------------
      integer :: id_psi,iertmp,iadr,iadx,inx
      integer :: isize,imul
      !-----------------------------

      call xplasma_common_ids(s,iertmp,id_psi=id_psi)

      if((iertmp.ne.0).or.(id_psi.eq.0)) then
         psimin = 0
         psimax = 0

      else

         iadr = s%dict(id_psi)%dindex

         iadx = s%profs(iadr)%gridIds(1)  ! get grid id
         iadx = s%dict(iadx)%dindex  ! get grid address
         inx = s%grids(iadx)%size    ! get grid size

         isize=s%profs(iadr)%size
         imul=isize/inx

         psimin = s%profs(iadr)%buf(1)
         psimax = s%profs(iadr)%buf(1+(inx-1)*imul)
      endif

    end subroutine xplasma_psi_range

    !---------------------------------------------------------
    subroutine xplasma_hermitize(s,id,ier)

      !  convert a Spline representation to a Hermite representation
      !  if the indicated profile is not a spline, silently do nothing.

      type (xplasma), pointer :: s
      integer, intent(in) :: id                   ! profile id
      integer, intent(out) :: ier                 ! status code on exit (0=OK)
      
      !------------------------
      integer :: iadr,irank,ispline,jspline(3),jsplndx(3)
      integer :: ii,jj,kk,mm,indx,inc,inc1
      integer, dimension(:,:), allocatable :: iderivs
      integer :: isize,idg1,idg2,idg3,inx1,inx2,inx3,ilim1
      real*8, dimension(:), allocatable :: xg1,xg2,xg3,vec2,vec3,buf
      real*8, dimension(:,:), allocatable :: dbuf
      !------------------------

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      irank = s%profs(iadr)%rank
      isize = s%profs(iadr)%size

      do
         idg1 = s%profs(iadr)%gridIds(1)
         call xplasma_grid_size(s,idg1,inx1,ier)
         if(ier.ne.0) exit
         allocate(xg1(inx1))
         call xplasma_grid(s,idg1,xg1,ier)
         if(ier.ne.0) exit

         if(irank.ge.2) then
            idg2 = s%profs(iadr)%gridIds(2)
            call xplasma_grid_size(s,idg2,inx2,ier)
            if(ier.ne.0) exit
            allocate(xg2(inx2))
            call xplasma_grid(s,idg2,xg2,ier)
            if(ier.ne.0) exit
         endif
         
         if(irank.ge.3) then
            idg3 = s%profs(iadr)%gridIds(3)
            call xplasma_grid_size(s,idg3,inx3,ier)
            if(ier.ne.0) exit
            allocate(xg3(inx3))
            call xplasma_grid(s,idg3,xg3,ier)
            if(ier.ne.0) exit
         endif

         if(irank.eq.1) then

            ! simple case of 1d profile

            ispline = s%profs(iadr)%kspline
            if(ispline.le.1) exit  ! NOT a spline

            ! evaluate derivatives at grid points

            allocate(vec2(inx1))
            call xplasma_eval_1dprofxs(s,id,xg1,vec2,ier,ideriv1=1)
            if(ier.ne.0) exit

            ! build buffer with values & derivatives (orig. values taken
            ! from old buffer)

            allocate(buf(isize))
            indx = 0
            do ii=1,inx1
               indx=indx+1
               buf(indx) = s%profs(iadr)%buf(indx)
               indx=indx+1
               buf(indx) = vec2(ii)
            enddo

            s%profs(iadr)%kspline = 1  ! Hermite now
            s%profs(iadr)%buf = buf    ! replace data buffer

         else

            allocate(vec2(inx1),vec3(inx1))

            ispline = s%profs(iadr)%kspline
            call spline_type_decode(ispline,jspline(1:irank),irank)
            if(minval(jspline(1:irank)).le.1) exit  ! NOT a spline

            inc=1
            indx=0
            do ii=1,irank
               if(jspline(ii).eq.2) then
                  inc=inc*2
                  indx=indx+1
                  jsplndx(indx)=ii
                  jspline(ii)=1  ! changing this dim to Hermite
               endif
            enddo

            !  calculate new spline type code now...

            if(indx.lt.irank) then
               call spline_type_ck_encode(s,'xplasma_hermitize','(id)', &
                    jspline,ispline,ier)
               if(ier.ne.0) exit
            else
               ispline=1
            endif

            inc1=inc    ! 2**N
            inc=inc1-1  ! 2**N - 1, where N = #spline dimensions

            allocate(iderivs(irank,inc)); iderivs = 0

            if(inc.eq.1) then
               indx=jsplndx(1)
               iderivs(indx,1)=1

            else if(inc.eq.3) then
               indx=jsplndx(1)
               iderivs(indx,1)=1    ! 1,0
               iderivs(indx,3)=1    ! 1,1  1st 1 set
               indx=jsplndx(2)
               iderivs(indx,2)=1    ! 0,1
               iderivs(indx,3)=1    ! 1,1  2nd 1 set

            else if(inc.eq.7) then
               !  derivate patterns: (1,0,0),(0,1,0),(0,0,1)
               !                     (1,1,0),(1,0,1),(0,1,1)  ... and (1,1,1)
               indx=jsplndx(1)
               iderivs(indx,1)=1    ! 1,0,0
               iderivs(indx,4)=1    ! 1,1,0
               iderivs(indx,5)=1    ! 1,0,1
               iderivs(indx,7)=1    ! 1,1,1

               indx=jsplndx(2)
               iderivs(indx,2)=1    ! 0,1,0
               iderivs(indx,4)=1    ! 1,1,0
               iderivs(indx,6)=1    ! 0,1,1
               iderivs(indx,7)=1    ! 1,1,1

               indx=jsplndx(3)
               iderivs(indx,3)=1    ! 0,0,1
               iderivs(indx,5)=1    ! 1,0,1
               iderivs(indx,6)=1    ! 0,1,1
               iderivs(indx,7)=1    ! 1,1,1

            else
               ier = 9999
               call xplasma_errmsg_append(s,'?xplasma_hermitize: spline rank exceeds 3d limit.')
               exit
            endif

            indx=-inc
            allocate(buf(isize))

            ilim1=1
            if(irank.eq.3) ilim1=inx3

            allocate(dbuf(inx1,inc))

            do kk=1,ilim1
               if(irank.eq.3) vec3 = xg3(kk)

               do jj=1,inx2
                  vec2 = xg2(jj)

                  do mm=1,inc
                     if(irank.eq.2) then
                        call xplasma_eval_2dprofxs(s,id, &
                             idg1,xg1, idg2,vec2, dbuf(1:inx1,mm), ier, &
                             ideriv1=iderivs(1,mm), ideriv2=iderivs(2,mm))
                             
                     else if(irank.eq.3) then
                        call xplasma_eval_3dprofxs(s,id, &
                             idg1,xg1, idg2,vec2, idg3,vec3, dbuf(1:inx1,mm), &
                             ier, &
                             ideriv1=iderivs(1,mm), &
                             ideriv2=iderivs(2,mm), &
                             ideriv3=iderivs(3,mm))

                     endif
                     if(ier.ne.0) exit
                  enddo
                  if(ier.ne.0) exit

                  do ii=1,inx1
                     indx=indx+inc1
                     buf(indx) = s%profs(iadr)%buf(indx)  ! copy data value
                     buf(indx+1:indx+inc) = dbuf(ii,1:inc)
                  enddo

               enddo
               if(ier.ne.0) exit
            enddo
            if(ier.ne.0) exit

            s%profs(iadr)%kspline = ispline  ! Hermite now
            s%profs(iadr)%buf = buf          ! replace data buffer

         endif

         exit
      enddo

      if(allocated(xg1)) deallocate(xg1)
      if(allocated(xg2)) deallocate(xg2)
      if(allocated(xg3)) deallocate(xg3)

      if(allocated(vec2)) deallocate(vec2)
      if(allocated(vec3)) deallocate(vec3)

      if(allocated(buf)) deallocate(buf)
      if(allocated(dbuf)) deallocate(dbuf)
      if(allocated(iderivs)) deallocate(iderivs)

    end subroutine xplasma_hermitize
    !---------------------------------------------------------
    subroutine xplasma_getProf_1data(s,id,zdata,ier, ccwflag1, zcoeff)
      
      !  return the original data associated with this profile

      type (xplasma), pointer :: s
      integer, intent(in) :: id                   ! profile id
      real*8, dimension(:), intent(out) :: zdata  ! data returned
      integer, intent(out) :: ier                 ! completion code 0=OK

      logical, intent(in), optional :: ccwflag1   ! CCW flag (default T)

      !  spline coefficients or nodal derivatives for Hermite function
      real*8, dimension(:), intent(out),  optional :: zcoeff

      !  **it is an error if zcoeff is present for a piecewise
      !    linear or zonal function

      !------------------------
      integer :: iadr,idg,iadx,inx,inxa,imul,isize,i,ix,ispline,jsign
      integer :: incr,iertmp
      character*128 msgbuf
      logical :: jccw,perio
      !------------------------

      zdata = 0
      if(present(zcoeff)) zcoeff = 0

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      if(s%profs(iadr)%rank.ne.1) then
         ier=512
         write(msgbuf,*) '  profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  profile rank: ',s%profs(iadr)%rank
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  passed array rank:  1 (xplasma_getProf_1data)'
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      ispline = s%profs(iadr)%kspline
      if(present(zcoeff).and.(ispline.le.0)) then
         ier=521
         write(msgbuf,*) '  profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  only Hermite & Spline functions have coefficients.'
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      idg = s%profs(iadr)%gridIds(1)  ! get grid id
      iadx = s%dict(idg)%dindex  ! get grid address
      inx = s%grids(iadx)%size    ! get grid size

      if(ispline.lt.0) then
         inxa=inx-1
      else
         inxa=inx
      endif

      if(size(zdata).ne.inxa) then
         ier=506
         write(msgbuf,*) '  profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  grid size:         ',inx
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  passed array size: ',size(zdata)
         call xplasma_errmsg_append(s,msgbuf)
         if(ispline.lt.0) call xplasma_errmsg_append(s, &
              '  NOTE: step function data, size = [grid size] - 1')
         return
      endif

      if(present(zcoeff)) then
         if(size(zcoeff).ne.inxa) then
            ier=506
            write(msgbuf,*) '  profile name: ',trim(s%dict(id)%name)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  grid size:         ',inx
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  passed "zcoeff" array size: ',size(zcoeff)
            call xplasma_errmsg_append(s,msgbuf)
            return
         endif
      endif

      jccw=.TRUE.
      jsign=1
      if(present(ccwflag1)) then
         if(.not.ccwflag1) then
            call xplasma_grid_info(s,idg,iertmp, perio=perio)
            if(perio) then
               jccw=ccwflag1
            endif
         endif
      endif
      if(.not.jccw) then
         if(ispline.eq.1) then
            jsign=-1
         endif
      endif

      isize=s%profs(iadr)%size
      imul=isize/inxa
      
      if(jccw) then
         ix=0
         incr=1
      else
         ix=inxa+1
         incr=-1
      endif
      do i=1,isize,imul
         ix=ix+incr
         zdata(ix)=s%profs(iadr)%buf(i)
         if(present(zcoeff)) then
            zcoeff(ix)=jsign*s%profs(iadr)%buf(i+1)
         endif
      enddo

    end subroutine xplasma_getProf_1data

    !---------------------------------------------------------
    subroutine xplasma_getProf_2data(s,id,zdata,ier, &
         ccwflag1,ccwflag2, zcoeff)
      
      !  return the original data associated with this profile

      type (xplasma), pointer :: s
      integer, intent(in) :: id                     ! profile id
      real*8, dimension(:,:), intent(out) :: zdata  ! data returned
      integer, intent(out) :: ier                   ! completion code 0=OK

      logical, intent(in), optional :: ccwflag1   ! dim. 1 CCW flag (default T)
      logical, intent(in), optional :: ccwflag2   ! dim. 2 CCW flag (default T)

      !  spline coefficients or nodal derivatives for Hermite function
      real*8, dimension(:,:,:), intent(out),  optional :: zcoeff

      !  **it is an error if zcoeff is present for a piecewise
      !    linear or zonal function; for other functions, the 1st dimension
      !    size must be 2**N - 1 where N is the number of dimensions along
      !    which the data is Spline or Hermite interpolated: 1 or 3, here.

      !------------------------
      integer :: iadr,iadx1,inx1,iadx2,inx2,imul,isize,i,ix1,ix2,jj
      integer :: i1,i2,iertmp,hspline(2),ksign(3),isize_exp
      integer :: idg1,idg2,incr1,incr2
      integer :: ispline,inx1a,inx2a
      character*128 msgbuf
      logical :: perio,jccw1,jccw2
      !------------------------

      zdata = 0
      if(present(zcoeff)) zcoeff = 0

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      if(s%profs(iadr)%rank.ne.2) then
         ier=512
         write(msgbuf,*) '  profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  profile rank: ',s%profs(iadr)%rank
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  passed array rank:  2 (xplasma_getProf_2data)'
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      idg1 = s%profs(iadr)%gridIds(1)  ! get grid id
      iadx1 = s%dict(idg1)%dindex  ! get grid address
      inx1 = s%grids(iadx1)%size    ! get grid size

      idg2 = s%profs(iadr)%gridIds(2)  ! get grid id
      iadx2 = s%dict(idg2)%dindex  ! get grid address
      inx2 = s%grids(iadx2)%size    ! get grid size

      call spline_type_decode(s%profs(iadr)%kspline,hspline,2)
      ispline = minval(hspline)

      if(present(zcoeff).and.(maxval(hspline).le.0)) then
         ier=521
         write(msgbuf,*) '  2d profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  only Hermite & Spline functions have coefficients.'
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      if(hspline(1).lt.0) then
         inx1a=inx1-1
      else
         inx1a=inx1
      endif

      if(hspline(2).lt.0) then
         inx2a=inx2-1
      else
         inx2a=inx2
      endif

      if((size(zdata,1).ne.inx1a).or.(size(zdata,2).ne.inx2a)) then
         ier=506
         write(msgbuf,*) '  2d profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  2d grid dimensions:    ',inx1,inx2
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  2d passed array sizes: ',size(zdata,1),size(zdata,2)
         if(ispline.lt.0) call xplasma_errmsg_append(s, &
              '  NOTE: step function data dims need size = [grid size] - 1')
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      if(present(zcoeff)) then
         isize_exp=1
         if(hspline(1).gt.0) isize_exp=isize_exp*2
         if(hspline(2).gt.0) isize_exp=isize_exp*2
         isize_exp=isize_exp-1

         if(size(zcoeff,1).ne.isize_exp) then
            ier=522
            write(msgbuf,*) '  2d profile name: ',trim(s%dict(id)%name)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  expecting size of 1st "zcoeff" dimension to be 2**N-1, where N=#of spline or hermite dimensions.'
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if((size(zcoeff,2).ne.inx1a).or.(size(zcoeff,3).ne.inx2a)) then
            ier=506
            write(msgbuf,*) '  size mismatch, 2d profile name: ',trim(s%dict(id)%name)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  2d grid dimensions:    ',inx1,inx2
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  3d passed "zcoeff" array dims 2 & 3: ',size(zcoeff,2),size(zcoeff,3)
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(ier.gt.0) return
      endif

      jccw1=.TRUE.
      if(present(ccwflag1)) then
         if(.not.ccwflag1) then
            call xplasma_grid_info(s,idg1,iertmp, perio=perio)
            if(perio) then
               jccw1=ccwflag1
            endif
         endif
      endif

      jccw2=.TRUE.
      if(present(ccwflag2)) then
         if(.not.ccwflag2) then
            call xplasma_grid_info(s,idg2,iertmp, perio=perio)
            if(perio) then
               jccw2=ccwflag2
            endif
         endif
      endif

      isize=s%profs(iadr)%size
      imul=isize/(inx1a*inx2a)

      if(present(zcoeff)) then
         call ksign_2d(jccw1,jccw2,hspline,isize_exp,ksign)
      endif
      
      ix2=1
      ix1=0
      do i=1,isize,imul
         if(ix1.eq.inx1a) then
            ix2=ix2+1
            ix1=1
         else
            ix1=ix1+1
         endif

         if(jccw1) then
            i1=ix1
         else
            i1=inx1a+1-ix1
         endif

         if(jccw2) then
            i2=ix2
         else
            i2=inx2a+1-ix2
         endif

         zdata(i1,i2)=s%profs(iadr)%buf(i)

         if(present(zcoeff)) then
            do jj=1,imul-1
               zcoeff(jj,i1,i2)=ksign(jj)*s%profs(iadr)%buf(i+jj)
            enddo
         endif
      enddo

    end subroutine xplasma_getProf_2data

    subroutine ksign_2d(jccw1,jccw2,hspline,isize_coeff,ksign)

      ! PRIVATE
      ! 2d sign factor for derivatives, depending on CCW reversal of
      ! periodic coordinates-- affects Hermite dimensions only

      logical, intent(in) :: jccw1,jccw2  ! CCW flags, dims 1 & 2
      integer, intent(in) :: hspline(2)   ! Hermite/Spline flags, dims 1 & 2
      integer, intent(in) :: isize_coeff  ! #of coefficients per node point
      !  1 if Spline/Hermite in 1 dim only; 3 if in 2 dims
      !  (could have computed this from hspline but callers already know)

      integer, intent(out) :: ksign(3)    ! output coefficient sign factors

      ! spline coefficients, involving 2nd derivatives and products of
      ! 2nd derivatives, always get a sign factor of +1

      !------------------------------
      ! local:
      integer :: jsign(2)
      !------------------------------

      jsign=1
      ksign=1
      if((.not.jccw1).and.(hspline(1).eq.1)) jsign(1)=-1
      if((.not.jccw2).and.(hspline(2).eq.1)) jsign(2)=-1
      if(isize_coeff.eq.1) then
         if(hspline(1).gt.0) ksign(1)=jsign(1)
         if(hspline(2).gt.0) ksign(1)=jsign(2)
      else
         ksign(1)=jsign(1)
         ksign(2)=jsign(2)
         ksign(3)=jsign(1)*jsign(2)
      endif

    end subroutine ksign_2d

    !---------------------------------------------------------
    subroutine xplasma_getProf_3data(s,id,zdata,ier, &
         ccwflag1,ccwflag2,ccwflag3, zcoeff)
      
      !  return the original data associated with this profile

      type (xplasma), pointer :: s
      integer, intent(in) :: id                       ! profile id
      real*8, dimension(:,:,:), intent(out) :: zdata  ! data returned
      integer, intent(out) :: ier                     ! completion code 0=OK

      logical, intent(in), optional :: ccwflag1   ! dim. 1 CCW flag (default T)
      logical, intent(in), optional :: ccwflag2   ! dim. 2 CCW flag (default T)
      logical, intent(in), optional :: ccwflag3   ! dim. 3 CCW flag (default T)

      !  spline coefficients or nodal derivatives for Hermite function
      real*8, dimension(:,:,:,:), intent(out),  optional :: zcoeff

      !  **it is an error if zcoeff is present for a piecewise
      !    linear or zonal function; for other functions, the 1st dimension
      !    size must be 2**N - 1 where N is the number of dimensions along
      !    which the data is Spline or Hermite interpolated: 1 3 or 7 here.

      !------------------------
      integer :: iadr,iadx1,inx1,iadx2,inx2,iadx3,inx3,imul,isize,i,ix1,ix2,ix3
      integer :: idg1,idg2,idg3,i1,i2,i3,iertmp
      integer :: isize_exp,jj,ksign(7)
      logical :: jccw1,jccw2,jccw3,perio
      integer :: ispline,hspline(3),inx1a,inx2a,inx3a
      character*128 msgbuf
      !------------------------

      zdata = 0
      if(present(zcoeff)) zcoeff = 0

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      if(s%profs(iadr)%rank.ne.3) then
         ier=512
         write(msgbuf,*) '  profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  profile rank: ',s%profs(iadr)%rank
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  passed array rank:  3 (xplasma_getProf_3data)'
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      idg1 = s%profs(iadr)%gridIds(1)  ! get grid id
      iadx1 = s%dict(idg1)%dindex  ! get grid address
      inx1 = s%grids(iadx1)%size    ! get grid size

      idg2 = s%profs(iadr)%gridIds(2)  ! get grid id
      iadx2 = s%dict(idg2)%dindex  ! get grid address
      inx2 = s%grids(iadx2)%size    ! get grid size

      idg3 = s%profs(iadr)%gridIds(3)  ! get grid id
      iadx3 = s%dict(idg3)%dindex  ! get grid address
      inx3 = s%grids(iadx3)%size    ! get grid size

      call spline_type_decode(s%profs(iadr)%kspline,hspline,3)
      ispline = minval(hspline)

      if(present(zcoeff).and.(maxval(hspline).le.0)) then
         ier=521
         write(msgbuf,*) '  3d profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  only Hermite & Spline functions have coefficients.'
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      if(hspline(1).lt.0) then
         inx1a=inx1-1
      else
         inx1a=inx1
      endif

      if(hspline(2).lt.0) then
         inx2a=inx2-1
      else
         inx2a=inx2
      endif

      if(hspline(3).lt.0) then
         inx3a=inx3-1
      else
         inx3a=inx3
      endif

      if((size(zdata,1).ne.inx1a).or.(size(zdata,2).ne.inx2a).or. &
           (size(zdata,3).ne.inx3a)) then
         ier=506
         write(msgbuf,*) '  3d profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  3d grid dimensions:    ',inx1,inx2,inx3
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  3d passed array sizes: ',size(zdata,1),size(zdata,2),size(zdata,3)
         call xplasma_errmsg_append(s,msgbuf)
         if(ispline.lt.0) call xplasma_errmsg_append(s, &
              '  NOTE: step function data dim needs size = [grid size] - 1')
         return
      endif

      if(present(zcoeff)) then
         isize_exp=1
         if(hspline(1).gt.0) isize_exp=isize_exp*2
         if(hspline(2).gt.0) isize_exp=isize_exp*2
         if(hspline(3).gt.0) isize_exp=isize_exp*2
         isize_exp=isize_exp-1

         if(size(zcoeff,1).ne.isize_exp) then
            ier=522
            write(msgbuf,*) '  3d profile name: ',trim(s%dict(id)%name)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  expecting size of 1st "zcoeff" dimension to be 2**N-1, where N=#of spline or hermite dimensions.'
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if((size(zcoeff,2).ne.inx1a).or.(size(zcoeff,3).ne.inx2a).or. &
              (size(zcoeff,4).ne.inx3a)) then
            ier=506
            write(msgbuf,*) '  size mismatch, 3d profile name: ',trim(s%dict(id)%name)
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  3d grid dimensions:    ',inx1,inx2,inx3
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  4d passed "zcoeff" array dims 2--4: ',size(zcoeff,2),size(zcoeff,3),size(zcoeff,4)
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(ier.gt.0) return
      endif

      jccw1=.TRUE.
      if(present(ccwflag1)) then
         if(.not.ccwflag1) then
            call xplasma_grid_info(s,idg1,iertmp, perio=perio)
            if(perio) then
               jccw1=ccwflag1
            endif
         endif
      endif

      jccw2=.TRUE.
      if(present(ccwflag2)) then
         if(.not.ccwflag2) then
            call xplasma_grid_info(s,idg2,iertmp, perio=perio)
            if(perio) then
               jccw2=ccwflag2
            endif
         endif
      endif

      jccw3=.TRUE.
      if(present(ccwflag3)) then
         if(.not.ccwflag3) then
            call xplasma_grid_info(s,idg3,iertmp, perio=perio)
            if(perio) then
               jccw3=ccwflag3
            endif
         endif
      endif

      isize=s%profs(iadr)%size
      imul=isize/(inx1a*inx2a*inx3a)
 
      if(present(zcoeff)) then
         call ksign_3d(jccw1,jccw2,jccw3,hspline,isize_exp,ksign)
      endif
     
      ix3=1
      ix2=1
      ix1=0
      do i=1,isize,imul
         if(ix1.eq.inx1a) then
            if(ix2.eq.inx2a) then
               ix3=ix3+1
               ix2=1
            else
               ix2=ix2+1
            endif
            ix1=1
         else
            ix1=ix1+1
         endif

         if(jccw1) then
            i1=ix1
         else
            i1=inx1a+1-ix1
         endif

         if(jccw2) then
            i2=ix2
         else
            i2=inx2a+1-ix2
         endif

         if(jccw3) then
            i3=ix3
         else
            i3=inx3a+1-ix3
         endif

         zdata(i1,i2,i3)=s%profs(iadr)%buf(i)

         if(present(zcoeff)) then
            do jj=1,imul-1
               zcoeff(jj,i1,i2,i3)=ksign(jj)*s%profs(iadr)%buf(i+jj)
            enddo
         endif
      enddo

    end subroutine xplasma_getProf_3data

    subroutine ksign_3d(jccw1,jccw2,jccw3,hspline,isize_coeff,ksign)

      ! PRIVATE
      ! 2d sign factor for derivatives, depending on CCW reversal of
      ! periodic coordinates-- affects Hermite dimensions only

      logical, intent(in) :: jccw1,jccw2,jccw3  ! CCW flags, dims 1 2 & 3
      integer, intent(in) :: hspline(3)   ! Hermite/Spline flags, dims 1 -- 3
      integer, intent(in) :: isize_coeff  ! #of coefficients per node point
      !  1 if Spline/Hermite in 1 dim only; 3 if in 2 dims; 7 if in 3 dims
      !  (could have computed this from hspline but callers already know)

      integer, intent(out) :: ksign(7)    ! output coefficient sign factors

      ! spline coefficients, involving 2nd derivatives and products of
      ! 2nd derivatives, always get a sign factor of +1

      !------------------------------
      ! local:
      integer :: jsign(3)
      !------------------------------

      jsign=1
      ksign=1
      if((.not.jccw1).and.(hspline(1).eq.1)) jsign(1)=-1
      if((.not.jccw2).and.(hspline(2).eq.1)) jsign(2)=-1
      if((.not.jccw3).and.(hspline(3).eq.1)) jsign(3)=-1

      if(isize_coeff.eq.1) then
         !  only 1 dim is Hermite/Spline
         if(hspline(1).gt.0) ksign(1)=jsign(1)
         if(hspline(2).gt.0) ksign(1)=jsign(2)
         if(hspline(3).gt.0) ksign(1)=jsign(3)

      else if(isize_coeff.eq.3) then
         !  2 dims Hermite or spline
         if(hspline(3).le.0) then
            ksign(1)=jsign(1)  ! dims 1 & 2
            ksign(2)=jsign(2)
         else if(hspline(2).le.0) then
            ksign(1)=jsign(1)  ! dims 1 & 3
            ksign(2)=jsign(3)
         else
            ksign(1)=jsign(2)  ! dims 2 & 3
            ksign(2)=jsign(3)
         endif
         ksign(3)=ksign(1)*ksign(2)

      else
         !  all 3 dims
         ksign(1:3)=jsign(1:3)
         ksign(4)=ksign(1)*ksign(2)
         ksign(5)=ksign(1)*ksign(3)
         ksign(6)=ksign(2)*ksign(3)
         ksign(7)=ksign(1)*ksign(2)*ksign(3)
      endif

    end subroutine ksign_3d

    !---------------------------------------------------------
    subroutine xplasma_diff_1dprof(s,id,zdata,diff_code,ier, &
         splineType, units, ccwflag1, coeffs)
      
      !  return information on differences between specified profile
      !  and the passed data

      type (xplasma), pointer :: s
      integer, intent(in) :: id                   ! profile id
      real*8, dimension(:), intent(in) :: zdata   ! data to be compared
      integer, intent(out) :: diff_code           ! difference code (see below)
      integer, intent(out) :: ier                 ! completion code 0=OK

      integer, intent(in), optional :: splineType ! interpolation type code
      character*(*), intent(in), optional :: units  ! units label
      logical, intent(in), optional :: ccwflag1   ! CCW flag (default T)

      real*8, dimension(:), intent(in), optional :: coeffs
      ! optional: compare coeffs (e.g. spline coeffs or Hermite derivatives).
      ! if this is present, it is compared to the stored data; differences
      ! are counted the same as differences in (zdata).
      !   This argument is ignored even if present, for non-spline non-Hermite
      !   data.

      !---------------------

      ! if an error occurs (e.g. id invalid), diff_code = -99 is returned.
      ! if NO differences are detected, diff_code = 0 is returned.

      ! if the data size is different, diff_code = -1 is returned
      ! if splineType is present and does not match the stored splineType,
      !    diff_code = -2 is returned.
      ! if the units label is present and does not match the stored label,
      !    diff_code = -3 is returned.

      ! if the coeffs array is present, but a size other than expected,
      !    diff_code = -4 is returned.

      ! if the above match, the passed data are compared point by point 
      ! with the data stored in the xplasma object; 

      ! tolerance is s%bdytol*max(abs(minval(zdata)),abs(maxval(zdata)))

      ! the value of diff_code returned is the number of points for which
      ! the absolute value of the difference exceeds the tolerance.

      !------------------------
      integer :: iadr,idg,iadx,inx,inxa,imul,isize,i,ix,ixc,ispline
      integer :: incr,iertmp,jsign
      character*128 msgbuf
      logical :: jccw,perio
      character*32 loc_units
      real*8 :: ztol
      !------------------------

      diff_code = -99

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      if(s%profs(iadr)%rank.ne.1) then
         ier=512
         write(msgbuf,*) '  profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  profile rank: ',s%profs(iadr)%rank
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  passed array rank:  1 (xplasma_diff_1dprof)'
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      ispline = s%profs(iadr)%kspline

      idg = s%profs(iadr)%gridIds(1)  ! get grid id
      iadx = s%dict(idg)%dindex  ! get grid address
      inx = s%grids(iadx)%size    ! get grid size

      if(ispline.lt.0) then
         inxa=inx-1
      else
         inxa=inx
      endif

      if(size(zdata).ne.inxa) then
         diff_code = -1
         return
      endif

      if(ispline.gt.0) then
         if(present(coeffs)) then
            if(size(coeffs).ne.inx) then
               diff_code = -4
               return
            endif
         endif
      endif

      if(present(splineType)) then
         if(ispline.ne.splineType) then
            diff_code = -2
            return
         endif
      endif

      if(present(units)) then
         call xplasma_prof_info(s,id,iertmp, units=loc_units)
         if(units.ne.loc_units) then
            diff_code = -3
            return
         endif
      endif

      diff_code = 0   ! no differences so far...

      ztol = max(abs(minval(zdata)),abs(maxval(zdata)))*s%bdytol

      jccw=.TRUE.
      if(present(ccwflag1)) then
         if(.not.ccwflag1) then
            call xplasma_grid_info(s,idg,iertmp, perio=perio)
            if(perio) then
               jccw=ccwflag1
            endif
         endif
      endif
            
      isize=s%profs(iadr)%size
      imul=isize/inxa

      ! imul=1 for step & pcwise linear; 2 for Hermite and Spline

      if(jccw) then
         ix=0
         incr=1
         jsign=1
      else
         ix=inxa+1
         incr=-1
         jsign=-1
      endif
      ixc=ix

      do i=1,isize,imul
         ix=ix+incr
         if(abs(zdata(ix)-s%profs(iadr)%buf(i)).gt.ztol) then
            diff_code = diff_code + 1
         endif
      enddo

      if(ispline.gt.0) then
         if(present(coeffs)) then

            ztol = max(abs(minval(coeffs)),abs(maxval(coeffs)))*s%bdytol

            ix=ixc
            do i=2,isize,imul
               ix=ix+incr
               if(abs(jsign*coeffs(ix)-s%profs(iadr)%buf(i)).gt.ztol) then
                  diff_code = diff_code + 1
               endif
            enddo

         endif
      endif

    end subroutine xplasma_diff_1dprof

    !---------------------------------------------------------
    subroutine xplasma_diff_2dprof(s,id,zdata,diff_code,ier, &
         splineType, units, &
         ccwflag1,ccwflag2, &
         coeffs)
      
      !  return information on differences between specified profile
      !  and the passed data

      type (xplasma), pointer :: s
      integer, intent(in) :: id                     ! profile id
      real*8, dimension(:,:), intent(in) :: zdata   ! data returned
      integer, intent(out) :: diff_code           ! difference code (see below)
      integer, intent(out) :: ier                   ! completion code 0=OK

      integer, intent(in), optional :: splineType ! interpolation type code
      character*(*), intent(in), optional :: units  ! units label
      logical, intent(in), optional :: ccwflag1   ! dim. 1 CCW flag (default T)
      logical, intent(in), optional :: ccwflag2   ! dim. 2 CCW flag (default T)

      real*8, intent(in), dimension(:,:,:), optional :: coeffs
      ! optional: compare coeffs (e.g. spline coeffs or Hermite derivatives).
      ! if this is present, it is compared to the stored data; differences
      ! are counted the same as differences in (zdata).
      !   This argument is ignored even if present, for non-spline non-Hermite
      !   data.

      !--------------------

      ! if an error occurs (e.g. id invalid), diff_code = -99 is returned.
      ! if NO differences are detected, diff_code = 0 is returned.

      ! if the data size is different, diff_code = -1 is returned
      ! if splineType is present and does not match the stored splineType,
      !    diff_code = -2 is returned.
      ! if the units label is present and does not match the stored label,
      !    diff_code = -3 is returned.

      ! if the coeffs array is present, but a size other than expected,
      !    diff_code = -4 is returned.

      ! if the above match, the passed data are compared point by point 
      ! with the data stored in the xplasma object; 

      ! tolerance is s%bdytol*max(abs(minval(zdata)),abs(maxval(zdata)))

      ! the value of diff_code returned is the number of points for which
      ! the absolute value of the difference exceeds the tolerance.

      !------------------------
      integer :: iadr,iadx1,inx1,iadx2,inx2,imul,isize,i,ix1,ix2
      integer :: i1,i2,iertmp
      integer :: idg1,idg2,incr1,incr2,jsign1,jsign2
      integer :: ispline,hspline(2),inx1a,inx2a
      integer :: isize_tot,isize_coeffs,isize_data
      character*128 msgbuf
      logical :: perio,jccw1,jccw2
      character*32 loc_units
      real*8 :: ztol,ztolc,zdiff
      !------------------------

      diff_code = -99

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      if(s%profs(iadr)%rank.ne.2) then
         ier=512
         write(msgbuf,*) '  profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  profile rank: ',s%profs(iadr)%rank
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  passed array rank:  2 (xplasma_diff_2dprof)'
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      idg1 = s%profs(iadr)%gridIds(1)  ! get grid id
      iadx1 = s%dict(idg1)%dindex  ! get grid address
      inx1 = s%grids(iadx1)%size    ! get grid size

      idg2 = s%profs(iadr)%gridIds(2)  ! get grid id
      iadx2 = s%dict(idg2)%dindex  ! get grid address
      inx2 = s%grids(iadx2)%size    ! get grid size

      ispline = s%profs(iadr)%kspline
      call spline_type_decode(ispline,hspline,2)

      isize_tot = 1
      isize_data = 1
      imul=1

      if(hspline(1).lt.0) then
         inx1a=inx1-1
         isize_tot=isize_tot*inx1a
         isize_data=isize_data*inx1a
      else if(hspline(1).eq.0) then
         inx1a=inx1
         isize_tot=isize_tot*inx1a
         isize_data=isize_data*inx1a
      else
         inx1a=inx1
         imul=imul*2
         isize_tot=isize_tot*inx1a*2
         isize_data=isize_data*inx1a
      endif

      if(hspline(2).lt.0) then
         inx2a=inx2-1
         isize_tot=isize_tot*inx2a
         isize_data=isize_data*inx2a
      else if(hspline(2).eq.0) then
         inx2a=inx2
         isize_tot=isize_tot*inx2a
         isize_data=isize_data*inx2a
      else
         inx2a=inx2
         imul=imul*2
         isize_tot=isize_tot*inx2a*2
         isize_data=isize_data*inx2a
      endif

      isize_coeffs = isize_tot - isize_data

      if((size(zdata,1).ne.inx1a).or.(size(zdata,2).ne.inx2a)) then
         diff_code = -1
         return
      endif

      if(present(splineType)) then
         if(ispline.ne.splineType) then
            diff_code = -2
            return
         endif
      endif

      if(present(units)) then
         call xplasma_prof_info(s,id,iertmp, units=loc_units)
         if(units.ne.loc_units) then
            diff_code = -3
            return
         endif
      endif

      if(imul.gt.1) then
         if(present(coeffs)) then
            if(size(coeffs,1).ne.(imul-1)) then
               diff_code = -4
               return
            endif
            if((size(coeffs,2).ne.inx1a).or.(size(coeffs,3).ne.inx2a)) then
               diff_code = -4
               return
            endif
         endif
      endif

      diff_code = 0   ! no differences so far...

      ztol = max(abs(minval(zdata)),abs(maxval(zdata)))*s%bdytol

      if(imul.gt.1) then
         if(present(coeffs)) then
            ztolc = max(abs(minval(coeffs)),abs(maxval(coeffs)))*s%bdytol
         endif
      endif

      jccw1=.TRUE.
      if(present(ccwflag1)) then
         if(.not.ccwflag1) then
            call xplasma_grid_info(s,idg1,iertmp, perio=perio)
            if(perio) then
               jccw1=ccwflag1
            endif
         endif
      endif

      jccw2=.TRUE.
      if(present(ccwflag2)) then
         if(.not.ccwflag2) then
            call xplasma_grid_info(s,idg2,iertmp, perio=perio)
            if(perio) then
               jccw2=ccwflag2
            endif
         endif
      endif

      jsign1=1
      jsign2=1
      if((.not.jccw1).and.(hspline(1).eq.1)) then
         jsign1=-1
      endif
      if((.not.jccw2).and.(hspline(2).eq.1)) then
         jsign2=-1
      endif

      isize=s%profs(iadr)%size
      imul=isize/(inx1a*inx2a)
      
      ix2=1
      ix1=0
      do i=1,isize,imul
         if(ix1.eq.inx1a) then
            ix2=ix2+1
            ix1=1
         else
            ix1=ix1+1
         endif

         if(jccw1) then
            i1=ix1
         else
            i1=inx1a+1-ix1
         endif

         if(jccw2) then
            i2=ix2
         else
            i2=inx2a+1-ix2
         endif

         if(abs(zdata(i1,i2)-s%profs(iadr)%buf(i)).gt.ztol) then
            diff_code = diff_code + 1
         endif

         if(imul.eq.2) then
            if(present(coeffs)) then
               if(hspline(1).gt.0) then
                  zdiff = coeffs(1,i1,i2)*jsign1 - s%profs(iadr)%buf(i+1)
                  if(abs(zdiff).gt.ztolc) then
                     diff_code=diff_code + 1
                  endif
               endif

               if(hspline(2).gt.0) then
                  zdiff = coeffs(1,i1,i2)*jsign2 - s%profs(iadr)%buf(i+1)
                  if(abs(zdiff).gt.ztolc) then
                     diff_code=diff_code + 1
                  endif
               endif
            endif

         else if(imul.gt.2) then
            if(present(coeffs)) then
               zdiff = coeffs(1,i1,i2)*jsign1 - s%profs(iadr)%buf(i+1)
               if(abs(zdiff).gt.ztolc) then
                  diff_code=diff_code + 1
               endif

               zdiff = coeffs(2,i1,i2)*jsign2 - s%profs(iadr)%buf(i+2)
               if(abs(zdiff).gt.ztolc) then
                  diff_code=diff_code + 1
               endif

               zdiff = coeffs(3,i1,i2)*jsign1*jsign2 - s%profs(iadr)%buf(i+3)
               if(abs(zdiff).gt.ztolc) then
                  diff_code=diff_code + 1
               endif
            endif
         endif

      enddo

    end subroutine xplasma_diff_2dprof

    !---------------------------------------------------------
    subroutine xplasma_diff_3dprof(s,id,zdata,diff_code,ier, &
         splineType, units, &
         ccwflag1,ccwflag2,ccwflag3)
      
      !  return information on differences between specified profile
      !  and the passed data

      type (xplasma), pointer :: s
      integer, intent(in) :: id                       ! profile id
      real*8, dimension(:,:,:), intent(in) :: zdata  ! data returned
      integer, intent(out) :: diff_code           ! difference code (see below)
      integer, intent(out) :: ier                     ! completion code 0=OK

      integer, intent(in), optional :: splineType ! interpolation type code
      character*(*), intent(in), optional :: units  ! units label
      logical, intent(in), optional :: ccwflag1   ! dim. 1 CCW flag (default T)
      logical, intent(in), optional :: ccwflag2   ! dim. 2 CCW flag (default T)
      logical, intent(in), optional :: ccwflag3   ! dim. 3 CCW flag (default T)

      ! if an error occurs (e.g. id invalid), diff_code = -99 is returned.
      ! if NO differences are detected, diff_code = 0 is returned.

      ! if the data size is different, diff_code = -1 is returned
      ! if splineType is present and does not match the stored splineType,
      !    diff_code = -2 is returned.
      ! if the units label is present and does not match the stored label,
      !    diff_code = -3 is returned.

      ! if the above match, the passed data are compared point by point 
      ! with the data stored in the xplasma object; 

      ! tolerance is s%bdytol*max(abs(minval(zdata)),abs(maxval(zdata)))

      ! the value of diff_code returned is the number of points for which
      ! the absolute value of the difference exceeds the tolerance.

      !------------------------
      integer :: iadr,iadx1,inx1,iadx2,inx2,iadx3,inx3,imul,isize,i,ix1,ix2,ix3
      integer :: idg1,idg2,idg3,i1,i2,i3,iertmp
      logical :: jccw1,jccw2,jccw3,perio
      integer :: ispline,hspline(3),inx1a,inx2a,inx3a
      character*128 msgbuf
      character*32 loc_units
      real*8 :: ztol
      !------------------------

      diff_code = -99

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      if(s%profs(iadr)%rank.ne.3) then
         ier=512
         write(msgbuf,*) '  profile name: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  profile rank: ',s%profs(iadr)%rank
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  passed array rank:  3 (xplasma_diff_3dprof)'
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      idg1 = s%profs(iadr)%gridIds(1)  ! get grid id
      iadx1 = s%dict(idg1)%dindex  ! get grid address
      inx1 = s%grids(iadx1)%size    ! get grid size

      idg2 = s%profs(iadr)%gridIds(2)  ! get grid id
      iadx2 = s%dict(idg2)%dindex  ! get grid address
      inx2 = s%grids(iadx2)%size    ! get grid size

      idg3 = s%profs(iadr)%gridIds(3)  ! get grid id
      iadx3 = s%dict(idg3)%dindex  ! get grid address
      inx3 = s%grids(iadx3)%size    ! get grid size

      ispline = s%profs(iadr)%kspline
      call spline_type_decode(ispline,hspline,3)

      if(hspline(1).lt.0) then
         inx1a=inx1-1
      else
         inx1a=inx1
      endif

      if(hspline(2).lt.0) then
         inx2a=inx2-1
      else
         inx2a=inx2
      endif

      if(hspline(3).lt.0) then
         inx3a=inx3-1
      else
         inx3a=inx3
      endif

      if((size(zdata,1).ne.inx1a).or.(size(zdata,2).ne.inx2a).or. &
           (size(zdata,3).ne.inx3a)) then
         diff_code = -1
         return
      endif

      if(present(splineType)) then
         if(ispline.ne.splineType) then
            diff_code = -2
            return
         endif
      endif

      if(present(units)) then
         call xplasma_prof_info(s,id,iertmp, units=loc_units)
         if(units.ne.loc_units) then
            diff_code = -3
            return
         endif
      endif

      diff_code = 0   ! no differences so far...

      ztol = max(abs(minval(zdata)),abs(maxval(zdata)))*s%bdytol

      jccw1=.TRUE.
      if(present(ccwflag1)) then
         if(.not.ccwflag1) then
            call xplasma_grid_info(s,idg1,iertmp, perio=perio)
            if(perio) then
               jccw1=ccwflag1
            endif
         endif
      endif

      jccw2=.TRUE.
      if(present(ccwflag2)) then
         if(.not.ccwflag2) then
            call xplasma_grid_info(s,idg2,iertmp, perio=perio)
            if(perio) then
               jccw2=ccwflag2
            endif
         endif
      endif

      jccw3=.TRUE.
      if(present(ccwflag3)) then
         if(.not.ccwflag3) then
            call xplasma_grid_info(s,idg3,iertmp, perio=perio)
            if(perio) then
               jccw3=ccwflag3
            endif
         endif
      endif

      isize=s%profs(iadr)%size
      imul=isize/(inx1a*inx2a*inx3a)
      
      ix3=1
      ix2=1
      ix1=0
      do i=1,isize,imul
         if(ix1.eq.inx1a) then
            if(ix2.eq.inx2a) then
               ix3=ix3+1
               ix2=1
            else
               ix2=ix2+1
            endif
            ix1=1
         else
            ix1=ix1+1
         endif

         if(jccw1) then
            i1=ix1
         else
            i1=inx1a+1-ix1
         endif

         if(jccw2) then
            i2=ix2
         else
            i2=inx2a+1-ix2
         endif

         if(jccw3) then
            i3=ix3
         else
            i3=inx3a+1-ix3
         endif

         if(abs(zdata(i1,i2,i3)-s%profs(iadr)%buf(i)).gt.ztol) then
            diff_code = diff_code + 1
         endif

      enddo

    end subroutine xplasma_diff_3dprof

    !----------------------------------------------
    subroutine xplasma_xpeval_init(xinfo)

      !  init an xpeval type -- nullify pointers, do NOT deallocate

      type (xpeval) :: xinfo

      xinfo%idx = 0
      xinfo%inum = 0
      xinfo%inx = 0
      xinfo%icoord = 0
      xinfo%ccwflag = .TRUE.

      nullify(xinfo%iderivs)
      nullify(xinfo%ix)
      nullify(xinfo%dxn)
      nullify(xinfo%hx)
      nullify(xinfo%hxi)

    end subroutine xplasma_xpeval_init

    !----------------------------------------------
    subroutine xpeval_free(xinfo)

      type(xpeval) :: xinfo

      !  deallocate & clean up xinfo

      xinfo%idx = 0
      xinfo%inum = 0
      xinfo%inx = 0
      xinfo%icoord = 0
      xinfo%ccwflag = .TRUE.

      if(associated(xinfo%iderivs)) deallocate(xinfo%iderivs)
      if(associated(xinfo%ix)) deallocate(xinfo%ix)
      if(associated(xinfo%dxn)) deallocate(xinfo%dxn)
      if(associated(xinfo%hx)) deallocate(xinfo%hx)
      if(associated(xinfo%hxi)) deallocate(xinfo%hxi)

    end subroutine xpeval_free

    !----------------------------------------------
    subroutine xplasma_x_lookup1(s,idx,x,xinfo,ier, &
         ccwflag,force_bounds,n_out_of_bounds)

      !  do x lookup -- scalar x value

      type (xplasma), pointer :: s
      
      integer, intent(in) :: idx                  ! x grid (grid) id
      real*8, intent(in) :: x                     ! single scalar x value
      type (xpeval) :: xinfo                      ! *** lookup results
      integer, intent(out) :: ier                 ! status code, 0=OK

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds
      
      !-------------------------------------------------
      integer :: iadx,iadc,inum,icoord,imode,iwarn
      real*8 :: xtmp(1)
      logical :: iforce,perio
      !-------------------------------------------------

      call xplasma_errck0_grid(s,idx,ier)
      if(ier.ne.0) return

      call xpeval_free(xinfo)

      xinfo%idx = idx
      xinfo%inum = 1
      inum = 1

      iadx = s%dict(idx)%dindex
      xinfo%inx = s%grids(iadx)%size

      icoord = s%grids(iadx)%coord
      xinfo%icoord = icoord
      iadc = s%dict(icoord)%dindex
      perio = s%coords(iadc)%periodic

      xinfo%ccwflag = .TRUE.
      if(present(ccwflag).and.perio) xinfo%ccwflag = ccwflag  ! defaults to T

      allocate(xinfo%iderivs(1)); xinfo%iderivs(1)=0

      !  verify that all xs elements are in range; for periodic
      !  grids normalize the range and perform CW->CCW reversal if necessary

      iforce=.FALSE.
      if(present(force_bounds)) iforce=force_bounds

      xtmp = x

      allocate(xinfo%ix(inum),xinfo%dxn(inum), &
           xinfo%hx(inum),xinfo%hxi(inum))
      imode=2
      call xplookup(s,inum,xtmp,xinfo%ccwflag, &
           s%grids(iadx)%size,s%grids(iadx)%xpkg,imode, &
           xinfo%ix,xinfo%dxn,xinfo%hx,xinfo%hxi,iwarn)

      if(present(n_out_of_bounds)) n_out_of_bounds = iwarn
      if(.not.iforce) then
         if(iwarn.gt.0) then
            ier=514
            call xpeval_free(xinfo)
         endif
      endif

    end subroutine xplasma_x_lookup1

    !----------------------------------------------
    subroutine xplasma_x_lookupn(s,idx,x,xinfo,ier, &
         ccwflag,force_bounds,n_out_of_bounds)

      !  do x lookup -- scalar x value

      type (xplasma), pointer :: s
      
      integer, intent(in) :: idx                  ! x grid (grid) id
      real*8, intent(in), dimension(:) :: x       ! vector of x values
      type (xpeval) :: xinfo                      ! *** lookup results
      integer, intent(out) :: ier                 ! status code, 0=OK

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds
      
      !-------------------------------------------------
      integer :: iadx,inum,icoord,iadc,imode,iwarn
      logical :: perio,iforce
      !-------------------------------------------------

      call xplasma_errck0_grid(s,idx,ier)
      if(ier.ne.0) return

      call xpeval_free(xinfo)

      xinfo%idx = idx
      xinfo%inum = size(x)
      inum = size(x)

      iadx = s%dict(idx)%dindex
      xinfo%inx = s%grids(iadx)%size

      icoord = s%grids(iadx)%coord
      xinfo%icoord = icoord
      iadc = s%dict(icoord)%dindex
      perio = s%coords(iadc)%periodic

      xinfo%ccwflag = .TRUE.
      if(present(ccwflag).and.perio) xinfo%ccwflag = ccwflag  ! defaults to T

      allocate(xinfo%iderivs(1)); xinfo%iderivs(1)=0

      !  verify that all xs elements are in range; for periodic
      !  grids normalize the range and perform CW->CCW reversal if necessary

      iforce=.FALSE.
      if(present(force_bounds)) iforce=force_bounds

      allocate(xinfo%ix(inum),xinfo%dxn(inum), &
           xinfo%hx(inum),xinfo%hxi(inum))
      imode=2
      call xplookup(s,inum,x,xinfo%ccwflag, &
           s%grids(iadx)%size,s%grids(iadx)%xpkg,imode, &
           xinfo%ix,xinfo%dxn,xinfo%hx,xinfo%hxi,iwarn)

      if(present(n_out_of_bounds)) n_out_of_bounds = iwarn
      if(.not.iforce) then
         if(iwarn.gt.0) then
            ier=514
            call xpeval_free(xinfo)
         endif
      endif

    end subroutine xplasma_x_lookupn

    !---------------------------------------------------------
    subroutine xplasma_xpeval_1dprof(s,id,xinfo,ans,ier, ideriv)

      !  evaluate a single 1d profile function f(x) at a vector of
      !  target points for which lookup has already been executed
      !  and stored in "xinfo".

      !  the x axis grid id of the profile and of the xinfo must match
      !  the size of the output and the no. of points in the xinfo must
      !  match.  

      !  xinfo was set up by a prior xplasma_x_lookup call.

      type(xplasma), pointer :: s
      integer, intent(in) :: id      ! 1d function to evaluate
      type (xpeval) :: xinfo         ! prepared interpolation information

      real*8, dimension(:), intent(out) :: ans   ! interpolation result
      integer, intent(out) :: ier    ! status code, 0=OK

      integer, intent(in), optional :: ideriv  ! derivative selection
      
      !-------------------------------
      integer :: idx,ids(2),ii,idxx,iadr,irank,iax1,iax2,iadx,icoord
      character*128 zmsg
      type(xpeval) :: xx(1)
      integer :: kk(maxRank+1),isign,jderiv,ispline,iertmp
      logical :: expand1,ccwflag1
      !-------------------------------

      !  check if id points to a profile object

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      jderiv=0
      if(present(ideriv)) jderiv=ideriv

      !  check if x grid matches

      iadr = s%dict(id)%dindex
      ids(1)=id
      ids(2)=s%profs(iadr)%profIds(1)

      do ii=2,1,-1
         if(ids(ii).eq.0) cycle

         iadr=s%dict(ids(ii))%dindex
         idx = s%profs(iadr)%gridIds(1)
         iadx = s%dict(idx)%dindex
         icoord = s%grids(iadx)%coord
         irank = s%profs(iadr)%rank
         if((icoord.eq.xinfo%icoord).and.(irank.eq.1)) exit
         if(ii.eq.1) exit
      enddo

      if(ii.eq.1) then
         if(irank.ne.1) then

            ier=519
            write(zmsg,*) ' ...profile function id is: ',trim(s%dict(id)%name)
            call xplasma_errmsg_append(s,zmsg)
            write(zmsg,*) '    profile rank is ',irank, &
                 ' but lookup information for only 1 x grid was provided.'
            call xplasma_errmsg_append(s,zmsg)
            return

         endif
      endif

      if(idx.ne.xinfo%idx) then

         ! equivalence still possible -- if so, no error

         iax1 = s%dict(idx)%dindex
         iax2 = xinfo%idx
         iax2 = s%dict(iax2)%dindex

         if(s%grid_equiv(iax1).ne.s%grid_equiv(iax2)) then

            ier=517
            call xplasma_errmsg_append(s,' ...profile function '// &
                 trim(s%dict(id)%name)//' defined vs. '// &
                 trim(s%dict(idx)%name))
            idxx=xinfo%idx
            call xplasma_errmsg_append(s,'    but xinfo was provided for '// &
                 trim(s%dict(idxx)%name))

            return
         endif
      endif

      !---------------------------
      !  check argument size

      if(size(ans).ne.xinfo%inum) then

         ier=518
         write(zmsg,*) ' ...xinfo has ',xinfo%inum, &
              ' points but answer vector has ',size(ans),' points.'
         call xplasma_errmsg_append(s,zmsg)
         return
         
      endif

      !---------------------------
      !  check derivative 

      ispline=s%profs(iadr)%kspline
      call chkderiv(s,jderiv,ispline,ids(ii),ier)
      if(ier.ne.0) return

      !---------------------------
      !  fill arguments for xplasma_evala; evaluate profile

      isign=1
      ccwflag1=xinfo%ccwflag
      if(.not.ccwflag1) then
         if((jderiv.eq.1).or.(jderiv.eq.3)) isign=-isign
      endif

      call xx_fill(xinfo,xx(1),jderiv,size(ans),expand1)
      kk(1)=iadr
      kk(2)=1

      call xplasma_evala(s,ans,kk,xx)

      if(isign.eq.-1) then
         ans=ans*isign
      endif

      call xx_cleanup(xx(1),expand1)

    end subroutine xplasma_xpeval_1dprof

    !---------------------------------------------------------
    subroutine xplasma_xpeval_2dprof(s,id,xinfo1,xinfo2,ans,ier, &
         ideriv1, ideriv2)

      type(xplasma), pointer :: s
      integer, intent(in) :: id      ! 2d function to evaluate
      type (xpeval) :: xinfo1,xinfo2 ! prepared interpolation information sets

      real*8, dimension(:), intent(out) :: ans   ! interpolation result
      integer, intent(out) :: ier    ! status code, 0=OK

      integer, intent(in), optional :: ideriv1  ! derivative selection, x1
      integer, intent(in), optional :: ideriv2  ! derivative selection, x2

      !-------------------------------
      integer :: idx1,idx2,idxx1,idxx2,iadr,irank,isize,iertmp,ids(2),ii
      integer :: iax1,iax2,iax3,icoord1,icoord2
      character*128 zmsg
      type(xpeval) :: xx(2)
      integer :: kk(maxRank+1),isign,jderiv1,jderiv2,ispline,hspline(2)
      logical :: expand1,expand2,ccwflag1,ccwflag2
      !-------------------------------

      !------------------------------
      !  check if id points to a profile object

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      jderiv1=0
      if(present(ideriv1)) jderiv1=ideriv1

      jderiv2=0
      if(present(ideriv2)) jderiv2=ideriv2

      !------------------------------
      !  check if x grids match & if rank is OK

      idxx1=xinfo1%idx
      idxx2=xinfo2%idx

      iadr = s%dict(id)%dindex
      ids(1)=id
      ids(2)=s%profs(iadr)%profIds(1)

      do ii=2,1,-1
         if(ids(ii).eq.0) cycle

         iadr=s%dict(ids(ii))%dindex
         irank = s%profs(iadr)%rank

         idx1 = s%profs(iadr)%gridIds(1)
         iax1 = s%dict(idx1)%dindex
         icoord1 = s%grids(iax1)%coord

         icoord2 = 0
         if(irank.gt.1) then
            idx2 = s%profs(iadr)%gridIds(2)
            iax2 = s%dict(idx2)%dindex
            icoord2 = s%grids(iax2)%coord
         endif

         if(irank.eq.2) then
            if((icoord1.eq.xinfo1%icoord).or.(icoord1.eq.xinfo2%icoord)) then
               if((icoord2.eq.xinfo1%icoord).or.(icoord2.eq.xinfo2%icoord)) &
                    exit
            endif
         endif
         if(ii.eq.1) exit
      enddo

      if(ii.eq.1) then
         if(irank.ne.2) then

            ier=519
            write(zmsg,*) ' ...profile function id is: ',trim(s%dict(id)%name)
            call xplasma_errmsg_append(s,zmsg)
            write(zmsg,*) '    profile rank is ',irank, &
                 ' but lookup information for 2 x axes were provided.'
            call xplasma_errmsg_append(s,zmsg)
            return

         endif
      endif

      if(  ((idx1.ne.xinfo1%idx).and.(idx1.ne.xinfo2%idx)).or. &
           ((idx2.ne.xinfo1%idx).and.(idx2.ne.xinfo2%idx)) ) then

         ! equivalence still possible -- if so, no error

         ier = 0

         iax1 = s%dict(idx1)%dindex

         iax2 = xinfo1%idx
         iax2 = s%dict(iax2)%dindex
         iax3 = xinfo2%idx
         iax3 = s%dict(iax3)%dindex

         if((s%grid_equiv(iax1).ne.s%grid_equiv(iax2)).and. &
              (s%grid_equiv(iax1).ne.s%grid_equiv(iax3))) ier = ier+1

         iax1 = s%dict(idx2)%dindex

         if((s%grid_equiv(iax1).ne.s%grid_equiv(iax2)).and. &
              (s%grid_equiv(iax1).ne.s%grid_equiv(iax3))) ier = ier+1

         if(ier.gt.0) then
            ier=517
            call xplasma_errmsg_append(s,' ...profile function '// &
                 trim(s%dict(id)%name)//' defined vs.: ')
            call xplasma_errmsg_append(s,'    '//trim(s%dict(idx1)%name)// &
                 ' and '//trim(s%dict(idx2)%name))
            idxx1=xinfo1%idx
            idxx2=xinfo2%idx
            call xplasma_errmsg_append(s,'    but xinfo was provided for '// &
                 trim(s%dict(idxx1)%name)//' and '//trim(s%dict(idxx2)%name))
            return
         endif

      endif

      isize=size(ans)
      if(  ((xinfo1%inum.ne.1).and.(xinfo1%inum.ne.isize)).or. &
           ((xinfo2%inum.ne.1).and.(xinfo2%inum.ne.isize)) ) then

         ier=518
         write(zmsg,*) ' ...answer vector has ',isize,' points.'
         call xplasma_errmsg_append(s,zmsg)
         write(zmsg,*) '    but xinfo packets have sizes: ', &
              xinfo1%inum,xinfo2%inum
         return

      endif

      !---------------------------
      !  check derivative 

      ispline=s%profs(iadr)%kspline
      call spline_type_decode(ispline,hspline,2)

      call chkderiv(s,jderiv1,hspline(1),ids(ii),ier)
      iertmp = ier

      call chkderiv(s,jderiv2,hspline(2),ids(ii),ier)
      ier=max(ier,iertmp)

      if(ier.ne.0) return

      !---------------------------
      !  fill arguments for xplasma_evala; evaluate profile

      isign=1
      ccwflag1=xinfo1%ccwflag
      ccwflag2=xinfo2%ccwflag
      if(.not.ccwflag1) then
         if((jderiv1.eq.1).or.(jderiv1.eq.3)) isign=-isign
      endif
      if(.not.ccwflag2) then
         if((jderiv2.eq.1).or.(jderiv2.eq.3)) isign=-isign
      endif

      if(icoord1.eq.xinfo1%icoord) then
         call xx_fill(xinfo1,xx(1),jderiv1,size(ans),expand1)
         if(xinfo2%inum.eq.1) then
            call xx_fill(xinfo2,xx(2),jderiv2,1,expand2)
         else
            call xx_fill(xinfo2,xx(2),jderiv2,size(ans),expand2)
         endif
      else
         call xx_fill(xinfo2,xx(1),jderiv2,size(ans),expand1)
         if(xinfo1%inum.eq.1) then
            call xx_fill(xinfo1,xx(2),jderiv1,1,expand2)
         else
            call xx_fill(xinfo1,xx(2),jderiv1,size(ans),expand2)
         endif
      endif

      kk(1)=iadr
      kk(2)=1
      kk(3)=2

      call xplasma_evala(s,ans,kk,xx)

      if(isign.eq.-1) then
         ans=ans*isign
      endif

      call xx_cleanup(xx(1),expand1)
      call xx_cleanup(xx(2),expand2)

    end subroutine xplasma_xpeval_2dprof

    !---------------------------------------------------------
    subroutine xplasma_xpeval_3dprof(s,id,xinfo1,xinfo2,xinfo3,ans,ier, &
         ideriv1, ideriv2, ideriv3)

      type(xplasma), pointer :: s
      integer, intent(in) :: id      ! 3d function to evaluate
      type (xpeval) :: xinfo1,xinfo2,xinfo3 
                                     ! prepared interpolation information sets

      real*8, dimension(:), intent(out) :: ans   ! interpolation result
      integer, intent(out) :: ier    ! status code, 0=OK

      integer, intent(in), optional :: ideriv1  ! derivative selection, x1
      integer, intent(in), optional :: ideriv2  ! derivative selection, x2
      integer, intent(in), optional :: ideriv3  ! derivative selection, x3
      
      !-------------------------------
      integer :: idx1,idx2,idx3,idxx1,idxx2,idxx3,iadr,irank,isize,iertmp
      integer :: ids(2),ii,iax1,iax2,iax3,iax4,icoord1,icoord2,icoord3
      character*128 zmsg
      type(xpeval) :: xx(3)
      integer :: kk(maxRank+1),isign,jderiv1,jderiv2,jderiv3,ispline,hspline(3)
      logical :: expand(3),ccwflag1,ccwflag2,ccwflag3
      integer :: ix1,ix2,ix3,isiz1,isiz2
      !-------------------------------

      !  check if id points to a profile object

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      jderiv1=0
      if(present(ideriv1)) jderiv1=ideriv1

      jderiv2=0
      if(present(ideriv2)) jderiv2=ideriv2

      jderiv3=0
      if(present(ideriv3)) jderiv3=ideriv3

      !------------------------------
      !  check if x grids match & if rank is OK

      iadr = s%dict(id)%dindex
      ids(1)=id
      ids(2)=s%profs(iadr)%profIds(1)

      idxx1=xinfo1%idx
      idxx2=xinfo2%idx
      idxx3=xinfo3%idx

      do ii=2,1,-1
         if(ids(ii).eq.0) cycle

         iadr=s%dict(ids(ii))%dindex
         irank = s%profs(iadr)%rank

         idx1 = s%profs(iadr)%gridIds(1)
         iax1 = s%dict(idx1)%dindex
         icoord1 = s%grids(iax1)%coord

         if(irank.gt.1) then
            idx2 = s%profs(iadr)%gridIds(2)
            iax2 = s%dict(idx2)%dindex
            icoord2 = s%grids(iax2)%coord
         endif

         if(irank.gt.2) then
            idx3 = s%profs(iadr)%gridIds(3)
            iax3 = s%dict(idx3)%dindex
            icoord3 = s%grids(iax3)%coord
         endif

         if(irank.eq.3) then
            if((icoord1.eq.xinfo1%icoord).or. &
                 (icoord1.eq.xinfo2%icoord).or. &
                 (icoord1.eq.xinfo3%icoord)) then
               if((icoord2.eq.xinfo1%icoord).or. &
                    (icoord2.eq.xinfo2%icoord).or. &
                    (icoord2.eq.xinfo3%icoord)) then
                  if((icoord3.eq.xinfo1%icoord).or. &
                       (icoord3.eq.xinfo2%icoord).or. &
                       (icoord3.eq.xinfo3%icoord)) then
                     exit
                  endif
               endif
            endif
         endif
         if(ii.eq.1) exit
      enddo

      if(ii.eq.1) then
         if(irank.ne.3) then

            ier=519
            write(zmsg,*) ' ...profile function id is: ',trim(s%dict(id)%name)
            call xplasma_errmsg_append(s,zmsg)
            write(zmsg,*) '    profile rank is ',irank, &
                 ' but lookup information for 3 x axes were provided.'
            call xplasma_errmsg_append(s,zmsg)
            return

         endif
      endif

      if(  ((idx1.ne.idxx1).and.(idx1.ne.idxx2).and.(idx1.ne.idxx3)).or. &
           ((idx2.ne.idxx1).and.(idx2.ne.idxx2).and.(idx2.ne.idxx3)).or. &
           ((idx3.ne.idxx1).and.(idx3.ne.idxx2).and.(idx3.ne.idxx3)) ) then

         ! equivalence still possible -- if so, no error

         ier = 0

         iax1 = s%dict(idx1)%dindex

         iax2 = xinfo1%idx
         iax2 = s%dict(iax2)%dindex
         iax3 = xinfo2%idx
         iax3 = s%dict(iax3)%dindex
         iax4 = xinfo3%idx
         iax4 = s%dict(iax4)%dindex

         if((s%grid_equiv(iax1).ne.s%grid_equiv(iax2)).and. &
              (s%grid_equiv(iax1).ne.s%grid_equiv(iax3)).and. &
              (s%grid_equiv(iax1).ne.s%grid_equiv(iax4))) ier = ier + 1

         iax1 = s%dict(idx2)%dindex

         if((s%grid_equiv(iax1).ne.s%grid_equiv(iax2)).and. &
              (s%grid_equiv(iax1).ne.s%grid_equiv(iax3)).and. &
              (s%grid_equiv(iax1).ne.s%grid_equiv(iax4))) ier = ier + 1

         iax1 = s%dict(idx3)%dindex

         if((s%grid_equiv(iax1).ne.s%grid_equiv(iax2)).and. &
              (s%grid_equiv(iax1).ne.s%grid_equiv(iax3)).and. &
              (s%grid_equiv(iax1).ne.s%grid_equiv(iax4))) ier = ier + 1

         if(ier.gt.0) then
            
            ier=517
            call xplasma_errmsg_append(s,' ...profile function '// &
                 trim(s%dict(id)%name)//' defined vs.: ')
            call xplasma_errmsg_append(s,'    '//trim(s%dict(idx1)%name)// &
                 ' and '//trim(s%dict(idx2)%name)//' and '// &
                 ' and '//trim(s%dict(idx3)%name))
            call xplasma_errmsg_append(s,'    but xinfo was provided for '// &
                 trim(s%dict(idxx1)%name)//' and '//trim(s%dict(idxx2)%name)//&
                 ' and '//trim(s%dict(idxx3)%name))
            return
         endif
      endif

      !---------------------------
      !  fill arguments for xplasma_evala; evaluate profile

      isize=size(ans)
      if(  ((xinfo1%inum.ne.1).and.(xinfo1%inum.ne.isize)).or. &
           ((xinfo2%inum.ne.1).and.(xinfo2%inum.ne.isize)).or. &
           ((xinfo3%inum.ne.1).and.(xinfo3%inum.ne.isize)) ) then

         ier=518
         write(zmsg,*) ' ...answer vector has ',isize,' points.'
         call xplasma_errmsg_append(s,zmsg)
         write(zmsg,*) '    but xinfo packets have sizes: ', &
              xinfo1%inum,xinfo2%inum,xinfo3%inum
         return
         
      endif

      !---------------------------
      !  check derivative 

      ispline=s%profs(iadr)%kspline
      call spline_type_decode(ispline,hspline,3)

      call chkderiv(s,jderiv1,hspline(1),id,ier)
      iertmp = ier

      call chkderiv(s,jderiv2,hspline(2),id,ier)
      iertmp=max(ier,iertmp)

      call chkderiv(s,jderiv3,hspline(3),id,ier)
      ier=max(ier,iertmp)

      if(ier.ne.0) return

      !---------------------------
      !  fill arguments for xplasma_evala; evaluate profile

      isiz2=1

      if(xinfo1%icoord.eq.icoord1) ix1=1
      if(xinfo1%icoord.eq.icoord2) then
         ix1=2
         isiz2=max(isiz2,xinfo1%inum)
      endif
      if(xinfo1%icoord.eq.icoord3) then
         ix1=3
         isiz2=max(isiz2,xinfo1%inum)
      endif

      if(xinfo2%icoord.eq.icoord1) ix2=1
      if(xinfo2%icoord.eq.icoord2) then
         ix2=2
         isiz2=max(isiz2,xinfo2%inum)
      endif
      if(xinfo2%icoord.eq.icoord3) then
         ix2=3
         isiz2=max(isiz2,xinfo2%inum)
      endif

      if(xinfo3%icoord.eq.icoord1) ix3=1
      if(xinfo3%icoord.eq.icoord2) then
         ix3=2
         isiz2=max(isiz2,xinfo3%inum)
      endif
      if(xinfo3%icoord.eq.icoord3) then
         ix3=3
         isiz2=max(isiz2,xinfo3%inum)
      endif
      
      if(isiz2.gt.1) isiz2=size(ans)

      if(ix1.eq.1) isiz1=size(ans)
      if(ix1.gt.1) isiz1=isiz2
      call xx_fill(xinfo1,xx(ix1),jderiv1,isiz1,expand(ix1))

      if(ix2.eq.1) isiz1=size(ans)
      if(ix2.gt.1) isiz1=isiz2
      call xx_fill(xinfo2,xx(ix2),jderiv2,isiz1,expand(ix2))

      if(ix3.eq.1) isiz1=size(ans)
      if(ix3.gt.1) isiz1=isiz2
      call xx_fill(xinfo3,xx(ix3),jderiv3,isiz1,expand(ix3))

      kk(1)=iadr
      kk(2)=1
      kk(3)=2
      kk(4)=3

      isign=1
      ccwflag1=xinfo1%ccwflag
      ccwflag2=xinfo2%ccwflag
      ccwflag3=xinfo3%ccwflag
      if(.not.ccwflag1) then
         if((jderiv1.eq.1).or.(jderiv1.eq.3)) isign=-isign
      endif
      if(.not.ccwflag2) then
         if((jderiv2.eq.1).or.(jderiv2.eq.3)) isign=-isign
      endif
      if(.not.ccwflag3) then
         if((jderiv3.eq.1).or.(jderiv3.eq.3)) isign=-isign
      endif

      call xplasma_evala(s,ans,kk,xx)

      if(isign.eq.-1) then
         ans=ans*isign
      endif

      call xx_cleanup(xx(1),expand(1))
      call xx_cleanup(xx(2),expand(2))
      call xx_cleanup(xx(3),expand(3))

    end subroutine xplasma_xpeval_3dprof

    !----------------------------------------------------
    subroutine xplasma_rhoth_inv(s,rvec,zvec,rhovec,thvec,nregion,ier)

      ! *public* but specialized -- rough inverse map for
      ! axisymmetric configuration
      !   ...used as initial guess; done in "kernel" for efficiency...

      type(xplasma), pointer :: s

      real*8, dimension(:), intent(in) :: rvec,zvec
      real*8, dimension(:), intent(out) :: rhovec,thvec
      integer, dimension(:), intent(out) :: nregion
      integer, intent(out) :: ier

      !--------------------------------------------
      integer :: id_R,iadR,ivec
      integer :: idrho,idth,inrho,inth,iarho,iath
      logical :: extrap
      integer :: id_Rx,iadRx,iadax
      integer :: idrhox,idthx,inrhox,inthx,iarhox,iathx
      integer :: id_atan0,id_atanx,id_atan00,id_atanx0,id_dist0,id_distx
      integer :: iin,iout,iad,iat,ia0,iatx,iadx
      real*8 :: atanv(size(rvec)),distv(size(rvec))
      real*8 :: raxis,zaxis,zdra,zdza
      !--------------------------------------------

      call xplasma_common_ids(s,ier, id_R=id_R)
      if(ier.ne.0) return

      if(.not.s%axisymm) then
         ier=107
         call xplasma_errmsg_append(s, &
              '  inverse map initial guess (xplasma_rhoth_inv):'// &
              '  axisymmetry is assumed.')
         return
      endif

      if(id_R.le.0) then
         ier=9999
         call xplasma_errmsg_append(s, &
              '  call to xplasma_rhoth_inv, but R & Z surfaces undefined.')
         return
      endif

      if((size(rvec).ne.size(zvec)).or.(size(rvec).ne.size(rhovec)).or. &
           (size(rvec).ne.size(thvec)).or.(size(rvec).ne.size(nregion))) then
         ier=9999
         call xplasma_errmsg_append(s, &
             '  argument vector sizes inconsistent in xplasma_rhoth_inv call.')
         return
      endif

      iadR = s%dict(id_R)%dindex

      idth=s%profs(iadR)%gridIds(1)
      idrho=s%profs(iadR)%gridIds(2)

      iath=s%dict(idth)%dindex
      iarho=s%dict(idrho)%dindex

      inth=s%grids(iath)%size
      inrho=s%grids(iarho)%size

      call xplasma_find_item(s,'__dist_axis',id_dist0,ier)
      if(ier.ne.0) return

      call xplasma_find_item(s,'__atan_axis',id_atan0,ier)
      if(ier.ne.0) return

      call xplasma_find_item(s,'__atan00_axis',id_atan00,ier)
      if(ier.ne.0) return

      extrap = (s%id_Rx .gt. 0)
      if(extrap) then

         id_Rx = s%id_Rx

         iadRx = s%dict(id_Rx)%dindex

         idthx=s%profs(iadRx)%gridIds(1)
         idrhox=s%profs(iadRx)%gridIds(2)

         iathx=s%dict(idthx)%dindex
         iarhox=s%dict(idrhox)%dindex

         inthx=s%grids(iathx)%size
         inrhox=s%grids(iarhox)%size

         call xplasma_find_item(s,'__distx_axis',id_distx,ier)
         if(ier.ne.0) return

         call xplasma_find_item(s,'__atanx_axis',id_atanx,ier)
         if(ier.ne.0) return

         call xplasma_find_item(s,'__atanx0_axis',id_atanx0,ier)
         if(ier.ne.0) return

      endif

      iadax=s%dict(s%id_axis)%dindex

      raxis=s%lists(iadax)%r8vals(1)
      zaxis=s%lists(iadax)%r8vals(2)

      do ivec=1,size(rvec)
         nregion(ivec)=0
         rhovec(ivec)=-1.0d0
         thvec(ivec)=0.0d0
         zdra=rvec(ivec)-raxis
         zdza=zvec(ivec)-zaxis
         distv(ivec)=sqrt(zdra*zdra+zdza*zdza)
         if(distv(ivec).gt.0.0d0) then
            atanv(ivec)=atan2(zdza,zdra)
         else
            atanv(ivec)=0
         endif
      enddo

      if(extrap) then
         iout=3
         iin=2
         iad=s%dict(id_distx)%dindex
         iat=s%dict(id_atanx)%dindex
         ia0=s%dict(id_atanx0)%dindex
         call eqi_dxa_search(iin,iout,size(rvec),distv,atanv, &
              inthx,s%grids(iathx)%xpkg,inrhox,s%grids(iarhox)%xpkg, &
              s%profs(iad)%buf,s%profs(iat)%buf,s%profs(ia0)%buf, &
              inthx,s%grids(iathx)%xpkg, &
              s%profs(iad)%buf((inrhox-1)*inthx+1:inrhox*inthx), &
              s%profs(iat)%buf((inrhox-1)*inthx+1:inrhox*inthx), &
              rhovec,thvec,nregion)
         iatx=iat
         iadx=iad

         iout=0
         iin=1
         iad=s%dict(id_dist0)%dindex
         iat=s%dict(id_atan0)%dindex
         ia0=s%dict(id_atan00)%dindex
         call eqi_dxa_search(iin,iout,size(rvec),distv,atanv, &
              inth,s%grids(iath)%xpkg,inrho,s%grids(iarho)%xpkg, &
              s%profs(iad)%buf,s%profs(iat)%buf,s%profs(ia0)%buf, &
              inthx,s%grids(iathx)%xpkg, &
              s%profs(iadx)%buf(1:inthx), &
              s%profs(iatx)%buf(1:inthx), &
              rhovec,thvec,nregion)

      else

         iout=3
         iin=1
         iad=s%dict(id_dist0)%dindex
         iat=s%dict(id_atan0)%dindex
         ia0=s%dict(id_atan00)%dindex
         call eqi_dxa_search(iin,iout,size(rvec),distv,atanv, &
              inth,s%grids(iath)%xpkg,inrho,s%grids(iarho)%xpkg, &
              s%profs(iad)%buf,s%profs(iat)%buf,s%profs(ia0)%buf, &
              inth,s%grids(iath)%xpkg, &
              s%profs(iad)%buf((inrho-1)*inth+1:inrho*inth), &
              s%profs(iat)%buf((inrho-1)*inth+1:inrho*inth), &
              rhovec,thvec,nregion)

      endif

    end subroutine xplasma_rhoth_inv

    !----------------------------------------------------
    !  data mgmt routine (private)  -- see xx_cleanup below...

    subroutine xx_fill(xin,xx,jderiv,isize,iexpand)

      type(xpeval) :: xin,xx
      integer, intent(in) :: jderiv,isize
      logical, intent(out) :: iexpand

      call xpeval_free(xx)

      xx%idx = xin%idx
      xx%inum= isize
      xx%inx = xin%inx
      xx%ccwflag = xin%ccwflag

      allocate(xx%iderivs(1)); xx%iderivs=jderiv

      if((xin%inum.eq.1).and.(xin%inum.lt.isize)) then

         iexpand = .TRUE.
         allocate(xx%ix(isize),xx%dxn(isize),xx%hx(isize),xx%hxi(isize))

         xx%ix = xin%ix(1)
         xx%dxn = xin%dxn(1)
         xx%hx = xin%hx(1)
         xx%hxi = xin%hxi(1)

      else

         iexpand = .FALSE.

         xx%ix => xin%ix
         xx%dxn => xin%dxn
         xx%hx => xin%hx
         xx%hxi => xin%hxi

      endif
    end subroutine xx_fill
         
    !----------------------------------------------------
    !  data mgmt routine (private)

    subroutine xx_cleanup(xx,iexpand)

      type(xpeval) :: xx
      logical :: iexpand

      if(iexpand) then
         call xpeval_free(xx)

      else
         deallocate(xx%iderivs)
         call xplasma_xpeval_init(xx)   ! nullify pointers; do NOT deallocate...
      endif

    end subroutine xx_cleanup

    !----------------------------------------------------
    !  error checking routine (private)

    subroutine chkderiv(s,ideriv,ispline,id,ier)

      type (xplasma), pointer :: s
      integer, intent(in) :: ideriv,ispline,id
      integer, intent(out) :: ier

      character*128 msgbuf
      character*3 zdlbl

      ier = 0
      if((ideriv.lt.0).or.(ideriv.gt.3)) then
         ier=513
         write(msgbuf,*) &
              '  xplasma_eval_prof -- derivative not in range [0:3]: ', &
              ideriv
         call xplasma_errmsg_append(s,msgbuf)
      else if(((ideriv.gt.1).and.(ispline.lt.2)).or. &
           ((ideriv.gt.0).and.(ispline.lt.0))) then
         ier=513
         if(ideriv.eq.1) zdlbl = '1st'
         if(ideriv.eq.2) zdlbl = '2nd'
         if(ideriv.eq.3) zdlbl = '3rd'
         if(ispline.eq.-1) then
            write(msgbuf,*) '  xplasma_eval_prof -- ',zdlbl,' derivative', &
                 ' unavailable with step function interpolation.'
            call xplasma_errmsg_append(s,msgbuf)
         else if(ispline.eq.0) then
            write(msgbuf,*) '  xplasma_eval_prof -- ',zdlbl,' derivative', &
                 ' unavailable with piecewise linear interpolation.'
            call xplasma_errmsg_append(s,msgbuf)
         else if(ispline.eq.1) then
            write(msgbuf,*) '  xplasma_eval_prof -- ',zdlbl,' derivative', &
                 ' unavailable with C1 Hermite interpolation.'
            call xplasma_errmsg_append(s,msgbuf)
         endif
         write(msgbuf,*) '    profile name is: ',trim(s%dict(id)%name)
         call xplasma_errmsg_append(s,msgbuf)
      endif
    end subroutine chkderiv

    !---------------------------------------------------------
    subroutine xplasma_eval_1dprofx(s,id,x,ans,ier, &
         ideriv1,ccwflag1,force_bounds,n_out_of_bounds)

      !  evaluate a single profile f(x) 1d profile at a single point

      type (xplasma), pointer :: s
      integer, intent(in) :: id                   ! profile id
      real*8, intent(in) :: x                     ! evaluation point
      real*8, intent(out) :: ans                  ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/dx; 2 for d2f/dx2; 3 for d3f/dx3
      !    2 & 3 available only for C2 spline fits

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !-----------------------------
      real*8 :: xv(1),ansv(1)
      integer :: jderiv1,ioutside
      logical :: mccwflag1,mforce
      !-----------------------------

      call xplasma_errck1_prof(s,id,ier)
      if(ier.ne.0) return

      jderiv1=0
      if(present(ideriv1)) jderiv1=ideriv1

      mccwflag1=.TRUE.
      if(present(ccwflag1)) mccwflag1=ccwflag1

      mforce=.FALSE.
      if(present(force_bounds)) mforce=force_bounds

      xv(1) = x
      call xplasma_eval_1dprofxs(s,id,xv,ansv,ier,jderiv1,mccwflag1, &
           mforce,ioutside)

      ans = ansv(1) 

      if(present(n_out_of_bounds)) n_out_of_bounds=ioutside

    end subroutine xplasma_eval_1dprofx

    !---------------------------------------------------------
    subroutine xplasma_eval_1dprofsx(s,ids,x,ans,ier, &
         ideriv1,ideriv1s,ccwflag1,force_bounds,n_out_of_bounds)

      !  evaluate a set of 1d profiles f(x) all at a single point x
      !  all profiles must be defined over the same coordinate, though
      !  not necessarily the same grid.

      type (xplasma), pointer :: s
      integer, dimension(:), intent(in) :: ids    ! profile id
      real*8, intent(in) :: x                     ! evaluation x point
      real*8, dimension(:), intent(out) :: ans    ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/dx; 2 for d2f/dx2; 3 for d3f/dx3
      !    2 & 3 available only for C2 spline fits

      integer, intent(in), dimension(:), optional :: ideriv1s
      !  as ideriv1, but separately specified for each evaluation

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !------------------------------
      real*8 :: xv(1),ansv(1,size(ids))
      !------------------------------
      integer :: jderiv1s(size(ids)),ioutside
      logical :: mccwflag1,mforce
      character*128 msgbuf
      !-----------------------------

      ier=0
      if(size(ids).ne.size(ans)) then
         ier=510
         write(msgbuf,*) '  size of profile id list: ',size(ids), &
              '; size of result vector: ',size(ans)
         call xplasma_errmsg_append(s,msgbuf)
      endif
      if(present(ideriv1s)) then
         if(size(ideriv1s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv1s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv1)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv1 and ideriv1s both present.')
         endif
      endif
      if(ier.ne.0) return

      jderiv1s=0
      if(present(ideriv1s)) then
         jderiv1s=ideriv1s
      else if(present(ideriv1)) then
         jderiv1s=ideriv1
      endif

      mccwflag1=.TRUE.
      if(present(ccwflag1)) mccwflag1=ccwflag1

      mforce=.FALSE.
      if(present(force_bounds)) mforce=force_bounds

      xv(1) = x
      call xplasma_eval_1dprofsxs(s,ids,xv,ansv,ier, &
           ideriv1s=jderiv1s,ccwflag1=mccwflag1, &
           force_bounds=mforce,n_out_of_bounds=ioutside)
      ans(1:size(ids))=ansv(1,1:size(ids))

      if(present(n_out_of_bounds)) n_out_of_bounds=ioutside

    end subroutine xplasma_eval_1dprofsx

    !---------------------------------------------------------
    subroutine xplasma_eval_1dprofxs(s,idi,xs,ans,ier, &
         ideriv1,ccwflag1,force_bounds,n_out_of_bounds)

      !  evaluate a single profile f(x) 1d profile at a vector of points

      type (xplasma), pointer :: s
      integer, intent(in) :: idi                  ! profile id
      real*8, dimension(:), intent(in) :: xs      ! evaluation vector
      real*8, dimension(:), intent(out) :: ans    ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/dx; 2 for d2f/dx2; 3 for d3f/dx3
      !    2 & 3 available only for C2 spline fits

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !-----------------------------
      integer :: jderiv1
      logical :: mccwflag1

      integer :: id,id2,iadr2,idx,icoord,isign
      logical :: perio1,iforce

      integer :: iadr,iadx,iadc,ispline,irank,imode,iwarn

      integer :: inaxes  ! number of distinct x axes in evaluation group

      !  here the zone lookup data will be gathered:
      type(xpeval), dimension(:), allocatable :: xx  ! (1:inaxes)

      !  here the rank and indices of zone lookup objects are stored:
      integer, dimension(:), allocatable :: kk     ! (maxRank+1)

      character*128 msgbuf
      !-----------------------------
      !  check if id points to a profile object

      call xplasma_errck1_prof(s,idi,ier)
      if(ier.ne.0) return

      id=idi
      iadr=s%dict(id)%dindex
      irank=s%profs(iadr)%rank

      if(irank.ne.1) then
         id2=s%profs(iadr)%profIds(1)
         if(id2.ne.0) then
            iadr2=s%dict(id2)%dindex
            if(s%profs(iadr2)%rank.eq.1) then

               id=id2
               iadr=iadr2
               irank=1

            endif
         endif
      endif

      !---------------------------
      !  check argument sizes

      if(size(xs).ne.size(ans)) then
         ier=510
         write(msgbuf,*) '  size of profile x vector: ',size(xs), &
              '; size of result vector: ',size(ans)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      ispline=s%profs(iadr)%kspline

      !---------------------------
      !  check profile rank

      if(irank.ne.1) then
         msgbuf = '  xplasma_eval_1dprof -- '//trim(s%dict(id)%name)// &
              ' is not a 1d "f(x)" profile.'
         call xplasma_errmsg_append(s,msgbuf)
         ier=512
         return
      endif

      !---------------------------
      !  check derivative request consistent with interpolation setup for
      !  profile

      jderiv1=0
      if(present(ideriv1)) jderiv1=ideriv1

      call chkderiv(s,jderiv1,ispline,id,ier)
      if(ier.ne.0) return

      !---------------------------
      !  check for periodic coordinate; check for coordinate reversal e.g.
      !  due to clockwise theta angle interpretation in caller

      idx =s%profs(iadr)%gridIds(1)
      iadx=s%dict(idx)%dindex
      icoord=s%grids(iadx)%coord     ! coordinate associated with grid
      iadc=s%dict(icoord)%dindex

      perio1=s%coords(iadc)%periodic ! periodicity attribute

      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1=ccwflag1

      isign=1
      if(.not.mccwflag1) then
         if((jderiv1.eq.1).or.(jderiv1.eq.3)) isign=-1
      endif
      !---------------------------
      !  allocate lookup arrays and evaluate

      inaxes = 1
      allocate(xx(inaxes))

      xx(1)%inum = size(xs)
      xx(1)%inx = s%grids(iadx)%size
      allocate(xx(1)%iderivs(1))
      xx(1)%iderivs(1)=jderiv1
      xx(1)%ccwflag=mccwflag1

      !  verify that all xs elements are in range; for periodic
      !  grids normalize the range and perform CW->CCW reversal if necessary

      iforce=.FALSE.
      if(present(force_bounds)) iforce=force_bounds

      !  OK -- perform lookup

      allocate(xx(1)%ix(size(xs)),xx(1)%dxn(size(xs)), &
           xx(1)%hx(size(xs)),xx(1)%hxi(size(xs)))

      imode=2
      call xplookup(s,size(xs),xs,mccwflag1, &
           xx(1)%inx,s%grids(iadx)%xpkg,imode, &
           xx(1)%ix,xx(1)%dxn,xx(1)%hx,xx(1)%hxi,iwarn)

      if(present(n_out_of_bounds)) n_out_of_bounds = iwarn
      if(.not.iforce) then
         if(iwarn.gt.0) then
            ier=514
            call xpeval_free(xx(1))
            deallocate(xx)
            ans=0
            return
         endif
      endif

      allocate(kk(maxRank+1))

      kk(1)=iadr
      kk(2)=1  ! rank 1, only 1 x grid in single 1d profile eval.

      !---------------------------
      !  evaluate interpolations

      call xplasma_evala(s,ans,kk,xx)

      ans = ans*isign

      call xpeval_free(xx(1))
      deallocate(kk,xx)

    end subroutine xplasma_eval_1dprofxs

    !---------------------------------------------------------
    subroutine xplasma_eval_1dprofsxs(s,ids,xs,ans,ier, &
         ideriv1,ideriv1s,ccwflag1,force_bounds,n_out_of_bounds)

      !  evaluate multiple profiles f(x) 1d at a vector of points

      type (xplasma), pointer :: s
      integer, dimension(:), intent(in) :: ids    ! profile ids
      real*8, dimension(:), intent(in) :: xs      ! evaluation vector
      real*8, dimension(:,:), intent(out) :: ans  ! result of evaluations
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  ans(1:size(xs),1) -- result of eval of profile ids(1)
      !  ans(1:size(xs),2) -- result of eval of profile ids(2) ...etc...

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/dx; 2 for d2f/dx2; 3 for d3f/dx3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv1s
      !  as ideriv1 but specified separately for each profile

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds
      !-----------------------------
      integer :: jderiv1s(size(ids))
      logical :: mccwflag1,iforce

      integer :: ids2(size(ids)),idq
      integer :: idxs(size(ids))  ! list all x axes (with duplicates)
      integer :: icoord,jcoord,i,j
      logical :: perio1

      integer :: iadr,iadx,iadc,ispline,irank,iwarn,imode,isigns(size(ids))

      integer :: ifound,istat
      integer :: inaxes  ! number of distinct x axes in evaluation group
      integer :: iaxes(size(ids)) ! distinct x axes (no duplicates) (1:inaxes)
      integer :: indax(size(ids)) ! index into iaxes

      !  here the zone lookup data will be gathered:
      type(xpeval), dimension(:), allocatable :: xx  ! (1:inaxes)

      !  here the rank and indices of zone lookup objects are stored:
      integer, dimension(:,:), allocatable :: kk  ! (maxRank+1,1)

      character*128 msgbuf
      !-----------------------------

      ier=0

      !---------------------------
      !  check argument sizes
      
      if((size(xs).ne.size(ans,1)).or.(size(ids).ne.size(ans,2))) then
         ier=510
         write(msgbuf,*) '  size of profile x vector: ',size(xs), &
              '; size of result vectors: ',size(ans,1)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  number of profil ids: ',size(ids), &
              '; number of result vectors: ',size(ans,2)
         call xplasma_errmsg_append(s,msgbuf)
      endif
      if(present(ideriv1s)) then
         if(size(ideriv1s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv1s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv1)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv1 and ideriv1s both present.')
         endif
      endif
      if(ier.ne.0) return

      jderiv1s=0
      if(present(ideriv1)) jderiv1s=ideriv1
      if(present(ideriv1s)) jderiv1s=ideriv1s

      !---------------------------
      !  pre-check profiles: possibility of use of associated profiles
      !  enough to check rank

      ids2=ids
      do i=1,size(ids)
         call xplasma_prof_info(s,ids(i),ier, gridId1=idxs(i))
         if(ier.ne.0) cycle

         iadr=s%dict(ids(i))%dindex
         if(s%profs(iadr)%rank.ne.1) then
            idq=s%profs(iadr)%profIds(1)
            if(idq.gt.0) then
               iadr=s%dict(idq)%dindex
               if(s%profs(iadr)%rank.eq.1) then
                  ids2(i)=idq
               endif
            endif
         endif
      enddo

      !---------------------------
      !  check all profiles -- all must be rank 1 and defined over the
      !  same coordinate.  Any requested derivative must be consistent
      !  with the interpolating function associated with each profile.

      icoord=-1
      inaxes=0

      do i=1,size(ids2)

         !--- error here if ids2(i) not pointing at a profile object

         call xplasma_prof_info(s,ids2(i),ier, gridId1=idxs(i))
         if(ier.ne.0) exit

         iadr=s%dict(ids2(i))%dindex
         iadx=s%dict(idxs(i))%dindex
         jcoord=s%grids(iadx)%coord

         !--- check independent coordinate consistency in profile set

         if(icoord.eq.-1) then
            icoord=jcoord
            inaxes=1
            iaxes(1)=idxs(i)  ! =idxs(1), i.e. this is the first one.
            indax(1)=1
         else
            if(icoord.ne.jcoord) then
               ier=511
               write(msgbuf,*) '  ',trim(s%dict(ids2(i-1))%name),' and ', &
                    trim(s%dict(ids2(i))%name)
               call xplasma_errmsg_append(s,msgbuf)
               msgbuf = '  do not have the same independent coordinate, so,'
               call xplasma_errmsg_append(s,msgbuf)
               msgbuf = '  cannot be evaluated in the same xplasma_eval_prof call.'
               call xplasma_errmsg_append(s,msgbuf)
            else
               ifound=0
               do j=1,inaxes
                  if(idxs(i).eq.iaxes(j)) then
                     ifound=j
                     exit
                  endif
               enddo
               if(ifound.eq.0) then
                  inaxes=inaxes+1
                  iaxes(inaxes)=idxs(i)
                  indax(i)=inaxes
               else
                  indax(i)=ifound
               endif
            endif
         endif
         if(ier.ne.0) exit

         !--- check rank

         irank=s%profs(iadr)%rank
         if(irank.ne.1) then
            msgbuf = '  xplasma_eval_1dprof -- '//trim(s%dict(ids2(i))%name)// &
                 ' is not a 1d "f(x)" profile.'
            call xplasma_errmsg_append(s,msgbuf)
            ier=512
         endif
         if(ier.ne.0) exit

         !---------------------------
         !  check derivative request consistent with interpolation setup for
         !  all profiles

         ispline=s%profs(iadr)%kspline
         iwarn=0

         call chkderiv(s,jderiv1s(i),ispline,ids2(i),ier)
         if(ier.ne.0) exit

      enddo
      if(ier.ne.0) return

      !---------------------------
      !  check for periodic coordinate; check for coordinate reversal e.g.
      !  due to clockwise theta angle interpretation in caller

      call xplasma_coord_isPeriodic(s,icoord,perio1,ier)

      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1=ccwflag1

      isigns=1
      if(.not.mccwflag1) then
         do j=1,size(ids2)
            if((jderiv1s(j).eq.1).or.(jderiv1s(j).eq.3)) isigns(j)=-1
         enddo
      endif

      !---------------------------
      !  allocate lookup arrays and evaluate

      allocate(xx(inaxes))
      do i=1,inaxes
         xx(i)%inum = size(xs)

         iadx = s%dict(iaxes(i))%dindex
         xx(i)%inx = s%grids(iadx)%size

         allocate(xx(i)%iderivs(size(ids2)))
         xx(i)%iderivs = jderiv1s

         !  verify that all xs elements are in range; for periodic
         !  grids normalize the range and perform CW->CCW reversal if necessary

         iforce=.FALSE.
         if(present(force_bounds)) iforce=force_bounds

         !  OK -- perform lookup

         allocate(xx(i)%ix(size(xs)),xx(i)%dxn(size(xs)), &
              xx(i)%hx(size(xs)),xx(i)%hxi(size(xs)))

         imode=2
         call xplookup(s,size(xs),xs,mccwflag1, &
              s%grids(iadx)%size,s%grids(iadx)%xpkg,imode, &
              xx(i)%ix,xx(i)%dxn,xx(i)%hx,xx(i)%hxi,iwarn)

         if(present(n_out_of_bounds)) n_out_of_bounds = iwarn
         if(.not.iforce) then
            if(iwarn.gt.0) then
               ier=514
               exit
            endif
         endif

      enddo

      if(ier.ne.0) then
         do j=1,i
            call xpeval_free(xx(j))
         enddo
         deallocate(xx)
         return
      endif

      !  set up eval table for each profile

      allocate(kk(maxRank+1,size(ids2)))
      do j=1,size(ids2)
         iadr=s%dict(ids2(j))%dindex
         kk(1,j)=iadr
         kk(2,j)=indax(j)
      enddo

      !---------------------------
      !  evaluate interpolations

      call xplasma_eval(s,ans,kk,xx)
      do j=1,size(ids2)
         if(isigns(j).ne.1) then
            ans(1:size(xs),j)=ans(1:size(xs),j)*isigns(j)
         endif
      enddo

      do j=1,inaxes
         call xpeval_free(xx(j))
      enddo
      deallocate(xx,kk)
         
    end subroutine xplasma_eval_1dprofsxs

    !---------------------------------------------------------
    subroutine xplasma_eval_2dprofx(s,idf,idx1,x1,idx2,x2,ans,ier, &
         ideriv1,ideriv2,ccwflag1,ccwflag2,force_bounds,n_out_of_bounds)

      !  evaluate a single profile f(x1,x2) 2d profile at a single point

      type (xplasma), pointer :: s
      integer, intent(in) :: idf                  ! profile id
      integer, intent(in) :: idx1                 ! id of x1 grid or coord.
      real*8, intent(in) :: x1                    ! evaluation point x1
      integer, intent(in) :: idx2                 ! id of x2 grid or coord.
      real*8, intent(in) :: x2                    ! evaluation point x2
      real*8, intent(out) :: ans                  ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/d[x1]; 2 for d2f/d[x1]2; 3 for d3f/d[x1]3
      !    2 & 3 available only for C2 spline fits

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv2    ! derivative control
      !    1 for df/d[x2]; 2 for d2f/d[x2]2; 3 for d3f/d[x2]3
      !    2 & 3 available only for C2 spline fits

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1  ! for x1
      logical, intent(in), optional :: ccwflag2  ! for x2
      !    (these are ignored except for theta and phi angle coordinates

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !-----------------------------
      real*8 :: x1v(1),x2v(1),ansv(1)
      integer :: jderiv1,jderiv2,ioutside
      logical :: mccwflag1,mccwflag2,mforce
      !-----------------------------

      call xplasma_errck1_prof(s,idf,ier)
      if(ier.ne.0) return

      jderiv1=0
      if(present(ideriv1)) jderiv1=ideriv1
      jderiv2=0
      if(present(ideriv2)) jderiv2=ideriv2

      mccwflag1=.TRUE.
      if(present(ccwflag1)) mccwflag1=ccwflag1
      mccwflag2=.TRUE.
      if(present(ccwflag2)) mccwflag2=ccwflag2

      mforce=.FALSE.
      if(present(force_bounds)) mforce=force_bounds

      x1v(1) = x1
      x2v(1) = x2
      call xplasma_eval_2dprofxs(s,idf,idx1,x1v,idx2,x2v,ansv,ier, &
           jderiv1,jderiv2,mccwflag1,mccwflag2, &
           mforce,ioutside)

      ans = ansv(1) 

      if(present(n_out_of_bounds)) n_out_of_bounds=ioutside

    end subroutine xplasma_eval_2dprofx

    !---------------------------------------------------------
    subroutine xplasma_eval_2dprofsx(s,ids,idx1,x1,idx2,x2,ans,ier, &
         ideriv1,ideriv1s,ideriv2,ideriv2s, &
         ccwflag1,ccwflag2,force_bounds,n_out_of_bounds)

      !  evaluate a set of 2d profiles f(x1,x2) all at a single point
      !  all profiles must be defined over the same coordinate2, though
      !  not necessarily the same grids.

      type (xplasma), pointer :: s
      integer, dimension(:), intent(in) :: ids    ! profile id
      integer :: idx1                             ! id of x1 grid or coord.
      real*8, intent(in) :: x1                    ! evaluation x1 point
      integer :: idx2                             ! id of x2 grid or coord.
      real*8, intent(in) :: x2                    ! evaluation x2 point
      real*8, dimension(:), intent(out) :: ans    ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/d[x1]; 2 for d2f/d[x1]2; 3 for d3f/d[x1]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv1s
      !  as ideriv1, 1 per evaluation profile id

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv2    ! derivative control
      !    1 for df/d[x2]; 2 for d2f/d[x2]2; 3 for d3f/d[x2]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv2s
      !  as ideriv2, 1 per evaluation profile id

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1  ! for x1
      logical, intent(in), optional :: ccwflag2  ! for x2
      !    (these are ignored except for theta and phi angle coordinates

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !------------------------------
      real*8 :: x1v(1),x2v(1),ansv(1,size(ids))
      !------------------------------
      integer :: jderiv1s(size(ids)),jderiv2s(size(ids)),ioutside
      logical :: mccwflag1,mccwflag2,mforce
      character*128 msgbuf
      !-----------------------------

      ier=0
      if(size(ids).ne.size(ans)) then
         ier=510
         write(msgbuf,*) '  size of profile id list: ',size(ids), &
              '; size of result vector: ',size(ans)
         call xplasma_errmsg_append(s,msgbuf)
      endif
      if(present(ideriv1s)) then
         if(size(ideriv1s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv1s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv1)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv1 and ideriv1s both present.')
         endif
      endif
      if(present(ideriv2s)) then
         if(size(ideriv2s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv2s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv2)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv2 and ideriv2s both present.')
         endif
      endif
      if(ier.ne.0) return

      jderiv1s=0
      if(present(ideriv1)) jderiv1s=ideriv1
      if(present(ideriv1s)) jderiv1s=ideriv1s
      jderiv2s=0
      if(present(ideriv2)) jderiv2s=ideriv2
      if(present(ideriv2s)) jderiv2s=ideriv2s

      mccwflag1=.TRUE.
      if(present(ccwflag1)) mccwflag1=ccwflag1
      mccwflag2=.TRUE.
      if(present(ccwflag2)) mccwflag2=ccwflag2

      mforce=.FALSE.
      if(present(force_bounds)) mforce=force_bounds

      x1v(1) = x1
      x2v(1) = x2
      call xplasma_eval_2dprofsxs(s,ids,idx1,x1v,idx2,x2v,ansv,ier, &
           ideriv1s=jderiv1s,ideriv2s=jderiv2s, &
           ccwflag1=mccwflag1,ccwflag2=mccwflag2, &
           force_bounds=mforce,n_out_of_bounds=ioutside)
      ans(1:size(ids))=ansv(1,1:size(ids))

      if(present(n_out_of_bounds)) n_out_of_bounds=ioutside

    end subroutine xplasma_eval_2dprofsx

    !---------------------------------------------------------
    subroutine xplasma_eval_2dprofxs(s,idi,idx1,x1s,idx2,x2s,ans,ier, &
         ideriv1,ideriv2, &
         ccwflag1,ccwflag2,force_bounds,n_out_of_bounds)

      !  evaluate a single profile f(x1,x2) 2d profile at a vector of points

      type (xplasma), pointer :: s
      integer, intent(in) :: idi                  ! profile id
      integer :: idx1                             ! id of x1 grid or coord.
      real*8, dimension(:), intent(in) :: x1s     ! evaluation vector x1
      integer :: idx2                             ! id of x2 grid or coord.
      real*8, dimension(:), intent(in) :: x2s     ! evaluation vector x2
      real*8, dimension(:), intent(out) :: ans    ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/d[x1]; 2 for d2f/d[x1]2; 3 for d3f/d[x1]3
      !    2 & 3 available only for C2 spline fits

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv2    ! derivative control
      !    1 for df/d[x2]; 2 for d2f/d[x2]2; 3 for d3f/d[x2]3
      !    2 & 3 available only for C2 spline fits

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1  ! for x1
      logical, intent(in), optional :: ccwflag2  ! for x2
      !    (these are ignored except for theta and phi angle coordinates

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !-----------------------------
      integer :: jderiv1,jderiv2
      logical :: mccwflag1,mccwflag2

      integer :: icoord1,icoord2
      logical :: perio1,perio2,iforce
      integer :: idf,idq,idx,jcoord1,jcoord2,ivec,isign,iertmp,i,istat

      integer :: iadr,iadx1,iadx2,iadc1,iadc2,ispline,irank,imode,iwarn
      integer :: hspline(2)

      integer :: inaxes  ! number of distinct x axes in evaluation group

      !  here the zone lookup data will be gathered:
      type(xpeval), dimension(:), allocatable :: xx  ! (1:inaxes)

      !  here the rank and indices of zone lookup objects are stored:
      integer, dimension(:), allocatable :: kk     ! (maxRank+1,1)

      character*128 msgbuf

      integer  :: k_tmp

      !-----------------------------
      !  check if id points to a profile object

      idf=idi

      call xplasma_errck1_prof(s,idf,ier)
      if(ier.ne.0) return

      !-----------------------------
      !  precheck if need to use associated profile

      call ck_coords(s,idx1,idx2,0,idf,istat)
      if(istat.gt.2) then
         iadr=s%dict(idf)%dindex
         idq=s%profs(iadr)%profIds(1)
         if(idq.gt.0) then
            call ck_coords(s,idx1,idx2,0,idq,istat)
            if(istat.le.2) then

               idf=idq  ! use assoc. profile

            endif
         endif
      endif

      !---------------------------
      !  check argument sizes

      if((size(x1s).ne.size(ans)).or.(size(x2s).ne.size(ans))) then
         ier=510
         write(msgbuf,*) '  size of result vector: ',size(ans)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x1 vector: ',size(x1s)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x2 vector: ',size(x2s)
         call xplasma_errmsg_append(s,msgbuf)
         return
      else
         ivec=size(x1s)
      endif

      iadr=s%dict(idf)%dindex
      irank=s%profs(iadr)%rank
      ispline=s%profs(iadr)%kspline
      call spline_type_decode(ispline,hspline,2)

      !---------------------------
      !  check profile rank

      if(irank.gt.2) then
         msgbuf = '  xplasma_eval_2dprof -- '//trim(s%dict(idf)%name)// &
              ' is not a 2d "f(x1,x2)" profile; rank > 2'
         call xplasma_errmsg_append(s,msgbuf)
         ier=512
         return
      endif

      !---------------------------
      !  check derivative request consistent with interpolation setup for
      !  profile

      jderiv1=0
      if(present(ideriv1)) jderiv1=ideriv1
      jderiv2=0
      if(present(ideriv2)) jderiv2=ideriv2

      !---------------------------
      !  check for periodic coordinate; check for coordinate reversal e.g.
      !  due to clockwise theta angle interpretation in caller

      if((idx1.ge.1).and.(idx1.le.s%nitems)) then
         k_tmp = s%dict(idx1)%dtype
      else
         k_tmp = 0
      endif
      if((idx1.lt.1).or.(idx1.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x1 grid/coord id=',idx1, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq. xplasma_gridType ) then
         iadx1=s%dict(idx1)%dindex
         icoord1=s%grids(iadx1)%coord     ! coordinate associated with grid x1
      else if( k_tmp .eq. xplasma_coordType) then
         icoord1=idx1
      else
         ier=105
         write(msgbuf,*) &
              '  x1 grid id points to "',trim(s%dict(idx1)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if((idx2.ge.1).and.(idx2.le.s%nitems)) then
         k_tmp = s%dict(idx2)%dtype
      else
         k_tmp = 0
      endif
      if((idx2.lt.1).or.(idx2.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x2 grid/coord id=',idx2, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq. xplasma_gridType) then
         iadx2=s%dict(idx2)%dindex
         icoord2=s%grids(iadx2)%coord     ! coordinate associated with grid x2
      else if( k_tmp .eq. xplasma_coordType) then
         icoord2=idx2
      else
         ier=105
         write(msgbuf,*) &
              '  x2 grid id points to "',trim(s%dict(idx2)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if(ier.ne.0) return

      if(icoord1.eq.icoord2) then
         ier=515
         write(msgbuf,*) '  profile name is: ', trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  coordinate is:   ', trim(s%dict(icoord1)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      iadc1=s%dict(icoord1)%dindex
      perio1=s%coords(iadc1)%periodic ! periodicity attribute
      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1=ccwflag1

      iadc2=s%dict(icoord2)%dindex
      perio2=s%coords(iadc2)%periodic ! periodicity attribute
      mccwflag2=.TRUE.
      if(perio2.and.present(ccwflag2)) mccwflag2=ccwflag2

      !  check for sign change of odd derivative evaluation-- for each
      !  periodic coordinate reversal...

      isign = 1
      if(.not.mccwflag1) then
         if((jderiv1.eq.1).or.(jderiv1.eq.3)) then
            isign=-isign
         endif
      endif
      if(.not.mccwflag2) then
         if((jderiv2.eq.1).or.(jderiv2.eq.3)) then
            isign=-isign
         endif
      endif
      !---------------------------
      !  allocate lookup arrays and evaluate
      !    the function being evaluated may be rank 1 or 2
      !    the dependent coordinate values provided must supply enough
      !    information to evaluate the function, information beyond this
      !    is allowed but ignored.

      iforce=.FALSE.
      if(present(force_bounds)) iforce=force_bounds
      if(present(n_out_of_bounds)) n_out_of_bounds = 0

      inaxes = irank
      allocate(xx(inaxes))

      idx=s%profs(iadr)%gridIds(1)
      iadx1=s%dict(idx)%dindex     ! profile coordinate #1
      jcoord1=s%grids(iadx1)%coord
      
      xx(1)%inum = ivec
      xx(1)%inx = s%grids(iadx1)%size
      allocate(xx(1)%iderivs(1))

      !  verify that all x1s or x2s elements are in range; for periodic
      !  grids normalize the range and perform CW->CCW reversal if necessary

      if(jcoord1.eq.icoord1) then

         xx(1)%iderivs=jderiv1
         xx(1)%ccwflag=mccwflag1

      else if(jcoord1.eq.icoord2) then

         xx(1)%iderivs=jderiv2
         xx(1)%ccwflag=mccwflag2

      else

         ier=511
         write(msgbuf,*) ' xplasma_eval_prof: 1st coordinate of profile "', &
              trim(s%dict(idf)%name),'"'
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  is: "',trim(s%dict(jcoord1)%name),'" but there was'
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  no such data provided in the subroutine call.'
         call xplasma_errmsg_append(s,msgbuf)

      endif

      ! check that requested derivative and interpolation method are
      ! consistent

      if(ier.eq.0) then
         call chkderiv(s,xx(1)%iderivs(1),hspline(1),idf,ier)
      endif

      if(irank.eq.1) then
         ! if derivative occurs along coordinate upon which fcn does not
         ! depend, the result is ZERO!

         if(jcoord1.eq.icoord1) then
            if(jderiv2.gt.0) isign=0
         else
            if(jderiv1.gt.0) isign=0
         endif
         
      else if(irank.eq.2) then

         idx=s%profs(iadr)%gridIds(2)
         iadx2=s%dict(idx)%dindex     ! profile coordinate #1
         jcoord2=s%grids(iadx2)%coord
      
         xx(2)%inum = ivec
         xx(2)%inx = s%grids(iadx2)%size
         allocate(xx(2)%iderivs(1))

         if(jcoord2.eq.icoord1) then

            xx(2)%iderivs=jderiv1
            xx(2)%ccwflag=mccwflag1

         else if(jcoord2.eq.icoord2) then

            xx(2)%iderivs=jderiv2
            xx(2)%ccwflag=mccwflag2

         else

            ier=511
            write(msgbuf,*) &
                 ' xplasma_eval_prof: 2nd coordinate of profile "', &
                 trim(s%dict(idf)%name),'"'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) &
                 '  is: "',trim(s%dict(jcoord2)%name),'" but there was'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  no such data provided in the subroutine call.'
            call xplasma_errmsg_append(s,msgbuf)

         endif

         ! check that requested derivative and interpolation method are
         ! consistent

         if(ier.eq.0) then
            call chkderiv(s,xx(2)%iderivs(1),hspline(2),idf,ier)
         endif

      endif
      if(ier.ne.0) return

      !  OK -- perform lookup

      imode=2

      allocate(xx(1)%ix(ivec),xx(1)%dxn(ivec), &
           xx(1)%hx(ivec),xx(1)%hxi(ivec))
      if(jcoord1.eq.icoord1) then
         call xplookup(s,ivec,x1s,xx(1)%ccwflag, &
              xx(1)%inx,s%grids(iadx1)%xpkg,imode, &
              xx(1)%ix,xx(1)%dxn,xx(1)%hx,xx(1)%hxi,iwarn)
      else if(ier.eq.0) then
         call xplookup(s,ivec,x2s,xx(1)%ccwflag, &
              xx(1)%inx,s%grids(iadx1)%xpkg,imode, &
              xx(1)%ix,xx(1)%dxn,xx(1)%hx,xx(1)%hxi,iwarn)
      endif
      if(present(n_out_of_bounds)) n_out_of_bounds=max(n_out_of_bounds,iwarn)
      if(.not.iforce) then
         if(iwarn.gt.0) then
            ier=514
         endif
      endif

      if(irank.eq.2) then

         allocate(xx(2)%ix(ivec),xx(2)%dxn(ivec), &
              xx(2)%hx(ivec),xx(2)%hxi(ivec))

         if(jcoord2.eq.icoord1) then
            call xplookup(s,ivec,x1s,xx(2)%ccwflag, &
                 xx(2)%inx,s%grids(iadx2)%xpkg,imode, &
                 xx(2)%ix,xx(2)%dxn,xx(2)%hx,xx(2)%hxi,iwarn)
         else if(ier.eq.0) then
            call xplookup(s,ivec,x2s,xx(2)%ccwflag, &
                 xx(2)%inx,s%grids(iadx2)%xpkg,imode, &
                 xx(2)%ix,xx(2)%dxn,xx(2)%hx,xx(2)%hxi,iwarn)
         endif
         if(present(n_out_of_bounds)) n_out_of_bounds=max(n_out_of_bounds,iwarn)
         if(.not.iforce) then
            if(iwarn.gt.0) then
               ier=514
            endif
         endif

      endif

      if(ier.ne.0) then

         call xpeval_free(xx(1))
         if(irank.eq.2) call xpeval_free(xx(2))
         deallocate(xx)
         ans=0
         return
      endif

      allocate(kk(maxRank+1))

      kk(1)=iadr
      kk(2)=1  ! rank 1, only 1 x grid in single 1d profile eval.
      if(irank.eq.2) then
         kk(3)=2
      endif

      !---------------------------
      !  evaluate interpolations

      call xplasma_evala(s,ans,kk,xx)

      ans = ans*isign

      do i=1,irank
         call xpeval_free(xx(i))
      enddo
      deallocate(kk,xx)

    end subroutine xplasma_eval_2dprofxs

    !---------------------------------------------------------
    subroutine xplasma_eval_2dprofsxs(s,ids,idx1,x1s,idx2,x2s,ans,ier, &
         ideriv1,ideriv1s,ideriv2,ideriv2s, &
         ccwflag1,ccwflag2,force_bounds,n_out_of_bounds)

      !  evaluate multiple profiles f(x1,x2) 2d at a vector of points

      type (xplasma), pointer :: s
      integer, dimension(:), intent(in) :: ids    ! profile ids
      integer :: idx1                             ! id of x1 grid or coord.
      real*8, dimension(:), intent(in) :: x1s     ! evaluation vector x1
      integer :: idx2                             ! id of x2 grid or coord.
      real*8, dimension(:), intent(in) :: x2s     ! evaluation vector x2
      real*8, dimension(:,:), intent(out) :: ans  ! result of evaluations
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  size(x1s)=size(x2s)=size(ans,1) expected...

      !  ans(1:size(x1s),1) -- result of eval of profile ids(1)
      !  ans(1:size(x1s),2) -- result of eval of profile ids(2) ...etc...

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/d[x1]; 2 for d2f/d[x1]2; 3 for d3f/d[x1]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv1s
      !  as ideriv1, 1 per evaluation profile id

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv2    ! derivative control
      !    1 for df/d[x2]; 2 for d2f/d[x2]2; 3 for d3f/d[x2]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv2s
      !  as ideriv2, 1 per evaluation profile id

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1  ! for x1
      logical, intent(in), optional :: ccwflag2  ! for x2
      !    (these are ignored except for theta and phi angle coordinates

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !-----------------------------
      integer :: jderiv1s(size(ids)),jderiv2s(size(ids))
      logical :: mccwflag1,mccwflag2,iforce

      integer :: idxs(2,size(ids))  ! list all x axes (with duplicates)
      !  if any rank 1 functions are evaluated with this routine, 0s are
      !  inserted...

      integer :: icoord1,icoord2,jcoord1,jcoord2,i,j,ivec
      logical :: perio1,perio2

      integer :: idq,ids2(size(ids))
      integer :: iadr,iadx1,iadx2,iadc1,iadc2,ispline,irank,iwarn,imode
      integer :: idxtmp(maxRank),iadx,iertmp
      integer :: hspline(2),indc

      integer :: ifound,istat
      integer :: inaxes  ! number of distinct x axes in evaluation group
      integer :: iaxes(2*size(ids)) ! distinct x axes (no duplicates)(1:inaxes)
      integer :: icoor(2*size(ids)) ! which coordinate vector 1 or 2
      integer :: indax(2,size(ids)) ! index into iaxes
      integer :: isigns(size(ids))  ! result sign or 0 factors for all evals

      !  here the zone lookup data will be gathered:
      type(xpeval), dimension(:), allocatable :: xx  ! (1:inaxes)

      !  here the rank and indices of zone lookup objects are stored:
      integer, dimension(:,:), allocatable :: kk  ! (maxRank+1,1)

      character*128 msgbuf
      integer :: k_tmp
      !-----------------------------

      ier=0

      !---------------------------
      !  check argument sizes
      
      if((size(x1s).ne.size(ans,1)).or.(size(x1s).ne.size(x2s)).or. &
           (size(ids).ne.size(ans,2))) then
         ier=510
         write(msgbuf,*) '  size of result vector: ',size(ans,1)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x1 vector: ',size(x1s)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x2 vector: ',size(x2s)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  number of profil ids: ',size(ids), &
              '; number of result vectors: ',size(ans,2)
         call xplasma_errmsg_append(s,msgbuf)
      else
         ivec=size(x1s)
      endif
      if(present(ideriv1s)) then
         if(size(ideriv1s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv1s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv1)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv1 and ideriv1s both present.')
         endif
      endif
      if(present(ideriv2s)) then
         if(size(ideriv2s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv2s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv2)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv2 and ideriv2s both present.')
         endif
      endif
      if(ier.ne.0) return

      jderiv1s=0
      if(present(ideriv1)) jderiv1s=ideriv1
      if(present(ideriv1s)) jderiv1s=ideriv1s
      jderiv2s=0
      if(present(ideriv2)) jderiv2s=ideriv2
      if(present(ideriv2s)) jderiv2s=ideriv2s

      !  will check later if 2nd or 3rd derivatives can be evaluated for *all*
      !  profiles...

      !---------------------------
      !  precheck for need to use associated profiles

      ids2 = ids
      do i=1,size(ids)
         if(s%dict(ids(i))%dtype.ne.xplasma_profType) cycle
         call ck_coords(s,idx1,idx2,0,ids(i),istat)
         if(istat.gt.2) then
            iadr=s%dict(ids(i))%dindex
            idq=s%profs(iadr)%profIds(1)
            if(idq.gt.0) then
               call ck_coords(s,idx1,idx2,0,idq,istat)
               if(istat.le.2) then

                  ids2(i) = idq  ! use assoc. profile

               endif
            endif
         endif
      enddo

      !---------------------------
      !  check for periodic coordinate; check for coordinate reversal e.g.
      !  due to clockwise theta angle interpretation in caller
      
      if((idx1.ge.1).and.(idx1.le.s%nitems)) then
         k_tmp = s%dict(idx1)%dtype
      else
         k_tmp = 0
      endif
      if((idx1.lt.1).or.(idx1.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x1 grid/coord id=',idx1, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq. xplasma_gridType) then
         iadx1=s%dict(idx1)%dindex
         icoord1=s%grids(iadx1)%coord     ! coordinate associated with grid x1
      else if( k_tmp .eq. xplasma_coordType) then
         icoord1=idx1
      else
         ier=105
         write(msgbuf,*) &
              '  x1 grid id points to "',trim(s%dict(idx1)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if((idx2.ge.1).and.(idx2.le.s%nitems)) then
         k_tmp = s%dict(idx2)%dtype
      else
         k_tmp = 0
      endif
      if((idx2.lt.1).or.(idx2.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x2 grid/coord id=',idx2, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq. xplasma_gridType) then
         iadx2=s%dict(idx2)%dindex
         icoord2=s%grids(iadx2)%coord     ! coordinate associated with grid x2
      else if ( k_tmp .eq. xplasma_coordType) then
         icoord2=idx2
      else
         ier=105
         write(msgbuf,*) &
              '  x2 grid id points to "',trim(s%dict(idx2)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if(ier.ne.0) return

      if(icoord1.eq.icoord2) then
         ier=515
         write(msgbuf,*) '  coordinate is:   ', trim(s%dict(icoord1)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      iadc1=s%dict(icoord1)%dindex
      perio1=s%coords(iadc1)%periodic ! periodicity attribute
      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1=ccwflag1

      iadc2=s%dict(icoord2)%dindex
      perio2=s%coords(iadc2)%periodic ! periodicity attribute
      mccwflag2=.TRUE.
      if(perio2.and.present(ccwflag2)) mccwflag2=ccwflag2

      !  check for sign change of odd derivative evaluation-- for each
      !  periodic coordinate reversal...

      isigns = 1
      if(.not.mccwflag1) then
         do i=1,size(ids2)
            if((jderiv1s(i).eq.1).or.(jderiv1s(i).eq.3)) then
               isigns(i)=-isigns(i)
            endif
         enddo
      endif
      if(.not.mccwflag2) then
         do i=1,size(ids2)
            if((jderiv2s(i).eq.1).or.(jderiv2s(i).eq.3)) then
               isigns(i)=-isigns(i)
            endif
         enddo
      endif

      !---------------------------
      !  check all profiles -- all must be rank 1 or 2 and defined over
      !  the coordinate values provided in the subroutine arguments.
      !  Any requested derivative must be consistent
      !  with the interpolating function associated with each profile.

      inaxes=0

      do i=1,size(ids2)

         !--- error here if ids2(i) not pointing at a profile object

         call xplasma_prof_info(s,ids2(i),ier, &
              gridId1=idxtmp(1), gridId2=idxtmp(2), gridId3=idxtmp(3))
         if(ier.ne.0) exit

         if(idxtmp(3).gt.0) then
            msgbuf = &
                 '  xplasma_eval_2dprof -- '//trim(s%dict(ids2(i))%name)// &
                 ' is not a 2d "f(x1,x2)" profile; rank > 2'
            call xplasma_errmsg_append(s,msgbuf)
            ier=512
            exit
         else
            if(idxtmp(2).eq.0) then
               irank=1
            else
               irank=2
            endif
         endif

         iadr=s%dict(ids2(i))%dindex
         iadx1=s%dict(idxtmp(1))%dindex
         jcoord1=s%grids(iadx1)%coord
         iertmp=0
         if((jcoord1.ne.icoord1).and.(jcoord1.ne.icoord2)) then
            iertmp=511
            write(msgbuf,*) '  xplasma_eval_prof: profile "', &
                 trim(s%dict(ids2(i))%name),'" depends on'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  coordinate "',trim(s%dict(jcoord1)%name),'" but'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) '  no such data was provided for evaluation.'
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(irank.eq.2) then
            iadx2=s%dict(idxtmp(2))%dindex
            jcoord2=s%grids(iadx2)%coord
            if((jcoord2.ne.icoord1).and.(jcoord2.ne.icoord2)) then
               iertmp=511
               write(msgbuf,*) '  xplasma_eval_prof: profile "', &
                    trim(s%dict(ids2(i))%name),'" depends on'
               call xplasma_errmsg_append(s,msgbuf)

               write(msgbuf,*) '  coordinate "',trim(s%dict(jcoord1)%name),'" but'
               call xplasma_errmsg_append(s,msgbuf)
               write(msgbuf,*) '  no such data was provided for evaluation.'
               call xplasma_errmsg_append(s,msgbuf)
            endif
         else
            iadx2=0
            jcoord2=0
            !  rank=1 -- if derivative is on other coordinate answer is zero...
            if(jcoord1.eq.icoord1) then
               if(jderiv2s(i).gt.0) isigns(i)=0
            else
               if(jderiv1s(i).gt.0) isigns(i)=0
            endif
         endif
         if(iertmp.ne.0) then
            ier=iertmp
            exit
         endif

         !--- check independent coordinate consistency in profile set

         idxs(1:2,i)=idxtmp(1:2)
         indax(1:2,i)=0

         if(inaxes.eq.0) then
            inaxes=irank
            iaxes(1)=idxtmp(1)
            indax(1,1)=1
            if(jcoord1.eq.icoord1) then
               icoor(1)=1
            else
               icoor(1)=2
            endif
            if(irank.eq.2) then
               iaxes(2)=idxtmp(2)
               indax(2,1)=2
               if(jcoord2.eq.icoord1) then
                  icoor(2)=1
               else
                  icoor(2)=2
               endif
            endif
         else
            ifound=0
            do j=1,inaxes
               if(idxs(1,i).eq.iaxes(j)) then
                  ifound=j
                  exit
               endif
            enddo
            if(ifound.eq.0) then
               inaxes=inaxes+1
               iaxes(inaxes)=idxs(1,i)
               indax(1,i)=inaxes
               if(jcoord1.eq.icoord1) then
                  icoor(inaxes)=1
               else
                  icoor(inaxes)=2
               endif
            else
               indax(1,i)=ifound
            endif
            if(irank.eq.2) then
               ifound=0
               do j=1,inaxes
                  if(idxs(2,i).eq.iaxes(j)) then
                     ifound=j
                     exit
                  endif
               enddo
               if(ifound.eq.0) then
                  inaxes=inaxes+1
                  iaxes(inaxes)=idxs(2,i)
                  indax(2,i)=inaxes
                  if(jcoord2.eq.icoord1) then
                     icoor(inaxes)=1
                  else
                     icoor(inaxes)=2
                  endif
               else
                  indax(2,i)=ifound
               endif
            endif
         endif
         if(ier.ne.0) exit

         !---------------------------
         !  check derivative request consistent with interpolation setup for
         !  all profiles

         ispline=s%profs(iadr)%kspline
         call spline_type_decode(ispline,hspline,2)

         indc=icoor(indax(1,i))
         if(indc.eq.1) then
            call chkderiv(s,jderiv1s(i),hspline(1),ids2(i),ier)
         else
            call chkderiv(s,jderiv2s(i),hspline(1),ids2(i),ier)
         endif

         if(irank.eq.2) then
            iertmp=ier
            indc=icoor(indax(2,i))
            if(indc.eq.1) then
               call chkderiv(s,jderiv1s(i),hspline(2),ids2(i),ier)
            else
               call chkderiv(s,jderiv2s(i),hspline(2),ids2(i),ier)
            endif
            ier=max(iertmp,ier)
         endif

      enddo
      if(ier.ne.0) return

      !---------------------------
      !  check for periodic coordinate; check for coordinate reversal e.g.
      !  due to clockwise theta angle interpretation in caller

      call xplasma_coord_isPeriodic(s,icoord1,perio1,ier)
      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1=ccwflag1

      call xplasma_coord_isPeriodic(s,icoord2,perio2,ier)
      mccwflag2=.TRUE.
      if(perio2.and.present(ccwflag2)) mccwflag2=ccwflag2

      !---------------------------
      !  allocate lookup arrays and evaluate

      allocate(xx(inaxes))
      if(present(n_out_of_bounds)) n_out_of_bounds = 0

      do i=1,inaxes
         xx(i)%inum = ivec

         iadx = s%dict(iaxes(i))%dindex
         xx(i)%inx = s%grids(iadx)%size

         allocate(xx(i)%iderivs(size(ids2)))
         if(icoor(i).eq.1) then
            xx(i)%iderivs = jderiv1s
            xx(i)%ccwflag = mccwflag1
         else
            xx(i)%iderivs = jderiv2s
            xx(i)%ccwflag = mccwflag2
         endif

         !  verify that all xs elements are in range; for periodic
         !  grids normalize the range and perform CW->CCW reversal if necessary

         iforce=.FALSE.
         if(present(force_bounds)) iforce=force_bounds

         !  OK -- perform lookup

         allocate(xx(i)%ix(ivec),xx(i)%dxn(ivec), &
              xx(i)%hx(ivec),xx(i)%hxi(ivec))

         imode=2
         if(icoor(i).eq.1) then
            call xplookup(s,ivec,x1s,xx(i)%ccwflag, &
                 s%grids(iadx)%size,s%grids(iadx)%xpkg,imode, &
                 xx(i)%ix,xx(i)%dxn,xx(i)%hx,xx(i)%hxi,iwarn)
         else
            call xplookup(s,ivec,x2s,xx(i)%ccwflag, &
                 s%grids(iadx)%size,s%grids(iadx)%xpkg,imode, &
                 xx(i)%ix,xx(i)%dxn,xx(i)%hx,xx(i)%hxi,iwarn)
         endif

         if(present(n_out_of_bounds)) n_out_of_bounds = &
              max(n_out_of_bounds,iwarn)
         if(.not.iforce) then
            if(iwarn.gt.0) then
               ier=514
               exit
            endif
         endif
      enddo

      if(ier.ne.0) then
         do j=1,i
            call xpeval_free(xx(j))
         enddo
         deallocate(xx)
         ans=0
         return
      endif

      !  set up eval table for each profile

      allocate(kk(maxRank+1,size(ids2)))
      do j=1,size(ids2)
         iadr=s%dict(ids2(j))%dindex
         kk(1,j)=iadr
         kk(2:3,j)=indax(1:2,j)
      enddo

      !---------------------------
      !  evaluate interpolations

      call xplasma_eval(s,ans,kk,xx)

      do j=1,size(ids2)
         ans(1:ivec,j)=ans(1:ivec,j)*isigns(j)
      enddo
      
      deallocate(kk)
      do j=1,inaxes
         call xpeval_free(xx(j))
      enddo
      deallocate(xx)
   
    end subroutine xplasma_eval_2dprofsxs

    !---------------------------------------------------------
    subroutine xplasma_RZeval_2d(s,ids,idx1,x1s,idx2,x2s,ans,ier, &
         ideriv1,ideriv1s,ideriv2,ideriv2s, &
         ccwflag1,ccwflag2,force_bounds,n_out_of_bounds)

      !  use xplasma_eval_2dprofsxs to
      !  evaluate multiple profiles f(x1,x2) 2d at a vector of points
      !    IF all the profiles are "R" and "Z",
      !    and IF there is a scrape-off region,
      !    split the evaluation so that points inside the plasma use
      !    the standard R & Z, and points outside use the extrapolated
      !    profiles.

      type (xplasma), pointer :: s
      integer, dimension(:), intent(in) :: ids    ! profile ids
      integer :: idx1                             ! id of x1 grid or coord.
      real*8, dimension(:), intent(in) :: x1s     ! evaluation vector x1
      integer :: idx2                             ! id of x2 grid or coord.
      real*8, dimension(:), intent(in) :: x2s     ! evaluation vector x2
      real*8, dimension(:,:), intent(out) :: ans  ! result of evaluations
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  size(x1s)=size(x2s)=size(ans,1) expected...

      !  ans(1:size(x1s),1) -- result of eval of profile ids(1)
      !  ans(1:size(x1s),2) -- result of eval of profile ids(2) ...etc...

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/d[x1]; 2 for d2f/d[x1]2; 3 for d3f/d[x1]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv1s
      !  as ideriv1, 1 per evaluation profile id

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv2    ! derivative control
      !    1 for df/d[x2]; 2 for d2f/d[x2]2; 3 for d3f/d[x2]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv2s
      !  as ideriv2, 1 per evaluation profile id

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1  ! for x1
      logical, intent(in), optional :: ccwflag2  ! for x2
      !    (these are ignored except for theta and phi angle coordinates

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !----------------------------
      integer :: ipro,i,ivec,jvec,icoord1,icoord2,iertmp
      integer :: i_out,idx1x,idx2x
      integer :: idsx(size(ids))
      real*8 :: zrho
      real*8, dimension(:), allocatable :: x1tmp,x2tmp
      real*8, dimension(:,:), allocatable :: anstmp
      logical :: inorm,irho1
      real*8, parameter :: ONE = 1.0d0

      integer :: k_tmp

      !----------------------------
      
      if(present(n_out_of_bounds)) n_out_of_bounds=0

      if(.not.s%scrapeoff) then

         inorm=.TRUE.

      else

         inorm=.FALSE.
         do ipro=1,size(ids)
            if((ids(ipro).ne.s%id_R).and.(ids(ipro).ne.s%id_Z)) then
               inorm=.TRUE.  ! not R nor Z
               exit
            else
               if(ids(ipro).eq.s%id_R) then
                  idsx(ipro)=s%id_Rx
               else
                  idsx(ipro)=s%id_Zx
               endif
            endif
         enddo

      endif

      !  report errors in 2dprfsxs...

      if((idx1.le.0).or.(idx2.le.0).or. &
           (idx1.gt.s%nitems).or.(idx2.gt.s%nitems)) inorm=.TRUE.
      
      if(.not.inorm) then
         k_tmp = s%dict(idx1)%dtype
         icoord1=0
         if(s%dict(idx1)%dtype.eq.xplasma_coordType) then
            icoord1=idx1
         else if( k_tmp .eq.xplasma_gridType) then
            call xplasma_grid_info(s,idx1,iertmp, coord=icoord1)
         endif

         icoord2=0
         k_tmp = s%dict(idx2)%dtype
         if(s%dict(idx2)%dtype.eq.xplasma_coordType) then
            icoord2=idx2
         else if( k_tmp .eq.xplasma_gridType) then
            call xplasma_grid_info(s,idx2,iertmp, coord=icoord2)
         endif

         if(icoord1.eq.icoord2) then
            inorm=.TRUE.
         else if( ((icoord1.eq.xplasma_rho_coord).or. &
              (icoord1.eq.xplasma_theta_coord)) .and. &
              ((icoord2.eq.xplasma_rho_coord).or. &
              (icoord2.eq.xplasma_theta_coord))) then
            irho1=(icoord1.eq.xplasma_rho_coord)
            if(irho1) then
               if(maxval(x1s).le.ONE) inorm=.TRUE.
               if(.not.inorm) then
                  idx1x=xplasma_rhox_coord
                  idx2x=xplasma_thx_coord
               endif
            else
               if(maxval(x2s).le.ONE) inorm=.TRUE.
               if(.not.inorm) then
                  idx1x=xplasma_thx_coord
                  idx2x=xplasma_rhox_coord
               endif
            endif
         else
            inorm=.TRUE.
         endif
      endif

      
      if(inorm) then

         call xplasma_eval_2dprofsxs(s,ids,idx1,x1s,idx2,x2s,ans,ier, &
              ideriv1,ideriv1s,ideriv2,ideriv2s, &
              ccwflag1,ccwflag2,force_bounds,n_out_of_bounds)

      else

         !  R & Z only eval -- allow bdy straddle

         ipro=min(size(ans,2),size(ids))
         ivec=min(size(ans,1),size(x1s),size(x2s))
         allocate(x1tmp(ivec),x2tmp(ivec),anstmp(ivec,ipro))

         !  eval all rho <= 1.000 pts

         jvec=0
         do i=1,ivec
            if(irho1) then
               zrho=x1s(i)
            else
               zrho=x2s(i)
            endif
            if(zrho.le.ONE) then

               jvec=jvec+1
               x1tmp(jvec)=x1s(i)
               x2tmp(jvec)=x2s(i)

            endif
         enddo

         
         if(jvec.gt.0) then
            call xplasma_eval_2dprofsxs(s,ids, &
                 idx1,x1tmp(1:jvec), &
                 idx2,x2tmp(1:jvec), &
                 anstmp(1:jvec,1:ipro),ier, &
                 ideriv1,ideriv1s,ideriv2,ideriv2s, &
                 ccwflag1,ccwflag2,force_bounds,i_out)
            if(present(n_out_of_bounds)) n_out_of_bounds =n_out_of_bounds+i_out

            if(ier.ne.0) then
               ans=0
               deallocate(x1tmp,x2tmp,anstmp)
               return
            endif

            jvec=0
            do i=1,ivec
               if(irho1) then
                  zrho=x1s(i)
               else
                  zrho=x2s(i)
               endif
               if(zrho.le.ONE) then
                  
                  jvec=jvec+1
                  ans(i,1:ipro)=anstmp(jvec,1:ipro)

               endif
            enddo
         endif

         !  eval all rho > 1.000 pts

         jvec=0
         do i=1,ivec
            if(irho1) then
               zrho=x1s(i)
            else
               zrho=x2s(i)
            endif
            if(zrho.gt.ONE) then

               jvec=jvec+1
               x1tmp(jvec)=x1s(i)
               x2tmp(jvec)=x2s(i)

            endif
         enddo

         
         if(jvec.gt.0) then
            call xplasma_eval_2dprofsxs(s,idsx, &    ! external R & Z eval
                 idx1x,x1tmp(1:jvec), &
                 idx2x,x2tmp(1:jvec), &
                 anstmp(1:jvec,1:ipro),ier, &
                 ideriv1,ideriv1s,ideriv2,ideriv2s, &
                 ccwflag1,ccwflag2,force_bounds,i_out)
            if(present(n_out_of_bounds)) n_out_of_bounds =n_out_of_bounds+i_out

            if(ier.ne.0) then
               ans=0
               deallocate(x1tmp,x2tmp,anstmp)
               return
            endif

            jvec=0
            do i=1,ivec
               if(irho1) then
                  zrho=x1s(i)
               else
                  zrho=x2s(i)
               endif
               if(zrho.gt.ONE) then
                  
                  jvec=jvec+1
                  ans(i,1:ipro)=anstmp(jvec,1:ipro)

               endif
            enddo
         endif

         deallocate(x1tmp,x2tmp,anstmp)
      endif
         
    end subroutine xplasma_RZeval_2d

    !---------------------------------------------------------
    subroutine xplasma_eval_3dprofx(s,idf,idx1,x1,idx2,x2,idx3,x3,ans,ier, &
         ideriv1,ideriv2,ideriv3,ccwflag1,ccwflag2,ccwflag3, &
         force_bounds,n_out_of_bounds)

      !  evaluate a single profile f(x1,x2,x3) 3d profile at a single point

      type (xplasma), pointer :: s
      integer, intent(in) :: idf                  ! profile id
      integer, intent(in) :: idx1                 ! id of x1 grid or coord.
      real*8, intent(in) :: x1                    ! evaluation point x1
      integer, intent(in) :: idx2                 ! id of x2 grid or coord.
      real*8, intent(in) :: x2                    ! evaluation point x2
      integer, intent(in) :: idx3                 ! id of x3 grid or coord.
      real*8, intent(in) :: x3                    ! evaluation point x3
      real*8, intent(out) :: ans                  ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/d[x1]; 2 for d2f/d[x1]2; 3 for d3f/d[x1]3
      !    2 & 3 available only for C2 spline fits

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv2    ! derivative control
      !    1 for df/d[x2]; 2 for d2f/d[x2]2; 3 for d3f/d[x2]3
      !    2 & 3 available only for C2 spline fits

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv3    ! derivative control
      !    1 for df/d[x3]; 2 for d2f/d[x3]2; 3 for d3f/d[x3]3
      !    2 & 3 available only for C2 spline fits

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1  ! for x1
      logical, intent(in), optional :: ccwflag2  ! for x2
      logical, intent(in), optional :: ccwflag3  ! for x3
      !    (these are ignored except for theta and phi angle coordinates

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !-----------------------------
      real*8 :: x1v(1),x2v(1),x3v(1),ansv(1)
      integer :: jderiv1,jderiv2,jderiv3,ioutside
      logical :: mccwflag1,mccwflag2,mccwflag3,mforce
      !-----------------------------

      call xplasma_errck1_prof(s,idf,ier)
      if(ier.ne.0) return

      jderiv1=0
      if(present(ideriv1)) jderiv1=ideriv1
      jderiv2=0
      if(present(ideriv2)) jderiv2=ideriv2
      jderiv3=0
      if(present(ideriv3)) jderiv3=ideriv3

      mccwflag1=.TRUE.
      if(present(ccwflag1)) mccwflag1=ccwflag1
      mccwflag2=.TRUE.
      if(present(ccwflag2)) mccwflag2=ccwflag2
      mccwflag3=.TRUE.
      if(present(ccwflag3)) mccwflag3=ccwflag3

      mforce=.FALSE.
      if(present(force_bounds)) mforce=force_bounds

      x1v(1) = x1
      x2v(1) = x2
      x3v(1) = x3
      call xplasma_eval_3dprofxs(s,idf,idx1,x1v,idx2,x2v,idx3,x3v,ansv,ier, &
           jderiv1,jderiv2,jderiv3,mccwflag1,mccwflag2,mccwflag3, &
           mforce,ioutside)

      ans = ansv(1) 

      if(present(n_out_of_bounds)) n_out_of_bounds=ioutside

    end subroutine xplasma_eval_3dprofx

    !---------------------------------------------------------
    subroutine xplasma_eval_3dprofsx(s,ids,idx1,x1,idx2,x2,idx3,x3,ans,ier, &
         ideriv1,ideriv1s,ideriv2,ideriv2s,ideriv3,ideriv3s, &
         ccwflag1,ccwflag2,ccwflag3, &
         force_bounds,n_out_of_bounds)

      !  evaluate a set of 3d profiles f(x1,x2,x3) all at a single point
      !  all profiles must be defined over the same coordinates, though
      !  not necessarily the same grids.

      type (xplasma), pointer :: s
      integer, dimension(:), intent(in) :: ids    ! profile id
      integer :: idx1                             ! id of x1 grid or coord.
      real*8, intent(in) :: x1                    ! evaluation x1 point
      integer :: idx2                             ! id of x2 grid or coord.
      real*8, intent(in) :: x2                    ! evaluation x2 point
      integer :: idx3                             ! id of x3 grid or coord.
      real*8, intent(in) :: x3                    ! evaluation x3 point
      real*8, dimension(:), intent(out) :: ans    ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/d[x1]; 2 for d2f/d[x1]2; 3 for d3f/d[x1]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv1s
      !  as ideriv1, but separately specified for each evaluation

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv2    ! derivative control
      !    1 for df/d[x2]; 2 for d2f/d[x2]2; 3 for d3f/d[x2]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv2s
      !  as ideriv2, but separately specified for each evaluation

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv3    ! derivative control
      !    1 for df/d[x3]; 2 for d2f/d[x3]2; 3 for d3f/d[x3]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv3s
      !  as ideriv3, but separately specified for each evaluation

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1  ! for x1
      logical, intent(in), optional :: ccwflag2  ! for x2
      logical, intent(in), optional :: ccwflag3  ! for x3
      !    (these are ignored except for theta and phi angle coordinates

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !------------------------------
      real*8 :: x1v(1),x2v(1),x3v(1),ansv(1,size(ids))
      !------------------------------
      integer :: jderiv1s(size(ids)),jderiv2s(size(ids)),jderiv3s(size(ids))
      integer :: ioutside
      logical :: mccwflag1,mccwflag2,mccwflag3,mforce
      character*128 msgbuf
      !-----------------------------

      ier=0
      if(size(ids).ne.size(ans)) then
         ier=510
         write(msgbuf,*) '  size of profile id list: ',size(ids), &
              '; size of result vector: ',size(ans)
         call xplasma_errmsg_append(s,msgbuf)
      endif
      if(present(ideriv1s)) then
         if(size(ideriv1s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv1s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv1)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv1 and ideriv1s both present.')
         endif
      endif
      if(present(ideriv2s)) then
         if(size(ideriv2s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv2s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv2)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv2 and ideriv2s both present.')
         endif
      endif
      if(present(ideriv3s)) then
         if(size(ideriv3s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv3s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv3)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv3 and ideriv3s both present.')
         endif
      endif
      if(ier.ne.0) return

      jderiv1s=0
      if(present(ideriv1)) jderiv1s=ideriv1
      if(present(ideriv1s)) jderiv1s=ideriv1s
      jderiv2s=0
      if(present(ideriv2)) jderiv2s=ideriv2
      if(present(ideriv2s)) jderiv2s=ideriv2s
      jderiv3s=0
      if(present(ideriv3)) jderiv3s=ideriv3
      if(present(ideriv3s)) jderiv3s=ideriv3s

      mccwflag1=.TRUE.
      if(present(ccwflag1)) mccwflag1=ccwflag1
      mccwflag2=.TRUE.
      if(present(ccwflag2)) mccwflag2=ccwflag2
      mccwflag3=.TRUE.
      if(present(ccwflag3)) mccwflag3=ccwflag3

      mforce=.FALSE.
      if(present(force_bounds)) mforce=force_bounds

      x1v(1) = x1
      x2v(1) = x2
      x3v(1) = x3
      call xplasma_eval_3dprofsxs(s,ids,idx1,x1v,idx2,x2v,idx3,x3v,ansv,ier, &
           ideriv1s=jderiv1s,ideriv2s=jderiv2s,ideriv3s=jderiv3s, &
           ccwflag1=mccwflag1,ccwflag2=mccwflag2,ccwflag3=mccwflag3, &
           force_bounds=mforce,n_out_of_bounds=ioutside)
      ans(1:size(ids))=ansv(1,1:size(ids))

      if(present(n_out_of_bounds)) n_out_of_bounds=ioutside

    end subroutine xplasma_eval_3dprofsx

    !---------------------------------------------------------
    subroutine xplasma_eval_3dprofxs(s,idi,idx1,x1s,idx2,x2s,idx3,x3s,ans,ier,&
         ideriv1,ideriv2,ideriv3,ccwflag1,ccwflag2,ccwflag3, &
         force_bounds,n_out_of_bounds)

      !  evaluate a single profile f(x1,x2) 2d profile at a vector of points

      type (xplasma), pointer :: s
      integer, intent(in) :: idi                  ! profile id
      integer :: idx1                             ! id of x1 grid or coord.
      real*8, dimension(:), intent(in) :: x1s     ! evaluation vector x1
      integer :: idx2                             ! id of x2 grid or coord.
      real*8, dimension(:), intent(in) :: x2s     ! evaluation vector x2
      integer :: idx3                             ! id of x3 grid or coord.
      real*8, dimension(:), intent(in) :: x3s     ! evaluation vector x3

      real*8, dimension(:), intent(out) :: ans    ! result of evaluation
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/d[x1]; 2 for d2f/d[x1]2; 3 for d3f/d[x1]3
      !    2 & 3 available only for C2 spline fits

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv2    ! derivative control
      !    1 for df/d[x2]; 2 for d2f/d[x2]2; 3 for d3f/d[x2]3
      !    2 & 3 available only for C2 spline fits

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv3    ! derivative control
      !    1 for df/d[x3]; 2 for d2f/d[x3]2; 3 for d3f/d[x3]3
      !    2 & 3 available only for C2 spline fits

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1  ! for x1
      logical, intent(in), optional :: ccwflag2  ! for x2
      logical, intent(in), optional :: ccwflag3  ! for x3
      !    (these are ignored except for theta and phi angle coordinates

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !-----------------------------
      integer :: jderiv1,jderiv2,jderiv3
      logical :: mccwflag1,mccwflag2,mccwflag3

      integer :: idf,idq
      integer :: icoord1,icoord2,icoord3
      logical :: perio1,perio2,perio3,iforce
      integer :: idx,jcoord1,jcoord2,jcoord3,ivec,isign,iertmp,i

      integer :: iadr,iadx1,iadx2,iadx3,iadc1,iadc2,iadc3
      integer :: ispline,irank,imode,iwarn,istat,hspline(3)

      integer :: inaxes  ! number of distinct x axes in evaluation group

      !  here the zone lookup data will be gathered:
      type(xpeval), dimension(:), allocatable :: xx  ! (1:inaxes)

      !  here the rank and indices of zone lookup objects are stored:
      integer, dimension(:), allocatable :: kk     ! (maxRank+1,1)

      character*128 msgbuf
      integer :: k_tmp

      !-----------------------------
      !  check if id points to a profile object

      idf=idi

      call xplasma_errck1_prof(s,idf,ier)
      if(ier.ne.0) return

      !-----------------------------
      !  precheck if need to use associated profile

      call ck_coords(s,idx1,idx2,idx3,idf,istat)
      if(istat.gt.2) then
         iadr=s%dict(idf)%dindex
         idq=s%profs(iadr)%profIds(1)
         if(idq.gt.0) then
            call ck_coords(s,idx1,idx2,idx3,idq,istat)
            if(istat.le.2) then

               idf=idq  ! use assoc. profile

            endif
         endif
      endif

      !---------------------------
      !  check argument sizes

      if((size(x1s).ne.size(ans)).or.(size(x2s).ne.size(ans)).or. &
           (size(x3s).ne.size(ans))) then
         ier=510
         write(msgbuf,*) '  size of result vector: ',size(ans)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x1 vector: ',size(x1s)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x2 vector: ',size(x2s)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x3 vector: ',size(x3s)
         call xplasma_errmsg_append(s,msgbuf)
         return
      else
         ivec=size(x1s)
      endif

      iadr=s%dict(idf)%dindex
      irank=s%profs(iadr)%rank
      ispline=s%profs(iadr)%kspline
      call spline_type_decode(ispline,hspline,3)

      !---------------------------
      !  check profile rank

      if(irank.gt.3) then
         msgbuf = '  xplasma_eval_3dprof -- '//trim(s%dict(idf)%name)// &
              ' is not a 3d "f(x1,x2,x3)" profile; rank > 3'
         call xplasma_errmsg_append(s,msgbuf)
         ier=512
         return
      endif

      !---------------------------
      !  check derivative request consistent with interpolation setup for
      !  profile

      jderiv1=0
      if(present(ideriv1)) jderiv1=ideriv1
      jderiv2=0
      if(present(ideriv2)) jderiv2=ideriv2
      jderiv3=0
      if(present(ideriv3)) jderiv3=ideriv3

      if(ier.ne.0) return

      !---------------------------
      !  check for periodic coordinate; check for coordinate reversal e.g.
      !  due to clockwise theta angle interpretation in caller

      if((idx1.ge.1).and.(idx1.le.s%nitems)) then
         k_tmp = s%dict(idx1)%dtype
      else
         k_tmp = 0
      endif
      if((idx1.lt.1).or.(idx1.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x1 grid/coord id=',idx1, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq.xplasma_gridType) then
         iadx1=s%dict(idx1)%dindex
         icoord1=s%grids(iadx1)%coord     ! coordinate associated with grid x1
      else if( k_tmp .eq.xplasma_coordType) then
         icoord1=idx1
      else
         ier=105
         write(msgbuf,*) &
              '  x1 grid id points to "',trim(s%dict(idx1)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if((idx2.ge.1).and.(idx2.le.s%nitems)) then
         k_tmp = s%dict(idx2)%dtype
      else
         k_tmp = 0
      endif
      if((idx2.lt.1).or.(idx2.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x2 grid/coord id=',idx2, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq.xplasma_gridType) then
         iadx2=s%dict(idx2)%dindex
         icoord2=s%grids(iadx2)%coord     ! coordinate associated with grid x2
      else if( k_tmp .eq.xplasma_coordType) then
         icoord2=idx2
      else
         ier=105
         write(msgbuf,*) &
              '  x2 grid id points to "',trim(s%dict(idx2)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if((idx3.ge.1).and.(idx3.le.s%nitems)) then
         k_tmp = s%dict(idx3)%dtype
      else
         k_tmp = 0
      endif
      if((idx3.lt.1).or.(idx3.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x3 grid/coord id=',idx3, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq.xplasma_gridType) then
         iadx3=s%dict(idx3)%dindex
         icoord3=s%grids(iadx3)%coord     ! coordinate associated with grid x3
      else if( k_tmp .eq.xplasma_coordType) then
         icoord3=idx3
      else
         ier=105
         write(msgbuf,*) &
              '  x3 grid id points to "',trim(s%dict(idx3)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if(ier.ne.0) return

      if((icoord1.eq.icoord2).or.(icoord1.eq.icoord3)) then
         ier=515
         write(msgbuf,*) '  profile name is: ', trim(s%dict(idf)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  coordinate 1 is:   ', trim(s%dict(icoord1)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  coordinate 2 is:   ', trim(s%dict(icoord2)%name)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  coordinate 3 is:   ', trim(s%dict(icoord3)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      iadc1=s%dict(icoord1)%dindex
      perio1=s%coords(iadc1)%periodic ! periodicity attribute
      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1=ccwflag1

      iadc2=s%dict(icoord2)%dindex
      perio2=s%coords(iadc2)%periodic ! periodicity attribute
      mccwflag2=.TRUE.
      if(perio2.and.present(ccwflag2)) mccwflag2=ccwflag2

      iadc3=s%dict(icoord3)%dindex
      perio3=s%coords(iadc3)%periodic ! periodicity attribute
      mccwflag3=.TRUE.
      if(perio3.and.present(ccwflag3)) mccwflag3=ccwflag3

      !  check for sign change of odd derivative evaluation-- for each
      !  periodic coordinate reversal...

      isign = 1
      if(.not.mccwflag1) then
         if((jderiv1.eq.1).or.(jderiv1.eq.3)) then
            isign=-isign
         endif
      endif
      if(.not.mccwflag2) then
         if((jderiv2.eq.1).or.(jderiv2.eq.3)) then
            isign=-isign
         endif
      endif
      if(.not.mccwflag3) then
         if((jderiv3.eq.1).or.(jderiv3.eq.3)) then
            isign=-isign
         endif
      endif
      !---------------------------
      !  allocate lookup arrays and evaluate
      !    the function being evaluated may be rank 1 or 2
      !    the dependent coordinate values provided must supply enough
      !    information to evaluate the function, information beyond this
      !    is allowed but ignored.

      iforce=.FALSE.
      if(present(force_bounds)) iforce=force_bounds
      if(present(n_out_of_bounds)) n_out_of_bounds = 0

      inaxes = irank
      allocate(xx(inaxes))

      idx=s%profs(iadr)%gridIds(1)
      iadx1=s%dict(idx)%dindex     ! profile coordinate #1
      jcoord1=s%grids(iadx1)%coord
      
      xx(1)%inum = ivec
      xx(1)%inx = s%grids(iadx1)%size
      allocate(xx(1)%iderivs(1))

      !  verify that all x1s or x2s elements are in range; for periodic
      !  grids normalize the range and perform CW->CCW reversal if necessary

      if(jcoord1.eq.icoord1) then

         xx(1)%iderivs=jderiv1
         xx(1)%ccwflag=mccwflag1

      else if(jcoord1.eq.icoord2) then

         xx(1)%iderivs=jderiv2
         xx(1)%ccwflag=mccwflag2

      else if(jcoord1.eq.icoord3) then

         xx(1)%iderivs=jderiv3
         xx(1)%ccwflag=mccwflag3

      else

         ier=511
         write(msgbuf,*) ' xplasma_eval_prof: 1st coordinate of profile "', &
              trim(s%dict(idf)%name),'"'
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  is: "',trim(s%dict(jcoord1)%name),'" but there was'
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  no such data provided in the subroutine call.'
         call xplasma_errmsg_append(s,msgbuf)

      endif

      if(ier.eq.0) then
         call chkderiv(s,xx(1)%iderivs(1),hspline(1),idf,ier)
      endif

      if(irank.eq.1) then
         ! if derivative occurs along coordinate upon which fcn does not
         ! depend, the result is ZERO!

         if(jcoord1.eq.icoord1) then
            if(max(jderiv2,jderiv3).gt.0) isign=0
         else if(jcoord1.eq.icoord2) then
            if(max(jderiv1,jderiv3).gt.0) isign=0
         else
            if(max(jderiv1,jderiv2).gt.0) isign=0
         endif
         
      else if(irank.ge.2) then

         idx=s%profs(iadr)%gridIds(2)
         iadx2=s%dict(idx)%dindex     ! profile coordinate #2
         jcoord2=s%grids(iadx2)%coord
      
         xx(2)%inum = ivec
         xx(2)%inx = s%grids(iadx2)%size
         allocate(xx(2)%iderivs(1))

         if(jcoord2.eq.icoord1) then

            xx(2)%iderivs=jderiv1
            xx(2)%ccwflag=mccwflag1

         else if(jcoord2.eq.icoord2) then

            xx(2)%iderivs=jderiv2
            xx(2)%ccwflag=mccwflag2

         else if(jcoord2.eq.icoord3) then

            xx(2)%iderivs=jderiv3
            xx(2)%ccwflag=mccwflag3

         else

            ier=511
            write(msgbuf,*) &
                 ' xplasma_eval_prof: 2nd coordinate of profile "', &
                 trim(s%dict(idf)%name),'"'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) &
                 '  is: "',trim(s%dict(jcoord2)%name),'" but there was'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) &
                 '  no such data provided in the subroutine call.'
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(ier.eq.0) then
            call chkderiv(s,xx(2)%iderivs(1),hspline(2),idf,ier)
         endif

         if(irank.eq.2) then
            ! rank3 call evaluates a 2d profile; figure out which passed
            ! coordinate does not affect the evaluation; set derivatives
            ! to zero along this coordinate if they are specified.

            if((icoord1.ne.jcoord1).and.(icoord1.ne.jcoord2)) then
               if(jderiv1.gt.0) isign=0
            endif
            if((icoord2.ne.jcoord1).and.(icoord2.ne.jcoord2)) then
               if(jderiv2.gt.0) isign=0
            endif
            if((icoord3.ne.jcoord1).and.(icoord3.ne.jcoord2)) then
               if(jderiv3.gt.0) isign=0
            endif

         else
            ! full rank3 profile evaluation...

            idx=s%profs(iadr)%gridIds(3)
            iadx3=s%dict(idx)%dindex     ! profile coordinate #3
            jcoord3=s%grids(iadx3)%coord
      
            xx(3)%inum = ivec
            xx(3)%inx = s%grids(iadx3)%size
            allocate(xx(3)%iderivs(1))

            if(jcoord3.eq.icoord1) then

               xx(3)%iderivs=jderiv1
               xx(3)%ccwflag=mccwflag1

            else if(jcoord3.eq.icoord2) then

               xx(3)%iderivs=jderiv2
               xx(3)%ccwflag=mccwflag2

            else if(jcoord3.eq.icoord3) then

               xx(3)%iderivs=jderiv3
               xx(3)%ccwflag=mccwflag3

            else

               ier=511
               write(msgbuf,*) &
                    ' xplasma_eval_prof: 2nd coordinate of profile "', &
                    trim(s%dict(idf)%name),'"'
               call xplasma_errmsg_append(s,msgbuf)
               write(msgbuf,*) &
                    '  is: "',trim(s%dict(jcoord2)%name),'" but there was'
               call xplasma_errmsg_append(s,msgbuf)
               write(msgbuf,*) &
                    '  no such data provided in the subroutine call.'
               call xplasma_errmsg_append(s,msgbuf)
            endif
            
            if(ier.eq.0) then
               call chkderiv(s,xx(3)%iderivs(1),hspline(3),idf,ier)
            endif

         endif
      endif

      if(ier.gt.0) then
         do i=1,irank
            call xpeval_free(xx(i))
         enddo
         deallocate(xx)
         ans=0
         return
      endif

      !  OK -- perform lookup

      imode=2

      allocate(xx(1)%ix(ivec),xx(1)%dxn(ivec), &
           xx(1)%hx(ivec),xx(1)%hxi(ivec))
      if(jcoord1.eq.icoord1) then
         call xplookup(s,ivec,x1s,xx(1)%ccwflag, &
              xx(1)%inx,s%grids(iadx1)%xpkg,imode, &
              xx(1)%ix,xx(1)%dxn,xx(1)%hx,xx(1)%hxi,iwarn)
      else if(jcoord1.eq.icoord2) then
         call xplookup(s,ivec,x2s,xx(1)%ccwflag, &
              xx(1)%inx,s%grids(iadx1)%xpkg,imode, &
              xx(1)%ix,xx(1)%dxn,xx(1)%hx,xx(1)%hxi,iwarn)
      else
         call xplookup(s,ivec,x3s,xx(1)%ccwflag, &
              xx(1)%inx,s%grids(iadx1)%xpkg,imode, &
              xx(1)%ix,xx(1)%dxn,xx(1)%hx,xx(1)%hxi,iwarn)
      endif
      if(present(n_out_of_bounds)) n_out_of_bounds = &
           max(n_out_of_bounds,iwarn)
      if(.not.iforce) then
         if(iwarn.gt.0) then
            ier=514
         endif
      endif

      if(irank.ge.2) then

         allocate(xx(2)%ix(ivec),xx(2)%dxn(ivec), &
              xx(2)%hx(ivec),xx(2)%hxi(ivec))
         if(jcoord2.eq.icoord1) then
            call xplookup(s,ivec,x1s,xx(2)%ccwflag, &
                 xx(2)%inx,s%grids(iadx2)%xpkg,imode, &
                 xx(2)%ix,xx(2)%dxn,xx(2)%hx,xx(2)%hxi,iwarn)
         else if(jcoord2.eq.icoord2) then
            call xplookup(s,ivec,x2s,xx(2)%ccwflag, &
                 xx(2)%inx,s%grids(iadx2)%xpkg,imode, &
                 xx(2)%ix,xx(2)%dxn,xx(2)%hx,xx(2)%hxi,iwarn)
         else
            call xplookup(s,ivec,x3s,xx(2)%ccwflag, &
                 xx(2)%inx,s%grids(iadx2)%xpkg,imode, &
                 xx(2)%ix,xx(2)%dxn,xx(2)%hx,xx(2)%hxi,iwarn)
         endif
         if(present(n_out_of_bounds)) n_out_of_bounds = &
              max(n_out_of_bounds,iwarn)
         if(.not.iforce) then
            if(iwarn.gt.0) then
               ier=514
            endif
         endif

      endif

      if(irank.eq.3) then

         allocate(xx(3)%ix(ivec),xx(3)%dxn(ivec), &
              xx(3)%hx(ivec),xx(3)%hxi(ivec))
         if(jcoord3.eq.icoord1) then
            call xplookup(s,ivec,x1s,xx(3)%ccwflag, &
                 xx(3)%inx,s%grids(iadx3)%xpkg,imode, &
                 xx(3)%ix,xx(3)%dxn,xx(3)%hx,xx(3)%hxi,iwarn)
         else if(jcoord3.eq.icoord2) then
            call xplookup(s,ivec,x2s,xx(3)%ccwflag, &
                 xx(3)%inx,s%grids(iadx3)%xpkg,imode, &
                 xx(3)%ix,xx(3)%dxn,xx(3)%hx,xx(3)%hxi,iwarn)
         else
            call xplookup(s,ivec,x3s,xx(3)%ccwflag, &
                 xx(3)%inx,s%grids(iadx3)%xpkg,imode, &
                 xx(3)%ix,xx(3)%dxn,xx(3)%hx,xx(3)%hxi,iwarn)
         endif
         if(present(n_out_of_bounds)) n_out_of_bounds = &
              max(n_out_of_bounds,iwarn)
         if(.not.iforce) then
            if(iwarn.gt.0) then
               ier=514
            endif
         endif

      endif

      if(ier.ne.0) then
         do i=1,irank
            call xpeval_free(xx(i))
         enddo
         deallocate(xx)
         ans=0
         return
      endif

      allocate(kk(maxRank+1))

      kk(1)=iadr
      kk(2)=1  ! rank 1, only 1 x grid in single 1d profile eval.
      if(irank.ge.2) then
         kk(3)=2
      endif
      if(irank.eq.3) then
         kk(4)=3
      endif

      !---------------------------
      !  evaluate interpolations

      call xplasma_evala(s,ans,kk,xx)

      ans = ans*isign

      do i=1,irank
         call xpeval_free(xx(i))
      enddo
      deallocate(kk,xx)

    end subroutine xplasma_eval_3dprofxs

    !---------------------------------------------------------
    subroutine xplasma_eval_3dprofsxs(s,ids,idx1,x1s,idx2,x2s,idx3,x3s, &
         ans,ier, &
         ideriv1,ideriv1s,ideriv2,ideriv2s,ideriv3,ideriv3s, &
         ccwflag1,ccwflag2,ccwflag3,force_bounds,n_out_of_bounds)

      !  evaluate multiple profiles f(x1,x2) 2d at a vector of points

      type (xplasma), pointer :: s
      integer, dimension(:), intent(in) :: ids    ! profile ids
      integer :: idx1                             ! id of x1 grid or coord.
      real*8, dimension(:), intent(in) :: x1s     ! evaluation vector x1
      integer :: idx2                             ! id of x2 grid or coord.
      real*8, dimension(:), intent(in) :: x2s     ! evaluation vector x2
      integer :: idx3                             ! id of x3 grid or coord.
      real*8, dimension(:), intent(in) :: x3s     ! evaluation vector x3
      real*8, dimension(:,:), intent(out) :: ans  ! result of evaluations
      integer, intent(out) :: ier                 ! completion code 0=OK

      !  size(x1s)=size(x2s)=size(x3s)=size(ans,1) expected...

      !  ans(1:size(x1s),1) -- result of eval of profile ids(1)
      !  ans(1:size(x1s),2) -- result of eval of profile ids(2) ...etc...

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv1    ! derivative control
      !    1 for df/d[x1]; 2 for d2f/d[x1]2; 3 for d3f/d[x1]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv1s
      !  as ideriv1, but separately specified for each evaluation

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv2    ! derivative control
      !    1 for df/d[x2]; 2 for d2f/d[x2]2; 3 for d3f/d[x2]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv2s
      !  as ideriv2, but separately specified for each evaluation

      !  default: 0= fcn value
      integer, intent(in), optional :: ideriv3    ! derivative control
      !    1 for df/d[x3]; 2 for d2f/d[x3]2; 3 for d3f/d[x3]3
      !    2 & 3 available only for C2 spline fits
      integer, intent(in), dimension(:), optional :: ideriv3s
      !  as ideriv3, but separately specified for each evaluation

      !  default: T -- CCW orientation, no CW->CCW transform
      logical, intent(in), optional :: ccwflag1  ! for x1
      logical, intent(in), optional :: ccwflag2  ! for x2
      logical, intent(in), optional :: ccwflag3  ! for x3
      !    (these are ignored except for theta and phi angle coordinates

      !  default: F -- T to force all points in bounds with min & max
      logical, intent(in), optional :: force_bounds

      !  optional output: return no. of target points out of range
      integer, intent(out), optional :: n_out_of_bounds

      !-----------------------------
      integer :: jderiv1s(size(ids)),jderiv2s(size(ids)),jderiv3s(size(ids))
      logical :: mccwflag1,mccwflag2,mccwflag3,iforce

      integer :: idq,ids2(size(ids))
      integer :: idxs(3,size(ids))  ! list all x axes (with duplicates)
      !  if any rank 1 or 2 functions are evaluated with this routine, 0s are
      !  inserted...

      integer :: icoord1,icoord2,icoord3,jcoord1,jcoord2,jcoord3
      integer :: i,j,k,ivec
      logical :: perio1,perio2,perio3

      integer :: iadr,iadx1,iadx2,iadx3,iadc1,iadc2,iadc3
      integer :: ispline,irank,iwarn,imode,hspline(3),indc
      integer :: idxtmp(maxRank),iadx,iertmp

      integer :: ifound,istat
      integer :: inaxes  ! number of distinct x axes in evaluation group
      integer :: iaxes(3*size(ids)) ! distinct x axes (no duplicates)(1:inaxes)
      integer :: icoor(3*size(ids)) ! which coordinate vector 1 or 2
      integer :: indax(3,size(ids)) ! index into iaxes
      integer :: isigns(size(ids))  ! result sign or 0 factors for all evals

      !  here the zone lookup data will be gathered:
      type(xpeval), dimension(:), allocatable :: xx  ! (1:inaxes)

      !  here the rank and indices of zone lookup objects are stored:
      integer, dimension(:,:), allocatable :: kk  ! (maxRank+1,1)

      character*128 msgbuf
      integer :: k_tmp

      !-----------------------------

      ier=0

      !---------------------------
      !  check argument sizes
      
      if((size(x1s).ne.size(ans,1)).or.(size(x1s).ne.size(x2s)).or. &
           (size(x2s).ne.size(x3s)).or.(size(ids).ne.size(ans,2))) then
         ier=510
         write(msgbuf,*) '  size of result vector: ',size(ans,1)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x1 vector: ',size(x1s)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x2 vector: ',size(x2s)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  size of x3 vector: ',size(x3s)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) '  number of profil ids: ',size(ids), &
              '; number of result vectors: ',size(ans,2)
         call xplasma_errmsg_append(s,msgbuf)
      else
         ivec=size(x1s)
      endif
      if(present(ideriv1s)) then
         if(size(ideriv1s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv1s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv1)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv1 and ideriv1s both present.')
         endif
      endif
      if(present(ideriv2s)) then
         if(size(ideriv2s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv2s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv2)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv2 and ideriv2s both present.')
         endif
      endif
      if(present(ideriv3s)) then
         if(size(ideriv3s).ne.size(ids)) then
            ier=510
            write(msgbuf,*) '  size of profile id list: ',size(ids), &
                 '; size of derivative control vector: ',size(ideriv3s)
            call xplasma_errmsg_append(s,msgbuf)
         endif
         if(present(ideriv3)) then
            ier=516
            call xplasma_errmsg_append(s,' ideriv3 and ideriv3s both present.')
         endif
      endif
      if(ier.ne.0) return

      jderiv1s=0
      if(present(ideriv1)) jderiv1s=ideriv1
      if(present(ideriv1s)) jderiv1s=ideriv1s
      jderiv2s=0
      if(present(ideriv2)) jderiv2s=ideriv2
      if(present(ideriv2s)) jderiv2s=ideriv2s
      jderiv3s=0
      if(present(ideriv3)) jderiv3s=ideriv3
      if(present(ideriv3s)) jderiv3s=ideriv3s

      !  will check later if 2nd or 3rd derivatives can be evaluated for *all*
      !  profiles...

      if(ier.ne.0) return

      !---------------------------
      !  precheck for need to use associated profiles

      ids2 = ids
      do i=1,size(ids)
         if(s%dict(ids(i))%dtype.ne.xplasma_profType) cycle
         call ck_coords(s,idx1,idx2,idx3,ids(i),istat)
         if(istat.gt.2) then
            iadr=s%dict(ids(i))%dindex
            idq=s%profs(iadr)%profIds(1)
            if(idq.gt.0) then
               call ck_coords(s,idx1,idx2,idx3,idq,istat)
               if(istat.le.2) then

                  ids2(i) = idq  ! use assoc. profile

               endif
            endif
         endif
      enddo

      !---------------------------
      !  check for periodic coordinate; check for coordinate reversal e.g.
      !  due to clockwise theta angle interpretation in caller

      if((idx1.ge.1).and.(idx1.le.s%nitems)) then
         k_tmp = s%dict(idx1)%dtype
      else
         k_tmp = 0
      endif
      if((idx1.lt.1).or.(idx1.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x1 grid/coord id=',idx1, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq.xplasma_gridType) then
         iadx1=s%dict(idx1)%dindex
         icoord1=s%grids(iadx1)%coord     ! coordinate associated with grid x1
      else if( k_tmp .eq.xplasma_coordType) then
         icoord1=idx1
      else
         ier=105
         write(msgbuf,*) &
              '  x1 grid id points to "',trim(s%dict(idx1)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if((idx2.ge.1).and.(idx2.le.s%nitems)) then
         k_tmp = s%dict(idx2)%dtype
      else
         k_tmp = 0
      endif
      if((idx2.lt.1).or.(idx2.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x2 grid/coord id=',idx2, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq.xplasma_gridType) then
         iadx2=s%dict(idx2)%dindex
         icoord2=s%grids(iadx2)%coord     ! coordinate associated with grid x2
      else if( k_tmp .eq.xplasma_coordType) then
         icoord2=idx2
      else
         ier=105
         write(msgbuf,*) &
              '  x2 grid id points to "',trim(s%dict(idx2)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if((idx3.ge.1).and.(idx3.le.s%nitems)) then
         k_tmp = s%dict(idx3)%dtype
      else
         k_tmp = 0
      endif
      if((idx3.lt.1).or.(idx3.gt.s%nitems)) then
         ier=105
         write(msgbuf,*) '  xplasma_eval_prof -- x3 grid/coord id=',idx3, &
              ' invalid.'
         call xplasma_errmsg_append(s,msgbuf)
      else if( k_tmp .eq.xplasma_gridType) then
         iadx3=s%dict(idx3)%dindex
         icoord3=s%grids(iadx3)%coord     ! coordinate associated with grid x3
      else if( k_tmp .eq.xplasma_coordType) then
         icoord3=idx3
      else
         ier=105
         write(msgbuf,*) &
              '  x3 grid id points to "',trim(s%dict(idx3)%name), &
              '" which is not a grid or coordinate object.'
         call xplasma_errmsg_append(s,msgbuf)
      endif

      if(ier.ne.0) return

      if((icoord1.eq.icoord2).or.(icoord1.eq.icoord3)) then
         ier=515
         write(msgbuf,*) '  coordinate is:   ', trim(s%dict(icoord1)%name)
         call xplasma_errmsg_append(s,msgbuf)
         return
      endif

      iadc1=s%dict(icoord1)%dindex
      perio1=s%coords(iadc1)%periodic ! periodicity attribute
      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1=ccwflag1

      iadc2=s%dict(icoord2)%dindex
      perio2=s%coords(iadc2)%periodic ! periodicity attribute
      mccwflag2=.TRUE.
      if(perio2.and.present(ccwflag2)) mccwflag2=ccwflag2

      iadc3=s%dict(icoord3)%dindex
      perio3=s%coords(iadc3)%periodic ! periodicity attribute
      mccwflag3=.TRUE.
      if(perio3.and.present(ccwflag3)) mccwflag3=ccwflag3

      !  check for sign change of odd derivative evaluation-- for each
      !  periodic coordinate reversal...

      isigns = 1
      if(.not.mccwflag1) then
         do i=1,size(ids2)
            if((jderiv1s(i).eq.1).or.(jderiv1s(i).eq.3)) then
               isigns(i)=-isigns(i)
            endif
         enddo
      endif
      if(.not.mccwflag2) then
         do i=1,size(ids2)
            if((jderiv2s(i).eq.1).or.(jderiv2s(i).eq.3)) then
               isigns(i)=-isigns(i)
            endif
         enddo
      endif
      if(.not.mccwflag3) then
         do i=1,size(ids2)
            if((jderiv3s(i).eq.1).or.(jderiv3s(i).eq.3)) then
               isigns(i)=-isigns(i)
            endif
         enddo
      endif

      !---------------------------
      !  check all profiles -- all must be rank 1 2 or 3 and defined over
      !  the coordinate values provided in the subroutine arguments.
      !  Any requested derivative must be consistent
      !  with the interpolating function associated with each profile.

      inaxes=0

      do i=1,size(ids2)

         !--- error here if ids2(i) not pointing at a profile object

         call xplasma_prof_info(s,ids2(i),ier, &
              gridId1=idxtmp(1), gridId2=idxtmp(2), gridId3=idxtmp(3))
         if(ier.ne.0) exit

         if(maxRank.gt.3) then
            if(idxtmp(min(maxRank,4)).gt.0) then
               msgbuf = &
                    '  xplasma_eval_3dprof -- '//trim(s%dict(ids2(i))%name)// &
                    ' is not a 3d "f(x1,x2,x3)" profile; rank > 3'
               call xplasma_errmsg_append(s,msgbuf)
               ier=512
               exit
            endif
         endif

         !  OK...
         if(idxtmp(2).eq.0) then
            irank=1
         else if(idxtmp(3).eq.0) then
            irank=2
         else
            irank=3
         endif

         iadr=s%dict(ids2(i))%dindex

         iadx1=s%dict(idxtmp(1))%dindex
         jcoord1=s%grids(iadx1)%coord
         iertmp=0
         if((jcoord1.ne.icoord1).and.(jcoord1.ne.icoord2).and. &
              (jcoord1.ne.icoord3)) then
            iertmp=511
            write(msgbuf,*) '  xplasma_eval_prof: profile "', &
                 trim(s%dict(ids2(i))%name),'" depends on'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) &
                 '  coordinate "',trim(s%dict(jcoord1)%name),'" but'
            call xplasma_errmsg_append(s,msgbuf)
            write(msgbuf,*) &
                 '  no such data was provided for evaluation.'
            call xplasma_errmsg_append(s,msgbuf)
         endif

         if(irank.ge.2) then
            iadx2=s%dict(idxtmp(2))%dindex
            jcoord2=s%grids(iadx2)%coord
            if((jcoord2.ne.icoord1).and.(jcoord2.ne.icoord2).and. &
               (jcoord2.ne.icoord3)) then
               iertmp=511
               write(msgbuf,*) '  xplasma_eval_prof: profile "', &
                    trim(s%dict(ids2(i))%name),'" depends on'
               call xplasma_errmsg_append(s,msgbuf)
               write(msgbuf,*) &
                    '  coordinate "',trim(s%dict(jcoord2)%name),'" but'
               call xplasma_errmsg_append(s,msgbuf)
               write(msgbuf,*) &
                    '  no such data was provided for evaluation.'
               call xplasma_errmsg_append(s,msgbuf)
            endif
            if(irank.eq.3) then
               iadx3=s%dict(idxtmp(3))%dindex
               jcoord3=s%grids(iadx3)%coord
               if((jcoord3.ne.icoord1).and.(jcoord3.ne.icoord2).and. &
                  (jcoord3.ne.icoord3)) then
                  iertmp=511
                  write(msgbuf,*) '  xplasma_eval_prof: profile "', &
                       trim(s%dict(ids2(i))%name),'" depends on'
                  call xplasma_errmsg_append(s,msgbuf)
                  write(msgbuf,*) &
                       '  coordinate "',trim(s%dict(jcoord3)%name),'" but'
                  call xplasma_errmsg_append(s,msgbuf)
                  write(msgbuf,*) &
                       '  no such data was provided for evaluation.'
                  call xplasma_errmsg_append(s,msgbuf)
               endif
            else
               !  rank 2...
               iadx3=0
               jcoord3=0
               if((icoord1.ne.jcoord1).and.(icoord1.ne.jcoord2)) then
                  if(jderiv1s(i).gt.0) isigns(i)=0
               endif
               if((icoord2.ne.jcoord1).and.(icoord2.ne.jcoord2)) then
                  if(jderiv2s(i).gt.0) isigns(i)=0
               endif
               if((icoord3.ne.jcoord1).and.(icoord3.ne.jcoord2)) then
                  if(jderiv3s(i).gt.0) isigns(i)=0
               endif
            endif
         else
            !  rank 1...
            iadx3=0
            jcoord3=0
            iadx2=0
            jcoord2=0
            !  rank=1 -- if derivative is on other coordinate answer is zero...
            if(jcoord1.eq.icoord1) then
               if(max(jderiv2s(i),jderiv3s(i)).gt.0) isigns(i)=0
            else if(jcoord1.eq.icoord2) then
               if(max(jderiv1s(i),jderiv3s(i)).gt.0) isigns(i)=0
            else
               if(max(jderiv1s(i),jderiv2s(i)).gt.0) isigns(i)=0
            endif
         endif
         if(iertmp.ne.0) then
            ier=iertmp
            exit
         endif

         !--- check independent coordinate consistency in profile set

         idxs(1:3,i)=idxtmp(1:3)
         indax(1:3,i)=0

         if(inaxes.eq.0) then
            inaxes=irank
            iaxes(1)=idxtmp(1)
            indax(1,1)=1
            if(jcoord1.eq.icoord1) then
               icoor(1)=1
            else if(jcoord1.eq.icoord2) then
               icoor(1)=2
            else
               icoor(1)=3
            endif
            if(irank.ge.2) then
               iaxes(2)=idxtmp(2)
               indax(2,1)=2
               if(jcoord2.eq.icoord1) then
                  icoor(2)=1
               else if(jcoord2.eq.icoord2) then
                  icoor(2)=2
               else
                  icoor(2)=3
               endif
               if(irank.eq.3) then
                  iaxes(3)=idxtmp(3)
                  indax(3,1)=3
                  if(jcoord3.eq.icoord1) then
                     icoor(3)=1
                  else if(jcoord3.eq.icoord2) then
                     icoor(3)=2
                  else
                     icoor(3)=3
                  endif
               endif
            endif
         else
            ifound=0
            do j=1,inaxes
               if(idxs(1,i).eq.iaxes(j)) then
                  ifound=j
                  exit
               endif
            enddo
            if(ifound.eq.0) then
               inaxes=inaxes+1
               iaxes(inaxes)=idxs(1,i)
               indax(1,i)=inaxes
               if(jcoord1.eq.icoord1) then
                  icoor(inaxes)=1
               else if(jcoord1.eq.icoord2) then
                  icoor(inaxes)=2
               else
                  icoor(inaxes)=3
               endif
            else
               indax(1,i)=ifound
            endif
            if(irank.ge.2) then
               ifound=0
               do j=1,inaxes
                  if(idxs(2,i).eq.iaxes(j)) then
                     ifound=j
                     exit
                  endif
               enddo
               if(ifound.eq.0) then
                  inaxes=inaxes+1
                  iaxes(inaxes)=idxs(2,i)
                  indax(2,i)=inaxes
                  if(jcoord2.eq.icoord1) then
                     icoor(inaxes)=1
                  else if(jcoord2.eq.icoord2) then
                     icoor(inaxes)=2
                  else
                     icoor(inaxes)=3
                  endif
               else
                  indax(2,i)=ifound
               endif
               if(irank.eq.3) then
                  ifound=0
                  do j=1,inaxes
                     if(idxs(3,i).eq.iaxes(j)) then
                        ifound=j
                        exit
                     endif
                  enddo
                  if(ifound.eq.0) then
                     inaxes=inaxes+1
                     iaxes(inaxes)=idxs(3,i)
                     indax(3,i)=inaxes
                     if(jcoord3.eq.icoord1) then
                        icoor(inaxes)=1
                     else if(jcoord3.eq.icoord2) then
                        icoor(inaxes)=2
                     else
                        icoor(inaxes)=3
                     endif
                  else
                     indax(3,i)=ifound
                  endif
               endif
            endif
         endif
         if(ier.ne.0) exit

         !---------------------------
         !  check derivative request consistent with interpolation setup for
         !  all profiles

         ispline=s%profs(iadr)%kspline
         call spline_type_decode(ispline,hspline,3)

         indc=icoor(indax(1,i))
         if(indc.eq.1) then
            call chkderiv(s,jderiv1s(i),hspline(1),ids2(i),ier)
         else if(indc.eq.2) then
            call chkderiv(s,jderiv2s(i),hspline(1),ids2(i),ier)
         else
            call chkderiv(s,jderiv3s(i),hspline(1),ids2(i),ier)
         endif

         if(irank.ge.2) then
            iertmp=ier
            indc=icoor(indax(2,i))
            if(indc.eq.1) then
               call chkderiv(s,jderiv1s(i),hspline(2),ids2(i),ier)
            else if(indc.eq.2) then
               call chkderiv(s,jderiv2s(i),hspline(2),ids2(i),ier)
            else
               call chkderiv(s,jderiv3s(i),hspline(2),ids2(i),ier)
            endif
            ier=max(iertmp,ier)
         endif

         if(irank.eq.3) then
            iertmp=ier
            indc=icoor(indax(3,i))
            if(indc.eq.1) then
               call chkderiv(s,jderiv1s(i),hspline(3),ids2(i),ier)
            else if(indc.eq.2) then
               call chkderiv(s,jderiv2s(i),hspline(3),ids2(i),ier)
            else
               call chkderiv(s,jderiv3s(i),hspline(3),ids2(i),ier)
            endif
            ier=max(iertmp,ier)
         endif

      enddo
      if(ier.ne.0) return

      !---------------------------
      !  check for periodic coordinate; check for coordinate reversal e.g.
      !  due to clockwise theta angle interpretation in caller

      call xplasma_coord_isPeriodic(s,icoord1,perio1,ier)
      mccwflag1=.TRUE.
      if(perio1.and.present(ccwflag1)) mccwflag1=ccwflag1

      call xplasma_coord_isPeriodic(s,icoord2,perio2,ier)
      mccwflag2=.TRUE.
      if(perio2.and.present(ccwflag2)) mccwflag2=ccwflag2

      call xplasma_coord_isPeriodic(s,icoord3,perio3,ier)
      mccwflag3=.TRUE.
      if(perio3.and.present(ccwflag3)) mccwflag3=ccwflag3

      !---------------------------
      !  allocate lookup arrays and evaluate

      if(present(n_out_of_bounds)) n_out_of_bounds=0

      allocate(xx(inaxes))
      do i=1,inaxes
         xx(i)%inum = ivec

         iadx = s%dict(iaxes(i))%dindex
         xx(i)%inx = s%grids(iadx)%size

         allocate(xx(i)%iderivs(size(ids2)))

         if(icoor(i).eq.1) then
            xx(i)%iderivs = jderiv1s
            xx(i)%ccwflag = mccwflag1
         else if(icoor(i).eq.2) then
            xx(i)%iderivs = jderiv2s
            xx(i)%ccwflag = mccwflag2
         else
            xx(i)%iderivs = jderiv3s
            xx(i)%ccwflag = mccwflag3
         endif

         !  verify that all xs elements are in range; for periodic
         !  grids normalize the range and perform CW->CCW reversal if necessary

         iforce=.FALSE.
         if(present(force_bounds)) iforce=force_bounds

         !  OK -- perform lookup

         allocate(xx(i)%ix(ivec),xx(i)%dxn(ivec), &
              xx(i)%hx(ivec),xx(i)%hxi(ivec))

         imode=2
         if(icoor(i).eq.1) then
            call xplookup(s,ivec,x1s,xx(i)%ccwflag, &
                 s%grids(iadx)%size,s%grids(iadx)%xpkg,imode, &
                 xx(i)%ix,xx(i)%dxn,xx(i)%hx,xx(i)%hxi,iwarn)
         else if(icoor(i).eq.2) then
            call xplookup(s,ivec,x2s,xx(i)%ccwflag, &
                 s%grids(iadx)%size,s%grids(iadx)%xpkg,imode, &
                 xx(i)%ix,xx(i)%dxn,xx(i)%hx,xx(i)%hxi,iwarn)
         else
            call xplookup(s,ivec,x3s,xx(i)%ccwflag, &
                 s%grids(iadx)%size,s%grids(iadx)%xpkg,imode, &
                 xx(i)%ix,xx(i)%dxn,xx(i)%hx,xx(i)%hxi,iwarn)
         endif

         if(present(n_out_of_bounds)) n_out_of_bounds = &
              max(n_out_of_bounds,iwarn)
         if(.not.iforce) then
            if(iwarn.gt.0) then
               ier=514
               exit
            endif
         endif

      enddo

      if(ier.ne.0) then
         do j=1,i
            call xpeval_free(xx(j))
         enddo
         deallocate(xx)
         ans=0
         return
      endif

      !  set up eval table for each profile

      allocate(kk(maxRank+1,size(ids2)))
      do j=1,size(ids2)
         iadr=s%dict(ids2(j))%dindex
         kk(1,j)=iadr
         kk(2:4,j)=indax(1:3,j)
      enddo

      !---------------------------
      !  evaluate interpolations

      call xplasma_eval(s,ans,kk,xx)

      do j=1,size(ids2)
         ans(1:ivec,j)=ans(1:ivec,j)*isigns(j)
      enddo
      
      deallocate(kk)
      do j=1,inaxes
         call xpeval_free(xx(j))
      enddo
      deallocate(xx)
   
    end subroutine xplasma_eval_3dprofsxs

    !-------------------------------------------
    !  private coord ck routine

    subroutine ck_coords(s,idx1,idx2,idx3,idf,istat)

      ! return status of set of coordinates w.r.t. a profile's set of
      ! coordinates:
      !   istat=1 means a match;
      !   istat=2 means idf set is a subset of {idx1,idx2,idx3}
      !   istat=3 means {idx1,idx2,idx3} is a subset of the idf set
      !   istat=4 means neither set contains the other

      ! the idx* could be grid Ids or coordinate Ids -- have to check
      ! some of the idx* may be zero; these are ignored.

      ! no assumption about the ordering of the set elements

      type (xplasma), pointer :: s

      integer, intent(in) :: idx1,idx2,idx3
      integer, intent(in) :: idf

      integer, intent(out) :: istat

      !-------------------
      !  local:
      integer :: icx(3),icf(3),inx,inf,innx,innf,iertmp,iadr,idx,i,j,iadx
      !-------------------

      icx=0; inx=0
      icf=0; inf=0

      if(idx1.ne.0) then
         inx=inx+1
         iadr=s%dict(idx1)%dindex
         if(s%dict(idx1)%dtype.eq.xplasma_gridType) then
            icx(inx)=s%grids(iadr)%coord
         else
            icx(inx)=idx1
         endif
      endif

      if(idx2.ne.0) then
         inx=inx+1
         iadr=s%dict(idx2)%dindex
         if(s%dict(idx2)%dtype.eq.xplasma_gridType) then
            icx(inx)=s%grids(iadr)%coord
         else
            icx(inx)=idx2
         endif
      endif

      if(idx3.ne.0) then
         inx=inx+1
         iadr=s%dict(idx3)%dindex
         if(s%dict(idx3)%dtype.eq.xplasma_gridType) then
            icx(inx)=s%grids(iadr)%coord
         else
            icx(inx)=idx3
         endif
      endif

      if(inx.eq.0) then
         istat=3  ! null set always a subset
         return
      endif

      iadr = s%dict(idf)%dindex

      do i=1,3
         idx=s%profs(iadr)%gridIds(i)
         if(idx.eq.0) exit

         inf=inf+1
         iadx=s%dict(idx)%dindex
         icf(inf)=s%grids(iadx)%coord
      enddo

      !  have the sets of coordinates; see how many of each set is in the other

      innx=0  ! # of x set in f set
      innf=0  ! # of f set in x set

      do i=1,inx
         do j=1,inf
            if(icx(i).eq.icf(j)) then
               innx=innx+1
               exit
            endif
         enddo
      enddo

      do i=1,inf
         do j=1,inx
            if(icf(i).eq.icx(j)) then
               innf=innf+1
               exit
            endif
         enddo
      enddo

      if((inx.eq.inf).and.(inx.eq.innx).and.(inf.eq.innf)) then
         istat=1  ! match

      else if((inf.lt.inx).and.(inf.eq.innf)) then
         istat=2  ! f set contained in x set

      else if((inx.lt.inf).and.(inx.eq.innx)) then
         istat=3  ! x set contained in f set

      else
         istat=4  ! neither is contained
      endif

    end subroutine ck_coords

    !-------------------------------------------

    subroutine xplasma_create_blackBox(s,iname,itype,id,ier, &
         iarray,r8array,iasize,r8asize,label,units, &
         ia_ptr,r8a_ptr)

      !  create a blackBox -- a structure which contains the following
      !    -- a scalar integer TYPE (user defined meaning)
      !    -- an integer array
      !    -- a real*8 array

      !  Xplasma provides a mechanism for storage and retrieval of
      !  Black boxes; they are included in Xplasma NetCDF state files.
      !  The only methods supported are:
      !    xplasma_create_blackBox -- create the item; store the data
      !    xplasma_blackBox_info   -- retrieve info on item (e.g. array sizes)
      !    xplasma_blackBox_retrieve -- retrieve black box item data

      !  Xplasma makes internal use of blackboxes, e.g. to define limiters
      !  or set up data for efficient numerical integration for flux surface
      !  averages.  But, there are no reserved TYPE codes-- the significance
      !  of the blackBox TYPE, like the data itself, is entirely up to the
      !  user.

      type (xplasma), pointer :: s
      character*(*), intent(in) :: iname  ! name of black box item
      integer, intent(in) ::  itype       ! Black Box type code (user defined).

      integer, intent(out) :: id          ! id of newly created black box item
      integer, intent(out) :: ier         ! completion code

      !  independently for the integer and floating point arrays,
      !  --either data must be provided, or, data sizes must be provided.

      integer, dimension(:), intent(in), optional :: iarray   ! integer data
      real*8, dimension(:), intent(in), optional :: r8array   ! real*8 data

      !  --but not both!
      integer, intent(in), optional :: iasize   ! size of integer data
      integer, intent(in), optional :: r8asize  ! size of real*8 data

      !  optional labels...
      character*(*), intent(in), optional :: label,units

      !  optional output: pointers to blackbox arrays (null if ierr.ne.0)

      integer, dimension(:), pointer, optional :: ia_ptr
      real*8, dimension(:), pointer, optional :: r8a_ptr

      !---------------------------------
      integer :: isize,r8size,iadr,iertmp,ipres
      !---------------------------------

      if(present(ia_ptr)) nullify(ia_ptr)
      if(present(r8a_ptr)) nullify(r8a_ptr)

      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_create_blackBox call)',ier)
      endif
      if(ier.ne.0) return

      ipres=0
      if(present(iarray)) ipres=ipres+1
      if(present(iasize)) ipres=ipres+1
      if(ipres.ne.1) then
         ier=616
         if(ipres.eq.0) then
            call xplasma_errmsg_append(s, &
                 ' Either an integer array (iarray) or an integer array size (iasize) must be provided.')
         else if(ipres.eq.2) then
            call xplasma_errmsg_append(s, &
                 ' Provide either an integer array (iarray) or integer array size (iasize) but not both.')
         endif
      endif

      ipres=0
      if(present(r8array)) ipres=ipres+1
      if(present(r8asize)) ipres=ipres+1
      if(ipres.ne.1) then
         ier=616
         if(ipres.eq.0) then
            call xplasma_errmsg_append(s, &
                 ' Either a real*8 array (r8array) or a real*8 array size (r8asize) must be provided.')
         else if(ipres.eq.2) then
            call xplasma_errmsg_append(s, &
                 ' Provide either a real*8 array (r8array) or real*8 array size (r8asize) but not both.')
         endif
      endif

      if(ier.ne.0) return

      id=0

      !  make the "black box" object...
      call xplasma_add_item(s,iname,xplasma_blkbxType,id,ier)
      if(ier.ne.0) return
      
      iadr=s%dict(id)%dindex

      !  store the data

      s%blkbxs(iadr)%type = itype

      if(present(iasize)) then
         isize=iasize
      else
         isize = size(iarray)
      endif

      if(present(r8asize)) then
         r8size=r8asize
      else
         r8size = size(r8array)
      endif

      allocate(s%blkbxs(iadr)%iwork(isize))
      if(present(iarray)) then
         s%blkbxs(iadr)%iwork = iarray
      else
         s%blkbxs(iadr)%iwork = 0
      endif

      allocate(s%blkbxs(iadr)%r8work(r8size))
      if(present(r8array)) then
         s%blkbxs(iadr)%r8work = r8array
      else
         s%blkbxs(iadr)%r8work = czero
      endif

      iertmp=0

      if(present(label)) then
         call xplasma_label_item(s,id,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_label_item(s,id,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      ier=iertmp

      if(ier.eq.0) then
         if(present(ia_ptr)) then
            ia_ptr => s%blkbxs(iadr)%iwork
         endif

         if(present(r8a_ptr)) then
            r8a_ptr => s%blkbxs(iadr)%r8work
         endif
      endif

    end subroutine xplasma_create_blackBox

    subroutine xplasma_blackBox_retrieve(s,id,ier, &
         itype,iarray,r8array,ia_ptr,r8a_ptr,name,label,units,author)

      !  retrieve blackBox data
      !    optional arguments; user controls what data is wanted.

      type (xplasma), pointer :: s
      integer, intent(in) :: id          ! id of black box item
      integer, intent(out) :: ier        ! completion code, 0=OK

      integer, intent(out), optional :: itype  ! "type" of black box
      !  meaning is user defined.

      !  if copies of arrays are wanted:
      integer, dimension(:), intent(out), optional :: iarray
      real*8, dimension(:), intent(out), optional :: r8array

      !  if pointers to arrays are wanted:
      integer, dimension(:), pointer, optional :: ia_ptr
      real*8, dimension(:), pointer, optional :: r8a_ptr

      !  optional labeling info...
      character*(*), intent(out), optional :: name     ! item name
      character*(*), intent(out), optional :: label    ! item label
      character*(*), intent(out), optional :: units    ! item units label
      character*(*), intent(out), optional :: author   ! author

      !------------------------------
      integer :: iadr,isize,r8size
      character*128 zmsg
      !------------------------------

      ier = 0
      if(present(itype)) itype=0

      call xplasma_errck1_blkbx(s,id,ier)
      if(ier.ne.0) return

      if(present(name)) call xplasma_get_item_info(s,id,ier,name=name)
      if(present(label)) call xplasma_get_item_info(s,id,ier,label=label)
      if(present(units)) call xplasma_get_item_info(s,id,ier,units=units)
      if(present(author)) call xplasma_get_item_info(s,id,ier,author=author)

      iadr = s%dict(id)%dindex

      if(present(itype)) itype = s%blkbxs(iadr)%type

      isize = size(s%blkbxs(iadr)%iwork)
      r8size = size(s%blkbxs(iadr)%r8work)

      if(present(iarray)) then
         if(size(iarray).ne.isize) then
            ier=612
            zmsg=' '
            write(zmsg,*) '  passed integer array (iarray) length: ',size(iarray)
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) '  length of black box integer array data: ',isize
            call xplasma_errmsg_append(s,zmsg)
            iarray = 0
         else
            iarray = s%blkbxs(iadr)%iwork
         endif
      endif

      if(present(ia_ptr)) then
         ia_ptr => s%blkbxs(iadr)%iwork
      endif

      if(present(r8array)) then
         if(size(r8array).ne.r8size) then
            ier=612
            zmsg=' '
            write(zmsg,*) '  passed real*8 array (r8array) length: ',size(r8array)
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) '  length of black box real*8 array data: ',r8size
            call xplasma_errmsg_append(s,zmsg)
            r8array = 0
         else
            r8array = s%blkbxs(iadr)%r8work
         endif
      endif

      if(present(r8a_ptr)) then
         r8a_ptr => s%blkbxs(iadr)%r8work
      endif

    end subroutine xplasma_blackBox_retrieve

    !========================================
    subroutine xplasma_blackBox_info(s,id,ier, &
         type,isize,r8size,name,label,units,author)
      !
      !  get the size of an existing list
      !
      type (xplasma), pointer :: s
      integer, intent(in) :: id           ! id of black box item
      integer, intent(out) :: ier         ! completion code

      integer, intent(out), optional :: type  ! "type" of black box
      ! generally, meaning of iitype is user defined...

      integer, intent(out), optional :: isize  ! number of integer words
      integer, intent(out), optional :: r8size ! number of floating pt words

      character*(*), intent(out), optional :: name
      character*(*), intent(out), optional :: label
      character*(*), intent(out), optional :: units
      character*(*), intent(out), optional :: author

      !---------------------
      integer :: iadr,iertmp
      !---------------------

      if(present(type)) type=0
      if(present(isize)) isize=0
      if(present(r8size)) r8size=0

      if(present(name)) name=' '
      if(present(label)) label=' '
      if(present(units)) units=' '
      if(present(author)) author=' '

      call xplasma_errck1_blkbx(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex

      if(present(type)) type = s%blkbxs(iadr)%type

      if(present(isize)) then
         if(allocated(s%blkbxs(iadr)%iwork)) then
            isize = size(s%blkbxs(iadr)%iwork)
         endif
      endif

      if(present(r8size)) then
         if(allocated(s%blkbxs(iadr)%r8work)) then
            r8size = size(s%blkbxs(iadr)%r8work)
         endif
      endif

      iertmp=0

      if(present(name)) then
         call xplasma_get_item_info(s,id,ier, name=name)
         iertmp=max(iertmp,ier)
      endif

      if(present(label)) then
         call xplasma_get_item_info(s,id,ier, label=label)
         iertmp=max(iertmp,ier)
      endif

      if(present(units)) then
         call xplasma_get_item_info(s,id,ier, units=units)
         iertmp=max(iertmp,ier)
      endif

      if(present(author)) then
         call xplasma_get_item_info(s,id,ier, author=author)
         iertmp=max(iertmp,ier)
      endif

      ier=iertmp

    end subroutine xplasma_blackBox_info

    !==========================================================
    subroutine xpmom_init(s,ier, rhoi,ikmom,ikmaxe)

      !  init Fourier Representation structure

      type (xplasma), pointer :: s
      integer, intent(out) :: ier             ! status code,  0=normal
 
      real*8,dimension(:),intent(in), optional :: rhoi ! rho values axis to bdy
      ! if omitted use "__RHO" grid

      integer, intent(in), optional :: ikmom  ! #moments can be set here
      !                    **allowed range of values** [2:64]
      !                    5 or 8 usually OK for conventional tokamaks
      !                    16 usually OK for STs

      integer, intent(in), optional :: ikmaxe ! max exponent n for /x**n
      ! scaling of moments in Fourier Spline (default: 1); range (1:2)

      !------------------------     
      type (xmom2d), pointer :: m
 
      integer :: kmom,ksurf,inx,ix,id_rho,iertmp,ikmax

      real*8, dimension(:), allocatable :: zx  ! radial grid

      real*8, parameter :: z0 = 0.0d0
      real*8, parameter :: z1 = 1.0d0

      !-----------------

      ikmax=1
      if(present(ikmaxe)) ikmax=min(2,max(1,ikmaxe))

      call xplasma_global_info(s,ier, kmom=kmom) ! # of moments
      if(ier.ne.0) return  ! this just checks that "s" was initialized...
 
      if(present(ikmom)) then
         s%kmom=max(2,min(64,ikmom))
         kmom = s%kmom
      endif

      if(present(rhoi)) then
         inx=size(rhoi)
         allocate(zx(inx))
         zx=rhoi

      else
         call xplasma_find_item(s,"__RHO",id_rho,ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xpmom_init: no grid argument AND no __RHO grid in xplasma object.')
            return
         endif
         call xplasma_grid_size(s,id_rho,inx,iertmp)
         allocate(zx(inx))
         call xplasma_grid(s,id_rho,zx,iertmp)

      endif

      if(abs(zx(1)).gt.1.0d-10) then
         call xplasma_errmsg_append(s, &
              ' ?xpmom_init: Fourier spline x grid must start at 0')
         ier=306
      endif

      if(abs(zx(inx)-z1).gt.1.0d-10) then
         call xplasma_errmsg_append(s, &
              ' ?xpmom_init: Fourier spline x grid must end at 1')
         ier=306
      endif

      do ix=2,inx
         if(zx(ix).le.zx(ix-1)) then
            ier=303
            call xplasma_errmsg_append(s, &
                 ' ?xpmom_init: passed Fourier spline x grid not strict ascending.')
            exit
         endif
      enddo
      if(ier.ne.0) return
      
      m => s%xm
      if(allocated(m%xrho)) call xpmom_free(s)

      m%kmom=kmom
      m%ksurf=inx
      m%kmaxe=ikmax

      ksurf=inx

      call xpmom_alloc(m)
 
      call r8genxpkg(ksurf,zx,m%xrho,0,0,0,z0,3,ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xmoments.f90: xmoments_init: unexpected error in genxpkg!')
         ier=9999
         return
      endif

      call xpmom_create_store(s)
 
    end subroutine xpmom_init

    subroutine xpmom_alloc(m)

      type(xmom2d), pointer :: m

      real*8, parameter :: z0 = 0.0d0
      real*8, parameter :: zn1 = -1.0d0

      integer :: kmom
      integer :: ksurf

      kmom = m%kmom
      ksurf = m%ksurf

      allocate(m%xrho(ksurf,4)); m%xrho=z0

      allocate(m%sxrmmc(4,ksurf,0:kmom),m%sxrmms(4,ksurf,0:kmom))
      m%sxrmmc=z0; m%sxrmms=z0
 
      allocate(m%sxzmmc(4,ksurf,0:kmom),m%sxzmms(4,ksurf,0:kmom))
      m%sxzmmc=z0; m%sxzmms=z0
 
    end subroutine xpmom_alloc

    subroutine xpmom_free(s)

      !  free (deallocate) Fourier Representation structure

      type (xplasma), pointer :: s
      integer :: id_xmom2d,ier

      type (xmom2d), pointer :: m
      
      m => s%xm

      call xpmom_ck_dealloc(m)

      call xplasma_find_item(s,'__XMOM2D',id_xmom2d,ier,nf_noerr=.TRUE.)
      if(id_xmom2d.gt.0) then
         call xplasma_remove_item(s,id_xmom2d,ier)
      endif

    end subroutine xpmom_free

    subroutine xpmom_ck_dealloc(m)

      type(xmom2d), pointer :: m

      if(allocated(m%xrho)) then
         m%ksurf=0
         m%kmom=0
         m%kmaxe=1

         deallocate(m%xrho)
         deallocate(m%sxrmmc,m%sxrmms,m%sxzmmc,m%sxzmms)
      endif

    end subroutine xpmom_ck_dealloc

    subroutine xpmom_create_store(s)

      !  save xmom2d object to black box storage

      type (xplasma), pointer :: s
      integer :: id_xmom2d,ier,isize,idata(12),ii
      real*8, dimension(:), pointer :: r8buf

      type (xmom2d), pointer :: m

      m => s%xm

      if(m%kmom.eq.0) return
      if(m%ksurf.eq.0) return

      idata(1)=m%kmom
      idata(2)=m%ksurf

      idata(3)=4*m%kmom+1
      do ii=4,7
         idata(ii)=idata(ii-1)+4*(m%kmom+1)*m%ksurf
      enddo
      do ii=8,10
         idata(ii)=idata(ii-1)+(m%kmom+1)
      enddo

      isize=idata(10)+m%kmom

      idata(11)=m%kmaxe
      idata(12)=0

      call xplasma_author_set(s,xplasma_xmhd,ier)

      call xplasma_create_blackBox(s,'__XMOM2D',66,id_xmom2d,ier, &
           iarray=idata, r8asize=isize, &
           label='equilibrium Fourier representation object', units='mixed')

      call xplasma_author_clear(s,xplasma_xmhd,ier)

      !  Black box exists; now use it...

      call xpmom_store(s)

      nullify(r8buf)

    end subroutine xpmom_create_store

    subroutine xpmom_store(s)

      type (xplasma), pointer :: s
      integer :: id_xmom2d,ier
      integer, dimension(:), pointer :: idata
      real*8, dimension(:), pointer :: r8buf

      type (xmom2d), pointer :: m

      m => s%xm

      if(m%kmom.eq.0) return
      if(m%ksurf.eq.0) return

      call xplasma_find_item(s,'__XMOM2D',id_xmom2d,ier,nf_noerr=.TRUE.)
      if(ier.ne.0) return
      if(id_xmom2d.eq.0) return

      call xplasma_blackBox_retrieve(s,id_xmom2d,ier, &
           ia_ptr=idata,r8a_ptr=r8buf)

      call xpblind_copy(4*m%ksurf,m%xrho,r8buf(1))

      call xpblind_copy(4*(m%kmom+1)*m%ksurf,m%sxrmmc,r8buf(idata(3)))
      call xpblind_copy(4*(m%kmom+1)*m%ksurf,m%sxrmms,r8buf(idata(4)))
      call xpblind_copy(4*(m%kmom+1)*m%ksurf,m%sxzmmc,r8buf(idata(5)))
      call xpblind_copy(4*(m%kmom+1)*m%ksurf,m%sxzmms,r8buf(idata(6)))

      r8buf(idata(7):idata(10)+m%kmom)=0

    end subroutine xpmom_store

    subroutine xpmom_fetch(s)

      type (xplasma), pointer :: s
      integer :: id_xmom2d,ier
      integer, dimension(:), pointer :: idata
      real*8, dimension(:), pointer :: r8buf

      type (xmom2d), pointer :: m

      call xplasma_find_item(s,'__XMOM2D',id_xmom2d,ier,nf_noerr=.TRUE.)
      if(ier.ne.0) return
      if(id_xmom2d.eq.0) return

      m => s%xm

      call xpmom_ck_dealloc(m)

      call xplasma_blackBox_retrieve(s,id_xmom2d,ier, &
           ia_ptr=idata,r8a_ptr=r8buf)

      s%kmom = idata(1)
      m%kmom = idata(1)
      m%ksurf= idata(2)

      if(size(idata).le.10) then
         m%kmaxe=idata(1)  ! in old datasets kmaxe=kmom
      else
         m%kmaxe=idata(11)
      endif

      call xpmom_alloc(m)

      call xpblind_copy(4*m%ksurf,r8buf(1),m%xrho)

      call xpblind_copy(4*(m%kmom+1)*m%ksurf,r8buf(idata(3)),m%sxrmmc)
      call xpblind_copy(4*(m%kmom+1)*m%ksurf,r8buf(idata(4)),m%sxrmms)
      call xpblind_copy(4*(m%kmom+1)*m%ksurf,r8buf(idata(5)),m%sxzmmc)
      call xpblind_copy(4*(m%kmom+1)*m%ksurf,r8buf(idata(6)),m%sxzmms)

    end subroutine xpmom_fetch

    subroutine xpmom_setup(s,rmmc,rmms,zmmc,zmms,ier)

      ! given the moments, set up Fourier Spline representation.
      ! Prior call to xpmom_init required to have defined x grid.

      type (xplasma), pointer :: s

      real*8, intent(in), dimension(0:,:) :: rmmc,rmms,zmmc,zmms
      !  moments coefficents (ksurf,0:kmom) Rcos Rsin Zcos Zsin

      integer, intent(out) :: ier   ! completion code, 0=OK

      !-------------------
      integer :: id_xmom2d,iertmp,ksurf,kmom,jsmoo,j,im,idum,ikmax,imexp
      integer, dimension(:), allocatable :: jmin
      real*8 :: zxbrk,zrtol  ! local copies
      type (xmom2d), pointer :: m
      real*8, parameter :: z0 = 0.0d0
      real*8, parameter :: ZERO=0.0d0
      real*8, parameter :: ZTENTH = 0.1d0
      real*8, parameter :: ZEPS3 = 1.0d-3
      real*8, parameter :: ZEPS4 = 1.0d-4
      real*8, parameter :: ZEPS5 = 1.0d-5
      real*8, parameter :: ZEPS6 = 1.0d-6
      !-----------------------------

      call xplasma_find_item(s,'__XMOM2D',id_xmom2d,ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' xpmom_setup: prior call to xpmom_init required to define grid')
         return
      endif

      !  __XMOM2D exists; assume that xplasma Fourier object has also been
      !  initialized (since __XMOM2D and this are set up together).

      m => s%xm

      ksurf=m%ksurf
      kmom=m%kmom
      ikmax=m%kmaxe

      ! check argument dimensions

      call ckm('rmmc',rmmc)
      call ckm('rmms',rmms)
      call ckm('zmmc',zmmc)
      call ckm('zmms',zmms)
      if(ier.ne.0) return

      ! 0'th moments
      m%sxrmmc(1,1:ksurf,0)=rmmc(0,1:ksurf)
      m%sxzmmc(1,1:ksurf,0)=zmmc(0,1:ksurf)
 
      ! higher moments:
 
      do im=1,kmom
         imexp=min(ikmax,im)
         m%sxrmmc(1,2:ksurf,im)=rmmc(im,2:ksurf)/m%xrho(2:ksurf,1)**imexp
         m%sxrmms(1,2:ksurf,im)=rmms(im,2:ksurf)/m%xrho(2:ksurf,1)**imexp
         m%sxzmmc(1,2:ksurf,im)=zmmc(im,2:ksurf)/m%xrho(2:ksurf,1)**imexp
         m%sxzmms(1,2:ksurf,im)=zmms(im,2:ksurf)/m%xrho(2:ksurf,1)**imexp
      enddo
 
      ! smoothing complete.
      do im=1,kmom
         m%sxrmmc(1,1,im)=(4*m%sxrmmc(1,2,im)-m%sxrmmc(1,3,im))/3  ! axial extrap.
         m%sxrmms(1,1,im)=(4*m%sxrmms(1,2,im)-m%sxrmms(1,3,im))/3  ! axial extrap.
         m%sxzmmc(1,1,im)=(4*m%sxzmmc(1,2,im)-m%sxzmmc(1,3,im))/3  ! axial extrap.
         m%sxzmms(1,1,im)=(4*m%sxzmms(1,2,im)-m%sxzmms(1,3,im))/3  ! axial extrap.
      enddo
 
      ! compute spline coefficients
 
      ier=0
      do im=0,kmom
         call r8cspline(m%xrho,ksurf,m%sxrmmc(1,1,im),0,z0,0,z0,m%xrho,ksurf, &
              idum,iertmp)
         if(iertmp.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xpmom_init: unexpected cspline error!')
            ier=9999
         endif

         call r8cspline(m%xrho,ksurf,m%sxrmms(1,1,im),0,z0,0,z0,m%xrho,ksurf, &
              idum,iertmp)
         if(iertmp.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xpmom_init: unexpected cspline error!')
            ier=9999
         endif

         call r8cspline(m%xrho,ksurf,m%sxzmmc(1,1,im),0,z0,0,z0,m%xrho,ksurf, &
              idum,iertmp)
         if(iertmp.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xpmom_init: unexpected cspline error!')
            ier=9999
         endif

         call r8cspline(m%xrho,ksurf,m%sxzmms(1,1,im),0,z0,0,z0,m%xrho,ksurf, &
              idum,iertmp)
         if(iertmp.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xpmom_init: unexpected cspline error!')
            ier=9999
         endif

      enddo

      call xpmom_store(s)

      contains

        subroutine ckm(str,arr)
          character*(*), intent(in) :: str   ! moments array label
          real*8, dimension(:,:) :: arr

          if(size(arr,1).ne.kmom+1) then
             ier=9999
             call xplasma_errmsg_append(s,' ?xpmom_setup: '// &
                  str//' 1st dimension size should be #moments + 1')
          endif

          if(size(arr,2).ne.ksurf) then
             ier=9999
             call xplasma_errmsg_append(s,' ?xpmom_setup: '// &
                  str//' 2nd dimension size should match radial grid size.')
          endif
        end subroutine ckm

    end subroutine xpmom_setup

    subroutine xpmom_protect(s,ier)

      ! protect the Fourier moments representation from an equilibrium 
      ! update.  This is useful if (e.g. as in eqi_fromgeqdsk) the moments
      ! representation is created prior to the standard spline representation.

      type(xplasma), pointer :: s
      integer, intent(out) :: ier

      !-----------------
      integer :: id_xmom2d
      !-----------------

      call xplasma_find_item(s,'__XMOM2D',id_xmom2d,ier)
      if(ier.ne.0) return

      s%dict(id_xmom2d)%author = xplasma_root

    end subroutine xpmom_protect

    subroutine xpmom_unprotect(s,ier)

      ! Return Fourier moments representation to normal status, i.e.
      ! "derived from equilibrium".

      type(xplasma), pointer :: s
      integer, intent(out) :: ier

      !-----------------
      integer :: id_xmom2d
      !-----------------

      call xplasma_find_item(s,'__XMOM2D',id_xmom2d,ier)
      if(ier.ne.0) return

      s%dict(id_xmom2d)%author = xplasma_xmhd

    end subroutine xpmom_unprotect

    subroutine xpmom_info(s,ier, kmom, ksurf, xrho, ikmax)

      !  get information on Fourier Spline module...

      type (xplasma), pointer :: s
      integer, intent(out) :: ier  ! completion code, 0=OK

      integer, intent(out), optional :: kmom  ! #moments
      integer, intent(out), optional :: ksurf ! #radial grid pts
      real*8, dimension(:), intent(out), optional :: xrho ! radial grid itself
      integer, intent(out), optional :: ikmax ! max exponent for scaling

      !------------------
      integer :: id_rho,inrho,id_R,id_Z,iflag,iertmp
      type (xmom2d), pointer :: m
      !------------------

      m => s%xm

      call xplasma_common_ids(s,ier,id_R=id_R,id_Z=id_Z)
      if(ier.ne.0) return

      if(min(id_R,id_Z).eq.0) then
         ier=608
         call xplasma_errmsg_append(s, &
              ' ?xpmom_info: define equilibrium first.')
         return
      endif

      if(present(kmom)) then
         kmom = m%kmom
         if(kmom.eq.0) kmom = s%kmom
      endif

      if(present(ksurf)) then
         ksurf = m%ksurf
         if(ksurf.eq.0) then
            call xplasma_find_item(s,'__RHO',id_rho,iertmp)
            if(iertmp.eq.0) then
               call xplasma_grid_size(s,id_rho,ksurf,iertmp)
            endif
         endif
      endif

      if(present(xrho)) then
         inrho = m%ksurf
         iflag = 0
         if(inrho.eq.0) then
            iflag = 1
            call xplasma_find_item(s,'__RHO',id_rho,iertmp)
            if(iertmp.eq.0) then
               call xplasma_grid_size(s,id_rho,inrho,iertmp)
            endif
         endif

         if(inrho.ne.size(xrho)) then
            ier=9999
            call xplasma_errmsg_append(s, &
                 ' ?xpmom_info: xrho grid size not correct.')
         else
            if(iflag.eq.0) then
               xrho = m%xrho(1:inrho,1)
            else
               call xplasma_grid(s,id_rho,xrho,iertmp)
            endif
         endif
      endif

      if(present(ikmax)) then
         ikmax = m%kmaxe
      endif

    end subroutine xpmom_info

    subroutine xpmom_get1(s,x,inorm,im,ier, rcmom,rsmom,zcmom,zsmom, ideriv)

      !  Fourier Spline Equilibrium (R,Z) data access...
      !  get a single moment (or 1 for each type Rcos Rsin Zcos Zsin)
      !  at the specified x grid, using the moment coeffcient spline.

      type (xplasma), pointer :: s

      real*8, intent(in), dimension(:) :: x  ! grid at which moments are needed
      integer, intent(in) :: inorm        ! 0: return moment; 1: moment/x**im
      integer, intent(in) :: im           ! moment number
      !  coefficient for Rcos(im*theta), Rsin(im*theta), etc...

      integer, intent(out) :: ier         ! status code,  0=normal

      real*8, dimension(:), intent(out), optional :: rcmom  ! Rcos moment out
      real*8, dimension(:), intent(out), optional :: rsmom  ! Rsin moment out
      real*8, dimension(:), intent(out), optional :: zcmom  ! Zcos moment out
      real*8, dimension(:), intent(out), optional :: zsmom  ! Zsin moment out
      
      !  optional: return 1st derivative instead of value
      integer, intent(in), optional :: ideriv  ! =1 for 1st derivative
      !  (default: ideriv=0, return value)

      !  The moments splines as stored internally are scaled by 1/x**j; the
      !  formula for flux surface reconstruction is:
      !
      !    R = R0 + sum(x**j*(Rj*cos(j*theta) + Rj~*sin(j*theta)))
      !    Z = Z0 + sum(x**j*(Zj~*cos(j*theta) + Zj*sin(j*theta)))
      !
      !  where {Rj,Rj~,Zj,Zj~} are the scaled moments.
      !  the unscaled moments (inorm=0) are e.g. x**j*Rj.
      !
      !  So e.g. "rcmom" value returned, im=0:
      !    ideriv=0, inorm=0:  R0
      !    ideriv=0, inorm=1:  R0
      !    ideriv=1, inorm=0:  R0'
      !    ideriv=1, inorm=1:  R0'
      !
      !  So e.g. "rcmom" value returned, j = im > 0:
      !    ideriv=0, inorm=0:  x**j*Rj
      !    ideriv=0, inorm=1:  Rj
      !    ideriv=1, inorm=0:  j*x**(j-1)*Rj + x**j*Rj'  (R1 + x*R1' for im=1)
      !    ideriv=1, inorm=1:  Rj'
      !
      !-----------------
      integer :: ix,id_xmom2d,kmom,ksurf,nx,ikmax,imexp,inx
      logical :: lderiv,sflag
      real*8 :: zxrel,zfac,xx(size(x))

      type (xmom2d), pointer :: m

      !-----------------

      ier=0

      nx=size(x)
      inx=nx

      if(present(rcmom)) then
         rcmom=0.0d0
         inx=min(inx,size(rcmom))
      endif
      if(present(rsmom)) then
         rsmom=0.0d0
         inx=min(inx,size(rsmom))
      endif
      if(present(zcmom)) then
         zcmom=0.0d0
         inx=min(inx,size(zcmom))
      endif
      if(present(zsmom)) then
         zsmom=0.0d0
         inx=min(inx,size(zsmom))
      endif

      if(inx.lt.nx) then
         ier=506
         call xplasma_errmsg_append(s, &
              ' xpmom_get1: moment output array size less than x array size')
         return
      endif

      !  check for derivative request

      lderiv=.FALSE.
      if(present(ideriv)) then
         lderiv = (ideriv.eq.1)
      endif

      !  check target x vector range...

      if(minval(x).lt.-1.0d-10) then
         call xplasma_errmsg_append(s,' x coordinate out of range low (.lt.0)')
         ier=514
      endif
      if(maxval(x).gt.1.0d0+1.0d-10) then
         call xplasma_errmsg_append(s,' x coordinate out of range high (.gt.1)')
         ier=514
      endif
      if(ier.gt.0) return

      do ix=1,nx
         xx(ix)=x(ix)
         if(xx(ix).lt.0.0d0) xx(ix)=0.0d0
         if(xx(ix).gt.1.0d0) xx(ix)=1.0d0
      enddo

      call xplasma_find_item(s,'__XMOM2D',id_xmom2d,ier,nf_noerr=.TRUE.)
      if(id_xmom2d.eq.0) call xpmom_getall(s,id_xmom2d,ier)
      if(ier.ne.0) return

      m => s%xm
      kmom = m%kmom
      if(kmom.eq.0) then
         call xpmom_getall(s,id_xmom2d,ier)
         if(ier.ne.0) return
         kmom = m%kmom
      endif

      ksurf= m%ksurf
      ikmax= m%kmaxe

      if((im.lt.0).or.(im.gt.kmom)) return

      imexp=min(im,ikmax)
      sflag = lderiv.AND.(imexp.gt.0).AND.(inorm.eq.0)

      if(present(rcmom)) then
         call xpmom_eval(s,im,inorm,imexp,lderiv,sflag, &
              xx,nx,rcmom, &
              ksurf,m%xrho,m%sxrmmc, ier)
      endif

      if(present(rsmom)) then
         call xpmom_eval(s,im,inorm,imexp,lderiv,sflag, &
              xx,nx,rsmom, &
              ksurf,m%xrho,m%sxrmms, ier)
      endif

      if(present(zcmom)) then
         call xpmom_eval(s,im,inorm,imexp,lderiv,sflag, &
              xx,nx,zcmom, &
              ksurf,m%xrho,m%sxzmmc, ier)
      endif

      if(present(zsmom)) then
         call xpmom_eval(s,im,inorm,imexp,lderiv,sflag, &
              xx,nx,zsmom, &
              ksurf,m%xrho,m%sxzmms, ier)
      endif

    end subroutine xpmom_get1

    subroutine xpmom_eval(s,im,inorm,imexp,lderiv,sflag, &
         xx,nx,rzmom, &
         ksurf,xrho,mspl, ier)

      ! eval routine for xpmom_get1

      type (xplasma), pointer :: s

      integer, intent(in) :: im     ! moment index
      integer, intent(in) :: inorm  ! =0: moment; =1: scaled moment as stored
      integer, intent(in) :: imexp  ! scaling power = min(<limit>,im)
      logical, intent(in) :: lderiv ! .TRUE. for 1st derivative
      logical, intent(in) :: sflag  ! .TRUE. if scaled derivative sum needed

      integer, intent(in) :: nx     ! xx array dimensino
      real*8, intent(in) :: xx(nx)  ! target evaluation values xx(...)
      real*8, dimension(:) :: rzmom  ! evaluation result (1:nx)

      integer, intent(in) :: ksurf  ! no. of x grid pts in Fourier spline
      real*8, dimension(:,:), intent(in) :: xrho ! spline x grid data
      real*8, dimension(:,:,:), intent(in) :: mspl ! spline data
      
      integer, intent(out) :: ier   ! status code returned, 0=normal

      !------------------------
      integer :: iwarn,imp1
      integer :: ictv(3),ictd(3),ict2(3)

      data ictv/1,0,0/  ! for value fetch
      data ictd/0,1,0/  ! for derivative fetch
      data ict2/1,1,0/  ! to fetch both...
      
      real*8, dimension(:,:), allocatable :: zbuf
      !------------------------
      imp1 = im+1

      if(sflag) then
         allocate(zbuf(nx,2))
         call r8spvec(ict2,nx,xx,nx,zbuf, &
              ksurf,xrho,mspl(1,1,imp1), iwarn,ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ? xpmom_get1: unexpected spline evaluation error.')
            ier=9999
            deallocate(zbuf)
            return
         endif

         ! note: sflag setting implies: imexp.gt.0, lderiv=.TRUE., inorm=0

         if(imexp.eq.1) then
            rzmom(1:nx)= zbuf(:,1) + xx*zbuf(:,2)
         else
            rzmom(1:nx)= (imexp*xx*(imexp-1))*zbuf(:,1) + (xx**imexp)*zbuf(:,2)
         endif

         deallocate(zbuf)

      else if(lderiv) then
         ! derivative evaluation

         call r8spvec(ictd,nx,xx,nx,rzmom, &
              ksurf,xrho,mspl(1,1,imp1), iwarn,ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ? xpmom_get1: unexpected spline evaluation error.')
            ier=9999
            deallocate(zbuf)
            return
         endif

      else
         ! value evaluation

         call r8spvec(ictv,nx,xx,nx,rzmom, &
              ksurf,xrho,mspl(1,1,imp1), iwarn,ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ? xpmom_get1: unexpected spline evaluation error.')
            ier=9999
            return
         endif
         if((inorm.eq.0).and.(imexp.gt.0)) then
            rzmom(1:nx)=rzmom(1:nx)*xx**imexp
         endif

      endif

    end subroutine xpmom_eval

    subroutine xpmom_getall(s,id,ier)

      !  build moments from equilibrium (if possible)

      type(xplasma), pointer :: s
      integer,intent(out) :: id   ! __XMOM2D id (returned)
      integer,intent(out) :: ier  ! status code, 0=normal (returned)

      !---------------------------------
      integer :: id_R,id_Z,id_rho,id_th,inrho,intk,kmom,iertmp
      integer :: ith,irho
      real*8, dimension(:), allocatable :: zrho,zth,zr,zz
      real*8, dimension(:,:), allocatable :: rmmc,rmms,zmmc,zmms
      real*8, dimension(:,:,:), allocatable :: zrmc,zzmc
      logical :: ilsym

      type (xpeval) :: th_intrp, rho_intrp
      !---------------------------------

      id=0

      !  check that an equilibrium and __RHO grid are available...

      call xplasma_common_ids(s,ier, id_R=id_R,id_Z=id_Z)
      if(ier.ne.0) return

      if(min(id_R,id_Z).eq.0) then
         call xplasma_errmsg_append(s, &
              ' ?xpmom_getall: no equilibrium, cannot produce Fourier spline')
         ier=608
         return
      endif

      call xplasma_find_item(s,'__RHO',id_rho,ier)
      if(ier.ne.0) return

      call xplasma_find_item(s,'__CHI',id_th,ier)
      if(ier.ne.0) return

      call xplasma_grid_size(s,id_rho,inrho,iertmp)
      allocate(zrho(inrho))
      call xplasma_grid(s,id_rho,zrho,iertmp)

      !  init structure to hold Fourier Spline; use __RHO rho grid...

      call xpmom_init(s,ier, rhoi=zrho)
      if(ier.ne.0) then
         deallocate(zrho)
         return
      endif

      kmom=s%kmom
      allocate(rmmc(0:kmom,inrho),rmms(0:kmom,inrho))
      allocate(zmmc(0:kmom,inrho),zmms(0:kmom,inrho))
      allocate(zrmc(3,0:kmom,2),zzmc(3,0:kmom,2)) ! for r8tkbmmc legacy routine

      call r8ntk(intk)
      allocate(zth(intk),zr(intk),zz(intk))

      do ith=1,intk
         zth(ith)=(ith-1)*C2PI/intk  ! even spacing, 2pi pt not duplicated.
      enddo

      ilsym=.FALSE.
      do
         call xplasma_x_lookup(s, id_th, zth, th_intrp, ier)
         if(ier.ne.0) exit

         do irho=1,inrho
            call xplasma_x_lookup(s, id_rho, zrho(irho), rho_intrp, ier)
            if(ier.ne.0) exit

            call xplasma_eval_prof(s, id_R, th_intrp, rho_intrp, zR, ier)
            if(ier.ne.0) exit

            call xplasma_eval_prof(s, id_Z, th_intrp, rho_intrp, zZ, ier)
            if(ier.ne.0) exit

            if(irho.eq.1) then
               rmmc(0:kmom,irho)=sum(zR)/intk
               zmmc(0:kmom,irho)=sum(zZ)/intk
               rmms(0:kmom,irho)=0.0d0
               zmms(0:kmom,irho)=0.0d0
            else
               call r8tkbmmc(ilsym,zrmc,zzmc,kmom,kmom,zR,zZ)
               rmmc(0:kmom,irho)=zrmc(3,0:kmom,1)
               rmms(1:kmom,irho)=zrmc(3,1:kmom,2); rmms(0,irho)=0.0d0
               zmmc(0:kmom,irho)=zzmc(3,0:kmom,1)
               zmms(1:kmom,irho)=zzmc(3,1:kmom,2); zmms(0,irho)=0.0d0
            endif

         enddo
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xpmom_getall: (R,Z) interpolation error.')
            exit
         endif

         call xpmom_setup(s,rmmc,rmms,zmmc,zmms,ier)

         call xplasma_find_item(s,'__XMOM2D',id,ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xpmom_getall: __XMOM2D should have been created.')
         endif
         exit
      enddo

      deallocate(zrho,zth,zr,zz,zrmc,zzmc)
      deallocate(rmmc,rmms,zmmc,zmms)
      call xpeval_free(rho_intrp)
      call xpeval_free(th_intrp)

    end subroutine xpmom_getall

    !==========================================================
    subroutine xplasma_create_integ(s,iname,x,id,ier, &
         icoordx, rhomin, nth_expand, icoordx2, x2, cache_enable)

      !  create a new "numerical integration control dataset" with specified 
      !  name
 
      type (xplasma), pointer :: s
      character*(*), intent(in) :: iname

      real*8, intent(in), target :: x(:) ! integration grid (strict ascending)
      !  this defines a sequence of points {x1,x2,x3,...} which are ranges
      !  over which definite numerical integrals will be performed:
      !  from x1 to x2, x2 to x3, etc.  If there are N points in the
      !  integration grid there will be N-1 integration results.
      
      integer, intent(out) :: id  ! id of integrator object
      integer, intent(out) :: ier ! status code (0=OK)

      integer, intent(in), optional :: icoordx
                                     ! id of coordinate object over which
                                     ! integrations are to be defined
      !  icoordx can also be a grid id, in which case the coordinate associated
      !  with the grid is the coordinate of integration
      !  DEFAULT: xplasma_rho_coord -- the only one supported at present!
      !    (dmc Aug 2006).

      real*8, intent(in), optional :: rhomin  ! minimum "rho"
      ! (optional control to adjust minimum "rho" to avoid axial
      ! singularity; default is generally OK

      integer, intent(in), optional :: nth_expand   ! optional -- subdivide
      ! poloidal angle grid into nth_expand subzones per base zone (default 1).

      integer, intent(in), optional :: icoordx2     ! optional -- specify
      ! 2nd integration coordinate-- DEFAULT: xplasma_theta_coord
      !   (only supported 2nd coordinate at present -- dmc Aug. 2006)
      
      real*8, intent(in), optional, target :: x2(:) ! 2nd coordinate grid
      
      ! NOTE also: icoordx and icoordx2 must be distinct.

      logical, intent(in), optional :: cache_enable  ! .TRUE. to enable 
      ! integrand evaluation cache-- for 1d integrands only; worthwhile if 
      ! several are to be evaluated on the same grid.

      !---------------------------
      ! DMC Feb 2006 -- only icoordx = (rho, radial flux coordinate)
      ! is accepted currently...
      !
      ! if x2 is present it must be poloidal angle theta
      !    (interface allows generalization at some future time)
      !---------------------------

      !---------------------------
      integer :: icoord,jcoordx,jcoordx2,iadr,i,isize,idloc
      integer :: iitype,irank,ith_expand
      real*8 :: zrhomin
      real*8, dimension(:), pointer :: xp
      logical :: icache
      !---------------------------

      ier = 0
      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         call xplasma_init(s,'(init due to xplasma_create_integ call)',ier)
      endif
      if(ier.ne.0) return

      id=0
      !  make a "black box" to hold numerical integration information...
      call xplasma_add_item(s,iname,xplasma_blkbxType,idloc,ier)
      if(ier.ne.0) return

      icoord=xplasma_rho_coord
      if(present(icoordx)) icoord=icoordx
      call xplasma_errck0_xinteg(s,icoord,x,jcoordx,ier)
      if(ier.ne.0) return
      
      if(jcoordx.ne.xplasma_rho_coord) then
         ier=601
         call xplasma_errmsg_append(s, &
              '  not yet supported as an integration coordinate: "'//trim(s%dict(icoordx)%name)//'".')
      endif
      if(ier.ne.0) return

      jcoordx2=0
      if(present(x2)) then
         icoord=xplasma_theta_coord
         if(present(icoordx2)) icoord=icoordx2
         call xplasma_errck0_xinteg(s,icoordx2,x2,jcoordx2,ier)
      endif
      if(ier.ne.0) return

      if(jcoordx2.gt.0) then
         if(jcoordx2.ne.xplasma_theta_coord) then
            ier=601
            call xplasma_errmsg_append(s, &
                 '  not yet supported as 2nd integration coordinate: "'//trim(s%dict(icoordx2)%name)//'".')
         endif
      endif

      ! {R,Z} must be defined...

      if(min(s%id_R,s%id_Z).eq.0) then
         ier=608
      endif
      if(ier.ne.0) return

      !  OK............

      id=idloc

      zrhomin=1.0d-7
      if(present(rhomin)) zrhomin=rhomin
      iitype=1

      icache=.FALSE.

      if(jcoordx2.eq.0) then
         xp => x
         if(present(cache_enable)) icache = cache_enable
      else
         xp => x2
      endif

      ith_expand=1
      if(present(nth_expand)) ith_expand=max(1,nth_expand)

      call gen_integ(s,id,iitype,jcoordx,jcoordx2,zrhomin,ith_expand,x,xp, &
           icache)

    end subroutine xplasma_create_integ

    subroutine gen_integ(s,id,iitype,jcoordx1,jcoordx2,zrhomin,ith_expand, &
         x1,x2,icache)

      ! private -- build integrator dataset (no error checking)

      type (xplasma), pointer :: s
      integer, intent(in) :: id      ! location of empty black box
      integer, intent(in) :: iitype  ! type of integrator dataset
      integer, intent(in) :: jcoordx1,jcoordx2  ! x coordinate codes
      real*8, intent(in) :: zrhomin  ! rhomin (singularity avoidance)
      integer, intent(in) :: ith_expand  ! theta grid multiplier
      real*8, intent(in), dimension(:) :: x1,x2
      logical, intent(in) :: icache

      ! the dataset contains integer sizes and addresses and a real*8 buffer
      ! with room for:
      !   integration grids
      !   merged spline & integration "break" grids
      !   volumes
      !   dV/drho

      !-------------------------------------
      integer :: iadr,ijtype,irank,isizep,isize,iertmp,inthvec,inthx
      integer :: idrho,idth,inrho,inth,idxr1,idxr2,inx1,inx2,ic,i,ii,j
      integer, dimension(:), pointer :: ip
      real*8, dimension(:), pointer :: r8p
      real*8, dimension(:), allocatable :: xrho,xth,xthx,x2c,xrho_brk,xth_brk
      integer :: inrho_brk,inth_brk
      integer :: inzones,insurfs,inzonth
      real*8 :: zfac
      real*8, parameter :: c1=1.0d0
      !-------------------------------------

      iadr=s%dict(id)%dindex

      ijtype=iitype
      irank=1
      if((ijtype.eq.1).and.(jcoordx2.gt.0)) then
         ijtype=2
      endif
      if(ijtype.eq.2) irank=2

      s%blkbxs(iadr)%type = iitype
      allocate(s%blkbxs(iadr)%iwork(50)); ip => s%blkbxs(iadr)%iwork

      ip = 0
      ip(1)=irank
      ip(2)=jcoordx1
      ip(3)=jcoordx2

      ip(9)=ith_expand
      ip(10)=s%eq_counter

      inx1=size(x1)
      inx2=size(x2)

      ip(11)=2
      ip(12)=inx1
      if(irank.eq.2) then
         allocate(x2c(inx2)); x2c=x2  ! rank 2 -- theta grid for integration
      else
         inx2=2
         allocate(x2c(inx2))  ! rank 1 -- integrate over all theta
      endif
      ip(13)=inx1+2
      ip(14)=inx2
      isize=inx1+inx2+1  ! so far...

      !  create integration grids with spline nodes inserted

      call xplasma_prof_info(s,s%id_R,iertmp, &
           gridId1=idxr1,gridId2=idxr2)

      call xplasma_grid_info(s,idxr1,iertmp, coord=ic)
      if(ic.eq.xplasma_rho_coord) then
         idrho=idxr1
         idth=idxr2
      else
         idrho=idxr2
         idth=idxr1
      endif

      call xplasma_grid_size(s,idrho,inrho,iertmp)
      call xplasma_grid_size(s,idth,inth,iertmp)

      inthx = 1 + (inth-1)*ith_expand

      allocate(xrho(inrho),xrho_brk(inrho+inx1))
      allocate(xth(inth),xthx(inthx),xth_brk(inthx+inx2))

      call xplasma_grid(s,idrho,xrho,iertmp)
      call xplasma_grid(s,idth,xth,iertmp)

      if(irank.eq.1) then
         x2c(1)=xth(1)
         x2c(2)=xth(inth)
      endif

      xthx(1)=xth(1)
      ii=1
      do i=2,inth
         xthx(ii+ith_expand)=xth(i)
         do j=ii+1,ii+ith_expand-1
            zfac=(j-ii)*c1/ith_expand
            xthx(j)=(c1-zfac)*xthx(ii) + zfac*xthx(ii+ith_expand)
         enddo
         ii=ii+ith_expand
      enddo

      !  insert nodes into integration grid; removed duplicates w/in tolerance

      call fluxav_brk_addsub(x1,inx1,xrho,inrho,xrho_brk,inrho_brk)
      call fluxav_brktol(xrho_brk,inrho_brk,s%bdytol)

      call fluxav_brk_addperio(x2c,inx2,xthx,inthx,xth_brk,inth_brk)
      call fluxav_brktol(xth_brk,inth_brk,s%bdytol)

      ip(21)= isize + 1
      ip(22)= inrho_brk

      inthvec = 10*(inth_brk-1) ! 10 evaluation points per theta break zone

      ip(23)= isize + inrho_brk + 1
      ip(24)= inthvec  ! expanded for vector integrating w/in theta breaks

      ip(25)= isize + inrho_brk + 1 + inthvec
      ip(26)= inthvec  ! poloidal integration weights here...

      isizep = isize
      isize = isizep + inrho_brk + 2*inthvec
      
      inzones=(inx1-1)*(inx2-1)
      insurfs=inx1*(inx2-1)

      if(irank.eq.1) then
         ip(30)=isize+1
         ip(31)=-inzones             ! dVol cache
         ip(32)=isize+inzones+1
         ip(33)=-insurfs             ! (dV/drho) cache
         ip(34)=isize+inzones+insurfs+1
         ip(35)=-inx1                ! Bmax(x) cache, 1/surface only.

         ip(38)=isize+inzones+insurfs+inx1+1

      else
         ! these items not cached for 2d output integrators...
         ip(38)=isize+1
      endif

      if(icache) then
         ip(39) = inthvec*insurfs*eqi_intg_neval

         ! ip(41:50) track cache status... start out as 0, empty
      else
         ip(39) = 0
      endif
      isize = ip(38)+ip(39)

      !  create real*8 block...

      allocate(s%blkbxs(iadr)%r8work(isize))

      s%blkbxs(iadr)%r8work=czero !MG added to please Intel compiler

      r8p => s%blkbxs(iadr)%r8work

      r8p(1)=zrhomin
      r8p(2:inx1+1) = x1
      r8p(inx1+2:inx1+inx2+1) = x2c

      r8p(ip(21):ip(21)+inrho_brk-1) = xrho_brk(1:inrho_brk)

      call xplasma_integ_mkvec(inth_brk,xth_brk, &
           r8p(ip(23):ip(23)+inthvec-1), &
           r8p(ip(25):ip(25)+inthvec-1))

      if(irank.eq.1) then
         r8p(ip(30):ip(30)+inzones-1) = 0  ! dVol not calculated yet
         r8p(ip(32):ip(32)+insurfs-1)= 0   ! dV/drho not calculated yet
         r8p(ip(34):ip(34)+inx1-1)= 0      ! Bmax not calculated yet
      endif

    end subroutine gen_integ

    subroutine regen_integ(s,id)

      ! re-gen integrator structure e.g. after equilibrium update

      type (xplasma), pointer :: s
      integer, intent(in) :: id

      !-------------------------
      integer :: iadr,iitype,jcoordx1,jcoordx2,icounter
      integer :: ilocx1,ilocx2,isizex1,isizex2,ith_expand
      real*8 :: zrhomin
      real*8, dimension(:), allocatable :: x1,x2
      integer, dimension(:), pointer :: ip
      real*8, dimension(:), pointer :: r8p
      logical :: icache
      !-------------------------

      iadr = s%dict(id)%dindex
      iitype = s%blkbxs(iadr)%type

      ip => s%blkbxs(iadr)%iwork
      r8p => s%blkbxs(iadr)%r8work

      ith_expand = ip(9)

      icounter = ip(10)
      if(icounter.eq.s%eq_counter) return

      jcoordx1 = ip(2)
      jcoordx2 = ip(3)

      ilocx1 = ip(11)
      isizex1 = ip(12)

      if(jcoordx2.eq.0) then
         ilocx2 = ilocx1
         isizex2 = 1
      else
         ilocx2 = ip(13)
         isizex2 = ip(14)
      endif

      zrhomin = r8p(1)

      allocate(x1(isizex1),x2(isizex2))
      x1 = r8p(ilocx1:ilocx1+isizex1-1)
      x2 = r8p(ilocx2:ilocx2+isizex2-1)

      icache = (ip(39).gt.0)

      call xpblkbx_free(s%blkbxs(iadr))
      nullify(ip,r8p)

      call gen_integ(s,id,iitype,jcoordx1,jcoordx2,zrhomin,ith_expand, &
           x1,x2,icache)

      deallocate(x1,x2)

    end subroutine regen_integ

    !==========================================================
    subroutine xplasma_augment_integ(s,iname,idi,x2,id,ier, icoordx2)

      !  create a new "numerical integration control dataset" with specified 
      !  name, by augmenting an existing set with a new dimension.
 
      type (xplasma), pointer :: s
      character*(*), intent(in) :: iname
      integer, intent(in) :: idi      ! id of existing integrator dataset

      real*8, intent(in) :: x2(:)  ! 2nd integration grid (strict ascending)
      
      integer, intent(out) :: id  ! id of integrator object
      integer, intent(out) :: ier ! status code (0=OK)

      integer, intent(in), optional :: icoordx2
                                      ! id of coordinate object over which
                                      ! integrations are to be defined
      !  icoordx2 can also be a grid id, in which case the coordinate 
      !  associated with the grid is the coordinate of integration
      !  DEFAULT: xplasma_theta_coord -- only one supported at present!

      ! (optional control to adjust minimum "rho" to avoid axial
      ! singularity; default is generally OK

      !---------------------------
      ! DMC Feb 2006 -- only icoordx = (rho, radial flux coordinate)
      ! is accepted currently; and x2, here, must be "theta", the poloidal
      ! angle flux coordinate.
      !    (interface allows generalization at some future time)
      !---------------------------

      !---------------------------
      integer :: icoord,jcoordx1,jcoordx2,idloc,ith_expand
      integer :: iitype,irank,iadi,iadr,isize,isizex1,ilx1
      real*8 :: zrhomin
      integer, dimension(:), pointer :: ipp
      real*8, dimension(:), pointer :: r8pp
      character*128 zmsg
      logical :: icache
      !---------------------------

      ier=0
      id=0
      call xplasma_errck0(s,idi,ier)
      if(ier.ne.0) return

      if(s%dict(idi)%dtype.ne.xplasma_blkbxType) then
         ier=605
         call xplasma_errmsg_append(s,'  id of item passed points to: "'// &
              trim(s%dict(idi)%name)//'".')
      else
         iadi=s%dict(idi)%dindex
         if(s%blkbxs(iadi)%type.ne.1) then
            ier=605
            call xplasma_errmsg_append(s, &
                 '  at present only "type 1" integrator dataset is compatible with xplasma_augment_integ(...')
            call xplasma_errmsg_append(s, &
                 '  dataset name is: "'//trim(s%dict(idi)%name)//'".')
         endif
      endif
      if(ier.gt.0) return

      iadr = s%dict(idi)%dindex
      ipp => s%blkbxs(iadr)%iwork
      r8pp => s%blkbxs(iadr)%r8work

      call xplasma_add_item(s,iname,xplasma_blkbxType,idloc,ier)
      if(ier.ne.0) return

      icoord=xplasma_theta_coord
      if(present(icoordx2)) icoord=icoordx2
      call xplasma_errck0_xinteg(s,icoord,x2,jcoordx2,ier)
      if(ier.ne.0) return
      
      if(jcoordx2.ne.xplasma_theta_coord) then
         ier=601
         call xplasma_errmsg_append(s, &
              '  not yet supported as 2nd integration coordinate: "'//trim(s%dict(icoordx2)%name)//'".')
      endif

      ! OK...

      id=idloc
      iitype=2
      
      zrhomin=r8pp(1)

      jcoordx1=ipp(2)
      ith_expand=ipp(9)
      ilx1=ipp(11)
      isizex1=ipp(12)

      icache=.FALSE.
      call gen_integ(s,id,iitype,jcoordx1,jcoordx2,zrhomin,ith_expand, &
           r8pp(ilx1:ilx1+isizex1-1),x2,icache)

    end subroutine xplasma_augment_integ

    !======================================================
    subroutine xplasma_integ_dva(s,idi,ier, dvol,dvdrho,bmaxa)

      !  store volume and dV/drho information

      type (xplasma), pointer :: s
      integer, intent(in) :: idi  ! set id (a "black box")
      integer, intent(out) :: ier ! status code 0=OK

      real*8, intent(in), optional :: dvol(:)   ! volume elements @zones

      real*8, intent(in), optional :: dvdrho(:) ! dV/drho elements @surfaces

      real*8, intent(in), optional :: bmaxa(:)  ! max(|B|) on surfaces

      !-----------------------------------
      integer :: iadr,inzons,insurfs,insurfb,iloc
      integer :: ii,ith,irho
      !-----------------------------------

      !------------
      ! error checks

      call xplasma_errck0(s,idi,ier)
      if(ier.ne.0) return

      if(s%dict(idi)%dtype.ne.xplasma_blkbxType) then
         ier=605
      else
         iadr = s%dict(idi)%dindex
         if((s%blkbxs(iadr)%type.ne.1).and. &
              (s%blkbxs(iadr)%type.ne.2)) then
            ier=605
         endif
      endif

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              '  passed id associated item name: '//trim(s%dict(idi)%name))
         return
      endif

      inzons = s%blkbxs(iadr)%iwork(31)
      insurfs = s%blkbxs(iadr)%iwork(33)
      insurfb = s%blkbxs(iadr)%iwork(35)

      if(present(dvol)) then
         iloc = s%blkbxs(iadr)%iwork(30)
         if(inzons.ge.0) then
            ier=613
         else
            inzons=-inzons
            if(size(dvol).ne.inzons) then
               ier=612
               call xplasma_errmsg_append(s, &
                    '   passed array DVOL of incorrect length')
            else
               s%blkbxs(iadr)%iwork(31)=inzons
               s%blkbxs(iadr)%r8work(iloc:iloc+inzons-1) = dvol
            endif
         endif
      endif
      if(ier.ne.0) return

      if(present(dvdrho)) then
         iloc = s%blkbxs(iadr)%iwork(32)
         if(insurfs.ge.0) then
            ier=613
         else
            insurfs=-insurfs
            if(size(dvdrho).ne.insurfs) then
               ier=612
               call xplasma_errmsg_append(s, &
                    '   passed array DVDRHO of incorrect length')
            else
               s%blkbxs(iadr)%iwork(33)=insurfs
               s%blkbxs(iadr)%r8work(iloc:iloc+insurfs-1) = dvdrho
            endif
         endif
      endif
      if(ier.ne.0) return

      if(present(bmaxa)) then
         iloc = s%blkbxs(iadr)%iwork(34)
         if(insurfb.ge.0) then
            ier=613
         else
            insurfb=-insurfb
            if(size(bmaxa).ne.insurfb) then
               ier=612
               call xplasma_errmsg_append(s, &
                    '   passed array BMAXA of incorrect length')
            else
               s%blkbxs(iadr)%iwork(35)=insurfb
               s%blkbxs(iadr)%r8work(iloc:iloc+insurfb-1) = bmaxa
            endif
         endif
      endif
      if(ier.ne.0) return

    end subroutine xplasma_integ_dva

    subroutine xplasma_integ_scache(s,idi,cache_init,ier)

      !  update cache status information

      type (xplasma), pointer :: s
      integer, intent(in) :: idi  ! set id (a "black box")
      logical, intent(in), dimension(:) :: cache_init
      integer, intent(out) :: ier ! status code 0=OK

      !------------------------------
      integer :: ii,iadr
      !------------------------------
      ! error checks

      call xplasma_errck0(s,idi,ier)
      if(ier.ne.0) return

      if(s%dict(idi)%dtype.ne.xplasma_blkbxType) then
         ier=605
      else
         iadr = s%dict(idi)%dindex
         if((s%blkbxs(iadr)%type.ne.1).and. &
              (s%blkbxs(iadr)%type.ne.2)) then
            ier=605
         endif
      endif

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              '  passed id associated item name: '//trim(s%dict(idi)%name))
         return
      endif

      if(min(s%id_R,s%id_Z).eq.0) then
         ier=608  ! no {R,Z} surfaces...
      endif
      if(ier.ne.0) return

      if(s%blkbxs(iadr)%iwork(39).eq.0) then
         ier=614
      endif
      if(ier.ne.0) return

      !------------

      do ii=1,min(eqi_intg_neval,size(cache_init))
         if(cache_init(ii)) then
            s%blkbxs(iadr)%iwork(40+ii)=1
         else
            s%blkbxs(iadr)%iwork(40+ii)=0
         endif
      enddo

    end subroutine xplasma_integ_scache

    !======================================================
    subroutine xplasma_integ_info(s,idi,ier, &
         n_rho_int,n_th_int,n_rho_brk, &
         rho_int,th_int,rho_brk, &
         dvol_avail,dvol,dvdrho_avail,dvdrho, &
         bmaxa_avail,bmaxa, &
         thvec_ptr, thwts_ptr, cache_ptr, &
         cache_init, rhomin)

      !  information on xplasma numeric integrator dataset

      type (xplasma), pointer :: s
      integer, intent(in) :: idi  ! set id (a "black box")
      integer, intent(out) :: ier ! status code 0=OK

      ! integration grids:
      integer, intent(out), optional :: n_rho_int
      real*8, intent(out), dimension(:), optional :: rho_int
      integer, intent(out), optional :: n_th_int
      real*8, intent(out), dimension(:), optional :: th_int

      ! rho grid with spline knot pts inserted
      integer, intent(out), optional :: n_rho_brk
      real*8, intent(out), dimension(:), optional :: rho_brk

      ! volume element info:
      logical, intent(out), optional :: dvol_avail
      real*8, intent(out), dimension(:), optional :: dvol

      ! surface dV/drho info:
      logical, intent(out), optional :: dvdrho_avail
      real*8, intent(out), dimension(:), optional :: dvdrho

      ! |B|max on surface info:
      logical, intent(out), optional :: bmaxa_avail
      real*8, intent(out), dimension(:), optional :: bmaxa

      ! pointers to work arrays

      real*8, dimension(:), pointer, optional :: thvec_ptr  ! theta eval vector
      real*8, dimension(:), pointer, optional :: thwts_ptr  ! integ. wts vector
      real*8, dimension(:), pointer, optional :: cache_ptr  ! data cache 

      real*8, intent(out), optional :: rhomin  !  singularity barrier @axis

      logical, intent(out), optional, dimension(:) :: cache_init
      !---------------------------------
      integer :: iadr,iarr,ilen,ii,inrhoi,inthi,ith,irho
      integer, dimension(:), pointer :: ip
      real*8, dimension(:), pointer :: r8p
      character*128 zmsg
      !---------------------------------

      if(present(n_rho_int)) n_rho_int = 0
      if(present(n_th_int)) n_th_int = 0
      if(present(n_rho_brk)) n_rho_brk = 0

      if(present(rho_int)) rho_int = 0
      if(present(th_int)) th_int = 0
      if(present(rho_brk)) rho_brk = 0

      if(present(dvol)) dvol = 0
      if(present(dvdrho)) dvdrho = 0
      if(present(bmaxa)) bmaxa = 0

      if(present(dvol_avail)) dvol_avail=.FALSE.
      if(present(dvdrho_avail)) dvdrho_avail=.FALSE.
      if(present(bmaxa_avail)) bmaxa_avail=.FALSE.

      if(present(thvec_ptr)) nullify(thvec_ptr)
      if(present(thwts_ptr)) nullify(thwts_ptr)
      if(present(cache_ptr)) nullify(cache_ptr)

      if(present(cache_init)) cache_init=.FALSE.

      if(present(rhomin)) rhomin=0

      !------------
      ! error checks

      call xplasma_errck0(s,idi,ier)
      if(ier.ne.0) return

      if(s%dict(idi)%dtype.ne.xplasma_blkbxType) then
         ier=605
      else
         iadr = s%dict(idi)%dindex
         if((s%blkbxs(iadr)%type.ne.1).and. &
              (s%blkbxs(iadr)%type.ne.2)) then
            ier=605
         endif
      endif

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              '  passed id associated item name: '//trim(s%dict(idi)%name))
         return
      endif

      if(min(s%id_R,s%id_Z).eq.0) then
         ier=608  ! no {R,Z} surfaces...
      endif
      if(ier.ne.0) return

      !---------------
      ! make sure dataset is up-to-date (equilibrium may have changed).

      call regen_integ(s,idi)

      !---------------

      iadr = s%dict(idi)%dindex
      ip => s%blkbxs(iadr)%iwork
      r8p => s%blkbxs(iadr)%r8work

      inrhoi = ip(12)
      if(present(n_rho_int)) then
         n_rho_int = inrhoi
      endif

      if(present(rho_int)) then
         iarr=ip(11)
         ilen=ip(12)
         if(size(rho_int).ne.ilen) then
            ier=612
            zmsg=' '
            write(zmsg,*) '   Passed array (rho_int) length: ',size(rho_int), &
                 ' Expected length: ',ilen
            call xplasma_errmsg_append(s,zmsg)
            call xplasma_errmsg_append(s, &
                 '  passed integrator dataset name: '//trim(s%dict(idi)%name))
         else
            rho_int = r8p(iarr:iarr+ilen-1)
         endif
      endif

      inthi=ip(14)
      if(present(n_th_int)) then
         n_th_int = inthi
      endif

      if(present(th_int)) then
         iarr=ip(13)
         ilen=ip(14)
         if(size(th_int).ne.ilen) then
            ier=612
            zmsg=' '
            write(zmsg,*) '   Passed array (th_int) length: ',size(th_int), &
                 ' Expected length: ',ilen
            call xplasma_errmsg_append(s,zmsg)
            call xplasma_errmsg_append(s, &
                 '  passed integrator dataset name: '//trim(s%dict(idi)%name))
         else
            th_int = r8p(iarr:iarr+ilen-1)
         endif
      endif

      if(present(n_rho_brk)) then
         n_rho_brk = ip(22)
      endif

      if(present(rho_brk)) then
         iarr=ip(21)
         ilen=ip(22)
         if(size(rho_brk).ne.ilen) then
            ier=612
            zmsg=' '
            write(zmsg,*) '   Passed array (rho_brk) length: ',size(rho_brk), &
                 ' Expected length: ',ilen
            call xplasma_errmsg_append(s,zmsg)
            call xplasma_errmsg_append(s, &
                 '  passed integrator dataset name: '//trim(s%dict(idi)%name))
         else
            rho_brk = r8p(iarr:iarr+ilen-1)
         endif
      endif

      if(present(dvol_avail)) then
         dvol_avail = ip(31).gt.0
      endif

      if(present(dvdrho_avail)) then
         dvdrho_avail = ip(33).gt.0
      endif

      if(present(bmaxa_avail)) then
         bmaxa_avail = ip(35).gt.0
      endif

      if(present(dvol)) then
         iarr=ip(30)
         ilen=abs(ip(31))
         if(size(dvol).ne.ilen) then
            ier=612
            zmsg=' '
            write(zmsg,*) '   Passed array (dvol) length: ',size(dvol), &
                 ' Expected length: ',ilen
            call xplasma_errmsg_append(s,zmsg)
            call xplasma_errmsg_append(s, &
                 '  passed integrator dataset name: '//trim(s%dict(idi)%name))
         else
            dvol = r8p(iarr:iarr+ilen-1)
         endif
      endif

      if(present(dvdrho)) then
         iarr=ip(32)
         ilen=abs(ip(33))
         if(size(dvdrho).ne.ilen) then
            ier=612
            zmsg=' '
            write(zmsg,*) '   Passed array (dvdrho) length: ',size(dvdrho), &
                 ' Expected length: ',ilen
            call xplasma_errmsg_append(s,zmsg)
            call xplasma_errmsg_append(s, &
                 '  passed integrator dataset name: '//trim(s%dict(idi)%name))
         else
            dvdrho = r8p(iarr:iarr+ilen-1)
         endif
      endif

      if(present(bmaxa)) then
         iarr=ip(34)
         ilen=abs(ip(35))
         if(size(bmaxa).ne.ilen) then
            ier=612
            zmsg=' '
            write(zmsg,*) '   Passed array (bmaxa) length: ',size(bmaxa), &
                 ' Expected length: ',ilen
            call xplasma_errmsg_append(s,zmsg)
            call xplasma_errmsg_append(s, &
                 '  passed integrator dataset name: '//trim(s%dict(idi)%name))
         else
            bmaxa = r8p(iarr:iarr+ilen-1)
         endif
      endif

      if(present(rhomin)) then
         rhomin=r8p(1)
      endif

      if(present(thvec_ptr)) then
         thvec_ptr => r8p(ip(23):ip(23)+ip(24)-1)
      endif

      if(present(thwts_ptr)) then
         thwts_ptr => r8p(ip(25):ip(25)+ip(26)-1)
      endif

      if(present(cache_ptr)) then
         if(ip(39).gt.0) then
            cache_ptr => r8p(ip(38):ip(38)+ip(39)-1)
         endif
      endif

      if(present(cache_init)) then
         if((ip(39).gt.0).and.(size(cache_init).ge.eqi_intg_neval)) then
            do ii=1,eqi_intg_neval
               cache_init(ii) = ip(40+ii).gt.0
            enddo
         endif
      endif
               
    end subroutine xplasma_integ_info

    !==========================================================
    !==========================================================
    !==========================================================
    !  private methods...
    !========================================

    subroutine xplookup(s,ivec,xvec,ccwflag, &
         inx,xpkg,imode, &
         ix,dxn,hx,hxi,iwarn)

      !  lookup points on a grid.  If the grid is evenly spaced do it here;
      !  if not, call pspline r8xlookup w/same arguments

      type(xplasma), pointer :: s

      integer, intent(in) :: ivec       ! evaluation vector size
      real*8, intent(in) :: xvec(ivec)  ! evaluation vector values
      logical, intent(in) :: ccwflag    ! periodic choord orientation flag
      ! (reversed if ccwflag=F)

      integer, intent(in) :: inx        ! grid to search on
      real*8, intent(in) :: xpkg(inx,4) ! grid xpkg (pspline r8genxpkg info)
      integer, intent(in) :: imode      ! argument for r8xlookup if needed

      integer, intent(out) :: ix(ivec)  ! grid cell indices for each xvec elem.
      real*8, intent(out) :: dxn(ivec)  ! normalized displacements w/in cell
      real*8, intent(out) :: hx(ivec)   ! grid cell widths
      real*8, intent(out) :: hxi(ivec)  ! inverse grid cell widths
      integer, intent(out) :: iwarn

      !---------------------------
      integer :: ilin,iper,i,inxm1
      real*8 :: zuse,zndx,zspan,zhav,zhavi
      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter :: ONE  = 1.0d0
      logical :: ireverse
      real*8, dimension(:), allocatable :: xuse
      !---------------------------

      iwarn=0

      if(xpkg(1,4).ne.ZERO) then
         ilin=0
      else
         ilin=1
      endif

      if(xpkg(2,4).ne.ZERO) then
         iper=1
      else
         iper=0
      endif

      ireverse=.FALSE.
      if((iper.eq.1).and.(.not.ccwflag)) ireverse=.TRUE.

      allocate(xuse(ivec))
      call xp_reverse(xvec,ireverse,xuse,xpkg(1,1),xpkg(inx,1))

      if(ilin.eq.0) then
         call r8xlookup(ivec,xuse,inx,xpkg,imode, &
              ix,dxn,hx,hxi,iwarn)

      else
         ! evenly spaced grid...
         inxm1=inx-1
         zspan=xpkg(inx,1)-xpkg(1,1)
         zhav=xpkg(inx,2)
         zhavi=xpkg(inx,3)

         if(iper.eq.0) then
            ! non-periodic grid, bounds checking applies
            do i=1,ivec
               zuse=xuse(i)
               if(zuse.lt.xpkg(1,1)) then
                  if((xpkg(1,1)-zuse).gt.s%bdytol) iwarn=iwarn+1
                  zuse=xpkg(1,1)
               else if(zuse.gt.xpkg(inx,1)) then
                  if((zuse-xpkg(inx,1)).gt.s%bdytol) iwarn=iwarn+1
                  zuse=xpkg(inx,1)
               endif
               zndx=ONE+zhavi*(zuse-xpkg(1,1))
               ix(i)=min(inxm1,int(zndx))
               dxn(i)=min(ONE,(zndx-ix(i)))
               hx(i)=zhav
               hxi(i)=zhavi
            enddo

         else
            ! periodic evenly spaced grid
            do i=1,ivec
               if((xuse(i).lt.xpkg(1,1)).or.(xuse(i).gt.xpkg(inx,1))) then
                  zuse=mod((xuse(i)-xpkg(1,1)),zspan)
                  if(zuse.lt.ZERO) zuse=zuse+zspan
                  zuse=zuse+xpkg(1,1)
                  zuse=max(xpkg(1,1),min(xpkg(inx,1),zuse))
               else
                  zuse=xuse(i)
               endif
               zndx=ONE+zhavi*(zuse-xpkg(1,1))
               ix(i)=min(inxm1,int(zndx))
               dxn(i)=min(ONE,(zndx-ix(i)))
               hx(i)=zhav
               hxi(i)=zhavi
            enddo

         endif
      endif
                     
    end subroutine xplookup

    subroutine xp_reverse(xin,ireverse,xout,x1,x2)

      ! private utility: if ireverse=F just copy; if ireverse=T, bring
      ! xin(:) into range x1 to x2, then reverse within this range

      real*8, dimension(:), intent(in) :: xin
      logical, intent(in) :: ireverse
      real*8, dimension(:), intent(out) :: xout

      real*8 :: x1,x2

      !------------------------
      integer :: i,isize
      real*8 :: zspan
      real*8, parameter :: ZERO=0.0d0
      !------------------------

      isize=size(xin)  ! size(xout) should be the same (not checked).

      if(.not.ireverse) then
         xout(1:isize) = xin(1:isize)
      else
         zspan=x2-x1
         do i=1,isize

            ! (a) bring in range

            if((xin(i).lt.x1).or.(xin(i).gt.x2)) then
               xout(i)=mod((xin(i)-x1),zspan)
               if(xout(i).lt.ZERO) xout(i)=xout(i)+zspan
               xout(i)=xout(i)+x1
            else
               xout(i)=xin(i)
            endif

            ! (b) reverse

            xout(i)=x1+(x2-xout(i))
            xout(i)=max(x1,min(x2,xout(i)))
         enddo
      endif

    end subroutine xp_reverse

    subroutine xplasma_evala(s,ans,kk,xx)

      !  profile evaluation utility -- all types & dimensionalities
      !  just one profile

      type (xplasma), pointer :: s
      real*8, intent(out), dimension(:) :: ans   ! answers returned
      integer, intent(in), dimension(:) :: kk    ! control packet
      type(xpeval), dimension(:), target :: xx   ! grid lookup info

      !-----------------------------
      integer :: inx,inprofs,j,ii
      integer :: iadr,ispline,irank,ideriv,ict(10),hspline(3)

      integer :: ideriv1,ideriv2,ideriv3
      integer :: inx1,inx2,inx3,ivec
      integer :: inxa(3),icoeff
      integer, dimension(:), pointer :: ix1,ix2,ix3
      real*8, dimension(:), pointer :: fspl
      real*8, dimension(:), pointer :: dxn1,dxn2,dxn3
      real*8, dimension(:), pointer :: hx1,hx2,hx3
      real*8, dimension(:), pointer :: hxi1,hxi2,hxi3
      logical :: eval2
      !-----------------------------

      !  evaluate profiles

      inx = size(ans)

      iadr = abs(kk(1))             ! profile data address
      ispline = s%profs(iadr)%kspline ! interpolation type
      irank = s%profs(iadr)%rank      ! profile rank

      fspl => s%profs(iadr)%buf

      eval2 = .TRUE.

      ideriv1 = xx(kk(2))%iderivs(1)
      inx1 = xx(kk(2))%inx
      ix1 => xx(kk(2))%ix
      dxn1=> xx(kk(2))%dxn
      hx1 => xx(kk(2))%hx
      hxi1=> xx(kk(2))%hxi

      if(irank.ge.2) then
         ideriv2 = xx(kk(3))%iderivs(1)
         inx2 = xx(kk(3))%inx
         ix2 => xx(kk(3))%ix
         dxn2=> xx(kk(3))%dxn
         hx2 => xx(kk(3))%hx
         hxi2=> xx(kk(3))%hxi
         if(xx(kk(3))%inum.gt.1) eval2=.FALSE.
      endif

      if(irank.ge.3) then
         ideriv3 = xx(kk(4))%iderivs(1)
         inx3 = xx(kk(4))%inx
         ix3 => xx(kk(4))%ix
         dxn3=> xx(kk(4))%dxn
         hx3 => xx(kk(4))%hx
         hxi3=> xx(kk(4))%hxi
         if(xx(kk(4))%inum.gt.1) eval2=.FALSE.
      endif

      ict = 0

      if(irank.eq.1) then
         !  eval rank 1 profile 

         if(ispline.eq.-1) then

            do ivec=1,inx
               ans(ivec)=fspl(ix1(ivec))
            enddo

         else if(ispline.eq.0) then
               
            !  1d piecewise linear interpolation

            call xplasma_build_ict1(ispline,ideriv1,ict)

            call r8pc1fcn(ict,inx,inx,ans, &
                 ix1,dxn1,hx1,hxi1, &
                 fspl,inx1)

         else if(ispline.eq.1) then
               
            !  1d cubic Hermite interpolation

            call xplasma_build_ict1(ispline,ideriv1,ict)

            call r8herm1fcn(ict,inx,inx,ans, &
                 ix1,dxn1,hx1,hxi1, &
                 fspl,inx1)

         else if(ispline.eq.2) then
               
            !  1d cubic Spline interpolation

            call xplasma_build_ict1(ispline,ideriv1,ict)

            call r8fvspline(ict,inx,inx,ans, &
                 ix1,dxn1,hx1,hxi1, &
                 fspl,inx1)

         endif

      else if(irank.eq.2) then
         !  eval rank 2 profile 

         if(ispline.eq.-1) then

            do ivec=1,inx
               ans(ivec)=fspl(ix1(ivec)+(ix2(ivec)-1)*(inx1-1))
            enddo

         else if(ispline.eq.0) then
               
            !  2d piecewise linear interpolation

            call xplasma_build_ict2(ispline,ideriv1,ideriv2,ict)

            if(.not.eval2) then
               call r8pc2fcn(ict,inx,inx,ans, &
                    ix1,ix2,dxn1,dxn2,hx1,hxi1,hx2,hxi2, &
                    fspl,inx1,inx2)
            else
               call r8pc2fcnx(ict,inx,inx,ans, &
                    ix1,ix2(1),dxn1,dxn2(1),hx1,hxi1,hx2(1),hxi2(1), &
                    fspl,inx1,inx2)
            endif

         else if(ispline.eq.1) then
               
            !  2d cubic Hermite interpolation

            call xplasma_build_ict2(ispline,ideriv1,ideriv2,ict)

            if(.not.eval2) then
               call r8herm2fcn(ict,inx,inx,ans, &
                    ix1,ix2,dxn1,dxn2,hx1,hxi1,hx2,hxi2, &
                    fspl,inx1,inx2)
            else
               call r8herm2fcnx(ict,inx,inx,ans, &
                    ix1,ix2(1),dxn1,dxn2(1),hx1,hxi1,hx2(1),hxi2(1), &
                    fspl,inx1,inx2)
            endif

         else if(ispline.eq.2) then
               
            !  2d cubic spline interpolation

            call xplasma_build_ict2(ispline,ideriv1,ideriv2,ict)

            if(.not.eval2) then
               call r8fvbicub(ict,inx,inx,ans, &
                    ix1,ix2,dxn1,dxn2,hx1,hxi1,hx2,hxi2, &
                    fspl,inx1,inx2)
            else
               call r8fvbicubx(ict,inx,inx,ans, &
                    ix1,ix2(1),dxn1,dxn2(1),hx1,hxi1,hx2(1),hxi2(1), &
                    fspl,inx1,inx2)
            endif

         else

            !  Hybrid
            !  build ict as for a spline

            call xplasma_build_ict2(2,ideriv1,ideriv2,ict)

            call spline_type_decode(ispline,hspline,2)
            icoeff=1
            inxa(1)=inx1
            inxa(2)=inx2
            do ii=1,2
               if(hspline(ii).eq.-1) then
                  inxa(ii)=inxa(ii)-1
               else if(hspline(ii).gt.0) then
                  icoeff=icoeff*2
               endif
            enddo

            call r8fvintrp2d(ict,inx,inx,ans, &
                 ix1,ix2,dxn1,dxn2,hx1,hxi1,hx2,hxi2, &
                 hspline,fspl,icoeff,inxa(1),inxa(2))
         endif

      else if(irank.eq.3) then
         !  eval rank 3 profile 

         if(ispline.eq.-1) then

            do ivec=1,inx
               ans(ivec)=fspl(ix1(ivec) +(ix2(ivec)-1)*(inx1-1) + &
                    (ix3(ivec)-1)*(inx2-1)*(inx1-1))
            enddo

         else if(ispline.eq.0) then
               
            !  3d piecewise linear interpolation

            call xplasma_build_ict3(ispline,ideriv1,ideriv2,ideriv3,ict)

            if(.not.eval2) then
               call r8pc3fcn(ict,inx,inx,ans, &
                    ix1,ix2,ix3,dxn1,dxn2,dxn3, &
                    hx1,hxi1,hx2,hxi2,hx3,hxi3, &
                    fspl,inx1,inx2,inx3)
            else
               call r8pc3fcnx(ict,inx,inx,ans, &
                    ix1,ix2(1),ix3(1),dxn1,dxn2(1),dxn3(1), &
                    hx1,hxi1,hx2(1),hxi2(1),hx3(1),hxi3(1), &
                    fspl,inx1,inx2,inx3)
            endif

         else if(ispline.eq.1) then
               
            !  3d cubic Hermite interpolation
            
            call xplasma_build_ict3(ispline,ideriv1,ideriv2,ideriv3,ict)

            if(.not.eval2) then
               call r8herm3fcn(ict,inx,inx,ans, &
                    ix1,ix2,ix3,dxn1,dxn2,dxn3, &
                    hx1,hxi1,hx2,hxi2,hx3,hxi3, &
                    fspl,inx1,inx2,inx3)
            else
               call r8herm3fcnx(ict,inx,inx,ans, &
                    ix1,ix2(1),ix3(1),dxn1,dxn2(1),dxn3(1), &
                    hx1,hxi1,hx2(1),hxi2(1),hx3(1),hxi3(1), &
                    fspl,inx1,inx2,inx3)
            endif

         else if(ispline.eq.2) then
               
            !  3d cubic Spline interpolation

            call xplasma_build_ict3(ispline,ideriv1,ideriv2,ideriv3,ict)

            if(.not.eval2) then
               call r8fvtricub(ict,inx,inx,ans, &
                    ix1,ix2,ix3,dxn1,dxn2,dxn3, &
                    hx1,hxi1,hx2,hxi2,hx3,hxi3, &
                    fspl,inx1,inx2,inx3)
            else
               call r8fvtricubx(ict,inx,inx,ans, &
                    ix1,ix2(1),ix3(1),dxn1,dxn2(1),dxn3(1), &
                    hx1,hxi1,hx2(1),hxi2(1),hx3(1),hxi3(1), &
                    fspl,inx1,inx2,inx3)
            endif


         else

            !  Hybrid
            !  build ict as for a spline

            call xplasma_build_ict3(2,ideriv1,ideriv2,ideriv3,ict)

            call spline_type_decode(ispline,hspline,3)
            icoeff=1
            inxa(1)=inx1
            inxa(2)=inx2
            inxa(3)=inx3
            do ii=1,3
               if(hspline(ii).eq.-1) then
                  inxa(ii)=inxa(ii)-1
               else if(hspline(ii).gt.0) then
                  icoeff=icoeff*2
               endif
            enddo

            call r8fvintrp3d(ict,inx,inx,ans, &
                 ix1,ix2,ix3,dxn1,dxn2,dxn3, &
                 hx1,hxi1,hx2,hxi2,hx3,hxi3, &
                 hspline,fspl,icoeff,inxa(1),inxa(2),inxa(3))
         endif

      endif
         
    end subroutine xplasma_evala

    subroutine xplasma_eval(s,ans,kk,xx)

      !  profile evaluation utility -- all types & dimensionalities

      type (xplasma), pointer :: s
      real*8, intent(out), dimension(:,:) :: ans   ! answers returned
      integer, intent(in), dimension(:,:) :: kk    ! control packet
      type(xpeval), dimension(:), intent(in) :: xx ! grid lookup info

      !-----------------------------
      integer :: inx,inprofs,i,j
      integer :: iadr,ispline,irank,ict(10)

      integer :: ideriv1,ideriv2,ideriv3
      integer :: inx1,inx2,inx3,ivec,ii,icoeff,inxa(3),hspline(3)
      real*8, dimension(:), pointer :: fspl
      integer, dimension(:), pointer :: ix1,ix2,ix3
      real*8, dimension(:), pointer :: dxn1,dxn2,dxn3
      real*8, dimension(:), pointer :: hx1,hx2,hx3
      real*8, dimension(:), pointer :: hxi1,hxi2,hxi3
      !-----------------------------

      !  evaluate profiles

      inx = size(ans,1)
      inprofs = size(ans,2)

      do i=1,inprofs
         iadr = abs(kk(1,i))             ! profile data address
         ispline = s%profs(iadr)%kspline ! interpolation type
         irank = s%profs(iadr)%rank      ! profile rank

         fspl => s%profs(iadr)%buf

         ideriv1 = xx(kk(2,i))%iderivs(i)
         inx1 = xx(kk(2,i))%inx
         ix1 => xx(kk(2,i))%ix
         dxn1=> xx(kk(2,i))%dxn
         hx1 => xx(kk(2,i))%hx
         hxi1=> xx(kk(2,i))%hxi

         if(irank.ge.2) then
            ideriv2 = xx(kk(3,i))%iderivs(i)
            inx2 = xx(kk(3,i))%inx
            ix2 => xx(kk(3,i))%ix
            dxn2=> xx(kk(3,i))%dxn
            hx2 => xx(kk(3,i))%hx
            hxi2=> xx(kk(3,i))%hxi
         endif

         if(irank.ge.3) then
            ideriv3 = xx(kk(4,i))%iderivs(i)
            inx3 = xx(kk(4,i))%inx
            ix3 => xx(kk(4,i))%ix
            dxn3=> xx(kk(4,i))%dxn
            hx3 => xx(kk(4,i))%hx
            hxi3=> xx(kk(4,i))%hxi
         endif

         ict = 0

         if(irank.eq.1) then
            !  eval rank 1 profile 

            if(ispline.eq.-1) then

               do ivec=1,inx
                  ans(ivec,i)=fspl(ix1(ivec))
               enddo

            else if(ispline.eq.0) then
               
               !  1d piecewise linear interpolation

               call xplasma_build_ict1(ispline,ideriv1,ict)

               call r8pc1fcn(ict,inx,inx,ans(1,i), &
                    ix1,dxn1,hx1,hxi1, &
                    fspl,inx1)

            else if(ispline.eq.1) then
               
               !  1d cubic Hermite interpolation

               call xplasma_build_ict1(ispline,ideriv1,ict)

               call r8herm1fcn(ict,inx,inx,ans(1,i), &
                    ix1,dxn1,hx1,hxi1, &
                    fspl,inx1)

            else if(ispline.eq.2) then
               
               !  1d cubic Spline interpolation

               call xplasma_build_ict1(ispline,ideriv1,ict)

               call r8fvspline(ict,inx,inx,ans(1,i), &
                    ix1,dxn1,hx1,hxi1, &
                    fspl,inx1)

            endif

         else if(irank.eq.2) then
            !  eval rank 2 profile 

            if(ispline.eq.-1) then

               do ivec=1,inx
                  ans(ivec,i)=fspl(ix1(ivec)+(ix2(ivec)-1)*(inx1-1))
               enddo

            else if(ispline.eq.0) then
               
               !  2d piecewise linear interpolation

               call xplasma_build_ict2(ispline,ideriv1,ideriv2,ict)

               call r8pc2fcn(ict,inx,inx,ans(1,i), &
                    ix1,ix2,dxn1,dxn2,hx1,hxi1,hx2,hxi2, &
                    fspl,inx1,inx2)

            else if(ispline.eq.1) then
               
               !  2d cubic Hermite interpolation

               call xplasma_build_ict2(ispline,ideriv1,ideriv2,ict)

               call r8herm2fcn(ict,inx,inx,ans(1,i), &
                    ix1,ix2,dxn1,dxn2,hx1,hxi1,hx2,hxi2, &
                    fspl,inx1,inx2)

            else if(ispline.eq.2) then
               
               !  2d cubic spline interpolation

               call xplasma_build_ict2(ispline,ideriv1,ideriv2,ict)

               call r8fvbicub(ict,inx,inx,ans(1,i), &
                    ix1,ix2,dxn1,dxn2,hx1,hxi1,hx2,hxi2, &
                    fspl,inx1,inx2)

            else

               !  Hybrid
               !  build ict as for a spline

               call xplasma_build_ict2(2,ideriv1,ideriv2,ict)

               call spline_type_decode(ispline,hspline,2)
               icoeff=1
               inxa(1)=inx1
               inxa(2)=inx2
               do ii=1,2
                  if(hspline(ii).eq.-1) then
                     inxa(ii)=inxa(ii)-1
                  else if(hspline(ii).gt.0) then
                     icoeff=icoeff*2
                  endif
               enddo

               call r8fvintrp2d(ict,inx,inx,ans(1,i), &
                    ix1,ix2,dxn1,dxn2,hx1,hxi1,hx2,hxi2, &
                    hspline,fspl,icoeff,inxa(1),inxa(2))
            endif

         else if(irank.eq.3) then
            !  eval rank 3 profile 

            if(ispline.eq.-1) then

               do ivec=1,inx
                  ans(ivec,i)=fspl(ix1(ivec) +(ix2(ivec)-1)*(inx1-1) + &
                       (ix3(ivec)-1)*(inx2-1)*(inx1-1))
               enddo

            else if(ispline.eq.0) then
               
               !  3d piecewise linear interpolation

               call xplasma_build_ict3(ispline,ideriv1,ideriv2,ideriv3,ict)

               call r8pc3fcn(ict,inx,inx,ans(1,i), &
                    ix1,ix2,ix3,dxn1,dxn2,dxn3, &
                    hx1,hxi1,hx2,hxi2,hx3,hxi3, &
                    fspl,inx1,inx2,inx3)

            else if(ispline.eq.1) then
               
               !  3d cubic Hermite interpolation

               call xplasma_build_ict3(ispline,ideriv1,ideriv2,ideriv3,ict)

               call r8herm3fcn(ict,inx,inx,ans(1,i), &
                    ix1,ix2,ix3,dxn1,dxn2,dxn3, &
                    hx1,hxi1,hx2,hxi2,hx3,hxi3, &
                    fspl,inx1,inx2,inx3)

            else if(ispline.eq.2) then
               
               !  3d cubic Spline interpolation

               call xplasma_build_ict3(ispline,ideriv1,ideriv2,ideriv3,ict)

               call r8fvtricub(ict,inx,inx,ans(1,i), &
                    ix1,ix2,ix3,dxn1,dxn2,dxn3, &
                    hx1,hxi1,hx2,hxi2,hx3,hxi3, &
                    fspl,inx1,inx2,inx3)

            else

               !  Hybrid
               !  build ict as for a spline

               call xplasma_build_ict3(2,ideriv1,ideriv2,ideriv3,ict)

               call spline_type_decode(ispline,hspline,3)
               icoeff=1
               inxa(1)=inx1
               inxa(2)=inx2
               inxa(3)=inx3
               do ii=1,3
                  if(hspline(ii).eq.-1) then
                     inxa(ii)=inxa(ii)-1
                  else if(hspline(ii).gt.0) then
                     icoeff=icoeff*2
                  endif
               enddo

               call r8fvintrp3d(ict,inx,inx,ans(1,i), &
                    ix1,ix2,ix3,dxn1,dxn2,dxn3, &
                    hx1,hxi1,hx2,hxi2,hx3,hxi3, &
                    hspline,fspl,icoeff,inxa(1),inxa(2),inxa(3))
            endif

         endif

      enddo
         
    end subroutine xplasma_eval

    subroutine xplasma_build_ict1(ispline,ideriv,ict)
      ! utility: build pspline ict control -- 1d eval

      integer, intent(in) :: ispline
      integer, intent(in) :: ideriv
      integer, dimension(:), intent(inout) :: ict

      if(ispline.le.1) then
         if(ideriv.eq.0) then
            ict(1)=1
         else
            ict(2)=1
         endif
      else
         if(ideriv.eq.0) then
            ict(1)=1
         else if(ideriv.eq.1) then
            ict(2)=1
         else if(ideriv.eq.2) then
            ict(3)=1
         else
            ict(1)=3
         endif
      endif

    end subroutine xplasma_build_ict1

    subroutine xplasma_build_ict2(ispline,ideriv1,ideriv2,ict)
      ! utility: build pspline ict control -- 1d eval

      integer, intent(in) :: ispline
      integer, intent(in) :: ideriv1,ideriv2
      integer, dimension(:), intent(inout) :: ict

      integer :: i1,i2,imark,isum,iii

      if(ispline.le.1) then
         if((ideriv1.eq.0).and.(ideriv2.eq.0)) then
            ict(1)=1
         else if((ideriv1.eq.1).and.(ideriv2.eq.0)) then
            ict(2)=1
         else if((ideriv1.eq.0).and.(ideriv2.eq.1)) then
            ict(3)=1
         else
            ict(4)=1
         endif
      else
         i1 = ideriv1
         i2 = ideriv2

         !  code borrowed from pspline/ezspline...

         !  (private utility for ezspline derivative2 subroutines)
         !  make ict(1:6) array
         !  for higher derivatives; d[i1+i2]f/dx[i1]dy[i2]
         !  expecting i1 & i2 in range [0:3] (NOT CHECKED)


         !  this generates the control argument needed by evbicub & similar
         !  routines...
         !----------------------

         isum = i1+i2
         ict(1)=isum

         imark=0

         if(isum.eq.0) then
            ict(1:6) = (/1, 0, 0, 0, 0, 0 /) ! seek f @ (p1, p2, p3)
         else if(isum.eq.1) then
            if(i1.eq.1) then
               ict(1:6) = (/0, 1, 0, 0, 0, 0 /) ! df/dx
            else
               ict(1:6) = (/0, 0, 1, 0, 0, 0 /) ! df/dy
            endif
         else if(isum.eq.2) then
            if(i1.eq.2) then
               ict(1:6) = (/0, 0, 0, 1, 0, 0 /) ! d2f/dx2
            else if(i2.eq.2) then
               ict(1:6) = (/0, 0, 0, 0, 1, 0 /) ! d2f/dy2
            else
               ict(1:6) = (/0, 0, 0, 0, 0, 1 /) ! d2f/dxdy
            endif
         else if(isum.eq.3) then
            if(i1.eq.3) then
               imark=2  ! fxxx
            else if(i1.eq.2) then
               imark=3  ! fxxy
            else if(i1.eq.1) then
               imark=4  ! fxyy
            else
               imark=5  ! fyyy
            endif
         else if(isum.eq.4) then
            if(i1.eq.3) then
               imark=2  ! fxxxy
            else if(i2.eq.3) then
               imark=4  ! fxyyy
            else
               imark=3  ! fxxyy
            endif
         else if(isum.eq.5) then
            if(i1.eq.3) then
               imark=2  ! fxxxyy
            else if(i2.eq.3) then
               imark=3  ! fxxyyy
            endif
         endif

         !  isum=6 --> fxxxyyy

         if(isum.gt.2) then
            do iii=2,6
               if(iii.eq.imark) then
                  ict(iii)=1
               else
                  ict(iii)=0
               endif
            enddo
         endif

      endif

    end subroutine xplasma_build_ict2

    subroutine xplasma_build_ict3(ispline,ideriv1,ideriv2,ideriv3,ict)
      ! utility: build pspline ict control -- 1d eval

      integer, intent(in) :: ispline
      integer, intent(in) :: ideriv1,ideriv2,ideriv3
      integer, dimension(:), intent(inout) :: ict

      integer :: i1,i2,i3,imark,isum,iii

      if(ispline.le.1) then
         !
         !  3d Hermite or pc lin derivative code
         !
         isum=ideriv1+ideriv2+ideriv3
         if(isum.eq.0) then
            ict(1)=1
         else if(isum.eq.3) then
            ict(8)=1
         else if(isum.eq.1) then
            if(ideriv1.eq.1) then
               ict(2)=1
            else if(ideriv2.eq.1) then
               ict(3)=1
            else
               ict(4)=1
            endif
         else
            if(ideriv3.eq.0) then
               ict(5)=1
            else if(ideriv2.eq.0) then
               ict(6)=1
            else
               ict(7)=1
            endif
         endif
      else
         !
         !  3d spline derivative code
         !  ... from ezspline_mod ezmake_ict3 originally...

         i1 = ideriv1
         i2 = ideriv2
         i3 = ideriv3

         isum = i1+i2+i3
         if(max(i1,i2,i3).eq.3) then
            isum=-isum
         endif
         ict(1)=isum

         imark=0

         if(isum.eq.0) then
            ict = (/1, 0, 0, 0, 0, 0, 0, 0, 0, 0 /) ! seek f @ (p1, p2, p3)

         else if(isum.eq.1) then
            !  1st derivatives
            if(i1.eq.1) then
               ict = (/0, 1, 0, 0, 0, 0, 0, 0, 0, 0 /)  ! df/dx
            else if(i2.eq.1) then
               ict = (/0, 0, 1, 0, 0, 0, 0, 0, 0, 0 /)  ! df/dy
            else
               ict = (/0, 0, 0, 1, 0, 0, 0, 0, 0, 0 /)  ! df/dz
            endif

         else if(isum.eq.2) then
            !  2nd derivatives-- legacy ordering; x-precedence ordering for all
            !  higher derivatives...

            if(i1.eq.2) then
               ict = (/0, 0, 0, 0, 1, 0, 0, 0, 0, 0 /)  ! d2f/dx2
            else if(i2.eq.2) then
               ict = (/0, 0, 0, 0, 0, 1, 0, 0, 0, 0 /)  ! d2f/dy2
            else if(i3.eq.2) then
               ict = (/0, 0, 0, 0, 0, 0, 1, 0, 0, 0 /)  ! d2f/dz2
            else if(i3.eq.0) then
               ict = (/0, 0, 0, 0, 0, 0, 0, 1, 0, 0 /)  ! d2f/dxdy
            else if(i2.eq.0) then
               ict = (/0, 0, 0, 0, 0, 0, 0, 0, 1, 0 /)  ! d2f/dxdz
            else
               ict = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 1 /)  ! d2f/dydz
            endif
        
         else if(isum.eq.3) then
            !  3rd derivative, continuous: max(i1,i2,i3)<3
            if(i1.eq.2) then
               if(i2.eq.1) then
                  imark=2     ! fxxy
               else
                  imark=3     ! fxxz
               endif
            else if(i1.eq.1) then
               if(i2.eq.2) then
                  imark=4     ! fxyy
               else if(i2.eq.1) then
                  imark=5     ! fxyz
               else
                  imark=6     ! fxzz
               endif
            else
               if(i2.eq.2) then
                  imark=7     ! fyyz
               else
                  imark=8     ! fyzz
               endif
            endif

         else if(isum.eq.-3) then
            !  3rd derivative
            if(i1.eq.3) then
               imark=2        ! fxxx
            else if(i2.eq.3) then
               imark=3        ! fyyy
            else if(i3.eq.3) then
               imark=4        ! fzzz
            endif
            
         else if(isum.eq.4) then
            !  4th derivative, continuous: max(i1,i2,i3)<3
            if(i1.eq.2) then
               if(i2.eq.2) then
                  imark=2     ! fxxyy
               else if(i2.eq.1) then
                  imark=3     ! fxxyz
               else
                  imark=4     ! fxxzz
               endif
            else if(i1.eq.1) then
               if(i2.eq.2) then
                  imark=5     ! fxyyz
               else
                  imark=6     ! fxyzz
               endif
            else
               imark=7        ! fyyzz
            endif

         else if(isum.eq.-4) then
            !  4th derivative
            if(i1.eq.3) then
               if(i2.eq.1) then
                  imark=2     ! fxxxy
               else
                  imark=3     ! fxxxz
               endif
            else if(i1.eq.1) then
               if(i2.eq.3) then
                  imark=4     ! fxyyy
               else
                  imark=5     ! fxzzz
               endif
            else
               if(i2.eq.3) then
                  imark=6     ! fyyyz
               else
                  imark=7     ! fyzzz
               endif
            endif

         else if(isum.eq.5) then
            !  5th derivative, continuous: max(i1,i2,i3)<3
            if(i3.eq.1) then
               imark=2     ! fxxyyz
            else if(i2.eq.1) then
               imark=3     ! fxxyzz
            else
               imark=4     ! fxyyzz
            endif

         else if(isum.eq.-5) then
            !  5th derivative
            if(i1.eq.3) then
               if(i2.eq.2) then
                  imark=2  ! fxxxyy
               else if(i2.eq.1) then
                  imark=3  ! fxxxyz
               else
                  imark=4  ! fxxxzz
               endif
            else if(i1.eq.2) then
               if(i2.eq.3) then
                  imark=5  ! fxxyyy
               else
                  imark=6  ! fxxzzz
               endif
            else if(i1.eq.1) then
               if(i2.eq.3) then
                  imark=7  ! fxyyyz
               else
                  imark=8  ! fxyzzz
               endif
            else
               if(i2.eq.3) then
                  imark=9  ! fyyyzz
               else
                  imark=10 ! fyyzzz
               endif
            endif

            !  isum=6 --> fxxyyzz  (i1=i2=i3=2)
         else if(isum.eq.-6) then
            !  6th derivative
            if(i1.eq.3) then
               if(i2.eq.3) then
                  imark=2  ! fxxxyyy
               else if(i2.eq.2) then
                  imark=3  ! fxxxyyz
               else if(i2.eq.1) then
                  imark=4  ! fxxxyzz
               else
                  imark=5  ! fxxxzzz
               endif
            else if(i1.eq.2) then
               if(i2.eq.3) then
                  imark=6  ! fxxyyyz
               else if(i2.eq.1) then
                  imark=7  ! fxxyzzz
               endif
            else if(i1.eq.1) then
               if(i2.eq.3) then
                  imark=8  ! fxyyyzz
               else
                  imark=9  ! fxyyzzz
               endif
            else
               imark=10    ! fyyyzzz
            endif

            !  isum=7 not possible
         else if(isum.eq.-7) then
            !  7th derivative
            if(i1.eq.3) then
               if(i2.eq.3) then
                  imark=2  ! fxxxyyyz
               else if(i2.eq.2) then
                  imark=3  ! fxxxyyzz
               else
                  imark=4  ! fxxxyzzz
               endif
            else if(i1.eq.2) then
               if(i2.eq.3) then
                  imark=5  ! fxxyyyzz
               else
                  imark=6  ! fxxyyzzz
               endif
            else
               imark=7     ! fxyyyzzz
            endif

            !  isum=8 not possible
         else if(isum.eq.-8) then
            !  8th derivative
            if(i3.eq.2) then 
               imark=2  ! fxxxyyyzz
            else if(i2.eq.2) then
               imark=3  ! fxxxyyzzz
            else
               imark=4  ! fxxyyyzzz
            endif

            !  isum=9 not possible
            !  isum=-9 --> fxxxyyyzzz
            
         endif

         if(abs(isum).gt.2) then
            do iii=2,10
               if(iii.eq.imark) then
                  ict(iii)=1
               else
                  ict(iii)=0
               endif
            enddo
         endif

      endif
    end subroutine xplasma_build_ict3

    subroutine xplasma_errck0(s,id,ier)

      !  perform common error checks for item access methods

      type (xplasma), pointer :: s
      integer, intent(in) :: id
      integer, intent(out) :: ier

      character*128 msgbuf

      ier = 0

      if(.not.associated(s)) then
         ier=101
      else if(s%nguard.eq.0) then
         ier=101
      endif
      if(ier.ne.0) return

      if((id.lt.1).or.(id.gt.s%nitems)) then
         write(msgbuf,*) 'id argument value out of range 1 to ',s%nitems, &
              ': ',id
         call xplasma_errmsg_append(s,msgbuf)
         ier=105 
         return
      endif

    end subroutine xplasma_errck0

    !========================================
    subroutine xplasma_errck0_list(s,id,ier)

      !  perform common error checks for list item access methods

      type (xplasma), pointer :: s
      integer, intent(in) :: id
      integer, intent(out) :: ier

      character*128 msgbuf

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      if(s%dict(id)%dtype.ne.xplasma_listType) then
         msgbuf = 'item name: "'//trim(s%dict(id)%name)//'".'
         call xplasma_errmsg_append(s,msgbuf)
         ier=200
      else
         ier=0
      endif

    end subroutine xplasma_errck0_list

    !========================================
    subroutine xplasma_errck0_coord(s,id,ier)

      !  perform common error checks for coord item access methods

      type (xplasma), pointer :: s
      integer, intent(in) :: id
      integer, intent(out) :: ier

      character*128 msgbuf

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      if(s%dict(id)%dtype.ne.xplasma_coordType) then
         msgbuf = 'item name: "'//trim(s%dict(id)%name)//'".'
         call xplasma_errmsg_append(s,msgbuf)
         ier=400
      else
         ier=0
      endif

    end subroutine xplasma_errck0_coord

    !========================================
    subroutine xplasma_errck0_grid(s,id,ier)

      !  perform common error checks for grid item access methods

      type (xplasma), pointer :: s
      integer, intent(in) :: id
      integer, intent(out) :: ier

      character*128 msgbuf

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      if(s%dict(id)%dtype.ne.xplasma_gridType) then
         msgbuf = 'item name: "'//trim(s%dict(id)%name)//'".'
         call xplasma_errmsg_append(s,msgbuf)
         ier=300
      else
         ier=0
      endif

    end subroutine xplasma_errck0_grid

    !========================================
    subroutine xplasma_errck0_prof(s,id,ier)

      !  perform common error checks for prof item access methods

      type (xplasma), pointer :: s
      integer, intent(in) :: id
      integer, intent(out) :: ier

      integer :: iadr

      character*128 msgbuf

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      if(s%dict(id)%dtype.ne.xplasma_profType) then
         msgbuf = 'item name: "'//trim(s%dict(id)%name)//'".'
         call xplasma_errmsg_append(s,msgbuf)
         ier=500
      else 
         ier=0
      endif

    end subroutine xplasma_errck0_prof

    subroutine xplasma_errck1_prof(s,id,ier)

      !  perform common error checks for prof item access methods

      type (xplasma), pointer :: s
      integer, intent(in) :: id
      integer, intent(out) :: ier

      integer :: iadr
      character*128 msgbuf

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      if(s%dict(id)%dtype.ne.xplasma_profType) then
         msgbuf = 'item name: "'//trim(s%dict(id)%name)//'".'
         call xplasma_errmsg_append(s,msgbuf)
         ier=500
      else 
         iadr = s%dict(id)%dindex
         if(s%profs(iadr)%rank.eq.0) then
            ier=503
         else
            ier=0
         endif
      endif

    end subroutine xplasma_errck1_prof

    subroutine xplasma_errck1_blkbx(s,id,ier)

      !  perform common error checks for blkbx item access methods

      type (xplasma), pointer :: s
      integer, intent(in) :: id
      integer, intent(out) :: ier

      integer :: iadr
      character*128 msgbuf

      call xplasma_errck0(s,id,ier)
      if(ier.ne.0) return

      if(s%dict(id)%dtype.ne.xplasma_blkbxType) then
         msgbuf = 'item name: "'//trim(s%dict(id)%name)//'".'
         call xplasma_errmsg_append(s,msgbuf)
         ier=606
      else 
         iadr = s%dict(id)%dindex
         if(s%blkbxs(iadr)%type.eq.0) then
            ier=607
         else
            ier=0
         endif
      endif

    end subroutine xplasma_errck1_blkbx

    subroutine xplasma_errck0_xinteg(s,icoordx,x,jcoordx,ier)

      !  check xplasma_create/augment_integ coordinate argument...

      type (xplasma), pointer :: s
      integer, intent(in) :: icoordx
      real*8, intent(in) :: x(:)
      integer, intent(out) :: jcoordx
      integer, intent(out) :: ier

      !------------------------
      integer :: i,iadr,isize
      real*8 :: xmin,xmax
      character*128 :: zmsg
      logical :: periodic

      integer :: k_tmp
      !------------------------

      jcoordx=0

      call xplasma_errck0(s,icoordx,ier)
      if(ier.ne.0) return

      k_tmp = s%dict(icoordx)%dtype
      if(s%dict(icoordx)%dtype.eq.xplasma_coordType) then
         jcoordx=icoordx
      else if( k_tmp .eq.xplasma_gridType) then
         iadr=s%dict(icoordx)%dindex
         jcoordx=s%grids(iadr)%coord
      else
         ier=600
         call xplasma_errmsg_append(s, &
              '  not a valid integration coordinate: "'//trim(s%dict(icoordx)%name)//'".')
      endif
      if(ier.ne.0) return

      isize=size(x)
      if(isize.lt.2) then
         ier=603
         if(isize.le.0) call xplasma_errmsg_append(s,'  zero length array was passed.')
         if(isize.le.1) call xplasma_errmsg_append(s,'  size of x array = 1.')
      else
         do i=1,isize-1
            if(x(i+1).le.x(i)) then
               ier=603
               call xplasma_errmsg_append(s,'  x points out of order.')
               exit
            endif
         enddo
      endif
      if(ier.ne.0) return

      call xplasma_coord_info(s,jcoordx,ier, &
           xmin=xmin,xmax=xmax,periodic=periodic)

      if(.not.periodic) then
         if((x(1).lt.xmin-s%bdytol).or.(x(size(x)).gt.xmax+s%bdytol)) then
            ier=609
            zmsg=' '
            write(zmsg,1001) 'coordinate range:',xmin,xmax
            write(zmsg,1001) 'integration range requested:',x(1),x(size(x))
         endif
1001     format(3x,a,2(1x,1pe13.6))
      endif

    end subroutine xplasma_errck0_xinteg

    subroutine xplasma_xpprof_lblgen(s,id,str)

      !  generate a label string showing a profile and its independent
      !  coordinates...

      type (xplasma), pointer :: s
      integer, intent(in) :: id
      character*(*), intent(out) :: str

      !-----------------------
      integer :: irank,iadr,ilen,jlen,ir,icoord,igrids(maxRank)
      !-----------------------

      iadr = s%dict(id)%dindex

      irank= s%profs(iadr)%rank
      igrids = s%profs(iadr)%gridIds

      str=' '
      
      ilen=len(trim(s%dict(id)%name))
      if(irank.eq.0) then
         str(1:ilen+2)=trim(s%dict(id)%name)//'()'
         return
      else
         str(1:ilen+1)=trim(s%dict(id)%name)//'('
         jlen=ilen+1
      endif

      do ir=1,irank
         iadr=s%dict(igrids(ir))%dindex
         icoord=s%grids(iadr)%coord
         ilen=len(trim(s%dict(icoord)%name))
         if(ir.eq.irank) then
            str(jlen+1:jlen+ilen+1)=trim(s%dict(icoord)%name)//')'
         else
            str(jlen+1:jlen+ilen+1)=trim(s%dict(icoord)%name)//','
            jlen=jlen+ilen+1
         endif
      enddo

    end subroutine xplasma_xpprof_lblgen

    !========================================
    subroutine xplasma_errck0_list_adr(s,id,nelems,iadr,ier)

      !  look up list address;
      !  perform common error checks for list item access methods

      type (xplasma), pointer :: s
      integer, intent(in) :: id
      integer, intent(in) :: nelems ! expected list size
      integer, intent(out) :: iadr  ! return list address
      integer, intent(out) :: ier

      integer :: inum
      character*128 msgbuf

      call xplasma_errck0_list(s,id,ier)
      if(ier.ne.0) return

      iadr = s%dict(id)%dindex
      inum = s%lists(iadr)%size

      if(inum.ne.nelems) then
         msgbuf = 'item name: "'//trim(s%dict(id)%name)//'".'
         call xplasma_errmsg_append(s,msgbuf)
         ier=204
         return
      endif
    end subroutine xplasma_errck0_list_adr

    !========================================
    subroutine xplasma_author_check(s,id,ier)

      !  check author id for write access prio

      type (xplasma), pointer :: s
      integer, intent(in) :: id    ! dictionary item in question
      integer, intent(out) :: ier  ! status code 0= access approved

      character*128 msgbuf

      integer :: k_tmp

      ier=0
  
      k_tmp = s%write_enable(s%cur_author)
      if((s%dict(id)%author.ne.s%authors(s%cur_author)).and. &
           (s%authors(s%cur_author).ne.xplasma_root)) then
         ier=51
         write(msgbuf,*) ' item name: "',trim(s%dict(id)%name), &
              '"; item author: ',trim(s%dict(id)%author)
         call xplasma_errmsg_append(s,msgbuf)
         write(msgbuf,*) ' name of author attempting modification: "', &
              trim(s%authors(s%cur_author))
         call xplasma_errmsg_append(s,msgbuf)
         if(s%write_enable(s%cur_author).ne.1) then
            write(msgbuf,*) ' note: write_enable attribute not set.'
            call xplasma_errmsg_append(s,msgbuf)
         endif

      else if( k_tmp .ne.1) then
         ier=50
         return
      endif

    end subroutine xplasma_author_check

    !========================================
    subroutine xplasma_checkType(itype,ier)

      !  check that type id code is in valid range...

      integer, intent(in) :: itype   ! data type code
      integer, intent(out) :: ier    ! status code returned, 0=OK

      if((itype.lt.xplasma_minType).or.(itype.gt.xplasma_maxType)) then
         ier=104
      else
         ier=0
      endif
    end subroutine xplasma_checkType

    !========================================
    subroutine xplasma_checkName(zname,ier)

      !  check that item name is legal

      character*(*), intent(in) :: zname
      integer, intent(out) :: ier

      call xplasma_checkName_sub(zname,ier)

      if(ier.ne.0) then
         if((ier.gt.2).or.(ier.lt.0)) then
            ier=9999
         else
            if(ier.eq.1) ier=102  ! illegal character
            if(ier.eq.2) ier=103  ! name too long
         endif
      endif

    end subroutine xplasma_checkName

    !========================================
    subroutine xplasma_expand(s,ineed,ier)

      ! expand xplasma container to make new room for new item...

      type (xplasma), pointer :: s
      integer, intent(in) :: ineed  ! minimum size needed
      integer, intent(out) :: ier   ! status code returned (always 0=OK)

      !---------------
      integer :: inold,innew,i,isize
      integer, dimension(:), pointer :: itmp
      type (xplasma_item), dimension(:), pointer :: xtmp
      !---------------

      ier=0

      inold=s%nmax
      if(ineed.le.inold) return

      innew=inold
      do
         innew=max(8,2*innew)
         if(innew.ge.ineed) exit
      enddo

      s%nmax = innew

      if(inold.gt.0) then
         itmp => s%order
         xtmp => s%dict
      endif

      allocate(s%order(innew),s%dict(innew))
      if(inold.gt.0) then
         s%order(1:inold)=itmp(1:inold)
         do i=1,inold
            call xpltest_alloc(s%dict(i))
            s%dict(i)%name = xtmp(i)%name
            s%dict(i)%author = xtmp(i)%author
            s%dict(i)%dtype = xtmp(i)%dtype
            s%dict(i)%dindex = xtmp(i)%dindex
            if(allocated(xtmp(i)%label)) then
               isize=size(xtmp(i)%label)
               allocate(s%dict(i)%label(isize))
               s%dict(i)%label = xtmp(i)%label
               deallocate(xtmp(i)%label)
            endif
            if(allocated(xtmp(i)%units)) then
               isize=size(xtmp(i)%units)
               allocate(s%dict(i)%units(isize))
               s%dict(i)%units = xtmp(i)%units
               deallocate(xtmp(i)%units)
            endif
         enddo
         deallocate(itmp,xtmp)
      endif
      do i=inold+1,innew
         call xpltest_alloc(s%dict(i))
         s%order(i)=0
         s%dict(i)%name=' '
         s%dict(i)%author=' '
         s%dict(i)%dtype=0
         s%dict(i)%dindex=0
      enddo

    end subroutine xplasma_expand

    subroutine xpltest_alloc(zitem)
      ! private test routine -- allocate and deallocate allocatable element
      ! (here for debugging reasons -- dmc)

      type (xplasma_item) :: zitem

      integer :: luntty = 6

      allocate(zitem%label(1))
      zitem%label(1)=zitem%name(1:1)

      allocate(zitem%units(1))
      zitem%units(1)=zitem%name(2:2)

      if(zitem%name.eq.'?12345') then
         !  (this will not happen...)
         write(luntty,*) " xplasma_obj debug: label & units: ", &
              zitem%label(1),zitem%units(1)
      endif

      deallocate(zitem%label)
      deallocate(zitem%units)

    end subroutine xpltest_alloc

    !========================================
    subroutine xplasma_locate(s,zname,iloc,iexact)

      ! find item position in ordered list; if exact item is in list
      ! then return order location (iloc) and iexact=1.
      !   if item is not in list return order location (iloc) of first
      ! item that follows (uppercase zname) in alphabetic order.
      !   if item is not in list and it is last in alphabetic order
      ! then return iloc = s%nitems + 1

      type (xplasma), pointer :: s
      character*(*), intent(in) :: zname  ! name of item sought
      integer, intent(out) :: iloc        ! location (as defined above)
      integer, intent(out) :: iexact      ! =1 for exact match

      !--------------------------

      character*32 znamu,znamt

      integer :: inum,ia1,ia2
      integer :: imin,imax

      !--------------------------

      iexact=0

      inum=s%nitems
      if(inum.eq.0) then
         iloc=1
         return
      endif

      znamu=zname
      call uupper(znamu)

      iloc = max(1,inum/2)
      imax = inum
      imin = 1

      !  binary search...

      do
         ia1=s%order(iloc)
         znamt = s%dict(ia1)%name
         if(znamu.eq.znamt) then
            iexact=1  ! exact match at this iloc
            exit
         else if(znamu.gt.znamt) then
            if(iloc.eq.inum) then
               iloc=iloc+1  ! comes after the last one...
               exit
            else
               imin=iloc+1
               iloc=max(imin,(iloc+imax)/2)
               cycle
            endif
         else
            if(iloc.eq.1) then
               exit  ! comes before 1st one...
            else
               ia2=s%order(iloc-1)
               znamt=s%dict(ia2)%name
               if(znamu.eq.znamt) then
                  iloc=iloc-1
                  iexact=1
                  exit
               else if(znamu.gt.znamt) then
                  exit  ! comes before iloc but after iloc predecessor
               else
                  imax=iloc-1
                  iloc=min(imax,(iloc+imin)/2)
                  cycle
               endif
            endif
         endif
      enddo

    end subroutine xplasma_locate

    !========================================
    subroutine xplasma_lists_expand(s,ineed,ier)

      ! expand list of lists

      type (xplasma), pointer :: s
      integer, intent(in) :: ineed
      integer, intent(out) :: ier

      !--------------------/local
      integer :: inold,innew,i
      type (xplist), dimension(:), allocatable :: xtmp
      !--------------------

      ier=0

      inold = s%nlists
      if(inold.ge.ineed) return

      innew = inold
      do
         innew=max(50,2*innew)
         if(innew.ge.ineed) exit
      enddo

      if(inold.gt.0) then
         allocate(xtmp(inold))
         do i=1,inold
            call xplist_fullcopy(s%lists(i),xtmp(i))
            call xplist_free(s%lists(i))
         enddo
         deallocate(s%list_freelist)
         deallocate(s%lists)
      endif

      allocate(s%lists(innew),s%list_freelist(innew))
      if(inold.gt.0) then
         do i=1,inold
            call xplist_fullcopy(xtmp(i),s%lists(i))
            call xplist_free(xtmp(i))
            s%list_freelist(innew-i+1)=0
         enddo
         deallocate(xtmp)
      endif

      do i=inold+1,innew
         call xplist_init(s%lists(i))
      enddo

      s%nlists=innew
      s%nlist_free=innew-inold
      do i=1,innew-inold
         s%list_freelist(i)=innew+1-i
      enddo

    end subroutine xplasma_lists_expand

    !========================================
    subroutine xplasma_coords_expand(s,ineed,ier)

      ! expand list of coords

      type (xplasma), pointer :: s
      integer, intent(in) :: ineed
      integer, intent(out) :: ier

      !--------------------/local
      integer :: inold,innew,i
      type (xpcoord), dimension(:), allocatable :: xtmp
      !--------------------

      ier=0

      inold = s%ncoords
      if(inold.ge.ineed) return

      innew = inold
      do
         innew=max(5,2*innew)
         if(innew.ge.ineed) exit
      enddo

      if(inold.gt.0) then
         allocate(xtmp(inold))
         do i=1,inold
            call xpcoord_fullcopy(s%coords(i),xtmp(i))
            call xpcoord_free(s%coords(i))
         enddo
         deallocate(s%coord_freelist)
         deallocate(s%coords)
      endif

      allocate(s%coords(innew),s%coord_freelist(innew))
      if(inold.gt.0) then
         do i=1,inold
            call xpcoord_fullcopy(xtmp(i),s%coords(i))
            call xpcoord_free(xtmp(i))
            s%coord_freelist(innew-i+1)=0
         enddo
         deallocate(xtmp)
      endif

      do i=inold+1,innew
         call xpcoord_init(s%coords(i))
      enddo

      s%ncoords=innew
      s%ncoord_free=innew-inold
      do i=1,innew-inold
         s%coord_freelist(i)=innew+1-i
      enddo

    end subroutine xplasma_coords_expand

    !========================================
    subroutine xplasma_grids_expand(s,ineed,ier)

      ! expand list of grids

      type (xplasma), pointer :: s
      integer, intent(in) :: ineed
      integer, intent(out) :: ier

      !--------------------/local
      integer :: inold,innew,i
      type (xpgrid), dimension(:), allocatable :: xtmp
      integer, dimension(:), allocatable :: itmp
      !--------------------

      ier=0

      inold = s%ngrids
      if(inold.ge.ineed) return

      innew = inold
      do
         innew=max(20,2*innew)
         if(innew.ge.ineed) exit
      enddo

      if(inold.gt.0) then
         allocate(xtmp(inold))
         allocate(itmp(inold))
         do i=1,inold
            itmp(i)=s%grid_equiv(i)
            call xpgrid_fullcopy(s%grids(i),xtmp(i))
            call xpgrid_free(s%grids(i))
         enddo
         deallocate(s%grid_freelist)
         deallocate(s%grid_equiv)
         deallocate(s%grids)
      endif

      allocate(s%grids(innew),s%grid_freelist(innew),s%grid_equiv(innew))
      if(inold.gt.0) then
         do i=1,inold
            call xpgrid_fullcopy(xtmp(i),s%grids(i))
            s%grid_equiv(i)=itmp(i)
            call xpgrid_free(xtmp(i))
            s%grid_freelist(innew-i+1)=0
         enddo
         deallocate(xtmp,itmp)
      endif

      do i=inold+1,innew
         call xpgrid_init(s%grids(i))
         s%grid_equiv(i)=0
      enddo

      s%ngrids=innew
      s%ngrid_free=innew-inold
      do i=1,innew-inold
         s%grid_freelist(i)=innew+1-i
      enddo

    end subroutine xplasma_grids_expand

    !========================================
    subroutine xplasma_blkbxs_expand(s,ineed,ier)

      ! expand list of blkbxs

      type (xplasma), pointer :: s
      integer, intent(in) :: ineed
      integer, intent(out) :: ier

      !--------------------/local
      integer :: inold,innew,i
      type (xpblkbx), dimension(:), allocatable :: xtmp
      !--------------------

      ier=0

      inold = s%nblkbxs
      if(inold.ge.ineed) return

      innew = inold
      do
         innew=max(50,2*innew)
         if(innew.ge.ineed) exit
      enddo

      if(inold.gt.0) then
         allocate(xtmp(inold))
         do i=1,inold
            call xpblkbx_fullcopy(s%blkbxs(i),xtmp(i))
            call xpblkbx_free(s%blkbxs(i))
         enddo
         deallocate(s%blkbx_freelist)
         deallocate(s%blkbxs)
      endif

      allocate(s%blkbxs(innew),s%blkbx_freelist(innew))
      if(inold.gt.0) then
         do i=1,inold
            call xpblkbx_fullcopy(xtmp(i),s%blkbxs(i))
            call xpblkbx_free(xtmp(i))
            s%blkbx_freelist(innew-i+1)=0
         enddo
         deallocate(xtmp)
      endif

      do i=inold+1,innew
         call xpblkbx_init(s%blkbxs(i))
      enddo

      s%nblkbxs=innew
      s%nblkbx_free=innew-inold
      do i=1,innew-inold
         s%blkbx_freelist(i)=innew+1-i
      enddo

    end subroutine xplasma_blkbxs_expand

    !========================================
    subroutine xplasma_profs_expand(s,ineed,ier)

      ! expand list of profs

      type (xplasma), pointer :: s
      integer, intent(in) :: ineed
      integer, intent(out) :: ier

      !--------------------/local
      integer :: inold,innew,i
      type (xpprof), dimension(:), allocatable :: xtmp
      !--------------------

      ier=0

      inold = s%nprofs
      if(inold.ge.ineed) return

      innew = inold
      do
         innew=max(200,2*innew)
         if(innew.ge.ineed) exit
      enddo

      if(inold.gt.0) then
         allocate(xtmp(inold))
         do i=1,inold
            call xpprof_fullcopy(s%profs(i),xtmp(i))
            call xpprof_free(s%profs(i))
         enddo
         deallocate(s%prof_freelist)
         deallocate(s%profs)
      endif

      allocate(s%profs(innew),s%prof_freelist(innew))
      if(inold.gt.0) then
         do i=1,inold
            call xpprof_fullcopy(xtmp(i),s%profs(i))
            call xpprof_free(xtmp(i))
            s%prof_freelist(innew-i+1)=0
         enddo
         deallocate(xtmp)
      endif

      do i=inold+1,innew
         call xpprof_init(s%profs(i))
      enddo

      s%nprofs=innew
      s%nprof_free=innew-inold
      do i=1,innew-inold
         s%prof_freelist(i)=innew+1-i
      enddo

    end subroutine xplasma_profs_expand

    !----------------------------------
    subroutine xplasma_authorStack_expand(s,ineed)

      !  expand authorstack storage if needed

      type (xplasma), pointer :: s
      integer, intent(in) :: ineed

      !--------------
      integer :: inold,innew
      character*32, dimension(:), pointer :: ctmp
      integer, dimension(:), pointer :: itmp
      !--------------

      if(s%max_author.ge.ineed) return

      inold = s%cur_author
      innew = s%max_author

      do
         innew = max(4,2*innew)
         if(innew.ge.ineed) exit
      enddo

      if(inold.gt.0) then
         ctmp => s%authors
         itmp => s%write_enable
      endif

      allocate(s%authors(innew),s%write_enable(innew))
      s%authors(inold+1:innew)=' '
      s%write_enable(inold+1:innew)=0
      if(inold.gt.0) then
         s%authors(1:inold)=ctmp(1:inold)
         s%write_enable(1:inold)=itmp(1:inold)
         deallocate(itmp,ctmp)
      endif

      s%max_author = innew

    end subroutine xplasma_authorStack_expand

    !----------------------------------
    subroutine xplist_init(d)
      ! initialize an xplist object

      type(xplist) :: d

      d%size = -1
      nullify(d%enames)
      nullify(d%chindx)
      nullify(d%r8vals)
      nullify(d%ivals)
      nullify(d%chbuf)

    end subroutine xplist_init

    !----------------------------------

    subroutine xplist_free(d)
  !.. xplist_free                      free memory associated with xplist object
  !                                    codesys/source/xplasma2/xplasma_obj
  !
  !
  !  Author                            <author>@<place>
  !  Version                           1.00,  00mmm2009
  !  Modifications
  !
  !  1.00  08Jan2010                   Jim.Conboy@ccfe.ac.uk
  !                                    BUG !
  !                                      "A DEALLOCATE statement cannot be executed for a
  !                                       pointer (1-th argument:d%enames) whose target 
  !                                       was not created by an ALLOCATE statement."
  !                                    1) Check associated status
  !                                    2) Check deallocate return status
  !                                    3) Not clear that this will work - ref Metcalf,
  !                                       'F95/2003 explained', p105 & footnote
  !                                    
  !        15Jan2010                   ?? BUT, if the pointers are allocated in
  !                                    xplasma_define_list, then they are created via
  !                                    pointer allocation, & so this error should not occur
  !                                    
  !                                    Please report any errors msgs from this routine
  !                                    
  !----_^---------------------------------=========================================|

      type(xplist) :: d

      integer      :: istat
      character(len=*), parameter  :: cSr = 'xplasma2/xplist_free'

    ! if(associated(d%enames)) deallocate(d%enames)
    ! if(associated(d%r8vals)) deallocate(d%r8vals)
    ! if(associated(d%ivals)) deallocate(d%ivals)
    ! if(associated(d%chindx)) deallocate(d%chindx)
    ! if(associated(d%chbuf)) deallocate(d%chbuf)
    !
      if(associated(d%enames))      then
    !    if( allocated(d%enames)) &
             deallocate(d%enames,stat=istat)
         if( istat /= 0 )  &
           call errorlog( cSr//' -E- enames :',istat)
         d%enames => null()         ! leaks ?? what leaks.. 
      endif
    !
      if(associated(d%r8vals))      then
    !    if( allocated(d%r8vals)) &
             deallocate(d%r8vals,stat=istat)
         if( istat /= 0 )  &
           call errorlog( cSr//' -E- r8vals :',istat)
         d%r8vals => null()         ! leaks ?? what leaks.. 
      endif
    !
      if(associated(d%ivals))      then
    !    if( allocated(d%ivals)) &
             deallocate(d%ivals,stat=istat)
         if( istat /= 0 )  &
           call errorlog( cSr//' -E- ivals :',istat)
         d%ivals => null()         ! leaks ?? what leaks.. 
      endif
    !
      if(associated(d%chindx))      then
    !    if( allocated(d%chindx)) &
             deallocate(d%chindx,stat=istat)
         if( istat /= 0 )  &
           call errorlog( cSr//' -E- chindx :',istat)
         d%chindx => null()         ! leaks ?? what leaks.. 
      endif
    !
      if(associated(d%chbuf))      then
    !    if( allocated(d%chbuf)) &
             deallocate(d%chbuf,stat=istat)
         if( istat /= 0 )  &
           call errorlog( cSr//' -E- chbuf :',istat)
         d%chbuf => null()         ! leaks ?? what leaks.. 
      endif
    !

      d%size = 0

    end subroutine xplist_free

    !----------------------------------

    subroutine xplist_fullcopy(d1,d2)
      ! copy d1->d2 -- new memory, not a pointer copy
      ! this is private method; d2 object initialization handled elsewhere.

      type(xplist) :: d1,d2
      
      integer :: i,isize,jsize,ilen

      isize = d1%size
      d2%size = isize

      if(associated(d1%enames)) then
         allocate(d2%enames(isize))
         d2%enames = d1%enames
      else
         nullify(d2%enames)
      endif

      if(associated(d1%r8vals)) then
         allocate(d2%r8vals(isize))
         d2%r8vals = d1%r8vals
      else
         nullify(d2%r8vals)
      endif

      if(associated(d1%ivals)) then
         allocate(d2%ivals(isize))
         d2%ivals  = d1%ivals
      else
         nullify(d2%ivals)
      endif

      if(associated(d1%chindx)) then
         allocate(d2%chindx(2,isize))
         d2%chindx = d1%chindx
         ilen = d2%chindx(2,isize)  ! last char address
         if(ilen.gt.0) then
            allocate(d2%chbuf(ilen))
            d2%chbuf = d1%chbuf
         else
            nullify(d2%chbuf)
         endif
      else
         nullify(d2%chindx)
         nullify(d2%chbuf)
      endif

    end subroutine xplist_fullcopy

    !-------------------------------------------------
    subroutine xplist_ncdefine(d,prefix,icdf,ier)
      ! (NetCDF) define xplist object

      type(xplist) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      ! mod dmc -- deal correctly with length<2 lists

      !-------------------------------------
      integer :: ielems,chtot
      character*32 :: enam2(2)
      character :: chdum(2)
      integer :: chindum(2,2)
      integer :: idum(2)
      real*8 :: r8dum(2)
      !-------------------------------------

      ier=0

      ielems = d%size

      chtot=d%chindx(2,ielems)
      if(chtot.gt.1) then
         call cdf_define(icdf,trim(prefix)//'%chbuf',d%chbuf)
      else
         chdum=" "
         if(chtot.eq.1) chdum(1)=d%chbuf(1)
         call cdf_define(icdf,trim(prefix)//'%chbuf',chdum)
      endif

      if(ielems.gt.1) then
         call cdf_define(icdf,trim(prefix)//'%enames',d%enames)
         call cdf_define(icdf,trim(prefix)//'%chindx',d%chindx)
         call cdf_define(icdf,trim(prefix)//'%ivals',d%ivals)
         call cdf_define(icdf,trim(prefix)//'%r8vals',d%r8vals)

      else
         enam2=' '
         chindum=0
         idum=0
         r8dum=0.0d0
         enam2(1:ielems)=d%enames(1:ielems)
         chindum(1:2,1:ielems)=d%chindx(1:2,1:ielems)
         idum(1:ielems)=d%ivals(1:ielems)
         r8dum(1:ielems)=d%r8vals(1:ielems)
         
         call cdf_define(icdf,trim(prefix)//'%enames',enam2)
         call cdf_define(icdf,trim(prefix)//'%chindx',chindum)
         call cdf_define(icdf,trim(prefix)//'%ivals',idum)
         call cdf_define(icdf,trim(prefix)//'%r8vals',r8dum)
      endif

    end subroutine xplist_ncdefine

    !-------------------------------------------------
    subroutine xplist_ncwrite(d,prefix,icdf,ier)
      ! (NetCDF) define xplist object

      type(xplist) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      ! mod dmc -- deal correctly with length<2 lists

      !-------------------------------------
      integer :: ielems,chtot
      character*32 :: enam2(2)
      character :: chdum(2)
      integer :: chindum(2,2)
      integer :: idum(2)
      real*8 :: r8dum(2)
      !-------------------------------------

      ier=0

      ielems = d%size

      chtot=d%chindx(2,ielems)
      if(chtot.gt.1) then
         call cdf_write(icdf,trim(prefix)//'%chbuf',d%chbuf)
      else
         chdum=" "
         if(chtot.eq.1) chdum(1)=d%chbuf(1)
         call cdf_write(icdf,trim(prefix)//'%chbuf',chdum)
      endif

      if(ielems.gt.1) then
         call cdf_write(icdf,trim(prefix)//'%enames',d%enames)
         call cdf_write(icdf,trim(prefix)//'%chindx',d%chindx)
         call cdf_write(icdf,trim(prefix)//'%ivals',d%ivals)
         call cdf_write(icdf,trim(prefix)//'%r8vals',d%r8vals)

      else
         enam2=' '
         chindum=0
         idum=0
         r8dum=0.0d0
         enam2(1:ielems)=d%enames(1:ielems)
         chindum(1:2,1:ielems)=d%chindx(1:2,1:ielems)
         idum(1:ielems)=d%ivals(1:ielems)
         r8dum(1:ielems)=d%r8vals(1:ielems)
         
         call cdf_write(icdf,trim(prefix)//'%enames',enam2)
         call cdf_write(icdf,trim(prefix)//'%chindx',chindum)
         call cdf_write(icdf,trim(prefix)//'%ivals',idum)
         call cdf_write(icdf,trim(prefix)//'%r8vals',r8dum)
      endif

    end subroutine xplist_ncwrite

    subroutine xplist_ncread(d,prefix,icdf,ier)
      ! (NetCDF) read xplist object

      type(xplist) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      !-------------------------------------
      integer :: chtot,ielems,idims(3),istat
      character*8 ztype
      character*32 :: enam2(2)
      character :: chdum(2)
      integer :: chindum(2,2)
      integer :: idum(2)
      real*8 :: r8dum(2)
      !-------------------------------------

      ier=0

      call cdf_inquire(icdf,trim(prefix)//'%ivals',idims,ztype,istat)
      if(istat.ne.0) ier=13
      if(ier.ne.0) return

      ielems=idims(1)
      if(ielems.eq.2) then
         ! list may actually be shorter...
         call cdf_read(icdf,trim(prefix)//'%enames',enam2,ier)
         if(enam2(2).eq.' ') ielems=1
         if(enam2(1).eq.' ') ielems=0
      endif

      d%size = ielems

      allocate(d%enames(ielems))
      allocate(d%chindx(2,ielems))
      allocate(d%ivals(ielems))
      allocate(d%r8vals(ielems))

      if(ielems.ge.2) then
         call cdf_read(icdf,trim(prefix)//'%enames',d%enames,ier)
         if(ier.ne.0) ier=13
         if(ier.eq.13)return

         call cdf_read(icdf,trim(prefix)//'%chindx',d%chindx,ier)
         if(ier.ne.0) ier=13
         if(ier.eq.13)return

         call cdf_read(icdf,trim(prefix)//'%r8vals',d%r8vals,ier)
         if(ier.ne.0) ier=13
         if(ier.eq.13)return

         call cdf_read(icdf,trim(prefix)//'%ivals',d%ivals,ier)
         if(ier.ne.0) ier=13
         if(ier.eq.13)return

      else
         call cdf_read(icdf,trim(prefix)//'%enames',enam2,ier)
         if(ier.ne.0) ier=13
         if(ier.eq.13)return

         call cdf_read(icdf,trim(prefix)//'%chindx',chindum,ier)
         if(ier.ne.0) ier=13
         if(ier.eq.13)return

         call cdf_read(icdf,trim(prefix)//'%r8vals',r8dum,ier)
         if(ier.ne.0) ier=13
         if(ier.eq.13)return

         call cdf_read(icdf,trim(prefix)//'%ivals',idum,ier)
         if(ier.ne.0) ier=13
         if(ier.eq.13)return

         d%enames(1:ielems)=enam2(1:ielems)
         d%chindx(1:2,1:ielems)=chindum(1:2,1:ielems)
         d%ivals(1:ielems)=idum(1:ielems)
         d%r8vals(1:ielems)=r8dum(1:ielems)
      endif

      chtot = d%chindx(2,ielems)

      if(chtot.gt.0) then
         allocate(d%chbuf(chtot))
         if(chtot.gt.1) then
            call cdf_read(icdf,trim(prefix)//'%chbuf',d%chbuf,ier)
            if(ier.ne.0) ier=13
            if(ier.eq.13)return
         else
            call cdf_read(icdf,trim(prefix)//'%chbuf',chdum,ier)
            if(ier.ne.0) ier=13
            if(ier.eq.13)return
            d%chbuf(1)=chdum(1)
         endif
      endif

    end subroutine xplist_ncread

    !----------------------------------

    subroutine xpcoord_init(d)
      ! initialize an xpcoord object

      type(xpcoord) :: d

      d%periodic = .FALSE.
      d%ngridc = 0
      d%ngrid_max = 0
      nullify(d%grid_ids)

    end subroutine xpcoord_init

    !----------------------------------
    
    subroutine xpcoord_fullcopy(d1,d2)
      !  copy xpcoord object d1->d2; copy memory not just pointer
      ! this is private method; d2 object initialization handled elsewhere.

      type(xpcoord) :: d1,d2

      d2%periodic = d1%periodic
      d2%ngridc = d1%ngridc
      d2%ngrid_max = d1%ngrid_max

      if(associated(d1%grid_ids)) then
         allocate(d2%grid_ids(d2%ngrid_max))
         d2%grid_ids = d1%grid_ids
      else
         nullify(d2%grid_ids)
      endif

    end subroutine xpcoord_fullcopy

    !----------------------------------

    subroutine xpcoord_free(d)
      ! free memory associated with xpcoord object

      type(xpcoord) :: d

      d%periodic = .FALSE.
      d%ngridc = 0
      d%ngrid_max = 0
      if(associated(d%grid_ids)) deallocate(d%grid_ids)

    end subroutine xpcoord_free

    !-------------------------------------------------
    subroutine xpcoord_ncdefine(d,prefix,icdf,ier)
      ! (NetCDF) define xpcoord object

      type(xpcoord) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      integer :: iperiodic = 0
      integer :: inum

      ier=0

      call cdf_define(icdf,trim(prefix)//'%periodic',iperiodic)

      inum = d%ngridc
      if(inum.gt.0) then
         call cdf_define(icdf,trim(prefix)//'%grid_ids',d%grid_ids(1:inum))
      endif

    end subroutine xpcoord_ncdefine

    subroutine xpcoord_ncwrite(d,prefix,icdf,ier)
      ! (NetCDF) write xpcoord object

      type(xpcoord) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      integer :: iperiodic,inum

      ier=0

      if(d%periodic) then
         iperiodic=1
      else
         iperiodic=0
      endif

      call cdf_write(icdf,trim(prefix)//'%periodic',iperiodic)

      inum = d%ngridc
      if(inum.gt.0) then
         call cdf_write(icdf,trim(prefix)//'%grid_ids',d%grid_ids(1:inum))
      endif

    end subroutine xpcoord_ncwrite

    subroutine xpcoord_ncread(d,prefix,icdf,ier)
      ! (NetCDF) read xpcoord object

      type(xpcoord) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      character*8 ztype
      integer :: iperiodic,idims(3),istat,inum

      ier=0

      call cdf_read(icdf,trim(prefix)//'%periodic',iperiodic,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return
      
      d%periodic = (iperiodic.eq.1)

      call cdf_inquire(icdf,trim(prefix)//'%grid_ids',idims,ztype,istat)
      if(istat.eq.0) then
         inum=idims(1)
         call xplasma_coord_expand(d,inum)
         call cdf_read(icdf,trim(prefix)//'%grid_ids',d%grid_ids(1:inum))
         d%ngridc=inum
      else
         d%ngridc=0
      endif

    end subroutine xpcoord_ncread

    !----------------------------------

    subroutine xpgrid_init(d)
      ! initialize an xpgrid object

      type(xpgrid) :: d

      d%size = 0
      if(allocated(d%xpkg)) deallocate(d%xpkg)
      d%coord = 0
      d%nrefs = 0

    end subroutine xpgrid_init

    !----------------------------------

    subroutine xpgrid_fullcopy(d1,d2)
      ! copy grid object d1->d2; full memory copy not just ptrs
      ! this is private method; d2 object initialization handled elsewhere.

      type(xpgrid) :: d1,d2

      d2%size = d1%size
      d2%coord = d1%coord
      d2%nrefs = d1%nrefs

      if(allocated(d1%xpkg)) then
         allocate(d2%xpkg(d2%size,4))
         d2%xpkg = d1%xpkg
      else
         if(allocated(d2%xpkg)) deallocate(d2%xpkg)
      endif

    end subroutine xpgrid_fullcopy

    !----------------------------------

    subroutine xpgrid_free(d)
      ! free memory associated with xpgrid object

      type(xpgrid) :: d

      if(allocated(d%xpkg)) deallocate(d%xpkg)
      d%size = 0
      d%coord = 0
      d%nrefs = 0

    end subroutine xpgrid_free

    !-------------------------------------------------
    subroutine xpgrid_ncdefine(d,prefix,icdf,ier)
      ! (NetCDF) define xpgrid object

      type(xpgrid) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      ier = 0

      call cdf_define(icdf,trim(prefix)//'%xpkg',d%xpkg)
      call cdf_define(icdf,trim(prefix)//'%coord',d%coord)
      call cdf_define(icdf,trim(prefix)//'%nrefs',d%nrefs)

    end subroutine xpgrid_ncdefine

    subroutine xpgrid_ncwrite(d,prefix,icdf,ier)
      ! (NetCDF) write xpgrid object

      type(xpgrid) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      ier=0

      call cdf_write(icdf,trim(prefix)//'%xpkg',d%xpkg)
      call cdf_write(icdf,trim(prefix)//'%coord',d%coord)
      call cdf_write(icdf,trim(prefix)//'%nrefs',d%nrefs)

    end subroutine xpgrid_ncwrite

    subroutine xpgrid_ncread(d,prefix,icdf,ier)
      ! (NetCDF) read xpgrid object

      type(xpgrid) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      integer :: idims(3),istat
      character*8 :: ztype

      ier=0

      call cdf_inquire(icdf,trim(prefix)//'%xpkg',idims,ztype,istat)
      d%size = idims(1)

      allocate(d%xpkg(d%size,4))
      d%xpkg=czero !MG, I added this line to please Intel compiler


      call cdf_read(icdf,trim(prefix)//'%xpkg',d%xpkg,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,trim(prefix)//'%coord',d%coord,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,trim(prefix)//'%nrefs',d%nrefs,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

    end subroutine xpgrid_ncread
    !----------------------------------

    subroutine xpprof_init(d)
      ! initialize an xpprof object

      type(xpprof) :: d

      d%rank = 0
      d%gridIds = 0
      d%profIds = 0
      d%kspline = 0
      d%prof_counter = 0
      d%size = 0
      if(allocated(d%buf)) deallocate(d%buf)

    end subroutine xpprof_init

    !----------------------------------

    subroutine xpprof_fullcopy(d1,d2)
      ! copy prof object d1->d2; full memory copy not just ptrs
      ! this is private method; d2 object initialization handled elsewhere.

      type(xpprof) :: d1,d2

      d2%rank = d1%rank
      d2%gridIds = d1%gridIds
      d2%profIds = d1%profIds
      d2%kspline = d1%kspline
      d2%prof_counter = d1%prof_counter

      d2%size = d1%size

      if(allocated(d1%buf)) then
         allocate(d2%buf(d2%size))
         d2%buf = d1%buf
      else
         if(allocated(d2%buf)) deallocate(d2%buf)
      endif

    end subroutine xpprof_fullcopy

    !----------------------------------

    subroutine xpprof_free(d)
      ! free memory associated with xpprof object

      type(xpprof) :: d

      d%rank = 0
      d%gridIds = 0
      d%profIds = 0
      d%kspline = 0
      d%prof_counter = 0
      d%size = 0
      if(allocated(d%buf)) deallocate(d%buf)

    end subroutine xpprof_free

    !----------------------------------
    subroutine xpprof_ncdefine(d,prefix,icdf,ier)
      !  Define xpprof object output records for NetCDF file

      type(xpprof) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      !----------------------------
      integer, dimension(:), pointer :: ibuf
      logical :: izero
      !----------------------------

      ier = 0

      allocate(ibuf(d%rank+6))  ! for cdf_define type def only...
      call cdf_define(icdf,trim(prefix)//'%ibuf',ibuf)
      deallocate(ibuf)

      if(d%size.gt.0) then
         izero = (maxval(abs(d%buf)) .eq. czero)
         if(izero) then
            d%size = -d%size  ! signals data of all ZEROs; do not write
            ! this sign flip MUST be undone in the subsequent xpprof_write call

         else
            call cdf_define(icdf,trim(prefix)//'%buf',d%buf)
         endif
      end if

    end subroutine xpprof_ncdefine

    !----------------------------------
    subroutine xpprof_ncwrite(d,prefix,icdf,ier)
      !  Write xpprof object output records for NetCDF file

      type(xpprof) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      !----------------------------
      integer, dimension(:), pointer :: ibuf
      !----------------------------

      ier = 0

      allocate(ibuf(d%rank+6))
      ibuf(1)=d%rank
      ibuf(2:d%rank+1)=d%gridIds(1:d%rank)
      ibuf(d%rank+2)=d%kspline
      ibuf(d%rank+3)=d%prof_counter
      ibuf(d%rank+4)=d%size
      ibuf(d%rank+5)=0  ! spare
      ibuf(d%rank+6)=d%profIds(1)
      call cdf_write(icdf,trim(prefix)//'%ibuf',ibuf)
      deallocate(ibuf)

      if(d%size.gt.0) then
         call cdf_write(icdf,trim(prefix)//'%buf',d%buf)
      else
         d%size = abs(d%size)  ! possibly undo sign flip from xpprof_ncdefine
      endif

    end subroutine xpprof_ncwrite

    !----------------------------------
    subroutine xpprof_ncread(d,prefix,icdf,ier)
      !  Read xpprof object records from NetCDF file

      type(xpprof) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      !----------------------------
      integer :: idims(3),istat
      integer, dimension(:), pointer :: ibuf
      character*8 :: ztype
      !----------------------------

      ier = 0

      call cdf_inquire(icdf,trim(prefix)//'%ibuf',idims,ztype,istat)
      if(istat.ne.0) then
         if(ier.ne.0) ier=13
         if(ier.ne.0) return
      endif

      allocate(ibuf(idims(1)))
      call cdf_read(icdf,trim(prefix)//'%ibuf',ibuf,ier)
      if(ier.ne.0) then
         ier=13
         deallocate(ibuf)
         return
      endif

      d%rank = ibuf(1)
      d%gridIds = 0
      d%gridIds(1:d%rank) = ibuf(2:d%rank+1)
      d%kspline = ibuf(d%rank+2)
      d%prof_counter = ibuf(d%rank+3)
      d%size = ibuf(d%rank+4)
      d%profIds(1)=ibuf(d%rank+6)
      deallocate(ibuf)

      if(d%size.gt.0) then
         allocate(d%buf(d%size))
         call cdf_read(icdf,trim(prefix)//'%buf',d%buf,ier)
         if(ier.ne.0) then
            ier=13
         endif

      else if(d%size.lt.0) then
         !  restore zeroes (without reading from file)

         d%size = -d%size
         allocate(d%buf(d%size))
         d%buf = czero

      endif

    end subroutine xpprof_ncread

    !----------------------------------

    subroutine xpblkbx_init(d)
      ! initialize an xpblkbx object

      type(xpblkbx) :: d

      d%type = 0
      if(allocated(d%iwork)) deallocate(d%iwork)
      if(allocated(d%r8work)) deallocate(d%r8work)

    end subroutine xpblkbx_init

    !----------------------------------
    ! private pack/unpack methods...

    subroutine xplasma_packmaster(s,iptr,r8ptr)
      !  pack all xplasma data into I and R8 streams...
      !  structure of this routine parallels xplasma_ncwrite(...)

      type (xplasma), pointer :: s
      integer, dimension(:), pointer :: iptr
      real*8, dimension(:), pointer :: r8ptr

      !--------------------
      integer :: iadr,r8adr,ilogical,ii,j,ktype,itchk
      integer :: luntty = 6
      character*32 :: tname
      !--------------------

      iadr = iptr(5)
      r8adr = iptr(6)

      if(xplasma_trace_name.ne.' ') then
         itchk=1
         tname=xplasma_trace_name
         call uupper(tname)
      else
         itchk=0
         tname=' '
      endif

      !  header...

      call pack_c1a(s%label,iadr,iptr) ! there must be a label...

      if(s%axisymm) then
         ilogical=1
      else
         ilogical=0
      endif
      call pack_iscalar(ilogical,iadr,iptr)

      if(s%scrapeoff) then
         ilogical=1
      else
         ilogical=0
      endif
      call pack_iscalar(ilogical,iadr,iptr)

      call pack_iscalar(s%bphi_ccw,iadr,iptr)
      call pack_iscalar(s%jphi_ccw,iadr,iptr)

      call pack_r8scalar(s%time,r8adr,r8ptr)

      !  author stack...

      call pack_iscalar(s%cur_author,iadr,iptr)
      do ii=1,s%cur_author
         call pack_str(s%authors(ii),iadr,iptr)
      enddo
      if(s%cur_author.gt.0) then
         call pack_ivec(s%write_enable(1:s%cur_author),iadr,iptr)
      endif

      !  book-keeping data for constituent objects...

      call pack_iscalar(s%nlists,iadr,iptr)
      call pack_iscalar(s%nlist_free,iadr,iptr)
      if(s%nlists.gt.0) then
         call pack_ivec(s%list_freelist,iadr,iptr)
      endif
    
      call pack_iscalar(s%ncoords,iadr,iptr)
      call pack_iscalar(s%ncoord_free,iadr,iptr)
      if(s%ncoords.gt.0) then
         call pack_ivec(s%coord_freelist,iadr,iptr)
      endif
   
      call pack_iscalar(s%ngrids,iadr,iptr)
      call pack_iscalar(s%ngrid_free,iadr,iptr)
      if(s%ngrids.gt.0) then
         call pack_ivec(s%grid_freelist,iadr,iptr)
         call pack_ivec(s%grid_equiv,iadr,iptr)
      endif

      call pack_iscalar(s%nprofs,iadr,iptr)
      call pack_iscalar(s%nprof_free,iadr,iptr)
      if(s%nprofs.gt.0) then
         call pack_ivec(s%prof_freelist,iadr,iptr)
      endif
  
      call pack_iscalar(s%nblkbxs,iadr,iptr)
      call pack_iscalar(s%nblkbx_free,iadr,iptr)
      if(s%nblkbxs.gt.0) then
         call pack_ivec(s%blkbx_freelist,iadr,iptr)
      endif

      !  mscellaneous scalars...

      call pack_iscalar(s%prof_counter,iadr,iptr)
      call pack_iscalar(s%eq_counter,iadr,iptr)

      call pack_r8scalar(s%bdytol,r8adr,r8ptr)
      call pack_r8scalar(s%ajac_maxVar,r8adr,r8ptr)

      call pack_iscalar(s%rzOrder,iadr,iptr)
      call pack_iscalar(s%kmom,iadr,iptr)

      !  the dictionary...

      call pack_iscalar(s%nitems,iadr,iptr)

      if(s%nitems.gt.0) then
         call pack_ivec(s%order(1:s%nitems),iadr,iptr)
         do ii=1,s%nitems
            call pack_str(s%dict(ii)%name,iadr,iptr)
            call pack_str(s%dict(ii)%author,iadr,iptr)
            if(allocated(s%dict(ii)%label)) then
               call pack_c1a(s%dict(ii)%label,iadr,iptr)
            else
               call pack_iscalar(0,iadr,iptr)
            endif
            if(allocated(s%dict(ii)%units)) then
               call pack_c1a(s%dict(ii)%units,iadr,iptr)
            else
               call pack_iscalar(0,iadr,iptr)
            endif
            call pack_iscalar(s%dict(ii)%dtype,iadr,iptr)
            call pack_iscalar(s%dict(ii)%dindex,iadr,iptr)
         enddo
      endif

      !  the items themselves...

      do ii=1,s%nitems
         j = s%dict(ii)%dindex
         ktype = s%dict(ii)%dtype

         if(ktype.eq.xplasma_listType) then
            call xplist_pack(s%lists(j),iadr,iptr,r8adr,r8ptr)

         else if(ktype.eq.xplasma_coordType) then
            call xpcoord_pack(s%coords(j),iadr,iptr,r8adr,r8ptr)

         else if(ktype.eq.xplasma_gridType) then
            call xpgrid_pack(s%grids(j),iadr,iptr,r8adr,r8ptr)

         else if(ktype.eq.xplasma_profType) then
            call xpprof_pack(s%profs(j),iadr,iptr,r8adr,r8ptr)

         else if(ktype.eq.xplasma_blkbxType) then
            call xpblkbx_pack(s%blkbxs(j),iadr,iptr,r8adr,r8ptr)

         endif

         if(itchk.gt.0) then
            if(s%dict(ii)%dindex.gt.0) then
               if(tname.eq.s%dict(ii)%name) then
                  call tchk(s,'xplasma_packmaster packed',tname)
               endif
            endif
         endif

      enddo

      !  store the final sizes...

      iptr(5)=iadr
      iptr(6)=r8adr

    end subroutine xplasma_packmaster

    subroutine xplist_pack(xx, iadr, iptr, r8adr, r8ptr)

      !  pack a single list

      type (xplist) :: xx
      integer :: iadr
      integer, dimension(:), pointer :: iptr
      integer :: r8adr
      real*8, dimension(:), pointer :: r8ptr

      !------------
      integer :: ielems,nelems,ichsize
      !------------

      nelems = xx%size
      call pack_iscalar(nelems, iadr, iptr)

      do ielems = 1,nelems
         call pack_str(xx%enames(ielems), iadr, iptr)
      enddo

      do ielems = 1,nelems
         call pack_iscalar(xx%chindx(1,ielems), iadr, iptr)
         call pack_iscalar(xx%chindx(2,ielems), iadr, iptr)
      enddo

      if(nelems.gt.0) then
         if(associated(xx%chbuf)) then
            ichsize = size(xx%chbuf)
         else
            ichsize = 0
         endif
         call pack_iscalar(ichsize, iadr, iptr)

         if(ichsize.gt.0) then
            call pack_c1a(xx%chbuf, iadr, iptr)
         endif

         call pack_ivec(xx%ivals(1:nelems), iadr, iptr)
         call pack_r8vec(xx%r8vals(1:nelems), r8adr, r8ptr)

      endif

    end subroutine xplist_pack

    subroutine xpcoord_pack(xx, iadr, iptr, r8adr, r8ptr)

      !  pack a single coord

      type (xpcoord) :: xx
      integer :: iadr
      integer, dimension(:), pointer :: iptr
      integer :: r8adr
      real*8, dimension(:), pointer :: r8ptr

      !------------------
      integer :: ilogical, iisize
      !------------------

      if(xx%periodic) then
         ilogical = 1
      else
         ilogical = 0
      endif

      call pack_iscalar(ilogical, iadr, iptr)

      call pack_iscalar(xx%ngridc, iadr, iptr)
      call pack_iscalar(xx%ngrid_max, iadr, iptr)

      if(xx%ngrid_max.gt.0) then
         iisize = size(xx%grid_ids)
         call pack_iscalar(iisize, iadr, iptr)
         call pack_ivec(xx%grid_ids, iadr, iptr)
      endif

    end subroutine xpcoord_pack

    subroutine xpgrid_pack(xx, iadr, iptr, r8adr, r8ptr)

      !  pack a single grid

      type (xpgrid) :: xx
      integer :: iadr
      integer, dimension(:), pointer :: iptr
      integer :: r8adr
      real*8, dimension(:), pointer :: r8ptr

      !---------------
      integer :: ii,isize
      !---------------

      isize = xx%size
      call pack_iscalar(xx%size, iadr, iptr)
      call pack_iscalar(xx%coord, iadr, iptr)
      call pack_iscalar(xx%nrefs, iadr, iptr)

      if(isize.gt.0) then
         do ii=1,4
            call pack_r8vec(xx%xpkg(1:isize,ii), r8adr, r8ptr)
         enddo
      endif
      
    end subroutine xpgrid_pack

    subroutine xpprof_pack(xx, iadr, iptr, r8adr, r8ptr)

      !  pack a single prof

      type (xpprof) :: xx
      integer :: iadr
      integer, dimension(:), pointer :: iptr
      integer :: r8adr
      real*8, dimension(:), pointer :: r8ptr

      call pack_iscalar(xx%rank, iadr, iptr)
      call pack_ivec(xx%gridIds, iadr, iptr)
      call pack_ivec(xx%profIds, iadr, iptr)
      call pack_iscalar(xx%kspline, iadr, iptr)
      call pack_iscalar(xx%prof_counter, iadr, iptr)
      call pack_iscalar(xx%size, iadr, iptr)

      if(xx%size.gt.0) then
         call pack_r8vec(xx%buf, r8adr, r8ptr)
      endif

    end subroutine xpprof_pack

    subroutine xpblkbx_pack(xx, iadr, iptr, r8adr, r8ptr)

      !  pack a single blkbx

      type (xpblkbx) :: xx
      integer :: iadr
      integer, dimension(:), pointer :: iptr
      integer :: r8adr
      real*8, dimension(:), pointer :: r8ptr

      integer :: izero = 0
      real*8 :: r8zero = 0.0d0

      call pack_iscalar(xx%type, iadr, iptr)

      if(allocated(xx%iwork)) then
         call pack_ivec(xx%iwork, iadr, iptr)
      else
         call pack_iscalar(izero, iadr, iptr)
      endif

      if(allocated(xx%r8work)) then
         call pack_r8vec(xx%r8work, r8adr, r8ptr)
      else
         call pack_r8scalar(r8zero, r8adr, r8ptr)
      endif

    end subroutine xpblkbx_pack

    subroutine xplasma_unpackmaster(s,iptr,r8ptr,ier)
      !  unpack all xplasma data from I and R8 streams into s elements
      !  s starts empty
      !  structure of this routine parallels xplasma_ncread(...)

      type (xplasma), pointer :: s
      integer, dimension(:) :: iptr
      real*8, dimension(:) :: r8ptr
      integer, intent(out) :: ier

      !--------------------
      integer :: iadr,r8adr,imaxa,r8maxa,ilogical,ii,j,ktype
      integer :: isize,r8size,ineed,itchk
      character*32 :: tname
      integer :: luntty = 6
      !--------------------

      imaxa = iptr(5)
      r8maxa = iptr(6)
      iadr = 6
      r8adr = 0

      if(xplasma_trace_name.ne.' ') then
         itchk=1
         tname=xplasma_trace_name
         call uupper(tname)
      else
         itchk=0
         tname=' '
      endif

      !  header...

      isize=iptr(iadr+1)
      allocate(s%label(isize))
      call unpack_c1a(s%label,iadr,iptr,ier)
      if(ier.ne.0) return

      call unpack_iscalar(ilogical,iadr,iptr,ier)
      s%axisymm = (ilogical.eq.1)
      if(ier.ne.0) return

      call unpack_iscalar(ilogical,iadr,iptr,ier)
      s%scrapeoff = (ilogical.eq.1)
      if(ier.ne.0) return

      call unpack_iscalar(s%bphi_ccw,iadr,iptr,ier)
      if(ier.ne.0) return
      call unpack_iscalar(s%jphi_ccw,iadr,iptr,ier)
      if(ier.ne.0) return

      call unpack_r8scalar(s%time,r8adr,r8ptr,ier)
      if(ier.ne.0) return

      !  author stack...

      call unpack_iscalar(ineed,iadr,iptr,ier)
      if(ier.ne.0) return

      call xplasma_authorStack_expand(s,ineed)
      s%cur_author=ineed

      do ii=1,s%cur_author
         call unpack_str(s%authors(ii),iadr,iptr,ier)
         if(ier.ne.0) return
      enddo
      if(s%cur_author.gt.0) then
         call unpack_ivec(s%write_enable(1:s%cur_author),iadr,iptr,ier)
      endif
      if(ier.ne.0) return

      !  book-keeping data for constituent objects...

      call unpack_iscalar(s%nlists,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%nlists.gt.0) then
         allocate(s%lists(s%nlists))
         do ii=1,s%nlists
            call xplist_init(s%lists(ii))
         enddo
         allocate(s%list_freelist(s%nlists))
      endif
      call unpack_iscalar(s%nlist_free,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%nlists.gt.0) then
         call unpack_ivec(s%list_freelist,iadr,iptr,ier)
      endif
      if(ier.ne.0) return

      call unpack_iscalar(s%ncoords,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%ncoords.gt.0) then
         allocate(s%coords(s%ncoords))
         do ii=1,s%ncoords
            call xpcoord_init(s%coords(ii))
         enddo
         allocate(s%coord_freelist(s%ncoords))
      endif
      call unpack_iscalar(s%ncoord_free,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%ncoords.gt.0) then
         call unpack_ivec(s%coord_freelist,iadr,iptr,ier)
      endif
      if(ier.ne.0) return

      call unpack_iscalar(s%ngrids,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%ngrids.gt.0) then
         allocate(s%grids(s%ngrids))
         do ii=1,s%ngrids
            call xpgrid_init(s%grids(ii))
         enddo
         allocate(s%grid_freelist(s%ngrids))
         allocate(s%grid_equiv(s%ngrids))
      endif
      call unpack_iscalar(s%ngrid_free,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%ngrids.gt.0) then
         call unpack_ivec(s%grid_freelist,iadr,iptr,ier)
         call unpack_ivec(s%grid_equiv,iadr,iptr,ier)
      endif
      if(ier.ne.0) return

      call unpack_iscalar(s%nprofs,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%nprofs.gt.0) then
         allocate(s%profs(s%nprofs))
         do ii=1,s%nprofs
            call xpprof_init(s%profs(ii))
         enddo
         allocate(s%prof_freelist(s%nprofs))
      endif
      call unpack_iscalar(s%nprof_free,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%nprofs.gt.0) then
         call unpack_ivec(s%prof_freelist,iadr,iptr,ier)
      endif
      if(ier.ne.0) return

      call unpack_iscalar(s%nblkbxs,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%nblkbxs.gt.0) then
         allocate(s%blkbxs(s%nblkbxs))
         do ii=1,s%nblkbxs
            call xpblkbx_init(s%blkbxs(ii))
         enddo
         allocate(s%blkbx_freelist(s%nblkbxs))
      endif
      call unpack_iscalar(s%nblkbx_free,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%nblkbxs.gt.0) then
         call unpack_ivec(s%blkbx_freelist,iadr,iptr,ier)
      endif
      if(ier.ne.0) return

      !  mscellaneous scalars...

      call unpack_iscalar(s%prof_counter,iadr,iptr,ier)
      if(ier.ne.0) return
      call unpack_iscalar(s%eq_counter,iadr,iptr,ier)
      if(ier.ne.0) return

      call unpack_r8scalar(s%bdytol,r8adr,r8ptr,ier)
      if(ier.ne.0) return
      call unpack_r8scalar(s%ajac_maxVar,r8adr,r8ptr,ier)
      if(ier.ne.0) return

      call unpack_iscalar(s%rzOrder,iadr,iptr,ier)
      if(ier.ne.0) return
      call unpack_iscalar(s%kmom,iadr,iptr,ier)
      if(ier.ne.0) return

      !  the dictionary...

      call unpack_iscalar(s%nitems,iadr,iptr,ier)
      if(ier.ne.0) return
      if(s%nitems.gt.0) then
         call xplasma_expand(s,s%nitems,ier)
         if(ier.ne.0) return
         call unpack_ivec(s%order(1:s%nitems),iadr,iptr,ier)
         if(ier.ne.0) return
         do ii=1,s%nitems
            call unpack_str(s%dict(ii)%name,iadr,iptr,ier)
            if(ier.ne.0) return
            call unpack_str(s%dict(ii)%author,iadr,iptr,ier)
            if(ier.ne.0) return

            isize=iptr(iadr+1)
            if(isize.gt.0) then
               allocate(s%dict(ii)%label(isize))
               call unpack_c1a(s%dict(ii)%label,iadr,iptr,ier)
            else
               iadr=iadr+1
            endif
            if(ier.ne.0) return

            isize=iptr(iadr+1)
            if(isize.gt.0) then
               allocate(s%dict(ii)%units(isize))
               call unpack_c1a(s%dict(ii)%units,iadr,iptr,ier)
            else
               iadr=iadr+1
            endif
            if(ier.ne.0) return

            call unpack_iscalar(s%dict(ii)%dtype,iadr,iptr,ier)
            if(ier.ne.0) return
            call unpack_iscalar(s%dict(ii)%dindex,iadr,iptr,ier)
            if(ier.ne.0) return

            if(itchk.gt.0) then
               if(s%dict(ii)%dindex.gt.0) then
                  if(tname.eq.s%dict(ii)%name) then
                     call tchk(s,'xplasma_unpackmaster unpacked',tname)
                  endif
               endif
            endif

         enddo
      endif

      !  the items themselves...

      do ii=1,s%nitems
         j = s%dict(ii)%dindex
         ktype = s%dict(ii)%dtype

         if(ktype.eq.xplasma_listType) then
            call xplist_unpack(s%lists(j),iadr,iptr,r8adr,r8ptr,ier)
            if(ier.ne.0) return

         else if(ktype.eq.xplasma_coordType) then
            call xpcoord_unpack(s%coords(j),iadr,iptr,r8adr,r8ptr,ier)
            if(ier.ne.0) return

         else if(ktype.eq.xplasma_gridType) then
            call xpgrid_unpack(s%grids(j),iadr,iptr,r8adr,r8ptr,ier)
            if(ier.ne.0) return

         else if(ktype.eq.xplasma_profType) then
            call xpprof_unpack(s%profs(j),iadr,iptr,r8adr,r8ptr,ier)
            if(ier.ne.0) return

         else if(ktype.eq.xplasma_blkbxType) then
            call xpblkbx_unpack(s%blkbxs(j),iadr,iptr,r8adr,r8ptr,ier)
            if(ier.ne.0) return

         endif
      enddo

      if(iadr.ne.imaxa) then
         call xplasma_errmsg_append(s, &
              '?xplasma_unpackmaster: iadr.ne.imaxa at end.')
         ier=9999
      endif

      if(r8adr.ne.r8maxa) then
         call xplasma_errmsg_append(s, &
              '?xplasma_unpackmaster: r8adr.ne.r8maxa at end.')
         ier=9999
      endif

    end subroutine xplasma_unpackmaster

    ! private routines to support pack/unpack

    subroutine xplist_unpack(xx, iadr, iptr, r8adr, r8ptr, ier)

      !  pack a single list

      type (xplist) :: xx
      integer :: iadr
      integer, dimension(:) :: iptr
      integer :: r8adr
      real*8, dimension(:) :: r8ptr
      integer :: ier

      !------------
      integer :: ielems,nelems,ichsize
      !------------

      call unpack_iscalar(nelems, iadr, iptr, ier)
      if(ier.ne.0) return

      xx%size = nelems

      allocate(xx%enames(nelems),xx%chindx(2,nelems))
      allocate(xx%ivals(nelems),xx%r8vals(nelems))

      do ielems = 1,nelems
         call unpack_str(xx%enames(ielems), iadr, iptr, ier)
         if(ier.ne.0) return
      enddo

      do ielems = 1,nelems
         call unpack_iscalar(xx%chindx(1,ielems), iadr, iptr, ier)
         if(ier.ne.0) return
         call unpack_iscalar(xx%chindx(2,ielems), iadr, iptr, ier)
         if(ier.ne.0) return
      enddo

      if(nelems.gt.0) then
         call unpack_iscalar(ichsize, iadr, iptr, ier)
         if(ier.ne.0) return

         if(ichsize.gt.0) then
            allocate(xx%chbuf(ichsize))
            call unpack_c1a(xx%chbuf, iadr, iptr, ier)
            if(ier.ne.0) return
         endif

         call unpack_ivec(xx%ivals(1:nelems), iadr, iptr, ier)
         if(ier.ne.0) return

         call unpack_r8vec(xx%r8vals(1:nelems), r8adr, r8ptr, ier)
         if(ier.ne.0) return

      endif

    end subroutine xplist_unpack

    subroutine xpcoord_unpack(xx, iadr, iptr, r8adr, r8ptr, ier)

      !  pack a single coord

      type (xpcoord) :: xx
      integer :: iadr
      integer, dimension(:) :: iptr
      integer :: r8adr
      real*8, dimension(:) :: r8ptr
      integer :: ier

      !------------------
      integer :: ilogical, iisize
      !------------------

      call unpack_iscalar(ilogical, iadr, iptr, ier)
      if(ier.ne.0) return
      xx%periodic = (ilogical.eq.1)

      call unpack_iscalar(xx%ngridc, iadr, iptr, ier)
      if(ier.ne.0) return

      call unpack_iscalar(xx%ngrid_max, iadr, iptr, ier)
      if(ier.ne.0) return

      if(xx%ngrid_max.gt.0) then
         call unpack_iscalar(iisize, iadr, iptr, ier)
         if(ier.ne.0) return

         allocate(xx%grid_ids(iisize))
         call unpack_ivec(xx%grid_ids, iadr, iptr, ier)
         if(ier.ne.0) return
      endif

    end subroutine xpcoord_unpack

    subroutine xpgrid_unpack(xx, iadr, iptr, r8adr, r8ptr, ier)

      !  pack a single grid

      type (xpgrid) :: xx
      integer :: iadr
      integer, dimension(:) :: iptr
      integer :: r8adr
      real*8, dimension(:) :: r8ptr
      integer :: ier

      !---------------
      integer :: ii,isize
      !---------------

      call unpack_iscalar(isize, iadr, iptr, ier)
      if(ier.ne.0) return
      xx%size = isize

      call unpack_iscalar(xx%coord, iadr, iptr, ier)
      if(ier.ne.0) return

      call unpack_iscalar(xx%nrefs, iadr, iptr, ier)
      if(ier.ne.0) return

      if(isize.gt.0) then
         allocate(xx%xpkg(isize,4))
         do ii=1,4
            call unpack_r8vec(xx%xpkg(1:isize,ii), r8adr, r8ptr, ier)
            if(ier.ne.0) return
         enddo
      endif

    end subroutine xpgrid_unpack

    subroutine xpprof_unpack(xx, iadr, iptr, r8adr, r8ptr, ier)

      !  pack a single prof

      type (xpprof) :: xx
      integer :: iadr
      integer, dimension(:) :: iptr
      integer :: r8adr
      real*8, dimension(:) :: r8ptr
      integer :: ier

      !----------------

      call unpack_iscalar(xx%rank, iadr, iptr, ier)
      if(ier.ne.0) return

      call unpack_ivec(xx%gridIds, iadr, iptr, ier)
      if(ier.ne.0) return
      call unpack_ivec(xx%profIds, iadr, iptr, ier)
      if(ier.ne.0) return

      call unpack_iscalar(xx%kspline, iadr, iptr, ier)
      if(ier.ne.0) return
      call unpack_iscalar(xx%prof_counter, iadr, iptr, ier)
      if(ier.ne.0) return
      call unpack_iscalar(xx%size, iadr, iptr, ier)
      if(ier.ne.0) return

      if(xx%size.gt.0) then

         allocate(xx%buf(xx%size))
         call unpack_r8vec(xx%buf, r8adr, r8ptr, ier)
         if(ier.ne.0) return

      endif

    end subroutine xpprof_unpack

    subroutine xpblkbx_unpack(xx, iadr, iptr, r8adr, r8ptr, ier)

      !  pack a single blkbx

      type (xpblkbx) :: xx
      integer :: iadr
      integer, dimension(:) :: iptr
      integer :: r8adr
      real*8, dimension(:) :: r8ptr
      integer :: ier

      integer :: isize
      real*8 :: zsize

      call unpack_iscalar(xx%type, iadr, iptr, ier)

      if(iadr.ge.size(iptr)) then
         write(ilun_dbg,*) ' ?iptr exceeded in xpblkbx_unpack'
         ier=1
         return
      endif

      if(r8adr.ge.size(r8ptr)) then
         write(ilun_dbg,*) ' ?r8ptr exceeded in xpblkbx_unpack'
         ier=1
         return
      endif

      isize = iptr(iadr+1)
      if(isize.gt.0) then
         allocate(xx%iwork(isize))
         call unpack_ivec(xx%iwork, iadr, iptr, ier)
         if(ier.ne.0) return
      endif

      zsize = r8ptr(r8adr+1)
      if((zsize.lt.0.0d0).or.(zsize.gt.2.0d9)) then
         write(ilun_dbg,*) ' ?r8ptr size value ',zsize,' sanity check failed.'
         ier=1
         return
      endif

      isize=zsize
      if(isize.gt.0) then
         allocate(xx%r8work(isize))
         call unpack_r8vec(xx%r8work, r8adr, r8ptr, ier)
         if(ier.ne.0) return
      endif

    end subroutine xpblkbx_unpack
 
    subroutine ipack_xpand(ineed,iptr)
      
      integer, intent(in) :: ineed    ! minimal buffer size needed
      integer, dimension(:), pointer :: iptr

      !----------
      integer :: isize,isize_new
      integer, dimension(:), pointer :: jptr
      !----------

      isize = size(iptr)
      if(ineed.le.isize) return

      !  expand to a larger buffer

      jptr => iptr

      isize_new = isize
      do
         isize_new = isize_new*2
         if(isize_new.ge.ineed) exit
      enddo

      allocate(iptr(isize_new))
      iptr(1:isize)=jptr(1:isize)
      deallocate(jptr)

      iptr(isize+1:isize_new)=0

    end subroutine ipack_xpand

    subroutine r8pack_xpand(ineed,r8ptr)
      
      integer, intent(in) :: ineed    ! minimal buffer size needed
      real*8, dimension(:), pointer :: r8ptr

      !----------
      integer :: isize,isize_new
      real*8, dimension(:), pointer :: r8tmp
      !----------

      isize = size(r8ptr)
      if(ineed.le.isize) return

      !  expand to a larger buffer

      r8tmp => r8ptr

      isize_new = isize
      do
         isize_new = isize_new*2
         if(isize_new.ge.ineed) exit
      enddo

      allocate(r8ptr(isize_new))
      r8ptr(1:isize)=r8tmp(1:isize)
      deallocate(r8tmp)

      r8ptr(isize+1:isize_new)=0

    end subroutine r8pack_xpand
      
    subroutine pack_c1a(c1a,iadr,iptr)

      !  pack var. length character array

      character, dimension(:), intent(in) :: c1a
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      integer, dimension(:), pointer :: iptr   ! expandable data array

      !----------------
      integer :: ic,c1size,ineed
      !----------------

      c1size = size(c1a)
      ineed = iadr + c1size + 1

      call ipack_xpand(ineed,iptr)

      iadr = iadr + 1
      iptr(iadr) = c1size
      do ic=1,c1size
         iptr(iadr+ic) = ichar(c1a(ic))
      enddo

      iadr = iadr + c1size

    end subroutine pack_c1a
      
    subroutine pack_str(str,iadr,iptr)

      !  pack character*(*) array

      character*(*), intent(in) :: str
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      integer, dimension(:), pointer :: iptr   ! expandable data array

      !----------------
      integer :: ic,c1size,ineed
      !----------------

      c1size = len(str)
      ineed = iadr + c1size + 1

      call ipack_xpand(ineed,iptr)

      iadr = iadr + 1
      iptr(iadr) = c1size
      do ic=1,c1size
         iptr(iadr+ic) = ichar(str(ic:ic))
      enddo

      iadr = iadr + c1size

    end subroutine pack_str
       
    subroutine pack_ivec(ivec,iadr,iptr)

      !  pack integer array

      integer, dimension(:), intent(in) :: ivec
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      integer, dimension(:), pointer :: iptr   ! expandable data array

      !----------------
      integer :: ii,iisize,ineed
      !----------------

      iisize = size(ivec)
      ineed = iadr + iisize + 1

      call ipack_xpand(ineed,iptr)

      iadr = iadr + 1
      iptr(iadr) = iisize
      do ii=1,iisize
         iptr(iadr+ii) = ivec(ii)
      enddo

      iadr = iadr + iisize

    end subroutine pack_ivec
       
    subroutine pack_r8vec(r8vec,iadr,r8ptr)

      !  pack real*8 array

      real*8, dimension(:), intent(in) :: r8vec
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      real*8, dimension(:), pointer :: r8ptr   ! expandable data array

      !----------------
      integer :: ii,iisize,ineed
      !----------------

      iisize = size(r8vec)
      ineed = iadr + iisize + 1

      call r8pack_xpand(ineed,r8ptr)

      iadr = iadr + 1
      r8ptr(iadr) = iisize + 0.1d0
      do ii=1,iisize
         r8ptr(iadr+ii) = r8vec(ii)
      enddo

      iadr = iadr + iisize

    end subroutine pack_r8vec
     
    subroutine pack_iscalar(intval,iadr,iptr)

      !  pack single integer

      integer, intent(in) :: intval
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      integer, dimension(:), pointer :: iptr   ! expandable data array

      !----------------
      integer :: ineed
      !----------------

      ineed = iadr  + 1

      call ipack_xpand(ineed,iptr)

      iadr = iadr + 1
      iptr(iadr) = intval

    end subroutine pack_iscalar
      
    subroutine pack_r8scalar(r8val,iadr,r8ptr)

      !  pack single real*8 number

      real*8, intent(in) :: r8val
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      real*8, dimension(:), pointer :: r8ptr   ! expandable data array

      !----------------
      integer :: ineed
      !----------------

      ineed = iadr  + 1

      call r8pack_xpand(ineed,r8ptr)

      iadr = iadr + 1
      r8ptr(iadr) = r8val

    end subroutine pack_r8scalar

    !-----------------

    subroutine unpack_c1a(c1a,iadr,iarray,ier)

      !  unpack var. length character array

      character, dimension(:), intent(out) :: c1a
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      integer, dimension(:), intent(in)  :: iarray   ! data array

      integer, intent(out) :: ier   ! status code 0=OK

      !-------------------
      integer :: ic,c1size
      !-------------------

      ier = 0

      c1size=size(c1a)

      iadr = iadr + 1
      if(iadr.gt.size(iarray)) then
         ier=1
         write(ilun_dbg,*) ' ?unpack_c1a: iarray exhausted.'
         return
      endif
      if(c1size.ne.iarray(iadr)) then
         write(ilun_dbg,*) ' ?unpack_c1a -- size mismatch, expected: ',iarray(iadr), &
              ' got: ',c1size
         ier = 1
         return
      endif

      do ic=1,c1size
         c1a(ic) = char(iarray(iadr+ic))
      enddo

      iadr = iadr + c1size

    end subroutine unpack_c1a

    subroutine unpack_str(str,iadr,iarray,ier)

      !  unpack character*(*) array

      character*(*), intent(out) :: str
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      integer, dimension(:), intent(in)  :: iarray   ! data array

      integer, intent(out) :: ier   ! status code 0=OK

      !-------------------
      integer :: ic,c1size
      !-------------------

      ier = 0

      c1size=len(str)

      iadr = iadr + 1
      if(iadr.gt.size(iarray)) then
         ier=1
         write(ilun_dbg,*) ' ?unpack_str: iarray exhausted.'
         return
      endif
      if(c1size.ne.iarray(iadr)) then
         write(ilun_dbg,*) ' ?unpack_str -- size mismatch, expected: ',iarray(iadr), &
              ' got: ',c1size
         ier = 1
         return
      endif

      do ic=1,c1size
         str(ic:ic) = char(iarray(iadr+ic))
      enddo

      iadr = iadr + c1size

    end subroutine unpack_str

    subroutine unpack_ivec(ivec,iadr,iarray,ier)

      !  unpack integer array

      integer, dimension(:), intent(out) :: ivec
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      integer, dimension(:), intent(in)  :: iarray   ! data array

      integer, intent(out) :: ier   ! status code 0=OK

      !-------------------
      integer :: ii,iisize
      !-------------------

      ier = 0

      iisize=size(ivec)

      iadr = iadr + 1
      if(iadr.gt.size(iarray)) then
         ier=1
         write(ilun_dbg,*) ' ?unpack_ivec: iarray exhausted.'
         return
      endif
      if(iisize.ne.iarray(iadr)) then
         write(ilun_dbg,*) ' ?unpack_ivec -- size mismatch, expected: ',iarray(iadr),&
              ' got: ',iisize
         ier = 1
         return
      endif

      do ii=1,iisize
         ivec(ii) = iarray(iadr+ii)
      enddo

      iadr = iadr + iisize

    end subroutine unpack_ivec

    subroutine unpack_r8vec(r8vec,iadr,r8array,ier)

      !  unpack real*8 array

      real*8, dimension(:), intent(out) :: r8vec
      integer, intent(inout) :: iadr
      !  on input, address of last word used;
      !  on output, address of last word used, incremented.

      real*8, dimension(:), intent(in)  :: r8array   ! data array

      integer, intent(out) :: ier   ! status code 0=OK

      !-------------------
      integer :: ii,iisize,iasize
      !-------------------

      ier = 0

      iisize=size(r8vec)

      iadr = iadr + 1
      if(iadr.gt.size(r8array)) then
         ier=1
         write(ilun_dbg,*) ' ?unpack_r8vec: r8array exhausted.'
         return
      endif
      iasize = r8array(iadr)
      if(iisize.ne.iasize) then
         write(ilun_dbg,*) ' ?unpack_r8vec -- size mismatch, expected: ',iasize,&
              ' got: ',iisize
         ier = 1
         return
      endif

      do ii=1,iisize
         r8vec(ii) = r8array(iadr+ii)
      enddo

      iadr = iadr + iisize

    end subroutine unpack_r8vec

    subroutine unpack_iscalar(intval,iadr,iarray,ier)

      ! unpack integer scalar

      integer, intent(out) :: intval
      integer, intent(inout) :: iadr   ! incremented
      integer, dimension(:), intent(in) :: iarray
      integer, intent(out) :: ier

      ier = 0

      iadr = iadr + 1
      if(iadr.gt.size(iarray)) then
         ier=1
         write(ilun_dbg,*) ' ?unpack_iscalar: iarray exhausted.'
         return
      endif

      intval = iarray(iadr)

    end subroutine unpack_iscalar

    subroutine unpack_r8scalar(r8val,iadr,r8array,ier)

      ! unpack real*8 scalar

      real*8, intent(out) :: r8val
      integer, intent(inout) :: iadr   ! incremented
      real*8, dimension(:), intent(in) :: r8array
      integer, intent(out) :: ier

      ier = 0

      iadr = iadr + 1
      if(iadr.gt.size(r8array)) then
         ier=1
         write(ilun_dbg,*) ' ?unpack_iscalar: r8array exhausted.'
         return
      endif

      r8val = r8array(iadr)

    end subroutine unpack_r8scalar
       
    !----------------------------------

    subroutine xpblkbx_fullcopy(d1,d2)
      ! copy integ object d1->d2; full memory copy not just ptrs
      ! this is private method; d2 object initialization handled elsewhere.

      type(xpblkbx) :: d1,d2

      d2%type = d1%type

      if(allocated(d1%iwork)) then
         allocate(d2%iwork(size(d1%iwork)))
         d2%iwork = d1%iwork
      else
         if(allocated(d2%iwork)) deallocate(d2%iwork)
      endif

      if(allocated(d1%r8work)) then
         allocate(d2%r8work(size(d1%r8work)))
         d2%r8work = d1%r8work
      else
         if(allocated(d2%r8work)) deallocate(d2%r8work)
      endif

    end subroutine xpblkbx_fullcopy

    !----------------------------------

    subroutine xpblkbx_free(d)
      ! free memory associated with xpblkbx object

      type(xpblkbx) :: d

      d%type = 0
      if(allocated(d%iwork)) deallocate(d%iwork)
      if(allocated(d%r8work)) deallocate(d%r8work)

    end subroutine xpblkbx_free

    !----------------------------------
    subroutine xpblkbx_ncdefine(d,prefix,icdf,ier)
      !  Define xpblkbx object output records for NetCDF file

      type(xpblkbx) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      !----------------------------

      ier = 0

      call cdf_define(icdf,trim(prefix)//'%type',d%type)

      if(d%type.gt.0) then
         call cdf_define(icdf,trim(prefix)//'%iwork',d%iwork)
         call cdf_define(icdf,trim(prefix)//'%r8work',d%r8work)
      end if

    end subroutine xpblkbx_ncdefine

    !----------------------------------
    subroutine xpblkbx_ncwrite(d,prefix,icdf,ier)
      !  Write xpblkbx object output records for NetCDF file

      type(xpblkbx) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      ier = 0

      call cdf_write(icdf,trim(prefix)//'%type',d%type)

      if(d%type.gt.0) then
         call cdf_write(icdf,trim(prefix)//'%iwork',d%iwork)
         call cdf_write(icdf,trim(prefix)//'%r8work',d%r8work)
      end if

    end subroutine xpblkbx_ncwrite

    !----------------------------------
    subroutine xpblkbx_ncread(d,prefix,icdf,ier)
      !  Read xpblkbx object records from NetCDF file

      type(xpblkbx) :: d
      character*(*), intent(in) :: prefix  ! naming prefix
      integer, intent(in) :: icdf   ! handle of open NetCDF file
      integer, intent(out) :: ier

      !----------------------------
      integer :: idims(3),istat
      character*8 :: ztype
      !----------------------------

      ier = 0

      call cdf_read(icdf,trim(prefix)//'%type',d%type,ier)
      if(ier.ne.0) then
         ier=13
         return
      endif

      if(d%type.gt.0) then
         call cdf_inquire(icdf,trim(prefix)//'%iwork',idims,ztype,istat)
         if(istat.ne.0) then
            if(ier.ne.0) ier=13
            if(ier.ne.0) return
         endif
         allocate(d%iwork(idims(1)))
         call cdf_read(icdf,trim(prefix)//'%iwork',d%iwork,ier)
         if(ier.ne.0) then
            ier=13
            return
         endif

         call cdf_inquire(icdf,trim(prefix)//'%r8work',idims,ztype,istat)
         if(istat.ne.0) then
            if(ier.ne.0) ier=13
            if(ier.ne.0) return
         endif
         allocate(d%r8work(idims(1)))
         call cdf_read(icdf,trim(prefix)//'%r8work',d%r8work,ier)
         if(ier.ne.0) then
            ier=13
            return
         endif
      endif

    end subroutine xpblkbx_ncread

    !========================================
    !========================================
    !========================================
    !========================================
    !========================================
    !========================================
    !========================================
    subroutine xplasma_ncdefine(s,icdf,ier)

      ! routine to define NetCDF output quantities

      type (xplasma), pointer :: s
      integer, intent(in) :: icdf   ! NetCDF i/o channel
      integer, intent(out) :: ier   ! status code on exit (0=OK)

      !-----------------------------------
      integer :: ilogical = 0
      integer :: i,j,lentot,lentotu,iauths
      character*32, dimension(:), allocatable :: zdnames,zauthors
      integer, dimension(:), allocatable :: zdtypes
      integer, dimension(:), allocatable :: zdindex
      character, dimension(:), allocatable :: zlabels,zunits
      integer, dimension(:), allocatable :: len_labels,len_units

      integer :: k_tmp
      !-----------------------------------

      ier = 0

      call cdf_define(icdf,'__Label',s%label)

      call cdf_define(icdf,'__Axisymm',ilogical)
      call cdf_define(icdf,'__Scrapeoff',ilogical)
      call cdf_define(icdf,'__Bphi_CCW',s%bphi_ccw)
      call cdf_define(icdf,'__Jphi_CCW',s%jphi_ccw)

      call cdf_define(icdf,'__Time',s%time)

      iauths=max(2,s%cur_author)
      allocate(zauthors(iauths))
      call cdf_define(icdf,'__author_stack',zauthors)
      call cdf_define(icdf,'__write_enable',s%write_enable(1:iauths))
      deallocate(zauthors)

      call cdf_define(icdf,'__nlists',s%nlists)
      call cdf_define(icdf,'__nlist_free',s%nlist_free)
      if(s%nlists.gt.0) then
         call cdf_define(icdf,'__list_freelist',s%list_freelist)
      endif

      call cdf_define(icdf,'__ncoords',s%ncoords)
      call cdf_define(icdf,'__ncoord_free',s%ncoord_free)
      if(s%ncoords.gt.0) then
         call cdf_define(icdf,'__coord_freelist',s%coord_freelist)
      endif

      call cdf_define(icdf,'__ngrids',s%ngrids)
      call cdf_define(icdf,'__ngrid_free',s%ngrid_free)
      if(s%ngrids.gt.0) then
         call cdf_define(icdf,'__grid_freelist',s%grid_freelist)
      endif

      call cdf_define(icdf,'__nprofs',s%nprofs)
      call cdf_define(icdf,'__nprof_free',s%nprof_free)
      if(s%nprofs.gt.0) then
         call cdf_define(icdf,'__prof_freelist',s%prof_freelist)
      endif

      call cdf_define(icdf,'__nblkbxs',s%nblkbxs)
      call cdf_define(icdf,'__nblkbx_free',s%nblkbx_free)
      if(s%nblkbxs.gt.0) then
         call cdf_define(icdf,'__blkbx_freelist',s%blkbx_freelist)
      endif

      call cdf_define(icdf,'__prof_counter',s%prof_counter)
      call cdf_define(icdf,'__eq_counter',s%eq_counter)
      call cdf_define(icdf,'__bdytol',s%bdytol)
      call cdf_define(icdf,'__ajac_maxVar',s%ajac_maxVar)
      call cdf_define(icdf,'__rzOrder',s%rzOrder)
      call cdf_define(icdf,'__kmom',s%kmom)

      if(s%nitems.gt.0) then

         !  package up dictionary into conventional type arrays
         !  and write them out...

         allocate(zdnames(s%nitems),zdtypes(s%nitems),zdindex(s%nitems))
         allocate(len_labels(s%nitems),len_units(s%nitems),zauthors(s%nitems))

         lentot=0
         do i=1,s%nitems
            if(allocated(s%dict(i)%label)) then
               len_labels(i)=size(s%dict(i)%label)
            else
               len_labels(i)=0
            endif
            lentot=lentot+len_labels(i)
         enddo
         if(lentot.gt.0) then
            allocate(zlabels(lentot))
         endif

         lentotu=0
         do i=1,s%nitems
            if(allocated(s%dict(i)%units)) then
               len_units(i)=size(s%dict(i)%units)
            else
               len_units(i)=0
            endif
            lentotu=lentotu+len_units(i)
         enddo
         if(lentotu.gt.0) then
            allocate(zunits(lentotu))
         endif

         call cdf_define(icdf,'__Order',s%order(1:s%nitems))
         call cdf_define(icdf,'__Dict_names',zdnames)
         call cdf_define(icdf,'__Dict_authors',zauthors)
         if(lentot.gt.0) then
            call cdf_define(icdf,'__Dict_labels',zlabels)
         endif
         call cdf_define(icdf,'__Dict_len_labels',len_labels)
         if(lentotu.gt.0) then
            call cdf_define(icdf,'__Dict_units',zunits)
         endif
         call cdf_define(icdf,'__Dict_len_units',len_units)
         call cdf_define(icdf,'__Dict_types',zdtypes)
         call cdf_define(icdf,'__Dict_indices',zdindex)

         deallocate(zdnames,zdtypes,zdindex,len_labels,len_units,zauthors)
         if(lentot.gt.0) then
            deallocate(zlabels)
         endif
         if(lentotu.gt.0) then
            deallocate(zunits)
         endif

         !  and define the individual items contents

         do i=1,s%nitems
            j=s%dict(i)%dindex
            
            k_tmp = s%dict(i)%dtype
            if(s%dict(i)%dtype.eq.xplasma_listType) then
               call xplist_ncdefine(s%lists(j),s%dict(i)%name,icdf,ier)

            else if( k_tmp .eq.xplasma_coordType) then
               call xpcoord_ncdefine(s%coords(j),s%dict(i)%name,icdf,ier)

            else if( k_tmp .eq.xplasma_gridType) then
               call xpgrid_ncdefine(s%grids(j),s%dict(i)%name,icdf,ier)

            else if( k_tmp .eq.xplasma_profType) then
               call xpprof_ncdefine(s%profs(j),s%dict(i)%name,icdf,ier)

            else if( k_tmp .eq.xplasma_blkbxType) then
               call xpblkbx_ncdefine(s%blkbxs(j),s%dict(i)%name,icdf,ier)

            endif
         enddo

      endif

    end subroutine xplasma_ncdefine

    !========================================
    subroutine xplasma_ncwrite(s,icdf,ier)

      ! routine to write NetCDF file containing entire xplasma object contents

      type (xplasma), pointer :: s
      integer, intent(in) :: icdf   ! NetCDF i/o channel
      integer, intent(out) :: ier   ! status code on exit (0=OK)

      !-----------------------------------
      integer :: ilogical
      integer :: i,j,lentot,lentotu,iauths,itchk
      character*32, dimension(:), allocatable :: zdnames,zauthors
      integer, dimension(:), allocatable :: zdtypes
      integer, dimension(:), allocatable :: zdindex
      character, dimension(:), allocatable :: zlabels,zunits
      integer, dimension(:), allocatable :: len_labels,len_units
      integer :: luntty=6
      character*32 :: tname
      integer :: k_tmp
      !-----------------------------------

      ier = 0

      if(xplasma_trace_name.ne.' ') then
         itchk=1
         tname=xplasma_trace_name
         call uupper(tname)
      else
         itchk=0
         tname=' '
      endif

      call cdf_write(icdf,'__Label',s%label)

      if(s%axisymm) then
         ilogical=1
      else
         ilogical=0
      endif
      call cdf_write(icdf,'__Axisymm',ilogical)

      if(s%scrapeoff) then
         ilogical=1
      else
         ilogical=0
      endif
      call cdf_write(icdf,'__Scrapeoff',ilogical)

      call cdf_write(icdf,'__Bphi_CCW',s%bphi_ccw)
      call cdf_write(icdf,'__Jphi_CCW',s%jphi_ccw)

      call cdf_write(icdf,'__Time',s%time)

      iauths=max(2,s%cur_author)
      allocate(zauthors(iauths))
      zauthors(1:s%cur_author) = s%authors(1:s%cur_author)
      if(s%cur_author.lt.2) zauthors(s%cur_author+1:2)=' '
      call cdf_write(icdf,'__author_stack',zauthors)
      call cdf_write(icdf,'__write_enable',s%write_enable(1:iauths))
      deallocate(zauthors)

      call cdf_write(icdf,'__nlists',s%nlists)
      call cdf_write(icdf,'__nlist_free',s%nlist_free)
      if(s%nlists.gt.0) then
         call cdf_write(icdf,'__list_freelist',s%list_freelist)
      endif

      call cdf_write(icdf,'__ncoords',s%ncoords)
      call cdf_write(icdf,'__ncoord_free',s%ncoord_free)
      if(s%ncoords.gt.0) then
         call cdf_write(icdf,'__coord_freelist',s%coord_freelist)
      endif

      call cdf_write(icdf,'__ngrids',s%ngrids)
      call cdf_write(icdf,'__ngrid_free',s%ngrid_free)
      if(s%ngrids.gt.0) then
         call cdf_write(icdf,'__grid_freelist',s%grid_freelist)
      endif

      call cdf_write(icdf,'__nprofs',s%nprofs)
      call cdf_write(icdf,'__nprof_free',s%nprof_free)
      if(s%nprofs.gt.0) then
         call cdf_write(icdf,'__prof_freelist',s%prof_freelist)
      endif

      call cdf_write(icdf,'__nblkbxs',s%nblkbxs)
      call cdf_write(icdf,'__nblkbx_free',s%nblkbx_free)
      if(s%nblkbxs.gt.0) then
         call cdf_write(icdf,'__blkbx_freelist',s%blkbx_freelist)
      endif

      call cdf_write(icdf,'__prof_counter',s%prof_counter)
      call cdf_write(icdf,'__eq_counter',s%eq_counter)
      call cdf_write(icdf,'__bdytol',s%bdytol)
      call cdf_write(icdf,'__ajac_maxVar',s%ajac_maxVar)
      call cdf_write(icdf,'__rzOrder',s%rzOrder)
      call cdf_write(icdf,'__kmom',s%kmom)

      if(s%nitems.gt.0) then

         !  here the structures to build the output file are allocated

         allocate(zdnames(s%nitems),zdtypes(s%nitems),zdindex(s%nitems))
         allocate(len_labels(s%nitems),len_units(s%nitems),zauthors(s%nitems))
         lentot=0
         do i=1,s%nitems
            if(allocated(s%dict(i)%label)) then
               len_labels(i)=size(s%dict(i)%label)
            else
               len_labels(i)=0
            endif
            lentot=lentot+len_labels(i)
         enddo
         if(lentot.gt.0) then
            allocate(zlabels(lentot))
         endif

         lentotu=0
         do i=1,s%nitems
            if(allocated(s%dict(i)%units)) then
               len_units(i)=size(s%dict(i)%units)
            else
               len_units(i)=0
            endif
            lentotu=lentotu+len_units(i)
         enddo
         if(lentotu.gt.0) then
            allocate(zunits(lentotu))
         endif

         lentot=0
         lentotu=0
         do i=1,s%nitems
            if(len_labels(i).gt.0) then
               zlabels(lentot+1:lentot+len_labels(i)) = &
                    s%dict(i)%label(1:len_labels(i))
               lentot = lentot + len_labels(i)
            endif
            if(len_units(i).gt.0) then
               zunits(lentotu+1:lentotu+len_units(i)) = &
                    s%dict(i)%units(1:len_units(i))
               lentotu = lentotu + len_units(i)
            endif
            zdnames(i) = s%dict(i)%name
            zauthors(i)= s%dict(i)%author
            zdtypes(i) = s%dict(i)%dtype
            zdindex(i) = s%dict(i)%dindex
            if(itchk.gt.0) then
               if(s%dict(i)%dindex.gt.0) then
                  if(tname.eq.s%dict(i)%name) then
                     call tchk(s,'written to NetCDF',tname)
                  endif
               endif
            endif
         enddo

         call cdf_write(icdf,'__Order',s%order(1:s%nitems))
         call cdf_write(icdf,'__Dict_names',zdnames)
         call cdf_write(icdf,'__Dict_authors',zauthors)
         call cdf_write(icdf,'__Dict_len_labels',len_labels)
         if(lentot.gt.0) then
            call cdf_write(icdf,'__Dict_labels',zlabels)
         endif
         call cdf_write(icdf,'__Dict_len_units',len_units)
         if(lentotu.gt.0) then
            call cdf_write(icdf,'__Dict_units',zunits)
         endif
         call cdf_write(icdf,'__Dict_types',zdtypes)
         call cdf_write(icdf,'__Dict_indices',zdindex)

         deallocate(zdnames,zdtypes,zdindex,len_labels,len_units)
         if(lentot.gt.0) then
            deallocate(zlabels)
         endif
         if(lentotu.gt.0) then
            deallocate(zunits)
         endif

         !  and write the individual items contents

         do i=1,s%nitems
            j=s%dict(i)%dindex

            k_tmp = s%dict(i)%dtype
            if( k_tmp .eq.xplasma_listType) then
               call xplist_ncwrite(s%lists(j),s%dict(i)%name,icdf,ier)

            else if( k_tmp .eq.xplasma_coordType) then
               call xpcoord_ncwrite(s%coords(j),s%dict(i)%name,icdf,ier)

            else if( k_tmp .eq.xplasma_gridType) then
               call xpgrid_ncwrite(s%grids(j),s%dict(i)%name,icdf,ier)

            else if( k_tmp .eq.xplasma_profType) then
               call xpprof_ncwrite(s%profs(j),s%dict(i)%name,icdf,ier)

            else if( k_tmp .eq.xplasma_blkbxType) then
               call xpblkbx_ncwrite(s%blkbxs(j),s%dict(i)%name,icdf,ier)

            endif
         enddo
      endif

    end subroutine xplasma_ncwrite

    !========================================
    subroutine xplasma_ncread(s,icdf,ier)

      ! routine to load entire xplasma object contents from NetCDF file

      type (xplasma), pointer :: s
      integer, intent(in) :: icdf   ! NetCDF i/o channel
      integer, intent(out) :: ier   ! status code on exit (0=OK)

      !-----------------------------------
      integer :: ilogical = 0
      integer :: i,j,lentot,lentotu,idims(3),istat,iauths,jauths
      integer :: itchk
      integer :: id_xmom2d
      character*32, dimension(:), allocatable :: zdnames,zauthors
      integer, dimension(:), allocatable :: zdtypes
      integer, dimension(:), allocatable :: zdindex
      character, dimension(:), allocatable :: zlabels,zunits
      integer, dimension(:), allocatable :: len_labels,len_units
      character*8 ztype
      character*32 tname

      integer :: k_tmp
      integer :: luntty = 6
      !-----------------------------------

      ier = 0

      if(xplasma_trace_name.ne.' ') then
         itchk=1
         tname=xplasma_trace_name
         call uupper(tname)
      else
         itchk=0
         tname=' '
      endif

      if(s%nguard.ne.0) call xplasma_freesub(s,ier,total=.TRUE.)
      if(ier.ne.0) return

      internal_counter = internal_counter + 1
      s%nguard=internal_counter

      if(xplasma_trace_name.ne.' ') then
         write(luntty,*) ' %xplasma_obj(xplasma_NCREAD): object init by read: '
         write(luntty,*) '  internal ID: ',s%nguard
      endif

      call cdf_inquire(icdf,'__Label',idims,ztype,istat)
      if(istat.eq.0) then
         lentot=max(1,idims(1))*max(1,idims(2))
         if(allocated(s%label)) deallocate(s%label)
         allocate(s%label(max(2,lentot)))
         if(lentot.eq.1) then
            s%label='??'
            ier=0
         else
            call cdf_read(icdf,'__Label',s%label,ier)
         endif
      endif
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,'__Axisymm',ilogical,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return
      s%axisymm = (ilogical.eq.1)

      call cdf_read(icdf,'__Scrapeoff',ilogical,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return
      s%scrapeoff = (ilogical.eq.1)

      call cdf_read(icdf,'__Bphi_CCW',s%bphi_ccw)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,'__Jphi_CCW',s%jphi_ccw)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,'__Time',s%time,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_inquire(icdf,'__author_stack',idims,ztype,istat)
      if(istat.eq.0) then
         iauths = idims(2)
         allocate(zauthors(iauths))
         call cdf_read(icdf,'__author_stack',zauthors)
         jauths=iauths
         do
            if(jauths.eq.1) exit
            if(zauthors(jauths).eq.' ') then
               jauths=jauths-1
            else
               exit
            endif
         enddo
         call xplasma_authorStack_expand(s,jauths)
         s%cur_author=jauths
         s%authors(1:s%cur_author)=zauthors(1:jauths)
         call cdf_read(icdf,'__write_enable',s%write_enable(1:iauths))
         deallocate(zauthors)
      endif
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,'__nlists',s%nlists,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return
      call cdf_read(icdf,'__nlist_free',s%nlist_free,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      if(allocated(s%lists)) deallocate(s%lists)
      if(allocated(s%list_freelist)) deallocate(s%list_freelist)
      if(s%nlists.gt.0) then
         allocate(s%lists(s%nlists))
         allocate(s%list_freelist(s%nlists))
         do i=1,s%nlists
            call xplist_init(s%lists(i))
         enddo
         call cdf_read(icdf,'__list_freelist',s%list_freelist,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return
      endif

      call cdf_read(icdf,'__ncoords',s%ncoords,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return
      call cdf_read(icdf,'__ncoord_free',s%ncoord_free,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      if(allocated(s%coords)) deallocate(s%coords)
      if(allocated(s%coord_freelist)) deallocate(s%coord_freelist)
      if(s%ncoords.gt.0) then
         allocate(s%coords(s%ncoords))
         allocate(s%coord_freelist(s%ncoords))
         do i=1,s%ncoords
            call xpcoord_init(s%coords(i))
         enddo
         call cdf_read(icdf,'__coord_freelist',s%coord_freelist,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return
      endif

      call cdf_read(icdf,'__ngrids',s%ngrids,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return
      call cdf_read(icdf,'__ngrid_free',s%ngrid_free,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      if(allocated(s%grids)) deallocate(s%grids)
      if(allocated(s%grid_freelist)) deallocate(s%grid_freelist)
      if(s%ngrids.gt.0) then
         allocate(s%grids(s%ngrids))
         allocate(s%grid_freelist(s%ngrids))
         do i=1,s%ngrids
            call xpgrid_init(s%grids(i))
         enddo
         call cdf_read(icdf,'__grid_freelist',s%grid_freelist,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return
      endif

      call cdf_read(icdf,'__nprofs',s%nprofs,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return
      call cdf_read(icdf,'__nprof_free',s%nprof_free,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      if(allocated(s%profs)) deallocate(s%profs)
      if(allocated(s%prof_freelist)) deallocate(s%prof_freelist)
      if(s%nprofs.gt.0) then
         allocate(s%profs(s%nprofs))
         allocate(s%prof_freelist(s%nprofs))
         do i=1,s%nprofs
            call xpprof_init(s%profs(i))
         enddo
         call cdf_read(icdf,'__prof_freelist',s%prof_freelist,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return
      endif

      call cdf_read(icdf,'__nblkbxs',s%nblkbxs,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return
      call cdf_read(icdf,'__nblkbx_free',s%nblkbx_free,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      if(allocated(s%blkbxs)) deallocate(s%blkbxs)
      if(allocated(s%blkbx_freelist)) deallocate(s%blkbx_freelist)
      if(s%nblkbxs.gt.0) then
         allocate(s%blkbxs(s%nblkbxs))
         allocate(s%blkbx_freelist(s%nblkbxs))
         do i=1,s%nblkbxs
            call xpblkbx_init(s%blkbxs(i))
         enddo
         call cdf_read(icdf,'__blkbx_freelist',s%blkbx_freelist,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return
      endif

      call cdf_read(icdf,'__prof_counter',s%prof_counter,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,'__eq_counter',s%eq_counter,ier)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,'__bdytol',s%bdytol)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,'__ajac_maxVar',s%ajac_maxVar)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,'__rzOrder',s%rzOrder)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_read(icdf,'__kmom',s%kmom)
      if(ier.ne.0) ier=13
      if(ier.ne.0) return

      call cdf_inquire(icdf,'__Order',idims,ztype,istat)
      if(istat.eq.0) then
         s%nitems = idims(1)
         call xplasma_expand(s,idims(1),ier)
         if(ier.ne.0) return
      else
         s%nitems = 0
      endif

      if(s%nitems.gt.0) then

         call cdf_read(icdf,'__Order',s%order(1:s%nitems),ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return

         !  here the structures to read the NetCDF file are allocated

         allocate(zdnames(s%nitems),zdtypes(s%nitems),zdindex(s%nitems))
         allocate(len_labels(s%nitems),len_units(s%nitems),zauthors(s%nitems))

         call cdf_read(icdf,'__Dict_names',zdnames,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return

         call cdf_read(icdf,'__Dict_authors',zauthors,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return

         call cdf_read(icdf,'__Dict_len_labels',len_labels,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return

         call cdf_read(icdf,'__Dict_len_units',len_units,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return

         call cdf_read(icdf,'__Dict_types',zdtypes,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return

         call cdf_read(icdf,'__Dict_indices',zdindex,ier)
         if(ier.ne.0) ier=13
         if(ier.ne.0) return

         lentot=0
         lentotu=0
         do i=1,s%nitems
            lentot=lentot+len_labels(i)
            lentotu=lentotu+len_units(i)
         enddo

         if(lentot.gt.0) then
            allocate(zlabels(lentot))
            call cdf_read(icdf,'__Dict_labels',zlabels,ier)
            if(ier.ne.0) ier=13
            if(ier.ne.0) return
         endif

         if(lentotu.gt.0) then
            allocate(zunits(lentotu))
            call cdf_read(icdf,'__Dict_units',zunits,ier)
            if(ier.ne.0) ier=13
            if(ier.ne.0) return
         endif

         lentot=0
         lentotu=0
         do i=1,s%nitems
            s%dict(i)%name = zdnames(i)
            s%dict(i)%author = zauthors(i)
            s%dict(i)%dtype = zdtypes(i)
            s%dict(i)%dindex = zdindex(i)
            if(itchk.gt.0) then
               if(s%dict(i)%dindex.gt.0) then
                  if(tname.eq.s%dict(i)%name) then
                     call tchk(s,'read from NetCDF',tname)
                  endif
               endif
            endif
            if(len_labels(i).gt.0) then
               allocate(s%dict(i)%label(len_labels(i)))
               s%dict(i)%label(1:len_labels(i)) = &
                    zlabels(lentot+1:lentot+len_labels(i))
               lentot = lentot + len_labels(i)
            endif
            if(len_units(i).gt.0) then
               allocate(s%dict(i)%units(len_units(i)))
               s%dict(i)%units(1:len_units(i)) = &
                    zunits(lentotu+1:lentotu+len_units(i))
               lentotu = lentotu + len_units(i)
            endif
         enddo
         deallocate(zdnames,zauthors,zdtypes,zdindex,len_labels,len_units)
         if(lentot.gt.0) then
            deallocate(zlabels)
         endif
         if(lentotu.gt.0) then
            deallocate(zunits)
         endif

         !  and read the individual items contents

         do i=1,s%nitems
            j=s%dict(i)%dindex

            k_tmp = s%dict(i)%dtype
            if(s%dict(i)%dtype.eq.xplasma_listType) then
               call xplist_ncread(s%lists(j),s%dict(i)%name,icdf,ier)
               if(ier.ne.0) return

            else if( k_tmp .eq.xplasma_coordType) then
               call xpcoord_ncread(s%coords(j),s%dict(i)%name,icdf,ier)
               if(ier.ne.0) return

            else if( k_tmp .eq.xplasma_gridType) then
               call xpgrid_ncread(s%grids(j),s%dict(i)%name,icdf,ier)
               if(ier.ne.0) return

            else if( k_tmp .eq.xplasma_profType) then
               call xpprof_ncread(s%profs(j),s%dict(i)%name,icdf,ier)
               if(ier.ne.0) return

            else if( k_tmp .eq.xplasma_blkbxType) then
               call xpblkbx_ncread(s%blkbxs(j),s%dict(i)%name,icdf,ier)
               if(ier.ne.0) return

            endif

         enddo

      endif

      call xplasma_internal_ids(s)

      call xplasma_find_item(s,'__XMOM2D',id_xmom2d,ier,nf_noerr=.TRUE.)
      if(id_xmom2d.ne.0) call xpmom_fetch(s)

      call xplasma_ncread_finish(s)

    end subroutine xplasma_ncread

    subroutine xplasma_ncread_finish(s)
      !  private routine, post-processing at end of read
      !  guild grid_equiv list -- detect grids with nearly identical values

      type(xplasma) :: s

      !----------------
      integer :: i,j,jequiv,isize,jsize,ipt,iequiv
      logical :: imatch
      real*8 :: x1,x2,ztol,ztol1,ztol2
      !----------------

      if(allocated(s%grid_equiv)) deallocate(s%grid_equiv)
      if(s%ngrids.gt.0) then
         allocate(s%grid_equiv(s%ngrids)); s%grid_equiv = 0

         iequiv=0
         do i=1,s%ngrids
            isize=s%grids(i)%size
            if(isize.gt.0) then
               if(i.eq.1) then
                  iequiv=iequiv+1
                  s%grid_equiv(i)=iequiv
               else
                  ! check for matching grids
                  do j=1,i-1
                     imatch=.TRUE.
                     jsize=s%grids(j)%size
                     if(isize.ne.jsize) then
                        imatch=.FALSE.
                     else
                        x1=s%grids(i)%xpkg(1,1)
                        x2=s%grids(i)%xpkg(isize,1)
                        ztol1=1.0d-12*max(abs(x1),abs(x2))

                        x1=s%grids(j)%xpkg(1,1)
                        x2=s%grids(j)%xpkg(isize,1)
                        ztol2=1.0d-12*max(abs(x1),abs(x2))

                        ztol=min(ztol1,ztol2)

                        do ipt=1,isize
                           x1=s%grids(i)%xpkg(ipt,1)
                           x2=s%grids(j)%xpkg(ipt,1)
                           if(abs(x2-x1).gt.ztol) then
                              imatch=.FALSE.
                              exit
                           endif
                        enddo
                     endif
                     if(imatch) then
                        jequiv=j
                        exit
                     endif
                  enddo
                  if(.not.imatch) then
                     ! no match: new equivalence number
                     iequiv=iequiv+1
                     s%grid_equiv(i)=iequiv
                  else
                     ! match: copy equivalence number
                     jequiv = s%grid_equiv(jequiv)
                     s%grid_equiv(i)=jequiv
                  endif
               endif
            endif
         enddo
      endif

    end subroutine xplasma_ncread_finish

    !------------------------
    subroutine xpeval_info(xp,ierr, &
         ix,dxn,hx,hxi)

      !  return select info on xpeval object

      type(xpeval) :: xp
      integer, intent(out) :: ierr   ! status code on exit 0=OK

      integer, dimension(:), pointer, optional :: ix
      real*8, dimension(:), pointer, optional :: dxn
      real*8, dimension(:), pointer, optional :: hx
      real*8, dimension(:), pointer, optional :: hxi

      if(xp%idx.eq.0) then
         ierr=1   ! object was not initialized

      else
         ierr=0

         if(present(ix)) then
            ix => xp%ix
         endif

         if(present(dxn)) then
            dxn => xp%dxn
         endif

         if(present(hx)) then
            hx => xp%hx
         endif

         if(present(hxi)) then
            hxi => xp%hxi
         endif
      endif
    end subroutine xpeval_info

    subroutine xpdebug(s)
      
      !  internal debug/display routine

      type(xplasma), pointer :: s

      call xpdbg_lists(s%nlists,s%nlist_free,s%list_freelist,s%lists)
      call xpdbg_coords(s%ncoords,s%ncoord_free,s%coord_freelist,s%coords)
      call xpdbg_grids(s%ngrids,s%ngrid_free,s%grid_freelist,s%grids)
      call xpdbg_profs(s%nprofs,s%nprof_free,s%prof_freelist,s%profs)
      call xpdbg_blkbxs(s%nblkbxs,s%nblkbx_free,s%blkbx_freelist,s%blkbxs)

    end subroutine xpdebug

    subroutine tchk(s,action,aname)

      ! print message if aname matches xplasma_trace_name

      type(xplasma), pointer :: s
      character*(*), intent(in) :: action
      character*(*), intent(in) :: aname

      character*32 :: anam1,anam2

      !-----------------
      integer :: luntty = 6
      !-----------------

      anam1=xplasma_trace_name
      anam2=aname

      call uupper(anam1)
      call uupper(anam2)

      if(anam1.eq.anam2) then
         write(luntty,*) ' '
         write(luntty,*) ' --> '//trim(action)//': '//trim(anam1)
         write(luntty,*) '      xplasma object ID: ',s%nguard
      endif

    end subroutine tchk

    subroutine xpdbg_lists(n,nfree,freelist,lists)

      integer, intent(in) :: n,nfree
      integer, intent(in) :: freelist(n)
      type (xplist), intent(in) :: lists(n)

      integer :: ilun

      ilun=6
      write(ilun,*) ' #lists = ',n,' #free = ',nfree

    end subroutine xpdbg_lists

    subroutine xpdbg_coords(n,nfree,freelist,coords)

      integer, intent(in) :: n,nfree
      integer, intent(in) :: freelist(n)
      type (xpcoord), intent(in) :: coords(n)

      integer :: ilun

      ilun=6
      write(ilun,*) ' #coords = ',n,' #free = ',nfree

    end subroutine xpdbg_coords

    subroutine xpdbg_grids(n,nfree,freelist,grids)

      integer, intent(in) :: n,nfree
      integer, intent(in) :: freelist(n)
      type (xpgrid), intent(in) :: grids(n)

      integer :: ilun

      ilun=6
      write(ilun,*) ' #grids = ',n,' #free = ',nfree

    end subroutine xpdbg_grids

    subroutine xpdbg_profs(n,nfree,freelist,profs)

      integer, intent(in) :: n,nfree
      integer, intent(in) :: freelist(n)
      type (xpprof), intent(in) :: profs(n)

      integer :: ilun

      ilun=6
      write(ilun,*) ' #profs = ',n,' #free = ',nfree

    end subroutine xpdbg_profs

    subroutine xpdbg_blkbxs(n,nfree,freelist,blkbxs)

      integer, intent(in) :: n,nfree
      integer, intent(in) :: freelist(n)
      type (xpblkbx), intent(in) :: blkbxs(n)

      integer :: ilun

      ilun=6
      write(ilun,*) ' #blkbxs = ',n,' #free = ',nfree

    end subroutine xpdbg_blkbxs

end module xplasma_obj
