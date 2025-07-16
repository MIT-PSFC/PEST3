subroutine eqi_fromGeqdsk(filename, ihermite, inew_device, &
     ns, nt1, i2pi, fbdy, cbdy, irz, Jopt, icur_psi, xbrk_cur, &
     xbrk_mom, rtol_mom, gs_maxerr, ier)
 
  !-----------------------------------------------------------------------
  ! read G-EQDSK file and and perform inverse map,
  !   i.e. derive {R(theta,x),Z(theta,x)} flux surfaces from Psi(R,Z)
  !   by contouring (interpolation) techniques
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! The following data from the G-EQDSK data object "geq" are used:
  !    grid dimensions: RLEFT_M, RDIM_M, ZMID_M, ZDIM_M, NW, NH
  !    limiter contour: LIMITR, RLIM_M, ZLIM_M
  !    total toroidal current: CURRENT_A
  !        (always used for sign of current; see icur_psi option 
  !        description in next section of comments).
  !    poloidal flux at axis and byd:  SIMAG_Wb__Rad, SIBRY_Wb__Rad
  !        (used, but adjusted based on analysis of Psi(R,Z))
  !    poloidal flux function Psi(R,Z):  PSIRZ_Wb__Rad -- flux surfaces.
  !    central toroidal field: BCENTR_T (sign only)
  !    plasma boundary contour:  NBBBS, RBBBS_M, ZBBBS_M
  !        (will be adjusted based on analysis of Psi(R,Z))
  !    plasma axis location:  RMAXIS_M, ZMAXIS_M
  !        (will be adjusted based on analysis of Psi(R,Z))
  !    toroidal field function:  FPOL_TM
  !    pressure profile:  PRES_Nt__M2
  ! 
  !-----------------------------------------------------------------------
  ! CONDITIONALLY USED:  if icur_psi = 0, QPSI is used for the q
  !    profile and |CURRENT_A| is used as a boundary condition for the
  !    enclosed toroidal current profile; if icur_psi = 1, QPSI and
  !    CURRENT_A values are discarded: the q profile and total enclosed
  !    current are computed from contour integrals over Psi(R,Z) instead.
  !
  !    (Advice: for EFIT data from e.g. TRANSP time slices, not based on a 
  !    free boundary GS solution, icur_psi = 0 is better; for EFIT data
  !    that is based on a good free boundary GS solution, icur_psi = 1 
  !    is better).
  !
  !    For either setting of icur_psi, xbrk_cur usually needs to be set
  !    to a point somewhat in from the boundary, to get a smooth dI/drho
  !    and <J.B> in the edge region.  Typical value: xbrk_cur=0.9; 
  !    constraint: 0.75 <= xbrk_cur <= 0.95 as of June 2008 (dmc).
  ! 
  !-----------------------------------------------------------------------
  ! NOT USED:  d/dpsi profiles FFPRIM_T2M2Rad__Wb and PPRIME_NtRad__M2Wb
  !            QPSI and CURRENT_A magnitude are not used, if icur_psi = 1.
  ! NB: PPRIME data is saved in xplasma profile "dPdPsi"
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! ORIGINAL COMMENTS:
  !
  ! dmc 10 Jan 2001
  ! based on Alex's i2mex_fromGeqdsk ...
  ! for general xplasma use:
  ! (1) choose boundary psi(bdy) from combination of data and user
  !     input, e.g. use the .99*psi(sep) surface for the boundary
  ! (2) form and use a normalized sqrt toroidal flux variable "xi":
  !     spline q(psi) and integrate dphi = q*dpsi; xi=1 @ psi=psi(bdy).
  ! (3) delineate the psi contours at the psi's corresponding to the
  !     set of selected "xi" surfaces.
  !
  ! dmc -- bugfix -- 9 Nov 2001 -- following correct Lang Lao convention
  !     on significance of signs of items in the G-EQDSK data
  !     "isnccwx" retained for old xplasma g-eqdsk files that did
  !     not follow the sign convention.
  ! (corresponsing bugfix for xplasma g-eqdsk file write is in
  ! eq_geqdsk.for)
  !     for all "real" g-eqdsk files, isnccwx=1 is to be expected.
  !
  !-----------------------------------------------------------------------
  ! CHANGES:
  !
  ! dmc June 2009
  ! --save "PSI0" = the offset applied to make PsiRZ = 0 at the mag. axis.
  !   so that PsiRZ + Psi0 restores the original EFIT Psi(R,Z).  This is
  !   written as element PSI0 in the PSI0_CONTAINER list in xplasma -- see
  !   xplasma_save_psi0 in xplasma_fromgeq module.  If the input is detected
  !   to be a g-eqdsk dataset from a core-only TRANSP run, Psi0 is is not
  !   available and so is not saved.
  !
  !-----------------------
  ! dmc March 2008
  ! 1.  Form theta adaptive grid (thadapt in eqi_geq_mod) to concentrate
  ! points near minima in curvature radius of contours, i.e. near the
  ! separatrix points beyond the boundary.
  ! 2.  If the boundary, as computed on this grid, starting from an EFIT
  ! specified boundary location, exceeds the curvature criterion (cbdy),
  ! it is accepted as the boundary for {R(theta,x),Z(theta,x)}.
  ! 3.  If the EFIT boundary minimum curvature is less than (cbdy), then,
  ! find a boundary shifted in from the EFIT boundary to the extent 
  ! necessary to match the (cbdy) normalized curvature limit PRECISELY.
  ! 4.  Taking the chosen boundary points on the (1:nt1) theta grid, find
  ! all interior points using a root finder, instead of the LSODE-based
  ! contour following method (for better accuracy/reproducibility).

  ! major related sources:
  !   eqi_contour.f90 -- LSODE-based contour follower, used for analysis of
  !   the boundary
  !   eq_geq_findrz.f90 -- root finder code (match chosen target Psi values
  !   along selected line segments in (R,Z) space).
 
  ! ***NOTE***
  !   there are many "bad" EFIT Psi(R,Z) datasets that will trip any number
  !   of error checks in this routine.  If an error is detected, Psi(R,Z) 
  !   should be examined: does it have well behaved nested Psi level contours
  !   inside the specified EFIT boundary ???  Do not assume that the error is
  !   in the code...
  !--------------------------------------

  use xplasma_definitions
  use eqi_rzbox_module

  use eqi_geq_mod  ! contains psi(R,Z) spline, (R0,Z0), and (Rpp,Zpp)
                   ! also contains G-EQDSK module & instance "geq".
 
  implicit none
  integer, parameter :: r8=selected_real_kind(12,100)
 
  ! ------------------ arguments ---------------------------------
 
  character*(*), intent(in) :: filename  ! Geqdsk filename or MDS+ path
  !                  ...or...EFIT_INDEX:filename(t=<time>) specification
  !                  see comments in source/old_xplasma/eqm_fromgeqdsk.f90

  logical, intent(in) :: ihermite  ! .TRUE. for Hermite {R,Z}(rho,theta)
  !  .FALSE. for bicubic spline {R,Z}(rho,theta)
  !  NOTE: in both cases the spline is computed; this flag causes the
  !  bicubic spline representation to be converted to bicubic Hermite.
  !  EITHER WAY, the representation is C2; the two representations are
  !  generally equivalent to machine precision, but there are technical
  !  reasons for making this control available.

  logical, intent(in) :: inew_device  ! .TRUE. for new device; limiter 
  !  definitions to be read in (otherwise assume they are already known).

  integer, intent(in) :: ns     ! no. of surfaces incl. mag. axis
  integer, intent(in) :: nt1    ! no. of theta points
  integer, intent(in) :: i2pi   ! =1: theta range -pi:pi; =2: range 0:2pi
  real(r8), intent(in) :: fbdy  ! psi(bdy) = fbdy * psi(separatrix) [***]

  real(r8), intent(in) :: cbdy  ! min. allowed relative curvature
  !  actually, min(0.2_R8,max(0.03_R8,cbdy)) is used

  integer, intent(in) :: irz    ! =1 to extend xplasma mapped region to (R,Z)

  integer, intent(in) :: jopt   ! =1 for equal arc poloidal angle coordinate
                                ! =0 for traditional VMEC/descur compressed
                                !    moments
                                ! =-1 to use eqi_contour output directly...

  integer, intent(in) :: icur_psi  ! =0: use EFIT scalar data for total current
                                ! and use the EFIT 1d q(psi) profile "QPSI".
                                ! =1: discard EFIT CURRENT_A and QPSI data;
                                ! recompute q(psi) and total current from
                                ! Psi(R,Z) data instead.

  real(r8), intent(in) :: xbrk_cur ! break point for quadratic interpolation
                                ! of I(rho) from interior to edge value
                                ! (whether taken from CURRENT_A or Psi(R,Z)).
                                ! This is necessary to get a smooth <J.B> in
                                ! the edge region.  Typical value: 0.9;
                                ! enforced range: 0.75 <= xbrk_cur <= 0.95.

  real(r8), intent(in) :: xbrk_mom ! x width of axial regularization region
  real(r8), intent(in) :: rtol_mom ! rel. tolerance for zeroing of moments

  real(r8), intent(out) :: gs_maxerr  ! max rel. GS error found
  !  (normalized measure of self-consistency of equilibrium)

  integer, intent(out) :: ier   ! completion code, 0=OK
 
  ! [***] more precisely...
  !
  !   psi(bdy) = psi(axis) + fbdy*(psi(separatrix) - psi(axis))
  !
  !   as Geqdsk psi(axis) can be non-zero
  !
  ! ------------------- local memory -----------------------------
 
  real(r8), dimension(:), allocatable :: p, q, g, psi, psi_geq, psis, &
       & the, rgrid, zgrid, xi, phi, psivec, qpsia, pprime
  integer :: inw,inh

  real(r8), dimension(:, :), allocatable :: r, z, psirz_geq, gserr_geq
  real(r8), dimension(:, :), allocatable :: rmmc, rmms, zmmc, zmms
  real(r8), dimension(:, :), allocatable :: psirz_sm, gpsi_geq, gpsi_wk
  real(r8), dimension(:, :), allocatable :: x2d,x2dq,wkint

  real(r8), dimension(:), allocatable :: csthtk,snthtk,xint,thint,wint
  logical, dimension(:), allocatable :: lsurf,lqpsi
  real(r8), dimension(:), allocatable :: icurt,xtmp

  real(r8) :: rmin, rmax, zmin, zmax, surf, rm, zm, rp, zp, &
       pmag, pi, twopi, zfac, zfacp, tol, ctol, ztol, zero, bdist, &
       rzdel, zxbrk_mom, zfac0, zfac_ok, zfac_fail, zqgen, zgtmp

  real(r8) :: wts(10),xwts(10)  ! gauss-kronrod-patterson quadrature data
 
  real(r8) :: dvdpsi,gpsi2r2i,mu0,xibrk,x1,x2,x3,xxx,aicur,bicur
  integer :: id4(4),ideriv4(4),nsbrk

  real(r8) :: psimin, psimax
  real(r8) :: ztime1,zrnew,zznew,zone,zrtol
 
  real(r8), dimension(ns,2) :: rzvec
  real(r8), dimension(ns) :: rpa,zpa,zero_vec,one_vec,xvec
 
  real(r8), dimension(1) :: wkvec1,wkvec2,wkvec3
  real(r8), dimension(1,4) :: wkvec4
  real(r8), dimension(1,2) :: wkvec5

  ! bdy cond arrays
 
  real(r8) bcr(nt1),bcz(nt1)
 
  integer ifail,jfail,icount,iertmp,imaxr,inumx
  logical lmask(max(ns,nt1))
 
  integer iok, i, ith, j, iter, itmax, imom, ismoo, iez
  integer idrv1,iknot,iwarn,im,ibc1,iiopt,ikmax,ins_test
  integer isnccwi, isnccwb, isnccwx, id_lim, isize_lim
 
  integer id_rho,id_chi,id_R,id_Z,id_psi,id_gpsi,id_phi,id_g,id_q,id_p,idprev
  integer id_Rgrid,id_Zgrid,id_gserr,id_curt,id_pp
 
  real(r8) fbdy_use,cbdy_use,cbdy_got
 
  real(r8) :: xbrkii,zbdyp,zbdym,zincr,zbd4(4),rbd4(4),zwp,zwm,zerr,zdrv1
  integer :: ibdyp,ibdym,idir,iorig

  character*120 msgbuf

  character*32 zauthor
  logical :: ienable,ifinish

  integer, parameter :: ilun_dbg=6
  logical, parameter :: dbg_print = .FALSE.

  integer :: inz,instep,inzmin,inzprev,ifac,ict,ict_ok,ict_fail
  real*8, dimension(:), allocatable :: psis0
  real*8 :: rbdy00,zbdy00,psibdy0,psibdy,psi_test,psi_ok,psi_avg,psibdy_use
  real*8 :: psi_avg_ok,psi_inc,psibdy_reset
  real*8 :: psi_avg_save,zfac_fail_save,zfac_ok_save
  real*8 :: zth,zth2,zth2p,diffmax,psitmp
  real*8 :: psidum1,psidum2
  integer :: ii,jj,jj1,jj2,istat,nsing

  real*8, dimension(:), allocatable :: zwk

  real*8, parameter :: limtol = 1.0d-4

  !--------------------
  !  PSI0 list...

  character*10 :: ztest
  real*8 :: psi0
  logical :: ipsi0

  ! ------------------- root finder passed routine ---------------
 
  external eqi_geq_findrzs  ! find (R,Z) pts along set of arrays to match psi
  external eqi_find_curv    ! find boundary with desired curvature
  external eqi_geq_psimin   ! find min(psi(R)) along line of fixed Z
 
  external eqm_bpsi         ! given psi(R,Z), g(psi), compute B(R,Z) vector
 
  ! ------------------- executable code --------------------------
 
  if(dbg_print) then
     write(ilun_dbg,*) ' ----------------> eqi_fromgeqdsk:'
     write(ilun_dbg,*) '  filename = "'//trim(filename)//'"'
     write(ilun_dbg,*) '  ihermite = ',ihermite
     write(ilun_dbg,*) '  inew_device = ',inew_device
     write(ilun_dbg,*) '  ns  = ',ns
     write(ilun_dbg,*) '  nt1 = ',nt1
     write(ilun_dbg,*) '  i2pi = ',i2pi
     write(ilun_dbg,*) '  fbdy = ',fbdy
     write(ilun_dbg,*) '  cbdy = ',cbdy
     write(ilun_dbg,*) '  irz = ',irz
     write(ilun_dbg,*) '  icur_psi = ',icur_psi
     write(ilun_dbg,*) '  xbrk_cur = ',xbrk_cur
     write(ilun_dbg,*) '  Jopt = ',jopt,' *** '
     write(ilun_dbg,*) '  xbrk_mom = ',xbrk_mom
     write(ilun_dbg,*) '  rtol_mom = ',xbrk_mom
     write(ilun_dbg,*) ' --------------------------------------- '
  endif

  ier = 0
  msgbuf = ' '

  gs_errmax = -1  ! default value (seen e.g. on an early error exit)
 
  call xplasma_author_query(sp,zauthor,ienable,ier)
  if(ier.ne.0) return

  if(zauthor.eq.xplasma_root) then
     call xplasma_author_set(sp,'eqi_fromgeqdsk',ier)
  else
     call xplasma_author_resume(sp,zauthor,ier)  ! ensure write enable
  endif
  if(ier.ne.0) return

  pi   =3.1415926535897931_R8
  twopi=6.2831853071795862_R8
  mu0  = 2*twopi*1.0e-7_R8
        
  zero=0.0_R8
  zone=1.0_R8
 
  iez  =0

  itmax=10

  xbrkii = max(0.75_R8,min(0.95_R8,xbrk_cur))
  if((icur_psi.lt.0).or.(icur_psi.gt.1)) then
     write(msgbuf,*) ' ?eqi_fromgeqdsk: invalid icur_psi switch: not 0 or 1.'
     call eqi_fromg_err(msgbuf)
     ier=101
     go to 1000
  endif

  ! ---> read Geqdsk data
 
  call geq_init(geq, filename, iok)
  call geq_error(iok)
  if(iok/=0) ier = 117
  if(iok/=0) go to 1000
 
  ! grid limits:

  rmin = geq%RLEFT_M
  rmax = geq%RLEFT_M + geq%RDIM_M
  zmin = geq%ZMID_M - geq%ZDIM_M/2.0_r8
  zmax = geq%ZMID_M + geq%ZDIM_M/2.0_r8

  if((rmax.le.rmin).or.(zmax.le.zmin)) then
     write(msgbuf,*) '  singular EFIT R or Z grid: Rmin,Rmax: ',rmin,rmax
     call eqi_fromg_err(msgbuf)
     write(msgbuf,*) '  singular EFIT R or Z grid: Zmin,Zmax: ',zmin,zmax
     call eqi_fromg_err(msgbuf)
     ier = 105
     go to 1000
  endif

  ! check limiter range:

  do ii=1,geq%Limitr
     if((geq%Rlim_m(ii).lt.(Rmin-limtol)).or. &
          (geq%Rlim_m(ii).gt.(Rmax+limtol))) then
        if(ier.eq.0) then
           write(msgbuf,*) &
                '  Limiter R point outside R grid range Rmin to Rmax: ', &
                Rmin,Rmax
           call eqi_fromg_err(msgbuf)
           write(msgbuf,*) &
                '  Limiter point #',ii,' R value: ',geq%Rlim_m(ii)
           call eqi_fromg_err(msgbuf)
        endif
        ier = ier + 1
     else
        geq%Rlim_m(ii) = max(Rmin, min(Rmax, geq%Rlim_m(ii)) )
     endif

     if((geq%Zlim_m(ii).lt.(Zmin-limtol)).or. &
          (geq%Zlim_m(ii).gt.(Zmax+limtol))) then
        if(ier.eq.0) then
           write(msgbuf,*) &
                '  Limiter Z point outside Z grid range Zmin to Zmax: ', &
                Zmin,Zmax
           call eqi_fromg_err(msgbuf)
           write(msgbuf,*) &
                '  Limiter point #',ii,' Z value: ',geq%Zlim_m(ii)
           call eqi_fromg_err(msgbuf)
        endif
        ier = ier + 1
     else
        geq%Zlim_m(ii) = max(Zmin, min(Zmax, geq%Zlim_m(ii)) )
     endif
  enddo

  if(ier.ne.0) then
     ier=104
     go to 1000
  endif

  ! ---> check signs
 
  if( geq%CURRENT_A .gt. 0.0_R8 ) then
 
     isnccwi = 1               ! toroidal current ccw looking down on tokamak
                               ! (g-eqdsk test)
 
  else
 
     isnccwi = -1              ! toroidal current clockwise
                               ! (g-eqdsk test)
 
  endif
 
  if( (geq%SIBRY_Wb__Rad - geq%SIMAG_Wb__Rad ) .gt. 0.0_R8 ) then
 
     isnccwx = 1               ! toroidal current ccw looking down on tokamak
                               ! (old xplasma g-eqdsk)
 
  else
 
     isnccwx = -1              ! toroidal current clockwise
                               ! (old xplasma g-eqdsk)

     call eqi_fromg_err(' %eqi_fromgeqdsk warning (fixup applied):')
     msgbuf=' '
     write(msgbuf,*) '  Psi(bdy) = ',geq%sibry_WB__Rad, &
         ' < Psi(0) = ',geq%SIMAG_Wb__Rad
     call eqi_fromg_err(msgbuf)
     call eqi_fromg_err('  (violates de facto EFIT G-eqdsk convention)')
 
  endif
 
  ! ---> and normalize psi definition:  0 on axis & increasing outward
 
  psimin=0.0_R8
  psimax=isnccwx*(geq%SIBRY_Wb__Rad-geq%SIMAG_Wb__Rad)
 
  inw = geq%NW   ! number of R pts = # of Psi pts in 1d profiles in EFIT data
  inh = geq%NH   ! number of Z pts

  allocate(psirz_geq(inw,inh),gserr_geq(inw,inh))
  allocate(psirz_sm(inw,inh))
  psirz_geq = isnccwx*(geq%psirz_Wb__Rad - geq%SIMAG_Wb__Rad)
 
  ! ---> form R & Z rectangular mesh & check data

  allocate(rgrid(inw))
  allocate(zgrid(inh))
 
  rgrid=rmin + (rmax - rmin)*(/ (real(i-1,r8)/real(inw-1,r8), i=1,inw) /)
  zgrid=zmin + (zmax - zmin)*(/ (real(i-1,r8)/real(inh-1,r8), i=1,inh) /)
 
  call psirz_minmax(inw,inh,rgrid,zgrid,psirz_geq, &
       geq%NBBBS,geq%RBBBS_M,geq%ZBBBS_M,psidum1,psidum2,nsing)

  if(nsing.gt.1) then
     call eqi_fromg_err( &
          ' ?eqi_fromgeqdsk: Psi(R,Z) screening failed, multiple local minima or maxima')
     ier=116
     go to 1000
  endif

  ztest(1:8) = geq%case_(1)
  ztest(9:10) = geq%case_(2)(1:2)

  ipsi0=.TRUE.   ! TRUE to save PSI0 offset
  if(ztest(1:6).eq.'TRXPL ') then

     ! a zero offset value will be saved if this appears to be a
     ! g-eqdsk data set extrapolated from a core only TRANSP equilibrium...

     if(ztest(7:10).ne.'(fb)') then
        ipsi0=.FALSE.
     endif
  endif

  if(ipsi0) then
     psi0 = geq%SIMAG_Wb__Rad
  else
     psi0 = 0.0_r8
  endif

  if (geq%BCENTR_T .gt. 0.0_R8 ) then
 
     isnccwb = 1               ! toroidal field ccw looking down on tokamak
 
  else
 
     isnccwb = -1              ! toroidal field clockwise
 
  endif

  ! ---> smooth Psi(R,Z) beyond the boundary only; use EFIT bdy for this...

  call psirz_smooth(inw,inh,rgrid,zgrid,psirz_geq,psirz_sm, &
       geq%NBBBS,geq%RBBBS_M,geq%ZBBBS_M)
 
  ! check GS error
 
  gserr_geq = 0
 
  bdist = 0.05_r8   ! screen out pts within bdist*(Rmax-Rmin) of bdy
 
  call geq_GSerror(isnccwx, &
       geq, inw, inh, rgrid, zgrid, bdist, gserr_geq, ier)
  if(ier.ne.0) then
     call eqi_fromg_err( &
          ' ?eqi_fromgeqdsk: geq_GSerror found no equilibrium data.')
     gserr_geq = 1
     ier=121
     go to 1000
  endif
 
  gs_errmax=maxval(abs(gserr_geq))  ! in eqi_geq_mod module
  gs_maxerr=gs_errmax               ! in return argument
 
  allocate(psi(ns), psi_geq(inw), p(ns), g(ns), pprime(ns), q(ns), the(nt1))
  allocate(phi(ns), xi(ns))
  allocate(r(nt1, ns), z(nt1, ns), gpsi_geq(nt1,ns))
  allocate(psis(nt1))
 
  ! ---> form psi(R,Z) spline
 
  iez = 1
  call ezspline_init(pspl, inw, inh, (/0,0/), (/0,0/), iok)
  call ezspline_error(iok)
  if(iok.eq.0) then
     pspl%x1 = rgrid
     pspl%x2 = zgrid
     !pspl%isHermite = 1 ! Akima Hermite
     call ezspline_setup(pspl, psirz_sm, iok)
     call ezspline_error(iok)
  endif

  deallocate(psirz_geq)  ! use psirz_sm now...

  if(iok.ne.0) then
     call eqi_fromg_err( &
          ' ?eqi_fromgeqdsk: Psi(R,Z) bicubic spline setup error.')
     ier = 106
     go to 1000
  endif

  ! ---> this will be the normalized sqrt(tor.flux) mesh
 
  xi = (/ (real(i-1,r8)/real(ns-1,r8), i=1,ns) /)
 
  ! ---> and the poloidal grid...

  if(i2pi.eq.1) then
     the = -pi + twopi * (/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /) ! -pi:pi
  else
     the = twopi * (/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /)  ! range 0:2pi
  endif
 
  ! ---> mag axis & point on G-EQDSK plasma boundary -- from EFIT
 
  rm = geq%RMAXIS_M
  zm = geq%ZMAXIS_M

  ! find point with max R value; also construct bdy sequence which does
  ! not close on the last point and which has no duplicate points

  inumx=geq%NBBBS
  if(allocated(rzwk)) deallocate(rzwk)
  allocate(rzwk(inumx,2))

  imaxr=1
  inumx=1
  rzwk(1,1)=geq%RBBBS_M(1)
  rzwk(1,2)=geq%ZBBBS_M(1)
  do ii=1,geq%NBBBS
     if((geq%RBBBS_M(ii).eq.rzwk(inumx,1)).and. &
          (geq%ZBBBS_M(ii).eq.rzwk(inumx,2))) cycle
     if((geq%RBBBS_M(ii).eq.rzwk(1,1)).and. &
          (geq%ZBBBS_M(ii).eq.rzwk(1,2))) cycle
     inumx=inumx+1
     rzwk(inumx,1)=geq%RBBBS_M(ii)
     rzwk(inumx,2)=geq%ZBBBS_M(ii)
     if(rzwk(inumx,1).gt.rzwk(imaxr,1)) then
        imaxr=inumx
     endif
  enddo

  rp = rzwk(imaxr,1)
  zp = rzwk(imaxr,2)
 
  rzdel = (rp-rm)/10

  ! --> try to adjust mag axis location (find true minimum of psi(R,Z) 
  !     interpolant)

  lmask(1)=.FALSE.
  tol=1.0e-9
  wkvec1=zm-rzdel
  wkvec2=zm+rzdel
  wkvec4(1,1)=rm
  wkvec4(1,2)=rzdel
  wkvec4(1,3)=tol
  wkvec4(1,4)=tol*(psimax-psimin)/rzdel
  
  call zriddery(1, lmask(1:1), wkvec1, wkvec2, tol, &
       wkvec4(1,4), eqi_geq_psimin, wkvec3, ifail, &
       1, wkvec4, 4, wkvec5, 2)
  jfail=abs(wkvec5(1,2))+0.5_R8

  if((ifail.ne.0).or.(jfail.ne.0)) then
     call eqi_fromg_err( &
          ' %eqi_fromgeqdsk: Psi_min search failed, using EFIT result')
     r0 = rm
     z0 = zm
  else
     z0=wkvec3(1)
     r0=wkvec5(1,1)
     if(max(abs(r0-rm),abs(z0-zm)).gt.0.005_R8*r0) then
        call eqi_fromg_err(' %eqi_fromgeqdsk: Psi_min axis adjustment:')
        msgbuf=' '
        write(msgbuf,'(5x,a,1pe12.5," -> ",1pe12.5)') 'R_axis: ',rm,r0
        call eqi_fromg_err(msgbuf)
        msgbuf=' '
        write(msgbuf,'(5x,a,1pe12.5," -> ",1pe12.5)') 'Z_axis: ',zm,z0
        call eqi_fromg_err(msgbuf)
     endif
     rm=r0
     zm=z0
  endif
          
  ! now find a new (EFIT) boundary point (rp,zp) with zp=z0 precisely
  ! initiate search

  if(zp.ne.z0) then
     if(zp.gt.z0) then
        zbdyp=zp                  ! zp > z0
        ibdyp=imaxr
        iorig=1

        zbdym=zbdyp
        ibdym=ibdyp
        idir=-1
        call fnxz(ibdym,zbdym)
        if(zbdym.gt.zp) then
           ibdym=ibdyp
           idir=1
           call fnxz(ibdym,zbdym)
        endif

     else
        zbdym=zp                 ! zp < z0
        ibdym=imaxr
        iorig=-1

        zbdyp=zbdym
        ibdyp=ibdym
        idir=1
        call fnxz(ibdyp,zbdyp)
        if(zbdyp.lt.zbdym) then
           ibdyp=ibdym
           idir=-1
           call fnxz(ibdyp,zbdyp)
        endif
     endif

     if(ier.ne.0) then
        call eqi_fromg_err('%eqi_fromgeqdsk: EFIT boundary singular.')
        ier=118
        go to 1000
     endif

     ict=0
     do while ((zbdyp-z0)*(zbdym-z0).gt.0.0_R8)
        ict = ict + 1
        if(ict.gt.inumx) then
           ier=119
           call eqi_fromg_err('%eqi_fromgeqdsk: Psimin outside EFIT boundary.')
           go to 1000
        endif

        call fnxz(ibdyp,zbdyp)
        call fnxz(ibdym,zbdym)
     enddo

     ! pts straddling z0...
     zbd4(2)=zbdym; rbd4(2)=rzwk(ibdym,1)
     zbd4(3)=zbdyp; rbd4(3)=rzwk(ibdyp,1)

     ! and their immediate neighbors...
     idir=idir*iorig
     call fnxz(ibdym,zbdym)
     zbd4(1)=zbdym; rbd4(1)=rzwk(ibdym,1)
     idir=-idir
     call fnxz(ibdyp,zbdyp)
     zbd4(4)=zbdyp; rbd4(4)=rzwk(ibdyp,1)

     ! Now: zbd4(1:4) is an ascending sequence with z0 btw zbd4(2) and zbd4(3).
     !    & rbd4(1:4) are the corresponding R values.

     zwp = (z0-zbd4(2))/(zbd4(3)-zbd4(2))
     zwm = 1.0_R8 - zwp

     rp = zwp*quadintrp(z0,zbd4(2:4),rbd4(2:4)) + &
          zwm*quadintrp(z0,zbd4(1:3),rbd4(1:3))
     zp = z0
  endif

  ! select a point slightly beyond the plasma boundary...
 
  if(rp.le.r0) then
     call eqi_fromg_err( &
          ' %eqi_fromgeqdsk:  outer boundary point Rp < axis R0 ')
     msgbuf=' '
     write(msgbuf,*) '  Rp = ',rp,' R0 = ',r0
     call eqi_fromg_err(msgbuf)
     ier = 120
     go to 1000
  endif
           
  zfac = min(1.025_R8,(rmax-r0)/(rp-r0))
  if(zfac.lt.1.0_R8) then
     call eqi_fromg_err( &
          ' %eqi_fromgeqdsk:  boundary point beyond mapped region.')
     msgbuf=' '
     write(msgbuf,*) '  rp = ',rp,' rmax = ',rmax
     call eqi_fromg_err(msgbuf)
     ier = 120
     go to 1000
  endif
 
  ! interpolate to get psi on axis
 
  call ezspline_interp(pspl, rm, zm, pmag, iok)
  call ezspline_error(iok)

  if(iok/=0) ier = 120
  if(iok/=0) go to 1000
 
  ! psi on separatrix (cleaned up EFIT bdy)
 
  allocate(psis0(inumx))

  call ezspline_interp(pspl, inumx, rzwk(1:inumx,1), rzwk(1:inumx,2), psis0, &
       iok)
  call ezspline_error(iok)

  if(iok/=0) ier = 120
  if(iok/=0) go to 1000

  psep=sum(psis0)/inumx

  deallocate(psis0)
 
  ! check against file values for psi on axis & separatrix
 
  fbdy_use = (psep-pmag)/(psimax-psimin)
  if((fbdy_use.lt.0.98_R8).or.(fbdy_use.gt.1.02_R8)) then
     call eqi_fromg_err(' %eqi_fromgeqdsk:  psi(R,Z) interpolation')
     call eqi_fromg_err('  yields psi(mag) and psi(sep) different')
     call eqi_fromg_err('  from scalar values in g-eqdsk file:')
     msgbuf=' '
     write(msgbuf,*) '  -->interpolation values:  ',pmag,psep
     call eqi_fromg_err(msgbuf)
     msgbuf=' '
     write(msgbuf,*) '  -->offset file values:    ',psimin,psimax
     call eqi_fromg_err(msgbuf)
     ier = 120
     go to 1000
  endif

  !  reset Psi(R,Z) so that actual minimum is zero

  psirz_sm = psirz_sm-pmag
  pspl%fspl(1,:,:) = pspl%fspl(1,:,:)-pmag
  psep=psep-pmag
  if(ipsi0) psi0 = psi0 + pmag
  pmag=0.0_R8
  psimax=psep
  psimin=pmag
 
  fbdy_use = min(1.0_R8, max(0.95_R8, fbdy))
 
  cbdy_use = max(0.03_R8, min(0.2_R8, cbdy))

  cbdy_chk = min(0.25_R8, 2.0_R8*cbdy_use)  ! criterion for "small curvature"
  ! used to identify regions needing finer th grid resolution near bdy
 
  psibdy=psep*fbdy_use  ! EFIT bdy, adjusted 
  psi_ok=0.9_r8*psep    ! should be a safe value well inside the bdy

  !  find a closed field line outside psibdy or as close to it as possible,
  !  relaxing the curvature requirement
  !     ...accept if found outside the EFIT bdy, or if its curvature ratio
  !     is less than half the constrained user specified value (cdby_adjust).

  if(allocated(thadapt)) deallocate(thadapt,radapt,zadapt,curvad,iadapt)
  inz=50*nt1  ! a safe upper limit on required size
  allocate(thadapt(inz,2),radapt(inz,2),zadapt(inz,2),curvad(inz,2))
  allocate(iadapt(inz))

  inzmin=15   ! min. no. of zones to cover each region of high curvature

  !  initial grid...

  isw=1
  nthad=2*(nt1-1) + 1  ! double the grid
  nthad0 = nthad       ! (store in eqi_geq_mod)

  do ii=1,nthad
     thadapt(ii,isw)=(ii-1)*twopi/(nthad-1)
  enddo
  
  zfac0 = 1.0011_R8
  zincr = 5.0e-4_R8

  zfac_ok = zero    ! have not yet seen a closed contour
  zfac_fail = zero  ! have not yet seen an open contour

  iter = 0
  ict_ok = 0        ! #of closed contours seen
  ict_fail = 0      ! #of non-closed (or otherwise unacceptable) contours seen

  psi_avg_ok = 1.000010_R8*psibdy
  psi_inc = 1.0e-6_R8*psibdy
  psibdy_reset = psibdy
  psi_avg_save = -psibdy

  do
     zfac = zfac0
     do
        iter=iter+1
        if(iter.gt.10*itmax) then
           call eqi_fromg_err( &
                ' ?eqi_fromgeqdsk:  constrained closed boundary search failure.')
           ier=129
           go to 1000
        endif

        zfac = zfac - zincr
        rbdy00 = r0 + zfac*(rp-r0)
        zbdy00 = z0 + zfac*(zp-z0)

        call eqi_contour(rbdy00, zbdy00, 0.1_r8*cbdy_use, &
             psep, cbdy_got, psi_avg, istat, iok)
 
        if(istat.ge.3) then
           call eqi_fromg_err( &
                ' ?eqi_fromgeqdsk: unexpected error in eqi_contour.')
           ier=129
           go to 1000

        else if(istat.eq.0) then
           ict_ok = ict_ok + 1
           zfac_ok=zfac
           exit
        else
           ! contour not closed or too kinked at psi=psi_avg
           zfac_fail = zfac
           ict_fail = ict_fail + 1  ! # of non-closed contours seen

           !  reduce psi_avg_ok to be inside open/kinked contours if necessary
           do while (psi_avg.lt.psi_avg_ok)
              psi_avg_ok = psi_avg_ok-psi_inc
           enddo

           psibdy_reset = min(psibdy_reset, psi_avg_ok-psi_inc)

           if(ict_ok.gt.0) then
              zincr=zincr/2         ! we are in a binary search now
           endif
        endif
        if(zfac.lt.0.75_R8) then
           ier=122
           call eqi_fromg_err( &
                ' ?eqi_fromgeqdsk:  closed boundary contour search failed.')
           go to 1000
        endif
     enddo

     if((cbdy_got.ge.0.5_r8*cbdy_use).and.(psi_avg.le.psi_avg_ok)) then
        if(ict_fail.eq.0) then
           ! no open contour seen; starting point of search was too far in
           zfac0 = zfac0 + zincr
        else
           if(psi_avg_save.gt.psi_avg_ok) then
              ! psi_avg_ok moved; prior surface OK now
              zfac_fail=zfac_fail_save
              zfac_ok=zfac_ok_save
           else
              psi_avg_save=psi_avg
              zfac_fail_save=zfac_fail
              zfac_ok_save=zfac_ok
           endif
           zincr = (zfac_fail-zfac_ok)/2  ! binary search btw known limits
           zfac0 = zfac + 2*zincr         ! 1*zincr gets subtracted off
        endif
        cycle
     else
        exit  ! have a usable reference surface now.
     endif
  enddo

  psibdy=psibdy_reset

  if(psi_avg.lt.psi_ok) then
     ! warning only
     call eqi_fromg_err(' %eqi_fromgeqdsk: found closed contour far from EFIT nominal boundary.')
  endif

  ! refine th grid, so that curvature can be measured on ray-psi intercepts
  ! without having to solve an ODE to follow each level-psi contour.  This
  ! is done because the regions of sharp curvature are highly localized in
  ! the geometrical angle space about (r0,z0).

  do 
     inzprev=nthad
     call nx_thadapt(nt1,inzmin)
     if(nthad.eq.inzprev) then
        isw=3-isw  ! no change from prior grid, so, reuse prior data rather
        ! than recompute
        exit
     endif

     !  tolerate detection of sharpened curvature; effect due to change of
     !  "nthad" is possible, so, pass 0.05_r8*cbdy_use instead of 0.1_r8*...

     call eqi_contour(rbdy00, zbdy00, 0.05_r8*cbdy_use, &
          psep, cbdy_got, psi_avg, istat, iok)

     ! iok must = 0 now
     if(iok.ne.0) then
        call eqi_fromg_err(' ?eqi_fromgeqdsk: algorithm failure (L581)')
        ier=138
        go to 1000
     endif
  enddo

  ! force bdy search points to be exactly on rays emanating from the
  ! origin (r0,z0) at the specified grid of angles atan2((Z-Z0),(R-R0)) 
  ! -- they are nearly so already

  if(allocated(rseek)) deallocate(rseek,zseek)
  allocate(rseek(-1:nthad+2),zseek(-1:nthad+2))
  allocate(zwk(-1:nthad+2))

  zth2=0.0_r8
  do ii=1,nthad
     zth2p=zth2
     zth2=atan2(zadapt(ii,isw)-z0,radapt(ii,isw)-r0)
     if(zth2.lt.zth2p-pi) zth2=zth2+twopi
     zwk(ii)=zth2
  enddo

  iter=0
  diffmax=1.0_r8
  do while(diffmax.gt.(twopi*1.0e-15_r8))
     iter=iter+1
     if(iter.gt.itmax) then
        call eqi_fromg_err( &
             ' ?eqi_fromgeqdsk:  search ray angle refinement failure.')
        ier=131
        go to 1000
     endif

     call load_seek

     diffmax=0.0_r8
     zth2=0.0_r8
     do ii=1,nthad
        zth=thadapt(ii,isw)
        if(zth.eq.zwk(ii)) then
           zth2=zth
           cycle
        else if(zth.lt.zwk(ii)) then
           zwp=(zth-zwk(ii-1))/(zwk(ii)-zwk(ii-1))
           zwm=1.0_r8-zwp
           radapt(ii,isw) = &
                zwp*quadintrp(zth,zwk(ii-1:ii+1),rseek(ii-1:ii+1)) + &
                zwm*quadintrp(zth,zwk(ii-2:ii),rseek(ii-2:ii))
           zadapt(ii,isw) = &
                zwp*quadintrp(zth,zwk(ii-1:ii+1),zseek(ii-1:ii+1)) + &
                zwm*quadintrp(zth,zwk(ii-2:ii),zseek(ii-2:ii))
        else
           zwm=(zwk(ii+1)-zth)/(zwk(ii+1)-zwk(ii))
           zwp=1.0_r8-zwm
           radapt(ii,isw) = &
                zwm*quadintrp(zth,zwk(ii-1:ii+1),rseek(ii-1:ii+1)) + &
                zwp*quadintrp(zth,zwk(ii-2:ii),rseek(ii-2:ii))
           zadapt(ii,isw) = &
                zwm*quadintrp(zth,zwk(ii-1:ii+1),zseek(ii-1:ii+1)) + &
                zwp*quadintrp(zth,zwk(ii:ii+2),zseek(ii:ii+2))
        endif
        ! verify
        zth2p=zth2
        zth2=atan2(zadapt(ii,isw)-z0,radapt(ii,isw)-r0)
        if(zth2.lt.zth2p-pi) zth2=zth2+twopi
        diffmax=max(diffmax,abs(zth2-zth))
        zwk(ii)=zth2
     enddo
  enddo

  deallocate(zwk)
  allocate(psis0(nthad))

  call ezspline_interp(pspl, nthad, radapt(:,isw), zadapt(:,isw), psis0, iok)
  call ezspline_error(iok)

  psibdy0=1.01_r8*minval(psis0)-0.01_r8*maxval(psis0)
               ! this Psi(R,Z) value corresponds safely to a closed surface
               ! within the boundary just found (psis0(:) values are all
               ! very nearly equal).
  deallocate(psis0)

  psi_test = min(psibdy0,psibdy)

  !  OK, eqi_geq_mod now contains contour safely beyond final boundary
  !     "zridderx" searches btw (r0,z0) and this can be safely carried out

  !---------------------
  !  find the boundary to be used...

  if(allocated(rzwk)) deallocate(rzwk)
  if(allocated(curvec)) deallocate(curvec)
  allocate(rzwk(nthad,2),curvec(nthad))

  call eqi_find_curvature(nthad,psi_test,0.8_r8,psep,rzwk,curvec,cbdy_got)

  ctol = 1.0e-12_r8
  ztol = ctol*cbdy_use

  if((psi_test.eq.psibdy).and.(cbdy_got.ge.cbdy_use)) then
     ! OK the EFIT bdy (adjusted by fbdy_use) curvature is OK -- use it

     psibdy_use = psi_test

  else if(abs(cbdy_got-cbdy_use).lt.ztol) then
     ! curvature is close enough

     psibdy_use = psi_test

  else if(cbdy_got.gt.cbdy_use) then

     call eqi_fromg_err( &
          ' ?eqi_fromgeqdsk: curvature-match boundary search algorithm failure!')
     ier = 120
     go to 1000

  else
     ! curvature is too sharp at psi=psi_test; find a psi value further
     ! in that gets the desired curvature

     lmask(1)=.FALSE.
     tol=ctol*psi_test

     call zriddery(1, lmask(1:1), psi_ok, psi_test, tol, ztol, &
          eqi_find_curv, psibdy_use, ifail, &
          1, cbdy_use, 1, zerr, 1)

     if((ifail.ne.0).or.(zerr.ne.0.0_r8)) then
        call eqi_fromg_err( &
             ' ?eqi_fromgeqdsk: boundary-curvature-match root finder failure!')
        ier = 120
        go to 1000
     endif

  endif

  !  pick off subset of points that are alligned with the the(1:nt1) grid
  deallocate(rseek,zseek)
  allocate(rseek(nt1),zseek(nt1))

  ii=1
  if(i2pi.eq.1) then
     ! contour covers [0:2pi]; the(1:nt1) convers [-pi:pi]
     jj=(nt1/2) + 1
  else
     ! both contour & the(1:nt1) cover [0:2pi]
     jj=1
  endif

  rseek(jj)=radapt(ii,isw)  ! safe surface beyond bdy
  zseek(jj)=zadapt(ii,isw)

  r(jj,ns)=rzwk(ii,1)       ! pt from psibdy_use bdy contour just found
  z(jj,ns)=rzwk(ii,2)
  jj=jj+1
  ifac=0

  do ii=2,nthad
     if(abs(the(jj)+ifac*twopi-thadapt(ii,isw)).lt.1.0e-14*twopi) then
        rseek(jj)=radapt(ii,isw)
        zseek(jj)=zadapt(ii,isw)
        r(jj,ns)=rzwk(ii,1)       ! pt from psibdy_use bdy contour just found
        z(jj,ns)=rzwk(ii,2)
        jj=jj+1
        if(jj.gt.nt1) then
           rseek(1)=rseek(nt1)
           zseek(1)=zseek(nt1)
           r(1,ns)=r(nt1,ns)
           z(1,ns)=z(nt1,ns)
           jj=2
           ifac=1
        endif
     endif
  enddo

  nthad = nt1
  isw = 1
  thadapt(1:nt1,isw) = the
  radapt(1:nt1,isw) = rseek
  zadapt(1:nt1,isw) = zseek

  deallocate(rseek,zseek)

  ! OK now ready to find rest of surfaces

  ! mag axis:
  r(1:nt1,1) = r0
  z(1:nt1,1) = z0

  if(allocated(rzwk)) deallocate(rzwk)
  if(allocated(curvec)) deallocate(curvec)
  allocate(rzwk(nthad,2))

  allocate(lsurf(ns)); lsurf=.FALSE.
  lsurf(1)=.TRUE.
  lsurf(ns)=.TRUE.

  allocate(x2d(nt1,ns),psivec(nt1))

  x2d(:,1)=0.0_R8
  x2d(:,ns)=1.0_R8

  !----------------------------------
  ! (equispaced) G-eqdsk Psi grid for q,P,g...

  psi_geq = psimin + (psimax-psimin)* &
       & (/ (real(i-1, r8)/real(inw-1,r8), i=1, inw) /)
  
  !----------------------------------
  ! reproduce q(psi) profile, by finding level contours on the equispaced
  ! Psi grid and performing the necessary integration

  if(icur_psi.eq.0) then
     ! use QPSI data from EFIT

     allocate(qpsia(inw))
     qpsia(1:inw) = geq%QPSI(1:inw)

  else
     ! reconstruct q profile from Psi(R,Z) contours

     allocate(qpsia(inw),lqpsi(inw),x2dq(nt1,inw),xtmp(nt1))

     lqpsi = .FALSE.
     x2dq(:,1) = 0.0_R8; lqpsi(1)=.TRUE.
     x2dq(:,inw) = 1.0_R8; lqpsi(inw)=.TRUE.

     ! bdy Qpsi recomputed here:
     call eqi_qpsi_regen(geq%FPOL_TM(inw), nt1, &
          r(1:nt1,ns), z(1:nt1,ns), zqgen)
     qpsia(inw) = zqgen

     ! axial q will be filled in later by an extrapolation

     inz=inw-1
     instep=inz
     do
        instep=instep/2
        if(instep.eq.0) exit

        do ii=1,inw,instep
           if(lqpsi(ii)) cycle

           ! find nearest known surfaces
        
           do jj=ii-1,1,-1
              if(lqpsi(jj)) then
                 jj1=jj
                 exit
              endif
           enddo

           do jj=ii+1,inw,1
              if(lqpsi(jj)) then
                 jj2=jj
                 exit
              endif
           enddo

           lmask(1:nt1-1)=.FALSE.
           lmask(nt1)=.TRUE.

           psivec = psi_geq(ii)

           tol = 1.0e-12_r8
           call zridderx(nt1,lmask, x2dq(:,jj1), x2dq(:,jj2), tol, tol*psep, &
                eqi_geq_findrzs, x2dq(:,ii), ifail, nt1, psivec, 1, rzwk, 2)

           if(ifail.ne.0) then
              call eqi_fromg_err(' %eqi_fromgeqdsk:  Psi root finder (zridderx) error.')
              msgbuf=' '
              write(msgbuf,*) ' ...ifail=',ifail,' (R,Z) pts not found.'
              call eqi_fromg_err(msgbuf)
              ier=120
              go to 1000

           else
              lqpsi(ii)=.TRUE.  ! succeeded with this surface
           endif

           rzwk(nt1,1)=rzwk(1,1)
           rzwk(nt1,2)=rzwk(1,2)

           call eqi_qpsi_regen(geq%FPOL_TM(ii), nt1, &
                rzwk(1:nt1,1), rzwk(1:nt1,2), zqgen)
           qpsia(ii) = zqgen
        enddo
     enddo

     ! axial extrapolation

     lmask(1:nt1-1)=.FALSE.
     lmask(nt1)=.TRUE.

     psivec = psi_geq(1) + 0.2_R8*(psi_geq(2)-psi_geq(1))

     tol = 1.0e-12_r8
     call zridderx(nt1,lmask, x2dq(:,1), x2dq(:,2), tol, tol*psep, &
          eqi_geq_findrzs, xtmp, ifail, nt1, psivec, 1, rzwk, 2)

     if(ifail.ne.0) then
        write(msgbuf,*) ' %zridderx psi contour error near axis; simplified q extrapolation there.'
        qpsia(1) = 2.0_R8*qpsia(2)-qpsia(3)
     else

        rzwk(nt1,1)=rzwk(1,1)
        rzwk(nt1,2)=rzwk(1,2)

        zgtmp = geq%FPOL_TM(1) + 0.2_R8*(geq%FPOL_TM(2)-geq%FPOL_TM(1))
        call eqi_qpsi_regen(zgtmp, nt1, rzwk(1:nt1,1), rzwk(1:nt1,2), zqgen)

        qpsia(1) = (5.0_R8*zqgen - qpsia(2))/4.0_R8
     endif
  endif

  !----------------------------------
  ! using the recomputed q(Psi),
  ! select the psi-surfaces to use
 
  call eqi_geq_psfind(psi_geq, qpsia, inw, xi, ns, psibdy_use/psimax, &
       psi, phi, q, ier)

  deallocate(qpsia)
  if(allocated(lqpsi)) deallocate(lqpsi,x2dq,xtmp)

  if(ier.ne.0) go to 1000
 
  ! get the 1d profiles interpolated to the selected surfaces
 
  call eqi_geq_intrp(psi_geq, geq%FPOL_TM, inw, psi, ns, g, ier)
  g=abs(g)         ! sign is saved in isnccwb
 
  call eqi_geq_intrp(psi_geq, geq%PRES_Nt__M2, inw, psi, ns, p, ier)
  if(ier.ne.0) go to 1000
 
  call eqi_geq_intrp(psi_geq, geq%PPRIME_NtRad__M2Wb, inw, psi, ns, pprime, &
       ier)
  if(ier.ne.0) go to 1000

  !----------------------------------
  ! axis and boundary surface are known; fill in interior surfaces

  inz=ns-1
  instep=inz
  do
     instep=instep/2
     if(instep.eq.0) exit

     do ii=1,ns,instep
        if(lsurf(ii)) cycle

        ! find nearest known surfaces
        
        do jj=ii-1,1,-1
           if(lsurf(jj)) then
              jj1=jj
              exit
           endif
        enddo

        do jj=ii+1,ns,1
           if(lsurf(jj)) then
              jj2=jj
              exit
           endif
        enddo

        lmask(1:nt1-1)=.FALSE.
        lmask(nt1)=.TRUE.

        psivec = psi(ii)

        tol = 1.0e-12_r8
        call zridderx(nt1,lmask, x2d(:,jj1), x2d(:,jj2), tol, tol*psep, &
             eqi_geq_findrzs, x2d(:,ii), ifail, nt1, psivec, 1, rzwk, 2)

        if(ifail.ne.0) then
           call eqi_fromg_err(' %eqi_fromgeqdsk:  root finder (zridderx) error.')
           msgbuf=' '
           write(msgbuf,*) ' ...ifail=',ifail,' (R,Z) pts not found.'
           call eqi_fromg_err(msgbuf)
           ier=120
           go to 1000

        else
           lsurf(ii)=.TRUE.  ! succeeded with this surface
        endif

        rzwk(nt1,1)=rzwk(1,1)
        rzwk(nt1,2)=rzwk(1,2)

        do jj=1,nt1
           r(jj,ii)=rzwk(jj,1)
           z(jj,ii)=rzwk(jj,2)
        enddo
  
     enddo
  enddo
 
  ! define rho grid for xplasma now...
 
  call ckgrid(xi,'__RHO',id_rho)
  if(id_rho.eq.0) then
     if(.not.inew_device) then
        ier=9999
        call xplasma_errmsg_append(sp, &
             ' ?eqi_fromgeqdsk: inew_device=.FALSE. __RHO mismatch.')
     else
        call xplasma_create_grid(sp,'__RHO',xplasma_rho_coord,xi,id_rho,ier)
     endif
     if(ier.ne.0) go to 1000
  endif

  ! *** nscrunch *** the theta parametrization to allow stable
  !  extrapolation beyond the boundary

  if(jopt.gt.-1) then

     call xplasma_global_info(sp,ier,kmom=imom)
     if(ier.ne.0) go to 1000

     if(jopt.eq.1) then
        if(imom.lt.16) then
           imom=16  ! for equal arc require this number of moments at minimum
           call xplasma_kmom_set(sp,imom,ier)
        endif
        ikmax=1
     else
        if(imom.gt.16) then
           imom=16  ! for scrunched moments, imom > 16 seems harmful...
           call xplasma_kmom_set(sp,imom,ier)
        endif
        ikmax=2
     endif

     allocate(rmmc(0:imom,ns),rmms(0:imom,ns),zmmc(0:imom,ns),zmms(0:imom,ns))

     ismoo = 0
     if(jopt.eq.0) then
        if(i2pi.eq.1) then
           iiopt=2
        else
           iiopt=0
        endif
        call scrunch_rzmom(r,z,nt1,xi,ns, rmmc,rmms,zmmc,zmms, imom,iiopt, ier)

        if(ier.ne.0) then
           msgbuf=' '
           write(msgbuf,*) ' ?eqi_fromgeqdsk:  scrunch_rz error, ier=',ier
           call eqi_fromg_err(msgbuf)
           go to 1000
        endif

        !  axial moments regularization; full or partial 
        !  regeneration of {R,Z} from moments.
        zxbrk_mom = xbrk_mom
        if(jopt.eq.1) then
           if(zxbrk_mom.lt.0.1d0) zxbrk_mom=0.1d0
        endif

        if(zxbrk_mom.gt.0.05d0) then
           !  axial moments regularization; regeneration of {R,Z} from moments.
           ismoo=2
        else
           !  moments trimming only; {R,Z} not modified.
           ismoo=1
        endif

        if(rtol_mom.lt.2.0e-6_R8) then
           zrtol=0.0_R8  ! some trimming done by scrunch_rzmom already...
        else
           zrtol=rtol_mom
        endif

        if(ismoo.ge.2) then

           call scrunch_cleanup(xi,ns, rmmc,rmms,zmmc,zmms, imom, &
                (ismoo.ge.2), zxbrk_mom, zrtol, ier)
           if(ier.ne.0) then
              call xplasma_errmsg_append(sp, &
                   ' ?eqi_fromgeqdsk: scrunch_cleanup error.')
              ier=9999
              go to 1000
           endif

           call xpmom_init(sp,ier, rhoi=xi(1:ns), ikmaxe=ikmax)
           call xpmom_setup(sp,rmmc,rmms,zmmc,zmms,ier)  ! save in xplasma
           if(ier.ne.0) then
              call xplasma_errmsg_append(sp, &
                   ' ?eqi_fromgeqdsk: scrunch_cleanup error.')
              ier=9999
              go to 1000
           endif
        endif

        if(ier.ne.0) go to 1000

        if(ismoo.ge.2) then
           !  reset surfaces (axis was moved slightly)
           r(1:nt1,1)=rmmc(0,1)
           z(1:nt1,1)=zmmc(0,1)
           allocate(csthtk(imom),snthtk(imom))
           do j=2,ns
              do ith=1,nt1
                 if(ith.lt.nt1) then
                    zrnew=rmmc(0,j)
                    zznew=zmmc(0,j)
                    call scrSinCos(the(ith),imom,snthtk,csthtk)
                    do im=1,imom
                       zrnew=zrnew+rmmc(im,j)*csthtk(im)+rmms(im,j)*snthtk(im)
                       zznew=zznew+zmmc(im,j)*csthtk(im)+zmms(im,j)*snthtk(im)
                    enddo

                    r(ith,j) = zrnew
                    z(ith,j) = zznew
                 else
                    r(nt1,j)=r(1,j)
                    z(nt1,j)=z(1,j)
                 endif
              enddo
           enddo
           deallocate(csthtk,snthtk)
        endif

     else if(jopt.eq.1) then

        ! rescrunch using equal arc length
        if(i2pi.eq.1) then
           iiopt=3
        else
           iiopt=1
        endif
        call scrunch_rzmom(r,z,nt1,xi,ns, rmmc,rmms,zmmc,zmms, imom,iiopt, ier)
        if(ier.ne.0) go to 1000

     endif
  endif
 
  ! ----------------------------------------------------------
  ! Ok, the inversion of the psi(R,Z) profile is complete; all
  ! G-EQDSK info can now be fed to XPLASMA
 
  call ckgrid(the,'__CHI',id_chi)
  if(id_chi.eq.0) then
     if(.not.inew_device) then
        ier=9999
        call xplasma_errmsg_append(sp, &
             ' ?eqi_fromgeqdsk: inew_device=.FALSE. __CHI mismatch.')
     else
        call xplasma_create_grid(sp,'__CHI',xplasma_theta_coord,the,id_chi,ier)
     endif
     if(ier.ne.0) go to 1000
  endif

  if(jopt.eq.0) then
     ! default BC (not a knot)
     ibc1=0
     bcr=0
     bcz=0
     ! bdy conditions in rho direction: not a knot on axis, 1st div. diff @edge
     ! ibc1=1
     ! bcr=(r(1:nt1,ns)-r(1:nt1,ns-1))/(xi(ns)-xi(ns-1))
     ! bcz=(z(1:nt1,ns)-z(1:nt1,ns-1))/(xi(ns)-xi(ns-1))
  else
     ! for equal arc: use default BC
     ibc1=0
     bcr=0
     bcz=0
  endif

  ! preserve moments representation during eq rebuild
  if(ismoo.ge.2) call xpmom_protect(sp,iertmp)

  call xplasma_rzmag(sp,id_chi,id_rho, &
       r,z,id_R,id_Z,ier, ibc0=0,ibc1=ibc1,Rbc1=bcr,Zbc1=bcz, &
       lhermite=ihermite)
  if(ier.ne.0) go to 1000
 
  if(ismoo.ge.2) call xpmom_unprotect(sp,iertmp)

  ! compute grad(psi) vs. (theta,rho)

  gpsi_geq(:,1) = ZERO
  ii = nt1-1
  allocate(gpsi_wk(ii,2))
  do i=2,ns
     call ezspline_gradient(pspl,ii,r(1:ii,i),z(1:ii,i),gpsi_wk,ier)
     if(ier.ne.0) go to 1000
     do j=1,ii
        gpsi_geq(j,i)=sqrt(gpsi_wk(j,1)**2+gpsi_wk(j,2)**2)
     enddo
     gpsi_geq(nt1,i)=gpsi_geq(1,i)
  enddo

  call xplasma_create_2dprof(sp,'gpsi_geq',id_chi,id_rho,gpsi_geq,id_gpsi,ier,&
       ispline=2,label='grad(Psi) from EFIT Psi(R,Z)',units='Wb/rad/m')
  deallocate(gpsi_wk)
  if(ier.ne.0) go to 1000

  ! Integration weights

  allocate(xint(10*ii),thint(10*ii),wint(10*ii),wkint(10*ii,4))

  call integ_wts(xwts,wts)
  do i=1,ii
     do j=1,10
        thint((i-1)*10 + j) = the(i) + xwts(j)*(the(i+1)-the(i))
        wint((i-1)*10 + j) = wts(j)*(the(i+1)-the(i))
     enddo
  enddo

  id4(1)=id_R;    ideriv4(1)=0
  id4(2)=id_gpsi; ideriv4(2)=0
  id4(3)=id_R;    ideriv4(3)=1
  id4(4)=id_Z;    ideriv4(4)=1

  allocate(icurt(ns))
  icurt(1)=0.0_R8
  if(icur_psi.eq.1) then
     !  q(psi) & total current from Psi(R,Z)
     call surfint(xi(ns))
     icurt(ns) = gpsi2r2i*dvdpsi/(twopi*mu0)
  else
     !  q(psi) & total current from EFIT file data
     icurt(ns)=abs(geq%current_a)
  endif

  ins_test = 1.0_R8/xbrkii + 0.5_R8

  if(ns.gt.ins_test) then
     nsbrk = ns - 2
     xibrk = xi(nsbrk)
  else
     nsbrk = ns - 1
  endif

  do j=2,nsbrk
     call surfint(xi(j))
     icurt(j)=gpsi2r2i*dvdpsi/(twopi*mu0)
     if(ns.gt.ins_test) then
        if(xi(j).ge.xbrkii) then
           xibrk=xi(j)
           nsbrk=j
           exit
        endif
     endif
  enddo

  if(ns.gt.ins_test) then
     ! parabolic piece fit to reach known total toroidal current
     x1=xibrk-xi(nsbrk-1)
     x2=xibrk-xi(nsbrk-2)
     ! parabolic fit to get slope d(icurt)/dx at xibrk:
     bicur = (x2*x2*(icurt(nsbrk-1)-icurt(nsbrk)) - &
          x1*x1*(icurt(nsbrk-2)-icurt(nsbrk))) / (x1*x2*(x1-x2))
     ! parabolic fit to match point and slope at xibrk and reach the
     ! boundary value icurt(ns):
     x3=xi(ns)-xibrk
     aicur = (icurt(ns)-icurt(nsbrk)-bicur*x3)/(x3*x3)
     do j=nsbrk+1,ns-1
        xxx = xi(j)-xibrk
        icurt(j)=icurt(nsbrk) + xxx*(bicur + xxx*aicur)
     enddo
  endif

  ! Boundary Conditions for f(xi) splines
 
  zero=0.0_R8
  iknot=0             ! not-a-knot BC
  idrv1=1             ! fix df/dxi BC, following arg gives value of derivative
 
  ! Splines...

  call xplasma_create_1dprof(sp,'ICURT',id_rho,icurt,id_curt,ier, &
       ispline=2,ibca=idrv1,zbca=ZERO,ibcb=iknot, &
       label='Enclosed toroidal current', units='A')
  if(ier.ne.0) go to 1000                   ! I_phi, dI/dxi = 0 on axis

  call xplasma_create_1dprof(sp,'PSI',id_rho,psi,id_psi,ier, &
       ispline=2,ibca=idrv1,zbca=ZERO,ibcb=iknot, &
       label='Poloidal Flux',units='Wb/rad')
  if(ier.ne.0) go to 1000                   ! psi, dpsi/dxi = 0 on axis
 
  call xplasma_create_1dprof(sp,'PHITOR',id_rho,phi,id_phi,ier, &
       ispline=2,ibca=idrv1,zbca=ZERO,ibcb=iknot, &
       label='Toroidal Flux',units='Wb')
  if(ier.ne.0) go to 1000                   ! phi, dphi/dxi = 0 on axis
 
  call xplasma_create_1dprof(sp,'G',id_rho,g,id_g,ier, &
       ispline=2,ibca=idrv1,zbca=ZERO,ibcb=iknot, &
       label='R*B_phi (toroidal field function)',units='T*m')
  if(ier.ne.0) go to 1000                   ! g, dg/dxi = 0 on axis, @edge
 
  call xplasma_create_1dprof(sp,'Q',id_rho,q,id_q,ier, &
       ispline=2,ibca=iknot,ibcb=iknot, &
       label='q profile',units=' ')
  if(ier.ne.0) go to 1000                   ! q, not-a-knot BCs
 
  ! Hermite pressure profile

  call xplasma_create_1dprof(sp,'P',id_rho,p,id_p,ier, &
       ispline=1,ibca=idrv1,zbca=ZERO,ibcb=iknot, &
       label='MHD equilibrium pressure',units='Pa')
  if(ier.ne.0) go to 1000                   ! pressure, dP/dxi=0 on axis

  call xplasma_set_pressure_id(sp,id_p,ier)
  if(ier.ne.0) go to 1000

  ! Linear dP/dpsi ... might be used to reconstruct pressure

  call xplasma_create_1dprof(sp,'dPdPsi',id_rho,pprime,id_pp,ier, &
       ispline=0,label='d[MHD eq pressure]/dPsi',units='Pa/(Wb/rad)')
 
  !  set the field and current orientations

  call xplasma_field_ccw_set(sp,ier,bphi_ccw=isnccwb,jphi_ccw=isnccwi)
  if(ier.ne.0) go to 1000
 
  ! ...OK, geom. and psi and g are set, compute B,BR,BZ vs. (xi,the)
  ! the field and current orientations are specified.
 
  call xplasma_eqcheck(sp,.FALSE.,icount,ier)
  if(ier.ne.0) go to 1000

  ! now the limiters and the R,Z grids (if required)
 
  if(irz.eq.1) then
 
     call cklim
     if(id_lim.eq.0) then
        isize_lim=geq%Limitr
        call xplasma_mklim_contour(sp, &
             geq%Rlim_m(1:isize_lim),geq%Zlim_m(1:isize_lim),iwarn,ier)

        if(ier.ne.0) go to 1000
     endif
 
     call ckgrid(Rgrid,'__RGRID',id_Rgrid)
     call ckgrid(Zgrid,'__ZGRID',id_Zgrid)
     if((id_Rgrid.eq.0).or.(id_Zgrid.eq.0)) then
        if(.not.inew_device) then
           ier=9999
           call xplasma_errmsg_append(sp, &
                ' ?eqi_fromgeqdsk: inew_device=.FALSE. __RGRID and/or __ZGRID  mismatch.')
        else
           call xplasma_create_RZgrid(sp,Rgrid,Zgrid,id_Rgrid,id_Zgrid,ier)
        endif
        if(ier.ne.0) go to 1000
     endif
 
     ! now:  psi(R,Z) with positive sign except psi(axis)=0.0
 
     idprev=id_psi
     call xplasma_create_2dprof(sp,'PSI_RZ', &
          id_Rgrid,id_Zgrid,psirz_sm,id_psi,ier, &
          ispline=1,assoc_id=idprev, &
          label='Poloidal Flux',units='Wb/rad')
     if(ier.ne.0) go to 1000
 
  ! and GS error
 
     call xplasma_create_2dprof(sp,'GS_ERROR', &
          id_Rgrid,id_Zgrid,gserr_geq,id_gserr,ier, &
          ispline=1, label='GS error estimate',units=' ')
     if(ier.ne.0) go to 1000
 
  ! converting psi(R,Z) -> BR,BZ takes into account isnccwi
  ! get B(R,Z) components
 
     call xplasma_brz(sp,eqm_bpsi,ier)
     if(ier.ne.0) go to 1000
 
     call eqm_geq_curtime(ztime1)
     call xplasma_time_set(sp,ztime1,ier)
     if(ier.ne.0) go to 1000

  endif

  if(ipsi0) then
     call xplasma_save_Psi0(sp,Psi0,ier)
     if(ier.ne.0) go to 1000
  endif

  ! ----------------------------------------------------------
 
1000 continue
 
  if(allocated(rgrid)) deallocate(rgrid)
  if(allocated(zgrid)) deallocate(zgrid)
  if(allocated(psi)) deallocate(psi, &
       psi_geq, psirz_sm, gserr_geq, p, g, q, the)
  if(allocated(pprime)) deallocate(pprime)
  if(allocated(phi)) deallocate(phi, xi)
  if(allocated(r)) deallocate(r, z, gpsi_geq)
  if(allocated(psis)) deallocate(psis)
  if(allocated(rmmc)) deallocate(rmmc,rmms,zmmc,zmms)
  if(allocated(x2d)) deallocate(x2d)
  if(allocated(lsurf)) deallocate(lsurf)
  if(allocated(wint)) deallocate(wint,thint)
  if(allocated(icurt)) deallocate(icurt)

  if(iez.eq.1) call ezspline_free(pspl, iok)

  if(zauthor.eq.xplasma_root) then
     call xplasma_author_clear(sp,'eqi_fromgeqdsk',iertmp)
  endif
  
  return
 
  contains

    subroutine load_seek

      ! extend contour indexing to ease handling of BCs

      rseek(1:nthad)=radapt(1:nthad,isw)
      zseek(1:nthad)=zadapt(1:nthad,isw)
      rseek(-1:0)=rseek(nthad-2:nthad-1)
      zseek(-1:0)=zseek(nthad-2:nthad-1)
      zwk(-1:0)=zwk(nthad-2:nthad-1)-twopi
      rseek(nthad+1:nthad+2)=rseek(2:3)
      zseek(nthad+1:nthad+2)=zseek(2:3)
      zwk(nthad+1:nthad+2)=zwk(2:3)+twopi

    end subroutine load_seek

    subroutine ckgrid(zdata,zname,id)

      real*8, dimension(:), intent(in) :: zdata
      character*(*), intent(in) :: zname
      integer, intent(out) :: id

      !-----------
      !  check if grid is already defined & identical to new grid

      integer :: i,isizo
      real*8 :: zgrido(size(zdata)),ztol

      real*8, parameter :: loctol = 1.0d-10
      !-----------

      call xplasma_find_item(sp,zname,id,iertmp,nf_noerr=.TRUE.)
      if(id.eq.0) return

      call xplasma_grid_size(sp,id,isizo,iertmp)
      if(iertmp.ne.0) then
         id=0
         return
      endif

      call xplasma_grid(sp,id,zgrido,iertmp)
      if(iertmp.ne.0) then
         id=0
         return
      endif

      ztol=loctol*max(abs(zgrido(1)),abs(zgrido(isizo)))

      do i=1,isizo
         if(abs(zgrido(i)-zdata(i)).gt.ztol) then
            id=0  ! value changed...
            exit
         endif
      enddo

    end subroutine ckgrid

    subroutine cklim

      !-------------
      !  check if limiter has changed...

      integer :: itype,isizo,i,iflag
      real*8, dimension(:), allocatable :: rlimo,zlimo
      !-------------

      iflag = 0

      call xplasma_find_item(sp,'__LIMITER',id_lim,iertmp,nf_noerr=.TRUE.)
      if(id_lim.eq.0) return

      call xplasma_lim_info(sp,iertmp,  &
           itype=itype)
      if(itype.ne.100) then
         id_lim=0
         return
      endif

      ! there is a contour limiter

      call xplasma_lim_info(sp,iertmp, npts=isizo)
      if(isizo.ne.geq%Limitr) then
         iflag = 1   ! size mismatch
      endif

      if(iflag.eq.0) then
         allocate(rlimo(isizo),zlimo(isizo))
         call xplasma_lim_info(sp,iertmp, rpts=rlimo, zpts=zlimo)

         do i=1,isizo
            if((rlimo(i).ne.geq%Rlim_m(i)).or.(zlimo(i).ne.geq%Zlim_m(i))) then
               iflag = 1   ! value mismatch
               exit
            endif
         enddo

         deallocate(rlimo,zlimo)
      endif

      if(iflag.eq.1) then
         msgbuf=' %eqi_fromgeqdsk: limiter contour appears to change; old limiter retained.'
         call eqi_fromg_err(msgbuf)
      endif

    end subroutine cklim
      
    subroutine fnxz(indx,zval)
      ! increment index in direction (idir); find next distinct
      ! boundary Z value; wrap index if reach either end of boundary array

      integer, intent(inout) :: indx
      real*8, intent(inout) :: zval

      indx = indx + idir
      if(indx.lt.1) indx=inumx
      if(indx.gt.inumx) indx=1

      zval = rzwk(indx,2)

    end subroutine fnxz

    REAL*8 function quadintrp(x,xx,ff)
      real*8, intent(in) :: x
      real*8, intent(in) :: xx(3),ff(3)

      ! form quadratic from 3 distinct pts 
      !   {(xx(1),ff(1)),(xx(2),ff(2)),(xx(3),ff(3))}
      ! evaluate at x and return result as function value

      real*8 :: xm,xp,fm,fp,aa,bb,denom,answer,xuse

      !--------------

      xuse=x-xx(2)

      xm = xx(1)-xx(2); xp = xx(3)-xx(2)
      fm = ff(1)-ff(2); fp = ff(3)-ff(2)

      denom=xp*xm*(xp-xm)
      aa=(fp*xm-fm*xp)/denom
      bb=(fm*xp*xp-fp*xm*xm)/denom

      answer = xuse*(aa*xuse + bb)

      quadintrp = answer + ff(2)

    end function quadintrp

    subroutine surfint(x)
      real*8, intent(in) :: x

      ! compute numerical integrations for: 
      !     dV/dpsi, <|grad(psi)|**2/R**2>, <1/R**2>

      !-------------------------
      integer :: ith
      real*8 :: dl,Rval,gpsival
      !-------------------------
      xint=x

      call xplasma_eval_prof(sp,id4,id_chi,thint,id_rho,xint,wkint,ier, &
           ideriv1s=ideriv4)

      dvdpsi=ZERO
      gpsi2r2i=ZERO

      do ith=1,size(thint)
         dl=sqrt(wkint(ith,3)**2+wkint(ith,4)**2)
         Rval=wkint(ith,1)
         gpsival=wkint(ith,2)

         dvdpsi = dvdpsi + dl*twopi*Rval*wint(ith)/gpsival
         gpsi2r2i = gpsi2r2i + dl*twopi*gpsival*wint(ith)/Rval
      enddo

      gpsi2r2i = gpsi2r2i/dvdpsi

    end subroutine surfint

end subroutine eqi_fromGeqdsk

subroutine psirz_smooth(nR,nZ,R,Z,psi,psism,nbdy,rbdy,zbdy)

  !  smooth Psi(R,Z) but only in the region beyond the plasma boundary
  !  idea is to reduce spline interpolation effects near the boundary
  !  for contouring...

  implicit NONE

  integer, intent(in) :: nR,nZ         ! array dimensions
  real*8, intent(in) :: R(nR),Z(nZ)    ! base grids

  real*8, intent(in) :: psi(nR,nZ)     ! original Psi
  real*8, intent(out) :: psism(nR,nZ)  ! smoothed Psi

  integer, intent(in) :: nbdy          ! #pts in boundary contour
  real*8, intent(in) :: rbdy(nbdy),zbdy(nbdy) ! bdy contour itself

  !--------------------
  !  local

  integer :: ibdy,i,j,inum
  real*8, dimension(:), allocatable :: rb,zb  ! bdy, cleaned up...

  real*8 :: rmin,rmax,zmin,zmax,wkmin,wkmax,rx,zx
  real*8 :: delr,delz,d

  integer :: iz1,iz2,ir1,ir2  ! index ranges which hit boundary
  integer :: im,ip,jm,jp

  !  index info for each bdy segment
  integer, dimension(:,:), allocatable :: ircr,izcr

  real*8, dimension(:,:), allocatable :: dsmoo
  real*8, dimension(:), allocatable :: zwk

  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE  = 1.0d0
  real*8, parameter :: THREE= 3.0d0
  real*8, parameter :: SIX  = 6.0d0

  !------------------------------------------------------------------------
  !  0. assume evenly spaced grid

  delr=R(2)-R(1)
  delz=Z(2)-Z(1)

  !--------------------
  !  1. clean up the boundary: remove duplicate points; make sure it closes
  !     at end.  Find global min & max of bdy

  ibdy = nbdy+1
  allocate(rb(ibdy),zb(ibdy))

  ibdy=1
  rb(1)=rbdy(1)
  zb(1)=zbdy(1)
  rmin=rb(1)
  rmax=rmin
  zmin=zb(1)
  zmax=zb(1)

  do i=2,nbdy
     rmin=min(rmin,rbdy(i))
     rmax=max(rmax,rbdy(i))
     zmin=min(zmin,zbdy(i))
     zmax=max(zmax,zbdy(i))
     if(max(abs(rbdy(i)-rbdy(i-1)),abs(zbdy(i)-zbdy(i-1))).gt.ZERO) then
        ibdy=ibdy+1
        rb(ibdy)=rbdy(i)
        zb(ibdy)=zbdy(i)
     endif
  enddo

  if((rb(ibdy).ne.rb(1)).or.(zb(ibdy).ne.zb(1))) then
     ibdy=ibdy+1
     rb(ibdy)=rb(1)
     zb(ibdy)=zb(1)
  endif

  !--------------------
  !  2. find global index range

  ir1=0
  do i=1,nR
     if(ir1.eq.0) then
        if(R(i).gt.rmin) ir1=i
     endif
     if(R(i).ge.rmax) then
        ir2=i-1
        exit
     endif
  enddo

  iz1=0
  do i=1,nZ
     if(iz1.eq.0) then
        if(Z(i).gt.zmin) iz1=i
     endif
     if(Z(i).ge.zmax) then
        iz2=i-1
        exit
     endif
  enddo

  !--------------------
  !  3. find intersecting segments for each horizontal and vertical line
     
  allocate(ircr(ibdy,nR),izcr(ibdy,nZ))

  ircr(1,:)=0
  izcr(1,:)=0

  do i=ir1,ir2
     do j=1,ibdy-1
        if((R(i)-rb(j))*(R(i)-rb(j+1)).le.ZERO) then
           if(rb(j).ne.rb(j+1)) then
              inum=ircr(1,i)
              inum=inum+1
              ircr(1,i)=inum
              ircr(inum+1,i)=j
           endif
        endif
     enddo
  enddo

  do i=iz1,iz2
     do j=1,ibdy-1
        if((Z(i)-zb(j))*(Z(i)-zb(j+1)).le.ZERO) then
           if(zb(j).ne.zb(j+1)) then
              inum=izcr(1,i)
              inum=inum+1
              izcr(1,i)=inum
              izcr(inum+1,i)=j
           endif
        endif
     enddo
  enddo

  !--------------------
  !  4. set up smoothing control based on distance to bdy

  allocate(dsmoo(nR,nZ)); dsmoo=ONE

  do i=ir1,ir2
     wkmin=Z(nZ)
     wkmax=Z(1)
     do j=1,ircr(1,i)
        inum=ircr(j+1,i)  ! bdy segment #
        zx = zb(inum) + &
             (R(i)-rb(inum))*(zb(inum+1)-zb(inum))/(rb(inum+1)-rb(inum))
        wkmin=min(wkmin,zx)
        wkmax=max(wkmax,zx)
     enddo
     do j=1,nZ
        if(z(j).lt.wkmin) then
           dsmoo(i,j)=min(dsmoo(i,j),(wkmin-z(j))/delz)
        else if(z(j).gt.wkmax) then
           dsmoo(i,j)=min(dsmoo(i,j),(z(j)-wkmax)/delz)
        else
           dsmoo(i,j)=ZERO
        endif
     enddo
  enddo

  do i=iz1,iz2
     wkmin=R(nR)
     wkmax=R(1)
     do j=1,izcr(1,i)
        inum=izcr(j+1,i)  ! bdy segment #
        rx = rb(inum) + &
             (Z(i)-zb(inum))*(rb(inum+1)-rb(inum))/(zb(inum+1)-zb(inum))
        wkmin=min(wkmin,rx)
        wkmax=max(wkmax,rx)
     enddo
     do j=1,nR
        if(r(j).lt.wkmin) then
           dsmoo(j,i)=min(dsmoo(j,i),(wkmin-r(j))/delr)
        else if(r(j).gt.wkmax) then
           dsmoo(j,i)=min(dsmoo(j,i),(r(j)-wkmax)/delr)
        else
           dsmoo(j,i)=ZERO
        endif
     enddo
  enddo

  !--------------------
  !  5.  apply smoothing -- hat function convolution
  !        d=dsmoo(i,j); 0 <= d <= 1
  !      f(i,j) -> (d/6)*f(i-1,j) + ((2+d)/3)*f(i,j) + (d/6)*f(i+1,j)
  !        and repeat in j dimension

  do j=1,nZ
     do i=1,nR
        im=max(1,i-1)
        ip=min(nR,i+1)
        d=dsmoo(i,j)
        if(d.gt.ZERO) then
           psism(i,j)=(d/SIX)*(psi(im,j)+psi(ip,j))+((THREE-d)/THREE)*psi(i,j)
        else
           psism(i,j)=psi(i,j)
        endif
     enddo
  enddo
        
  allocate(zwk(nZ))
  do i=1,nR
     do j=1,nZ
        zwk(j)=psism(i,j)
     enddo
     do j=1,nZ
        jm=max(1,j-1)
        jp=min(nZ,j+1)
        d=dsmoo(i,j)
        if(d.gt.ZERO) then
           psism(i,j)=(d/SIX)*(zwk(jm)+zwk(jp))+((THREE-d)/THREE)*zwk(j)
        else
           psism(i,j)=zwk(j)
        endif
     enddo
  enddo

  deallocate(zwk)
  deallocate(rb,zb)
  deallocate(ircr,izcr)
  deallocate(dsmoo)

!dbg  allocate(zwk(2*max(nR,nZ))); zwk=0.0d0

!dbg  call r8_grf3f2(R,psi,psism,Z,zwk,nR,nZ,nR, &
!dbg       'R','Psi','Psi(smoothed)','Z', &
!dbg       'm','Wb/rad','m', &
!dbg       'eqi_fromgeqdsk:psirz_smooth debug plot','xplasma2',0)

!dbg  call r8_grafx1(rbdy,zbdy,nbdy,'m','m','EFIT est. Psi(a) contour', &
!dbg       'eqi_fromgeqdsk:psirz_smooth debug plot','xplasma2')

!dbg  deallocate(zwk)

end subroutine psirz_smooth

subroutine psirz_minmax(nR,nZ,R,Z,psi, &
     nbdy,rbdy,zbdy,psimin,psimax,nsing)

  !  Find min & max Psi(R,Z) inside the passed boundary contour
  !  This uses the same method as psirz_smooth to find what is inside
  !  and what is outside; but instead of smoothing it just computes the
  !  min and max...

  !  also count the number of "singularities" = local maxima or minima
  !  that are fully enclosed w/in the boundary

  implicit NONE

  integer, intent(in) :: nR,nZ         ! array dimensions
  real*8, intent(in) :: R(nR),Z(nZ)    ! base grids

  real*8, intent(in) :: psi(nR,nZ)     ! original Psi

  integer, intent(in) :: nbdy          ! #pts in boundary contour
  real*8, intent(in) :: rbdy(nbdy),zbdy(nbdy) ! bdy contour itself

  real*8, intent(out) :: psimin,psimax   ! results...
  integer, intent(out) :: nsing        ! singularity count

  !--------------------
  !  local

  integer :: ibdy,i,j,inum
  real*8, dimension(:), allocatable :: rb,zb  ! bdy, cleaned up...

  real*8 :: wkmin,wkmax,rx,zx,rmin,rmax,zmin,zmax

  integer :: iz1,iz2,ir1,ir2  ! index ranges which hit boundary
  integer :: im,ip,jm,jp,ict

  !  index info for each bdy segment
  integer, dimension(:,:), allocatable :: ircr,izcr

  real*8, dimension(:,:), allocatable :: dsmoo

  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE  = 1.0d0
  real*8, parameter :: THREE= 3.0d0
  real*8, parameter :: SIX  = 6.0d0

  !------------------------------------------------------------------------
  !  0. set up evenly spaced R and Z grids

  psimin=0.0d0
  psimax=0.0d0
  nsing=0
  
  !--------------------
  !  1. clean up the boundary: remove duplicate points; make sure it closes
  !     at end.  Find global min & max of bdy

  ibdy = nbdy+1
  allocate(rb(ibdy),zb(ibdy))

  ibdy=1
  rb(1)=rbdy(1)
  zb(1)=zbdy(1)
  rmin=rb(1)
  rmax=rmin
  zmin=zb(1)
  zmax=zb(1)

  do i=2,nbdy
     rmin=min(rmin,rbdy(i))
     rmax=max(rmax,rbdy(i))
     zmin=min(zmin,zbdy(i))
     zmax=max(zmax,zbdy(i))
     if(max(abs(rbdy(i)-rbdy(i-1)),abs(zbdy(i)-zbdy(i-1))).gt.ZERO) then
        ibdy=ibdy+1
        rb(ibdy)=rbdy(i)
        zb(ibdy)=zbdy(i)
     endif
  enddo

  if((rb(ibdy).ne.rb(1)).or.(zb(ibdy).ne.zb(1))) then
     ibdy=ibdy+1
     rb(ibdy)=rb(1)
     zb(ibdy)=zb(1)
  endif

  !--------------------
  !  2. find global index range

  ir1=0
  do i=1,nR
     if(ir1.eq.0) then
        if(R(i).gt.rmin) ir1=i
     endif
     if(R(i).ge.rmax) then
        ir2=i-1
        exit
     endif
  enddo

  iz1=0
  do i=1,nZ
     if(iz1.eq.0) then
        if(Z(i).gt.zmin) iz1=i
     endif
     if(Z(i).ge.zmax) then
        iz2=i-1
        exit
     endif
  enddo

  !--------------------
  !  3. find intersecting segments for each horizontal and vertical line
     
  allocate(ircr(ibdy,nR),izcr(ibdy,nZ))

  ircr(1,:)=0
  izcr(1,:)=0

  do i=ir1,ir2
     do j=1,ibdy-1
        if((R(i)-rb(j))*(R(i)-rb(j+1)).le.ZERO) then
           if(rb(j).ne.rb(j+1)) then
              inum=ircr(1,i)
              inum=inum+1
              ircr(1,i)=inum
              ircr(inum+1,i)=j
           endif
        endif
     enddo
  enddo

  do i=iz1,iz2
     do j=1,ibdy-1
        if((Z(i)-zb(j))*(Z(i)-zb(j+1)).le.ZERO) then
           if(zb(j).ne.zb(j+1)) then
              inum=izcr(1,i)
              inum=inum+1
              izcr(1,i)=inum
              izcr(inum+1,i)=j
           endif
        endif
     enddo
  enddo

  !--------------------
  !  4. set up smoothing control based on distance to bdy

  allocate(dsmoo(nR,nZ)); dsmoo=ZERO

  do i=ir1,ir2
     wkmin=Z(nZ)
     wkmax=Z(1)
     do j=1,ircr(1,i)
        inum=ircr(j+1,i)  ! bdy segment #
        zx = zb(inum) + &
             (R(i)-rb(inum))*(zb(inum+1)-zb(inum))/(rb(inum+1)-rb(inum))
        wkmin=min(wkmin,zx)
        wkmax=max(wkmax,zx)
     enddo
     do j=1,nZ
        if(z(j).lt.wkmin) then
           continue
        else if(z(j).gt.wkmax) then
           continue
        else
           dsmoo(i,j)=ONE
        endif
     enddo
  enddo

  do i=iz1,iz2
     wkmin=R(nR)
     wkmax=R(1)
     do j=1,izcr(1,i)
        inum=izcr(j+1,i)  ! bdy segment #
        rx = rb(inum) + &
             (Z(i)-zb(inum))*(rb(inum+1)-rb(inum))/(zb(inum+1)-zb(inum))
        wkmin=min(wkmin,rx)
        wkmax=max(wkmax,rx)
     enddo
     do j=1,nR
        if(r(j).lt.wkmin) then
           continue
        else if(r(j).gt.wkmax) then
           continue
        else
           dsmoo(j,i)=ONE
        endif
     enddo
  enddo

  !--------------------
  !  5.  apply smoothing -- hat function convolution
  !        d=dsmoo(i,j); 0 <= d <= 1
  !      f(i,j) -> (d/6)*f(i-1,j) + ((2+d)/3)*f(i,j) + (d/6)*f(i+1,j)
  !        and repeat in j dimension

  ict=0

  do j=1,nZ
     do i=1,nR
        im=max(1,i-1)
        ip=min(nR,i+1)
        if(dsmoo(i,j).gt.ZERO) then
           if(ict.eq.0) then
              psimin=psi(i,j)
              psimax=psi(i,j)
           else
              psimin=min(psimin,psi(i,j))
              psimax=max(psimax,psi(i,j))
              call cksing
           endif
           ict = ict + 1
        endif
     enddo
  enddo
        
  deallocate(rb,zb)
  deallocate(ircr,izcr)
  deallocate(dsmoo)

CONTAINS

  subroutine cksing

    integer :: ii,jj,ilt,igt
    real*8 :: dmin

    if((i.eq.1).or.(i.eq.nR)) return
    if((j.eq.1).or.(j.eq.nZ)) return

    dmin=min(dsmoo(i-1,j),dsmoo(i,j),dsmoo(i+1,j))

    ilt=0
    igt=0

    if(psi(i-1,j).lt.psi(i,j)) ilt=ilt+1
    if(psi(i+1,j).lt.psi(i,j)) ilt=ilt+1

    if(psi(i-1,j).gt.psi(i,j)) igt=igt+1
    if(psi(i+1,j).gt.psi(i,j)) igt=igt+1

    do jj=j-1,j+1,2
       do ii=i-1,i+1
          dmin=min(dmin,dsmoo(ii,jj))
          if(psi(ii,jj).lt.psi(i,j)) ilt=ilt+1
          if(psi(ii,jj).gt.psi(i,j)) igt=igt+1
       enddo
    enddo

    if(dmin.le.ZERO) return

    ! OK all neighboring points are inside boundary
    !   so, eligible for singularity check.
    !     singularity means: less than, or greater than, all neighbors

    if(ilt.eq.8) nsing = nsing + 1
    if(igt.eq.8) nsing = nsing + 1

  end subroutine cksing
end subroutine psirz_minmax
