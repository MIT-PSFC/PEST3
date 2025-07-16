module xplasma_sol

  !  xplasma "scrape off layer" module -- xplasma_sol
  !  At present this contains:
  !     -- routines for defining a limiter (2 methods provided)
  !     -- a routine for checking consistency of MHD equilibrium with the
  !        given limiter definition.

  use xplasma_obj
  use xplasma_ctran
  use eqi_rzbox_module  ! for sp ptr
  implicit NONE

  private

  public :: xplasma_mklim_contour
  public :: xplasma_mklim_special
  public :: xplasma_lim_info
  public :: xplasma_limcheck
  public :: xplasma_limcon,xplasma_bdycon
  public :: xplasma_lim_rzminmax
  public :: xplasma_lim_distance, xplasma_lim_dist1,xplasma_lim_distn
  public :: xplasma_fxcep

  interface xplasma_lim_distance
     module procedure xplasma_lim_dist1
     module procedure xplasma_lim_distn
  end interface

  ! (default) surface to use to test for intersection with
  ! mechanical limiter.  The consistency requirement of the equilibrium 
  ! with a limiter specification is that the equilibrium boundary not 
  ! intersect a limiter.  In practice we test a surface offset somewhat
  ! inside the boundary; this parameter identifies what surface is the
  ! default.  See "xplasma_limcheck"...

  real*8, parameter, public :: xplasma_rho_limcheck = 0.95d0

  contains
    subroutine xplasma_mklim_contour(s,rpts,zpts,iwarn,ier)
 
      ! define an axisymmetric limiter/wall from a closed (R,Z) contour

      type (xplasma), pointer :: s
      real*8,intent(in), dimension(:) :: Rpts,Zpts  ! the contour
      !  Both arrays must be of the same length;

      !  Note: must be a closed curve with Rpts(1)=Rpts(npts) and 
      !  Zpts(1)=Zpts(npts).  The curve must not cross itself (this is tested).

      integer, intent(out) :: iwarn  ! warning flag: if limiter is redefined.
      integer, intent(out) :: ier    ! completion status code: 0=OK

      !----------------------------------------------
      integer :: id,itype,iertmp,isize,id_R
      integer :: idata(7),ilocmmx,npts
      real*8, dimension(:), allocatable :: r8data
      !----------------------------------------------

      ier=0
      iwarn=0

      npts=size(rpts)
      if(npts.ne.size(zpts)) then
         ier=510
         iwarn=510
         call xplasma_errmsg_append(s,&
              'limiter contour Rpts(:) and Zpts(:) array lengths unequal.')
         return
      endif

      !  set warning flag if something named "__LIMITER" already exists.

      call xplasma_find_item(s,'__LIMITER',id,iertmp, nf_noerr=.TRUE.)
      if(id.ne.0) iwarn=1

      itype=100

      idata = 0
      idata(1)=1
      idata(4)=Npts

      isize= 6*Npts
      allocate(r8data(isize+4))
      ilocmmx=isize+1
      idata(7)=ilocmmx

      ! call adapted f77-xplasma routine to check contour
      sp => s
      call eqi_cbdy_pts(r8data,Npts,Rpts,Zpts,ier)
      if(ier.ne.0) then
         deallocate(r8data)
         return
      endif

      r8data(ilocmmx)=minval(Rpts)
      r8data(ilocmmx+1)=maxval(Rpts)
      r8data(ilocmmx+2)=minval(Zpts)
      r8data(ilocmmx+3)=maxval(Zpts)

      !----------------------------------------------
      !  create limiter dataset "black box"...

      call xplasma_author_set(s,xplasma_root,iertmp)

      call xplasma_create_blackbox(s,'__LIMITER',itype,id,ier, &
           iarray=idata, r8array=r8data, &
           label='axisymmetric contour limiter data')

      call xplasma_author_clear(s,xplasma_root,iertmp)

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' xplasma_mklim_contour: "__LIMITER" create failure.')
      endif

      deallocate(r8data)

      call xplasma_common_ids(s,iertmp, id_R)
      if(id_R.gt.0) then
         call xplasma_limcheck(s,ier)
      endif

    end subroutine xplasma_mklim_contour

    !==========================================================

    subroutine xplasma_mklim_special(s,iwarn,ier, &
         dist, Rmin,Rmax, Zmin,Zmax, &
         nlines, Rl, Zl, thl, &
         ncircs, Rc, Zc, rad)

      ! define an axisymmetric limiter/wall by various means other than
      ! an (R,Z) contour, as described in the optional arguments

      type (xplasma), pointer :: s

      integer, intent(out) :: iwarn  ! warning flag: if limiter is redefined.
      integer, intent(out) :: ier    ! completion status code: 0=OK

      real*8, intent(in), optional :: dist   ! if >=0, limiter at/outside
      !  plasma boundary.  if >0, limiter is located (dist) meters from the
      !  plasma boundary all the way around.  Default: -1.0, i.e. no such 
      !  limiter.

      real*8, intent(in),  optional :: Rmin,Rmax, Zmin,Zmax
      !  Rectangular "box" limiter.  It is expected that the plasma boundary
      !  will not cross this box boundary.

      !--------------------------------------------
      !  "Legacy" limiter specifications...

      integer, intent(in), optional :: nlines   ! # of "line limiters"
      real*8, intent(in), dimension(:), optional :: Rl,Zl,thl

      integer, intent(in), optional :: ncircs   ! # of "circle limiters"
      real*8, intent(in), dimension(:), optional :: Rc,Zc,rad

      !  set up a TRANSP-style machine boundary specification.  This consists
      !  two lists:  a list of circles (R,Z) & radius, and a list of infinite
      !  lines (R,Z) and inclination angle, thl (degrees).

      !  for each circle, the plasma is excluded from the space inside
      !  the circle (sign=-1.) or is excluded from the space outside
      !  (sign=+1.).  The sign is determined by whether the plasma magnetic
      !  axis point (R0,Z0) is inside or outside the circle.

      !  for each infinite line, a start point and orientation angle are
      !  given.  the angle (thl) is in degrees, with convention:

      !    thl(i)=0.0d0 -- horizontal orientation
      !    thl(i)=45.0d0 -- line goes up & to right (+R,+Z) & 
      !                     down to left (-R,-Z)
      !    thl(i)=90.0d0 -- vertical orientation

      !  **CAUTION TRANSP USERS** in TRANSP limiters theta=45.0d is up and to
      !  the left (-R,+Z). -- this is a different convention ***

      !  each line partitions space into two regions, from one of which the
      !  plasma is excluded.  Using (R0,Z0) to test, a sign value is assigned
      !  which determines which half plane is excluded.

      !  the result:  excluded space = union of all excluded spaces of
      !                 individual circle/line "limiters".

      !               included space = all space - excluded space

      !  the plasma boundary must lie entirely inside the included space.
      !  this routine test for this-- if this test fails the error flag is
      !  set.

      !  plasma boundary should never intersect a circle limiter; this is
      !  tested and reported as an error of the plasma 0.95 surface crosses.

      !-------------------------------
      real*8 :: zR0,zZ0
      real*8 :: zRmin_plasma,zRmax_plasma,zZmin_plasma,zZmax_plasma
      real*8 :: zRmin_use,zRmax_use,zZmin_use,zZmax_use,zdist_use
      real*8, parameter :: zrho_bdy = 1.0d0
      integer :: icircs,ilines,ilin,ilocmmx
      real*8, dimension(:), allocatable :: zZl,zRl,zthl

      real*8, parameter :: zsmall = 0.01d0 ! 1cm
      real*8, parameter :: zbig = 100.00d0 ! 100 meters

      integer :: idata(7),ir8size,itype,id,iertmp,id_R,i,j
      real*8, dimension(:), allocatable :: r8data

      integer, parameter :: isrch=300
      real*8 :: rsrch(isrch),zsrch(isrch),zswrk(isrch),zphwrk(isrch)
      real*8 :: zdist(isrch),zdum1(isrch),zdum2(isrch),zdum3(isrch)
      real*8 :: zR0_srch,zZ0_srch,zRmin_srch,zRmax_srch,zZmin_srch,zZmax_srch
      real*8 :: zdistRmin,zdistRmax,zdistZmin,zdistZmax
      real*8 :: zincr
      !-------------------------------

      ier=0
      iwarn=0

      !  set warning flag if something named "__LIMITER" already exists.

      call xplasma_find_item(s,'__LIMITER',id,iertmp, nf_noerr=.TRUE.)
      if(id.ne.0) iwarn=1

      call xplasma_mag_axis(s,ier, zR0, zZ0)
      if(ier.ne.0) then
         ier=3001
         call xplasma_errmsg_append(s, &
              ' ?xplasma_mklim_special: could not get (R0,Z0).')
         return
      endif

      call xplasma_RZminmax(s, zrho_bdy, ier, &
           Rmin=zRmin_plasma, Rmax=zRmax_plasma, &
           Zmin=zZmin_plasma, Zmax=zZmax_plasma)
      if(ier.ne.0) then
         ier=3001
         call xplasma_errmsg_append(s, &
              ' ?xplasma_mklim_special: could not get R & Z min & max at plasma bdy.')
         return
      endif

      zR0_srch=(zRmin_plasma+zRmax_plasma)/2
      zZ0_srch=(zZmin_plasma+zZmax_plasma)/2

      !  check consistency of optional arguments
      if(present(ncircs)) then
         if(ncircs.lt.0) then
            ier=3002
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_mklim_special: ncircs < 0.')
         else if(ncircs.gt.0) then
            if(.not.(present(Rc).and.present(Zc).and.present(rad))) then
               ier=3002
               call xplasma_errmsg_append(s, &
                    ' ?xplasma_mklim_special: if ncircs is present, then')
               call xplasma_errmsg_append(s, &
                    '   Rc, Zc, and rad arrays must also be present.')
            else
               if(min(size(Rc),size(Zc),size(rad)).lt.ncircs) then
                  ier=3002
                  call xplasma_errmsg_append(s, &
                       ' ?xplasma_mklim_special: if ncircs is present, then')
                  call xplasma_errmsg_append(s, &
                       '   Rc, Zc, and rad array sizes must all be >= ncircs.')
               endif
            endif
         endif
      else
         if(present(Rc).or.present(Zc).or.present(rad)) then
            ier=3002
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_mklim_special: if ncircs is absent, then')
            call xplasma_errmsg_append(s, &
                 '  array arguments Rc, Zc, and rad must also all be absent.')
         endif
      endif

      if(present(nlines)) then
         if(nlines.lt.0) then
            ier=3002
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_mklim_special: nlines < 0.')
         else if(nlines.gt.0) then
            if(.not.(present(Rl).and.present(Zl).and.present(thl))) then
               ier=3002
               call xplasma_errmsg_append(s, &
                    ' ?xplasma_mklim_special: if nlines is present, then')
               call xplasma_errmsg_append(s, &
                    '   Rl, Zl, and thl arrays must also be present.')
            else
               if(min(size(Rl),size(Zl),size(thl)).lt.nlines) then
                  ier=3002
                  call xplasma_errmsg_append(s, &
                       ' ?xplasma_mklim_special: if nlines is present, then')
                  call xplasma_errmsg_append(s, &
                       '   Rl, Zl, and thl array sizes must all be >= nlines.')
               endif
            endif
         endif
      else
         if(present(Rl).or.present(Zl).or.present(thl)) then
            ier=3002
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_mklim_special: if nlines is absent, then')
            call xplasma_errmsg_append(s, &
                 '  array arguments Rl, Zl, and thl must also all be absent.')
         endif
      endif
      if(ier.ne.0) return

      !------------------------------------------------------------
      !  set default values; get optional inputs as provided

      !  dist < 0 means: use circles & lines limiters only...

      zdist_use=-1.0d0
      if(present(dist)) zdist_use = dist
      if(zdist_use.lt.0.0d0) zdist_use=-1.0d0

      !  provided "sanity limiters" to guarantee a bounded space
      !  in case user's lines and circles do not enclose plasma...
      !     these are expressed here in terms of R & Z min & max but
      !     will be turned into line limiters below...

      if(present(Rmin)) then
         zRmin_use=Rmin
      else
         zRmin_use=zsmall ! 1cm from machine axis
      endif

      if(present(Rmax)) then
         zRmax_use=Rmax
      else
         zRmax_use=zRmax_plasma + zbig
      endif

      if(present(Zmin)) then
         zZmin_use=Zmin
      else
         zZmin_use=zZmin_plasma - zbig
      endif

      if(present(Zmax)) then
         zZmax_use=Zmax
      else
         zZmax_use=zZmax_plasma + zbig
      endif

      ! circle limiters...
      icircs = 0
      if(present(ncircs)) icircs=ncircs

      ! line limiters
      ilines = 4
      if(present(nlines)) ilines=nlines+4

      allocate(zZl(ilines),zRl(ilines),zthl(ilines))
      if(ilines.gt.4) then
         zZl(1:nlines)=Zl(1:nlines)
         zRl(1:nlines)=Rl(1:nlines)
         zthl(1:nlines)=thl(1:nlines)
      endif

      ! box limiters -- horizontal and vertical lines...

      ilin = ilines-3
      zZl(ilin)=0.0d0
      zRl(ilin)=zRmin_use
      zthl(ilin)=90.0d0

      ilin = ilines-2
      zZl(ilin)=0.0d0
      zRl(ilin)=zRmax_use
      zthl(ilin)=90.0d0

      ilin = ilines-1
      zZl(ilin)=zZmin_use
      zRl(ilin)=0.0d0
      zthl(ilin)=0.0d0

      ilin = ilines
      zZl(ilin)=zZmax_use
      zRl(ilin)=0.0d0
      zthl(ilin)=0.0d0

      !--------------------------------------------
      !  OK -- set up limiter data

      idata(1)=1                ! address of lines data
      idata(2)=1 + 6*ilines     ! address of circles data
      idata(3)=idata(2)+4*icircs           ! address of distance data
      idata(4)=ilines           ! #lines
      idata(5)=icircs           ! #circs, can be zero
      idata(6)=0                ! set =1 for plasma bdy as limiter

      ir8size = 6*ilines + 4*icircs + 1

      ilocmmx = ir8size+1
      idata(7)=ilocmmx

      itype = 101

      allocate(r8data(ir8size+4)); r8data = 0

      r8data(idata(3)) = -1.0d0  ! no distance information yet

      call eqi_tbdy_lines(r8data(idata(1)),ilines,zRl,zZl,zthl,zR0,zZ0)
      if(icircs.gt.0) then
         call eqi_tbdy_circs(r8data(idata(2)),icircs,Rc,Zc,rad,zR0,zZ0)
      endif

      !----------------------------------------------
      !  find approximate min and max of vacuum region--
      !  scan in vicinity of sample plasma; cast a wide net though...

      zincr = zR0_srch/1000

      zRmax_srch=3*zR0_srch
      zRmin_srch=zRmax_srch/isrch

      zZmin_srch=zZ0_srch-3*zR0_srch
      zZmax_srch=zZ0_srch+3*zR0_srch

      do i=1,isrch
         rsrch(i)=((isrch-i)*zRmin_srch + (i-1)*zRmax_srch)/(isrch-1)
         zsrch(i)=((isrch-i)*zZmin_srch + (i-1)*zZmax_srch)/(isrch-1)
         zphwrk(i)=0
      enddo

      !  reuse these variables...

      zRmin_srch=1.0d30
      zRmax_srch=-zRmin_srch

      zZmin_srch=1.0d30
      zZmax_srch=-zZmin_srch

      do j=1,isrch
         zswrk=zsrch(j)  ! one z value for rsrch(1:isrch)
         call eqi_rzdist(isrch,rsrch,zswrk,zphwrk,zdist, &
              itype,idata,r8data, &
              0,zdum1,zdum2,zdum3,ier)
         do i=1,isrch
            if(zdist(i).lt.0.0d0) then
               ! point is inside all limiters...

               if(rsrch(i).lt.zRmin_srch) then
                  zRmin_srch=rsrch(i)
                  zdistRmin=-zdist(i)
               else if(rsrch(i).eq.zRmin_srch) then
                  zdistRmin=max(zdistRmin,-zdist(i))
               endif

               if(rsrch(i).gt.zRmax_srch) then
                  zRmax_srch=rsrch(i)
                  zdistRmax=-zdist(i)
               else if(rsrch(i).eq.zRmax_srch) then
                  zdistRmax=max(zdistRmax,-zdist(i))
               endif

               if(zswrk(i).lt.zZmin_srch) then
                  zZmin_srch=zswrk(i)
                  zdistZmin=-zdist(i)
               else if(zswrk(i).eq.zZmin_srch) then
                  zdistZmin=max(zdistZmin,-zdist(i))
               endif

               if(zswrk(i).gt.zZmax_srch) then
                  zZmax_srch=zswrk(i)
                  zdistZmax=-zdist(i)
               else if(zswrk(i).eq.zZmax_srch) then
                  zdistZmax=max(zdistZmax,-zdist(i))
               endif

            endif
         enddo
      enddo

      r8data(ilocmmx)=max(zRmin_use,(zRmin_srch-zdistRmin-zincr))
      r8data(ilocmmx+1)=min(zRmax_use,(zRmax_srch+zdistRmax+zincr))
      r8data(ilocmmx+2)=max(zZmin_use,(zZmin_srch-zdistZmin-zincr))
      r8data(ilocmmx+3)=min(zZmax_use,(zZmax_srch+zdistZmax+zincr))

      !----------------------------------------------
      !  insert distance information

      r8data(idata(3)) = zdist_use

      !----------------------------------------------
      !  create limiter dataset "black box"...

      call xplasma_author_set(s,xplasma_root,iertmp)

      call xplasma_create_blackbox(s,'__LIMITER',itype,id,ier, &
           iarray=idata, r8array=r8data, &
           label='axisymmetrix limiter data')

      call xplasma_author_clear(s,xplasma_root,iertmp)

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' xplasma_mklim_special: "__LIMITER" create failure.')
      endif

      deallocate(r8data)

      call xplasma_common_ids(s,iertmp, id_R)
      if(id_R.gt.0) then
         call xplasma_limcheck(s,ier)
      endif

    end subroutine xplasma_mklim_special

    !==========================================================
    subroutine xplasma_lim_info(s,ier, &
         itype, npts,rpts,zpts, &
         dist, &
         nlines, Rl, Zl, thl, &
         ncircs, Rc, Zc, rad)

      ! retrieve limiter information

      type(xplasma), pointer :: s
      integer, intent(out) :: ier

      integer, intent(out), optional :: itype
      !  itype=100 means piecewise linear closed contour (axisymmetric limiter)
      !  itype=101 means list of circles & infinite lines

      integer, intent(out), optional :: npts ! #pts in contour (itype=100 only)

      real*8, dimension(:), intent(out), optional :: rpts,zpts
      !  the closed contour; rpts(npts)=rpts(1) and zpts(npts)=zpts(1) (m).

      real*8, intent(out), optional :: dist  ! itype=101 limiters only...
      !  if .ge.0.0d0, the plasma boundary + dist is considered a limiter
      !  if .lt.0.0d0-- no plasma boundary based limiter.

      integer, intent(out), optional :: nlines  ! #lines in itype=101 limiter
      real*8, dimension(:), intent(out), optional :: Rl,Zl,thl
      !  point (m) through which line passes, and line orientation in DEGREES
      !  e.g. 45.0 means up and to the right...

      integer, intent(out), optional :: ncircs ! #circles in itype=101 limiter
      real*8, dimension(:), intent(out), optional :: Rc,Zc,rad
      !  center location and radius of each circle (m)

      !-------------------------------------------------------
      integer :: idlim,ilines,icircs,jtype,iadrl,iadrc,iadist,isize_min
      integer :: iadr,iadz,i
      character*128 zmsg
      integer, dimension(:), pointer :: iwork
      real*8, dimension(:), pointer :: r8work
      real*8 :: zargr,zargz,zarg3
      !-------------------------------------------------------

      if(present(itype)) itype=0

      if(present(npts)) npts=0
      if(present(rpts)) rpts=0
      if(present(zpts)) zpts=0

      if(present(dist)) dist=0

      if(present(nlines)) nlines=0
      if(present(Rl)) Rl=0
      if(present(Zl)) Zl=0
      if(present(thl)) thl=0

      if(present(nlines)) nlines=0
      if(present(Rc)) Rc=0
      if(present(Zc)) Zc=0
      if(present(rad)) rad=0

      call xplasma_find_item(s,'__LIMITER',idlim,ier,nf_noerr=.TRUE.)
      if(ier.ne.0) return

      if(idlim.eq.0) then
         ier=3004          ! no __LIMITER object found.
         return
      endif

      !----------------------------------
      !  get limiter information

      call xplasma_blackBox_retrieve(s,idlim,ier, &
           itype=jtype, ia_ptr=iwork, r8a_ptr=r8work)

      if(present(itype)) itype=jtype

      iadrl =iwork(1)
      ilines=iwork(4)

      iadist=iwork(3)

      iadrc =iwork(2)
      icircs=iwork(5)

      !----------------------------------

      if(jtype.eq.100) then

         if(present(npts)) npts=ilines

         isize_min=ilines
         if(present(rpts)) isize_min=min(isize_min,size(rpts))
         if(present(zpts)) isize_min=min(isize_min,size(zpts))

         if(isize_min.lt.ilines) then
            ier=510
            zmsg=' '
            write(zmsg,*) ' ?xplasm_lim_info:  there are ',ilines,' limiter data points.'
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) '  but the output vector size provided is just ',isize_min,'.'
            call xplasma_errmsg_append(s,zmsg)

         else

            iadr=-5
            iadz=-4
            do i=1,ilines-1
               iadr = iadr + 6
               iadz = iadz + 6
               if(present(rpts)) rpts(i)=r8work(iadr)
               if(present(zpts)) zpts(i)=r8work(iadz)
            enddo
            if(present(rpts)) rpts(ilines)=rpts(1)
            if(present(zpts)) zpts(ilines)=zpts(1)
         endif

      else if(jtype.eq.101) then

         if(present(nlines)) nlines=ilines
         if(present(ncircs)) ncircs=icircs

         if(present(dist)) dist = r8work(iadist)

         isize_min=ilines
         if(present(Rl)) isize_min=min(isize_min,size(Rl))
         if(present(Zl)) isize_min=min(isize_min,size(Zl))
         if(present(thl)) isize_min=min(isize_min,size(thl))

         if(isize_min.lt.ilines) then
            ier=510
            zmsg=' '
            write(zmsg,*) ' ?xplasm_lim_info:  there are ',ilines,' limiter lines.'
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) '  but the output vector size provided is just ',isize_min,'.'
            call xplasma_errmsg_append(s,zmsg)

         else
            iadr=-5
            do i=1,ilines
               iadr=iadr+6
               call linparse(r8work(iadr:iadr+5),zargR,zargZ,zarg3)
               if(present(Rl)) Rl(i)=zargR
               if(present(Zl)) Zl(i)=zargZ
               if(present(thl)) thl(i)=zarg3
            enddo

         endif

         isize_min=icircs
         if(present(Rc)) isize_min=min(isize_min,size(Rc))
         if(present(Zc)) isize_min=min(isize_min,size(Zc))
         if(present(rad)) isize_min=min(isize_min,size(rad))

         if(isize_min.lt.icircs) then
            ier=510
            zmsg=' '
            write(zmsg,*) ' ?xplasm_lim_info:  there are ',icircs,' limiter circles.'
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) '  but the output vector size provided is just ',isize_min,'.'
            call xplasma_errmsg_append(s,zmsg)

         else
            iadr=-3
            do i=1,icircs
               iadr=iadr+4
               call cirparse(r8work(iadr:iadr+3),zargR,zargZ,zarg3)
               if(present(Rc)) Rc(i)=zargR
               if(present(Zc)) Zc(i)=zargZ
               if(present(rad)) rad(i)=zarg3
            enddo

         endif

      endif

      nullify(iwork,r8work)

      contains
        subroutine linparse(zdata,Rpt,Zpt,ang_degrees)
          real*8,intent(in) :: zdata(6)
          real*8,intent(out) :: Rpt,Zpt,ang_degrees ! reference pt & angle

          real*8 zang_rad

          Rpt=zdata(1)
          Zpt=zdata(2)
          zang_rad=atan2(-zdata(3),zdata(4))
          ang_degrees=360.00d0*zang_rad/6.2831853071795862D+00

        end subroutine linparse

        subroutine cirparse(zdata,Rc,Zc,radius)
          real*8,intent(in) :: zdata(4)
          real*8,intent(out) :: Rc,Zc,radius  ! center & radius of circle

          Rc=zdata(1)
          Zc=zdata(2)
          radius=zdata(3)

        end subroutine cirparse

    end subroutine xplasma_lim_info
    !==========================================================

    subroutine xplasma_lim_rzminmax(s,zRmin,zRmax,zZmin,zZmax,ier, &
         itype)

      !  get the Rmin,Rmax,Zmin,Zmax of the vacuum region
      !  enclosed by the mechanical limiter
      !    optional: get a limiter type code also.

      type (xplasma), pointer :: s

      real*8, intent(out) :: zRmin,zRmax,zZmin,zZmax  ! R & Z min & max (m)

      integer, intent(out) :: ier

      integer, intent(out), optional :: itype  ! limiter representation
      ! within xplasma: 100 for piecewise contour; 101 for "circles and
      ! lines" possibly augmented by a fixed distance from the plasma
      ! boundary.

      !------------------------------------------
      integer :: idlim,ilocmmx,jtype
      integer, dimension(:), pointer :: iwork
      real*8, dimension(:), pointer :: r8work
      !------------------------------------------

      zRmin = 0; zRmax = 0; zZmin = 0; zZmax = 0

      call xplasma_find_item(s,'__LIMITER',idlim,ier,nf_noerr=.TRUE.)
      if(ier.ne.0) return

      if(idlim.eq.0) then
         ier=3004          ! no __LIMITER object found.
         return
      endif

      !----------------------------------
      !  get limiter information

      call xplasma_blackBox_retrieve(s,idlim,ier, &
           itype=jtype, ia_ptr=iwork, r8a_ptr=r8work)

      if(present(itype)) itype=jtype

      ilocmmx=iwork(7)

      zRmin=r8work(ilocmmx)
      zRmax=r8work(ilocmmx+1)
      zZmin=r8work(ilocmmx+2)
      zZmax=r8work(ilocmmx+3)

      nullify(iwork,r8work)

    end subroutine xplasma_lim_rzminmax

    !==========================================================

    subroutine xplasma_limcheck(s,ier, rho_check)

      !  check if plasma boundary intersects mechanical limiter as 
      !  specified...

      type (xplasma), pointer :: s
      integer, intent(out) :: ier

      real*8, intent(in), optional :: rho_check

      !----------------------------------
      integer :: idlim  ! id of limiter object (must exist)
      real*8, parameter :: zrhomin = 0.9d0
      real*8, parameter :: zrhomax = 1.0d0
      real*8 :: zrho    ! rho of surface to use for testing
      integer :: id_R,id_Z,id2(2),ith

      integer, parameter :: intest=201
      real*8 :: zrho_test(intest),zth_test(intest),zRZ(intest,2),zdist(intest)
      real*8 :: zphiang(intest),zdum1(intest),zdum2(intest),zdum3(intest)
      real*8 :: zdist_max
      integer :: idist_max

      real*8, parameter :: C2PI = 6.2831853071795862D+00

      integer :: itype
      integer, dimension(:), pointer :: iwork
      real*8, dimension(:), pointer :: r8work
      character*128 zmsg
      !----------------------------------

      ier = 0

      zrho = xplasma_rho_limcheck   ! the default
      if(present(rho_check)) then
         zrho = max(zrhomin,min(zrhomax,rho_check))
      endif

      call xplasma_find_item(s,'__LIMITER',idlim,ier,nf_noerr=.TRUE.)
      if(ier.ne.0) return

      if(idlim.eq.0) then
         ier=3004          ! no __LIMITER object found.
         return
      endif

      call xplasma_common_ids(s,ier, id_R=id_R,id_Z=id_Z)
      if(ier.ne.0) return

      if(min(id_R,id_Z).eq.0) then
         ier=3005          ! no flux surfaces found.
         return
      endif

      !----------------------------------
      !  OK the needed data is available.  Get the specified flux surface...
      !  Direct R & Z profile evaluation is done, although a coordinate
      !  transform routine could just as well have been used...

      id2(1)=id_R
      id2(2)=id_Z

      do ith=1,intest
         zrho_test(ith)=zrho
         zth_test(ith)=(ith-1)*C2PI/(intest-1)
         zphiang(ith)=0
      enddo

      call xplasma_RZeval_2d(s,id2, &
           xplasma_theta_coord,zth_test, xplasma_rho_coord,zrho_test, &
           zRZ,ier)

      if(ier.ne.0) return

      !----------------------------------
      !  get limiter information

      call xplasma_blackBox_retrieve(s,idlim,ier, &
           itype=itype, ia_ptr=iwork, r8a_ptr=r8work)

      if(ier.ne.0) then
         nullify(iwork,r8work)
         return
      endif

      !  treat as pure mechanical limiter, i.e. ignore plasma boundary
      !  based limiter specification for now

      !  this routine checks mechanical limiter components only, in 
      !  calculating the distance...

      call eqi_rzdist(intest,zRZ(1:intest,1),zRZ(1:intest,2),zphiang, &
           zdist,itype,iwork,r8work, &
           0,zdum1,zdum2,zdum3,ier)

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_limcheck: error in eqi_rzdist call!')
         nullify(iwork,r8work)
         return
      endif

      idist_max=0
      zdist_max=-1.0d30

      do ith=1,intest
         if(zdist(ith).gt.zdist_max) then
            zdist_max=zdist(ith)
            idist_max=ith
         endif
      enddo

      if(zdist_max.gt.0.0) then
         ier=3003
         zmsg=' '
         write(zmsg,*) '  at rho=',zrho,' theta=',zth_test(idist_max)
         call xplasma_errmsg_append(s,zmsg)
         zmsg=' '
         write(zmsg,*) '  location R=',zRZ(idist_max,1),' Z=',zRZ(idist_max,2)
         call xplasma_errmsg_append(s,zmsg)
         zmsg=' '
         write(zmsg,*) '  plasma interior point outside limiter by ', &
              zdist_max,' meters.'
         call xplasma_errmsg_append(s,zmsg)
      endif

      nullify(iwork,r8work)

    end subroutine xplasma_limcheck

    subroutine xplasma_bdycon(s,rbdy,zbdy,ier)

      !  return a closed contour around the plasma boundary, similar to
      !  xplasma_limcon...

      type (xplasma), pointer :: s

      real*8, dimension(:), intent(out) :: rbdy,zbdy  ! bdy contour
      !  the sizes of these vectors must be the same.

      integer, intent(out) :: ier

      !-------------
      integer :: id_thetax,id_R,ith,inth,isize,indx
      real*8 :: zindx,zfac
      real*8, dimension(:),allocatable :: thetax,rho_use,th_use
      real*8, parameter :: ONE=1.0d0
      real*8, parameter :: ZERO=0.0d0
      !-------------

      call xplasma_find_item(s,'__THETAX',id_thetax,ier,nf_noerr=.TRUE.)
      if(ier.ne.0) return

      isize=size(rbdy)
      if(isize.ne.size(zbdy)) then
         ier=510
         call xplasma_errmsg_append(s, &
              ' ?xplasma_bdycon: size of "rbdy" and "zbdy" vectors must match.')
         return
      endif

      if(id_thetax.eq.0) then
         call xplasma_find_item(s,'R',id_R,ier)
         if(ier.ne.0) return

         call xplasma_prof_info(s,id_R,ier, gridId1=id_thetax)
         if(ier.ne.0) return
      endif

      call xplasma_grid_size(s,id_thetax,inth,ier)
      if(ier.ne.0) return

      allocate(thetax(inth),rho_use(isize),th_use(isize))

      call xplasma_grid(s,id_thetax,thetax,ier)
      if(ier.ne.0) return

      !  here the thetax grid is used, to take advantage of possible
      !  uneven spacing that might compensate for distortion of theta
      !  with respect to arclength along the plasma boundary

      rho_use=ONE
      th_use(1)=thetax(1)
      th_use(isize)=thetax(inth)
      do ith=2,isize-1
         zindx=ONE+(ith-1)*(inth-ONE)/(isize-1)
         indx=zindx
         indx=max(1,min(inth-1,indx))
         zfac=zindx-indx
         zfac=max(ZERO,min(ONE,zfac))
         th_use(ith)=thetax(indx)+zfac*(thetax(indx+1)-thetax(indx))
      enddo

      call xplasma_ctrans(s,.TRUE.,ier, &
           rho_in=rho_use, theta_in=th_use, r_out=rbdy, z_out=zbdy)

      deallocate(thetax,rho_use,th_use)

    end subroutine xplasma_bdycon

    subroutine xplasma_limcon(s,rlim,zlim,inum_got,ier, &
         tol,maptype)

      !  return a closed contour of length size(rlim)=size(zlim) or less
      !  describing the limiter.

      !  if tol.le.0.0d0 -- the length of the description is size(rlim) exactly
      !  if tol.gt.0.0d0 -- length may be shortened to remove colinear points.
      !      also, if the limiter was originally given in contour form, and,
      !      size(rlim) is sufficient, the original specification will be
      !      returned.

      type (xplasma), pointer :: s

      real*8, dimension(:), intent(out) :: rlim,zlim  ! limiter contour
      !  the sizes of these vectors must be the same.
      !  unless the original specification is being returned, the code
      !  requires a minimum size of 30 points.

      integer, intent(out) :: inum_got     ! contour length returned
      !  if inum_got.lt.size(rlim), rlim(inum_got+1:size(rlim))=0 on exit.
      !  similarly for zlim.

      integer, intent(out) :: ier  ! status code, 0=OK

      real*8, intent(in), optional :: tol  ! tolerance specification
      !  if zero or negative, no shortening of the limiter description is
      !  attempted.  If .gt. 0.0, adjacent segments that are colinear to within
      !  tol*[Rmax] are replaced by a single segment.  DEFAULT: s%bdytol

      integer, intent(in), optional :: maptype  ! mapping option
      !  1=slowest but most exact; 2=intermediate; 3=bilinear interpolation
      !  2 or 3 recommended.  DEFAULT: 2

      !----------------------------------
      integer :: idlim  ! id of limiter object (must exist)

      real*8, parameter :: C2PI = 6.2831853071795862D+00

      logical :: icopy
      integer :: itype,ipts,imaptype,isize,i,ibase
      integer, dimension(:), pointer :: iwork
      real*8, dimension(:), pointer :: r8work
      real*8 :: ztol,zbdy_tol,zrmax
      character*128 zmsg

      real*8, dimension(:), allocatable :: ra,za,da,phia, rb,zb,db,phib
      real*8, dimension(:), allocatable :: rx,zx,phix
      real*8 :: lperp(2),lnorm,dtest
      integer, dimension(:), allocatable :: istat
      !----------------------------------

      ier = 0

      inum_got = 0
      rlim = 0
      zlim = 0

      call xplasma_find_item(s,'__LIMITER',idlim,ier,nf_noerr=.TRUE.)
      if(ier.ne.0) return

      if(idlim.eq.0) then
         ier=3004          ! no __LIMITER object found.
         return
      endif

      if(size(rlim).ne.size(zlim)) then
         ier=510
         call xplasma_errmsg_append(s,' in xplasma_limcon:')
         call xplasma_errmsg_append(s, &
              ' sizes of "Rlim" and "Zlim" vector arguments do not match.')
         return
      endif

      !----------------------------------
      !  handle optional arguments; set defaults...

      call xs_max_maptype(s,imaptype,ier)
      if(ier.ne.0) return

      if(present(maptype)) then
         imaptype=max(1,min(imaptype,maptype))
      else
         imaptype=2
      endif

      call xplasma_global_info(s,ier,bdytol=zbdy_tol)
      ztol=zbdy_tol
      if(present(tol)) ztol=tol
      ztol=min(1.0d-5,ztol)

      !----------------------------------
      !  get limiter information

      call xplasma_blackBox_retrieve(s,idlim,ier, &
           itype=itype, ia_ptr=iwork, r8a_ptr=r8work)

      if(ier.ne.0) then
         nullify(iwork,r8work)
         return
      endif

      !----------------------------------
      !  see if we can just copy the original contour...

      icopy=.False.
      if(itype.eq.100) then
         ipts=iwork(4)
         if(size(rlim).eq.ipts) then
            icopy=.TRUE.
         else if((ztol.gt.0.0d0).and.(size(rlim).ge.ipts)) then
            icopy=.TRUE.
         endif
      endif

      if(icopy) then

         inum_got = ipts
         call xplasma_lim_info(s,ier, rpts=rlim(1:ipts), zpts=zlim(1:ipts))
         return

      endif

      !----------------------------------
      !  OK -- need to probe for limiter location...

      if(size(rlim).lt.30) then
         ier=510
         call xplasma_errmsg_append(s,' in xplasma_limcon:')
         call xplasma_errmsg_append(s, &
              ' sizes of "Rlim" and "Zlim" vector arguments too small, need at least 30 pts.')
         return
      endif

      !----------------------------------

      isize=size(rlim)

      allocate(ra(isize),za(isize),da(isize),phia(isize))
      allocate(rb(isize),zb(isize),db(isize),phib(isize))
      allocate(rx(isize),zx(isize),phix(isize),istat(isize))

      sp => s

      do
         call eqi_bprobe(isize,2,imaptype,ra,za,da,rb,zb,db,ier)
         if(ier.ne.0) exit

         phia=0
         phib=0

         call xplasma_fxcep(s,2, &
              ra,za,phia, rb,zb,phib, rx,zx,phix, istat,ier, &
              maptype=imaptype, tol=zbdy_tol, da_input=da, db_input=db)
         if(ier.ne.0) exit

         if(ztol.le.0.0d0) then

            inum_got=isize
            rlim=rx
            zlim=zx
            ier=0
            exit

         endif

         ! remove colinear points...

         ibase=1
         inum_got=1

         zrmax=maxval(rx)

         rlim(1)=rx(1)
         zlim(1)=zx(1)

         do i=3,isize

            ! colinearity test of point (i-1) w.r.t. (ibase) and (i)

            lperp(1)=-(zx(ibase+1)-zx(ibase))
            lperp(2)=(rx(ibase+1)-rx(ibase))
            lnorm=sqrt(lperp(1)**2+lperp(2)**2)
            lperp=lperp/lnorm

            dtest=lperp(1)*(rx(i)-rx(ibase))+ lperp(2)*(zx(i)-zx(ibase))

            if(abs(dtest).gt.ztol*zrmax) then

               ! not sufficiently colinear; add to sequence & update base

               ibase = i-1
               inum_got = inum_got + 1
               rlim(inum_got)=rx(ibase)
               zlim(inum_got)=zx(ibase)

            endif
         enddo

         !  add in last point; matches 1st point

         inum_got = inum_got + 1
         rlim(inum_got) = rlim(1)
         zlim(inum_got) = zlim(1)
         ier = 0
         exit
      enddo

      deallocate(ra,za,phia,rb,zb,phib,da,db,rx,zx,phix)
      nullify(iwork,r8work)

    end subroutine xplasma_limcon

    subroutine xplasma_lim_dist1(s,r,z,d,ier, maptype, rlim1,zlim1)

      !  scalar version of xplasma_lim_distn -- axisymmetry assumed.

      type (xplasma), pointer :: s
      real*8, intent(in) :: R,Z      ! input location (m)
      real*8, intent(out) :: d       ! output distance from nearest point
      !                                on limiter or wall (m).

      integer, intent(out) :: ier

      integer, intent(in), optional :: maptype  ! see xplasma_lim_distn
      real*8, intent(out), optional :: rlim1,zlim1  ! nearest limiter point

      !-----------------------

      real*8 :: rvec(1),zvec(1),dist(1),rlim(1),zlim(1)

      rvec = R
      zvec = Z

      call xplasma_lim_distn(s,rvec,zvec,dist,ier, &
           maptype=maptype, rlim=rlim, zlim=zlim)

      d = dist(1)
      if(present(rlim1)) rlim1=rlim(1)
      if(present(zlim1)) zlim1=zlim(1)

    end subroutine xplasma_lim_dist1

    subroutine xplasma_lim_distn(s,rvec,zvec,dist,ier, &
         maptype,rlim,zlim)

      !  find distance from each element of a vector of (R,Z) pairs
      !  to the nearest point on a limiter.

      !  optionally return the location of that point on the limiter
      !  in (R,Z).

      !  the routine assumes axisymmetry.  Optional phi arguments may
      !  be added someday.

      type (xplasma), pointer :: s
      real*8, intent(in), dimension(:) :: rvec,zvec  ! (R,Z) input vector

      real*8, intent(out), dimension(:) :: dist      ! distance values returned
      integer, intent(out) :: ier

      integer, intent(in), optional :: maptype       ! distance map option
      ! if it is required to compute the distance from the (R,Z) point to 
      ! the plasma boundary, this specifies the option (for xplasma_bdfind).
      !   Default is the circle-fit option (2).  (1) is slower but more
      ! precise; (3) is less accurate but would be faster if a bilinear
      ! distance map has already been computed.

      real*8, intent(out), dimension(:), optional :: rlim,zlim  ! locations
      !                      of nearest contact points on limiter, for each
      !                      input location.
      !-----------------------------------------------------------
      real*8 :: phitmp1(size(rvec)),phitmp2(size(rvec))
      real*8 :: rhotmp(size(rvec)),thtmp(size(rvec)),dtmp(size(rvec))
      integer :: idlim,id_R,id_Z,itype,isize,isize_min,isize_max,ictrz,ipdist
      real*8 :: plasma_lim,delR,delZ,dnorm
      integer :: imap_type,j,jvec,ids(4),idchis(4)
      integer, dimension(:), pointer :: iwork
      real*8, dimension(:), pointer :: r8work
      real*8, dimension(:,:), allocatable :: zrzbuf
      !-----------------------------------------------------------

      call xplasma_find_item(s,'__LIMITER',idlim,ier,nf_noerr=.TRUE.)
      if(ier.ne.0) return

      if(idlim.eq.0) then
         ier=3004          ! no __LIMITER object found.
         return
      endif

      call xplasma_common_ids(s,ier, id_R=id_R,id_Z=id_Z)
      if(ier.ne.0) return

      if(min(id_R,id_Z).eq.0) then
         ier=3005          ! no flux surfaces found.
         return
      endif

      !----------------------------------
      !  OK data is available...

      isize=size(rvec)
      isize_min=isize
      isize_max=isize

      isize_min=min(isize_min,size(zvec))
      isize_min=min(isize_min,size(dist))

      isize_max=max(isize_max,size(zvec))
      isize_max=max(isize_max,size(dist))

      ictrz=0

      if(present(rlim)) then
         ictrz=ictrz+1
         isize_min=min(isize_min,size(rlim))
         isize_max=max(isize_max,size(rlim))
      endif

      if(present(zlim)) then
         ictrz=ictrz+1
         isize_min=min(isize_min,size(zlim))
         isize_max=max(isize_max,size(zlim))
      endif

      if((isize_min.ne.isize_max).or.(ictrz.eq.1)) then
         ier=510
         if(isize_min.ne.isize_max) call xplasma_errmsg_append(s, &
              ' input/output vector sizes not consistent in xplasma_lim_distance call.')
         if(ictrz.eq.1) call xplasma_errmsg_append(s, &
              ' if one of the optional vectors {rlim,zlim} are present, both must be present.')
         return
      endif

      !----------------------------------
      !  get limiter information

      call xplasma_blackBox_retrieve(s,idlim,ier, &
           itype=itype, ia_ptr=iwork, r8a_ptr=r8work)

      if(ier.ne.0) then
         nullify(iwork,r8work)
         return
      endif

      !  first: treat as pure mechanical limiter, i.e. ignore plasma boundary
      !  based limiter specification for now

      !  this routine checks mechanical limiter components only, in 
      !  calculating the distance...

      phitmp1=0
      if(ictrz.gt.0) then
         call eqi_rzdist(isize,rvec,zvec,phitmp1, &
              dist,itype,iwork,r8work, &
              1,rlim,zlim,phitmp2,ier)
      else
         call eqi_rzdist(isize,rvec,zvec,phitmp1, &
              dist,itype,iwork,r8work, &
              0,phitmp2,phitmp2,phitmp2,ier)
      endif

      !  see if a limiter based on plasma boundary location is specified

      plasma_lim = -1.0d0

      if(itype.eq.101) then
         plasma_lim = r8work(iwork(3))
      endif

      if(plasma_lim.lt.0.0d0) return

      !------------------------
      !  a plasma boundary based limiter IS specified...

      imap_type=2
      if(present(maptype)) imap_type=max(1,min(3,maptype))
      if(imap_type.eq.3) call xs_max_maptype(s,imap_type,ier)

      call xplasma_bdfind(s,isize,rvec,zvec,ier, &
           maptype=maptype, theta_out=thtmp, dist=dtmp)

      dtmp = dtmp - plasma_lim  ! distance to limiter, not bdy...

      ! points j for which dtmp(j) > dist(j) [from mechanical limiters]
      ! are situated such that the plasma based limiter is the one closest
      ! to the plasma, i.e. the one to count.

      if(ictrz.eq.0) then
         do j=1,isize
            if(dtmp(j).gt.dist(j)) then
               dist(j)=dtmp(j)
            endif
         enddo
         jvec=0

      else
         ! rlim & zlim are requested
         ! compress the vector... need to find the (R,Z) location which
         ! is "plasma_lim" out on a normal vector from the plasma boundary...

         jvec=0
         do j=1,isize

            if(dtmp(j).gt.dist(j)) then
               ! defer update of dist(j) to loop below...
               jvec=jvec+1
               rhotmp(jvec)=1.0d0
               thtmp(jvec)=thtmp(j)
            endif
         enddo

      endif

      if(jvec.eq.0) return

      ids(1)=id_R
      ids(2)=id_Z
      ids(3)=id_R
      ids(4)=id_Z
      idchis(1:2)=0
      idchis(3:4)=1

      allocate(zrzbuf(jvec,4))

      call xplasma_eval_prof(s,ids, &
           xplasma_theta_coord,thtmp(1:jvec),xplasma_rho_coord,rhotmp(1:jvec),&
           zrzbuf,ier, &
           ideriv1s=idchis)

      jvec=0
      do j=1,isize
         if(dtmp(j).gt.dist(j)) then
            jvec=jvec+1
            dnorm=sqrt(zrzbuf(jvec,3)**2+zrzbuf(jvec,4)**2)
            delr=plasma_lim*zrzbuf(jvec,4)/dnorm
            delz=-plasma_lim*zrzbuf(jvec,3)/dnorm

            dist(j)=dtmp(j)
            rlim(j)=zrzbuf(jvec,1)+delr
            zlim(j)=zrzbuf(jvec,2)+delz
         endif
      enddo

      nullify(iwork,r8work)

    end subroutine xplasma_lim_distn
      
    subroutine xplasma_fxcep(s,ibdy_type, &
         ra,za,phia, rb,zb,phib, rx,zx,phix, istat, ier, &
         maptype, tol, rhosurf, da_input, db_input, da_output, db_output)

      ! given a set of segments, find their intersection with either
      ! the plasma boundary or the limiters

      type (xplasma), pointer :: s

      integer, intent(in) :: ibdy_type  ! =1: plasma bdy; =2: limiter

      real*8, dimension(:), intent(in) :: ra,za,phia  ! endpoints A
      real*8, dimension(:), intent(in) :: rb,zb,phib  ! endpoints B

      real*8, dimension(:), intent(out) :: rx,zx,phix ! intersections found
      integer, dimension(:), intent(out) :: istat     ! status summary

      !  istat(j)=2 -- both endpoints outside bdy
      !  istat(j)=1 -- 1st endpoint outside, 2nd endpoint inside
      !  istat(j)=-1 -- 1st endpoint inside, 2nd endpoint outside
      !  istat(j)=-2 -- both endpoints inside bdy
      !    (rx(j),zx(j)) contains bdy intercept if abs(istat(j))=1.

      integer, intent(out) :: ier   ! status code (0=OK)

      integer, intent(in), optional :: maptype  ! control for type of map
      !  to use when searching for intercept
      !  1=very accurate & slow, 2=quite accurate & faster,
      !  3=bilinear map, very fast, accurate if sufficient resolution provided
      !      (usually 2 or 3 are sufficient).  DEFAULT=2

      real*8, intent(in), optional :: tol       ! search tolerance (relative)
      !      default = s%bdytol

      real*8, intent(in), optional :: rhosurf   ! for ibdy_type=1, set this
      !      to a value other than 1, to search for a surface other than the
      !      plasma boundary

      real*8, intent(in), dimension(:), optional :: da_input,db_input
      ! distances at end points-- can be provided if already known.

      real*8, intent(out), dimension(:), optional :: da_output,db_output
      ! distances at end points-- optional output.

      !-------------------------------------------
      real*8 :: zbdy_tol,ztol
      integer :: imap_type,isize,isize_min,isize_max,iin
      real*8, dimension(:), allocatable :: zda,zdb
      real*8 :: zbdy_info
      !-------------------------------------------

      call xplasma_global_info(s,ier, bdytol=zbdy_tol)
      if(ier.ne.0) return

      ztol=zbdy_tol
      if(present(tol)) ztol=min(1.0d-4,tol)

      call xs_max_maptype(s,imap_type,ier)
      if(ier.ne.0) return

      if(present(maptype)) then
         imap_type=max(1,min(imap_type,maptype))
      else
         imap_type=2
      endif

      !-------------
      if(ibdy_type.eq.1) then
         if(present(rhosurf)) then
            if(rhosurf.le.0.0d0) then
               ier=9999
               call xplasma_errmsg_append(s, &
                    ' xplasma_sol:xplasma_fxcep, optional argument "rhosurf" if present must be >0.')
               return
            endif
            zbdy_info=rhosurf  ! any surface
         else
            zbdy_info=1.0d0  ! boundary surface
         endif
      else
         zbdy_info = -1.0d0  ! i.e. for ibdy_type=2 (limiter)
      endif
      !-------------
      !  check vector sizes-- all must match

      isize=size(ra)
      isize_min=isize
      isize_max=isize

      call schk(za); call schk(phia)
      call schk(rb); call schk(zb); call schk(phib)
      call schk(rx); call schk(zx); call schk(phix)
      isize_min=min(isize_min,size(istat))
      isize_max=max(isize_max,size(istat))

      if(present(da_output)) call schk(da_output)
      if(present(db_output)) call schk(db_output)

      iin=0
      if(present(da_input)) then
         call schk(da_input)
         iin=iin+1
      endif
      if(present(db_input)) then
         call schk(db_input)
         iin=iin+1
      endif

      if(isize_min.lt.isize_max) then
         ier=510
         call xplasma_errmsg_append(s, &
              '  vector sizes inconsistent in call to xplama_fxcep')
      endif

      if(iin.eq.1) then
         ier=510
         call xplasma_errmsg_append(s, &
              '  in either of {da_input,db_input} are present, *both* must be present.')
      endif

      if(ier.ne.0) return

      !  OK.................

      allocate(zda(isize),zdb(isize))
      if(iin.gt.0) then

         iin=1
         zda = da_input
         zdb = db_input

      endif

      sp => s

      call eqi_fxcep(isize,zbdy_info,imap_type,iin, &
           ra,za,phia,zda, rb,zb,phib,zdb, ztol, &
           rx,zx,phix,istat, ier)

      if(ier.ne.0) then

         if(present(da_output)) da_output = zda
         if(present(db_output)) db_output = zdb

      endif

      contains

        subroutine schk(arr)
          real*8, dimension(:) :: arr

          isize_min=min(isize_min,size(arr))
          isize_max=max(isize_max,size(arr))
        end subroutine schk

    end subroutine xplasma_fxcep

    subroutine xs_max_maptype(s,maptype,ier)
      !  ** PRIVATE **

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

    end subroutine xs_max_maptype

end module xplasma_sol
