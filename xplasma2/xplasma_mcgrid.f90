module xplasma_mcgrid

  !  xplasma module for "MCgrid" objects-- these are objects defined
  !  over the irregular MC grid of NUBEAM, etc.  At the xplasma kernel
  !  level these are implemented as blackbox objects.

  !  Comments from original f77 implementation:
  !  ------------------------------------------
  !  module for implementation of a set of "MC" grids within the xplasma
  !  environment.  These are semi-irregular grids that partition an
  !  (x,theta,phi) space into zones of roughly equal volume, by setting
  !  up a series of "zone rows" in the x direction, with a variable number
  !  of "theta zones" in each zone row, proportional to the value of x
  !
  !  there are two variants:  updown symmetric and updown asymmetric.
  !    NTH0 = no. of theta zones from 0 to pi in the zone row touching
  !    the magnetic axis (this is HALF the number of zones in the first
  !    row, for updown asymmetric cases).
  !
  !            0 to pi                               0 to 2pi
  !            updown symmetric:                     asymmetric:
  !  1st x row:  NTH0 zones                            2*NTH0 zones
  !  2nd x row:  2*NTH0 zones                          2*2*NTH0 zones
  !  3rd x row:  3*NTH0 zones                          3*2*NTH0 zones
  !                ...                                   ...
  !  Nth x row: N*NTH0 zones                           N*2*NTH0 zones
  !
  !     (the above figures are for axisymmetry, nphi=1; if there
  !     are more than one toroidal (phi) zones, all numbers are
  !     multiplied by the factor nphi)
  !
  !  although x=1 is the normal boundary of the space, the grid is
  !  allowed to extend beyond x=1 into xplasma's "extrapolated" polar
  !  region.  But, zone volumes are only evaluated for the region
  !  inside the boundary.

  use xplasma_obj
  use xplasma_flxint, only: xplasma_2d_zonint
  implicit NONE

  private

  public :: xplasma_mcgrid_find
  public :: xplasma_mcgrid_findProfs
  public :: xplasma_mcgrid_define
  public :: xplasma_mcgrid_info
  public :: xplasma_mcgrid_profInfo
  public :: xplasma_mcgrid_volumes
  public :: xplasma_mcgrid_putobj
  public :: xplasma_mcgrid_getobj

  real*8, parameter :: ONE = 1.0d0
  real*8, parameter :: ZERO= 0.0d0
  real*8, parameter :: EPS7= 1.0d-7
  real*8, parameter :: CPI=  3.1415926535897931D+00
  real*8, parameter :: C2PI= 6.2831853071795862D+00
  integer, parameter, public :: xplasma_bbgtype=17  ! MC grid blackbox type
  integer, parameter, public :: xplasma_bbftype=18  ! MC profile blackbox type

  contains

    subroutine xplasma_mcgrid_find(s,ierr, author_only, &
         id_mcgrid1,num_mcgrids,id_mcgrids)

      !  find MCgrid grid definitions that are available in the current
      !  xplasma-- usually there will be only one we think.

      type(xplasma), pointer :: s
      integer, intent(out) :: ierr          ! status code returned (0=OK)

      character*(*), intent(in), optional :: author_only
      !  restrict search to MCgrid objects written by this author or code

      integer, intent(out), optional :: id_mcgrid1  ! MCgrid id returned, if...
      !  if exactly one exists.  If none exist, id_mcgrid1=0 is returned; if
      !  >1 exists,id_mcgrid1=-<the number of MCgrids> is returned.

      integer, intent(out), optional :: num_mcgrids ! number of MCgrids found.

      integer, dimension(:), optional :: id_mcgrids ! MCgrid ids returned.

      !-----------------------
      integer :: inum,ids(1)
      !-----------------------

      ierr=0
      if(present(id_mcgrid1)) id_mcgrid1=0
      if(present(num_mcgrids)) num_mcgrids=0
      if(present(id_mcgrids)) id_mcgrids=0

      call xplasma_find_blkbxs(s,ierr, author_only, itype=xplasma_bbgtype, &
           num_blkbxs=inum)

      if(ierr.ne.0) return
      if(inum.eq.0) return

      if(present(num_mcgrids)) num_mcgrids=inum

      if(present(id_mcgrids)) then
         call xplasma_find_blkbxs(s,ierr, author_only, itype=xplasma_bbgtype, &
              id_blkbxs=id_mcgrids)

         if(present(id_mcgrid1)) then
            if(inum.eq.1) then
               id_mcgrid1=id_mcgrids(1)
            else
               id_mcgrid1=-inum
            endif
         endif

      else
         if(present(id_mcgrid1)) then
            if(inum.eq.1) then
               call xplasma_find_blkbxs(s,ierr, author_only, itype=xplasma_bbgtype, &
                    id_blkbxs=ids)
               id_mcgrid1=ids(1)
            else
               id_mcgrid1=-inum
            endif
         endif
      endif

    end subroutine xplasma_mcgrid_find

    subroutine xplasma_mcgrid_findProfs(s,id_mcgrid,ierr, &
         author_only, num_profs, id_profs)

      ! find profiles defined over given MCgrid

      type(xplasma), pointer :: s
      integer, intent(in) :: id_mcgrid      ! MC grid id
      integer, intent(out) :: ierr          ! status code returned (0=OK)

      character*(*), intent(in), optional :: author_only
      !  restrict search to MCgrid profiles written by this author or code

      integer, intent(out), optional :: num_profs ! number of profiles found.

      integer, dimension(:), optional :: id_profs ! profile ids returned.

      !-----------------------
      integer :: isize,inum,ii,jj,iertmp
      integer, dimension(:), allocatable :: ids
      integer, dimension(:), pointer :: iadata
      !-----------------------

      ierr=0
      if(present(num_profs)) num_profs=0
      if(present(id_profs)) then
         id_profs=0
         isize=size(id_profs)
      endif

      call xplasma_find_blkbxs(s,ierr,author_only, itype=xplasma_bbftype, &
           num_blkbxs=inum)

      if(ierr.ne.0) return
      if(inum.eq.0) return

      allocate(ids(inum))

      call xplasma_find_blkbxs(s,iertmp,author_only, itype=xplasma_bbftype, &
           id_blkbxs=ids)

      ii=0
      do jj=1,inum

         call xplasma_blackBox_retrieve(s,ids(jj),iertmp, ia_ptr=iadata)

         if(iadata(1).ne.id_mcgrid) cycle

         ii=ii+1
         if(present(id_profs)) then
            if(isize.lt.ii) then
               ierr=61
            else
               id_profs(ii)=ids(jj)
            endif
         endif

      enddo

      if(present(num_profs)) num_profs = ii

    end subroutine xplasma_mcgrid_findProfs

    subroutine xplasma_mcgrid_define(s,mcname,inth0,inrow,id,ierr, &
         nphi,udsym,label,nrow_ext_in,nrow_ext_out)

      type(xplasma), pointer :: s
      character*(*), intent(in) :: mcname   ! name of MC grid

      integer, intent(in) :: inth0          ! no. of zones covering [0,pi]
                                            ! in 1st row at axis

      integer, intent(in) :: inrow          ! no. of rows covering plasma

      integer, intent(out) :: id            ! MC grid id returned (0 if error)
      integer, intent(out) :: ierr          ! status code returned (0=OK)

      integer, intent(in), optional :: nphi ! number of Phi=const zones
      !  default:1 (axisymmetry)

      logical, intent(in), optional :: udsym   ! updown symmetry flag
      !  default:.FALSE. (do not assume updown symmetric MHD equilibrium).

      character*(*), intent(in), optional :: label  ! optional label for MCgrid

      integer, intent(in), optional :: nrow_ext_in   ! no. of external rows
      !  requested (default: 0, meaning: use internal calculation).  If non-
      !  zero value is given, inrow >= nrow_ext_in >= 3 will be enforced.

      integer, intent(out), optional :: nrow_ext_out ! actual no. of external
      !  zone rows created.

      !---------------------------------------------------
      !  define a 2d irregular "MC" polar grid
      !  fewer zones near the axis, more zones towards the edge, each
      !  zone has roughly the same volume...
      !-----------------------------------------
      
      integer :: inznbmr,inznbmri,inphi,idum,iudsym,idi,ir,ifac,ifbznsi,isize
      character*120 zmsg

      integer, dimension(:), allocatable :: idata,nthzns,nthzsum
      real*8, dimension(:), allocatable :: zx,zvols
      real*8 :: zxmin
      logical :: axisymm

      !-----------------------------------------
 
      id = 0
      call xplasma_global_info(s,ierr,axisymm=axisymm)
      if(ierr.ne.0) return

      if(.not.axisymm) then
         ierr=107
         call xplasma_errmsg_append(s,' ?xplasma_mcgrid_define: axisymmetry still required.')
         return
      endif

      ierr = 0
      if(inrow.le.0) then
         zmsg=' '
         write(zmsg,*) ' ?xplasma mcgrid_define:  inrow=',inrow
         call xplasma_errmsg_append(s,trim(zmsg))
         zmsg=' '
         write(zmsg,*) '  inrow.gt.0 required!'
         call xplasma_errmsg_append(s,trim(zmsg))
         ierr = 9999
      endif

      inznbmri = inrow
      idum=0
      if(present(nrow_ext_in)) idum=max(3,min(inznbmri,nrow_ext_in))
      if(idum.eq.0) idum = max(3,inznbmri/3)

      inznbmr = inznbmri+idum
      if(present(nrow_ext_out)) nrow_ext_out = idum
 
      if(inth0.le.0) then
         zmsg=' '
         write(zmsg,*) ' ?xplasma_mcgrid_define:  inth0=',inth0
         call xplasma_errmsg_append(s,trim(zmsg))
         zmsg=' '
         write(zmsg,*) '  inth0.gt.0 required!'
         call xplasma_errmsg_append(s,trim(zmsg))
         ierr = 9999
      endif

      inphi=1
      if(present(nphi)) inphi=max(1,inphi)
 
      if(max(inznbmr,inth0,inphi).gt.1000) then
         zmsg=' '
         write(zmsg,*) ' ?xplasma_mcgrid_define: sanity check:'
         call xplasma_errmsg_append(s,trim(zmsg))
         zmsg=' '
         write(zmsg,*) '  max(inznbmr,inth0,inphi) = ',max(inznbmr,inth0,inphi)
         call xplasma_errmsg_append(s,trim(zmsg))
         zmsg=' '
         write(zmsg,*) '  inznbmr = ',inznbmr
         call xplasma_errmsg_append(s,trim(zmsg))
         zmsg=' '
         write(zmsg,*) '  inth0 = ',inth0
         call xplasma_errmsg_append(s,trim(zmsg))
         zmsg=' '
         write(zmsg,*) '  inphi = ',inphi
         call xplasma_errmsg_append(s,trim(zmsg))
         ierr = 9999
      endif
 
      iudsym=2
      if(present(udsym)) then
         if(udsym) then
            iudsym=1
         else
            iudsym=2
         endif
      endif

      if (ierr.ne.0) return

      !--------------------------------------------
      ! error checks complete.

      ! generate evenly spaced extrapolated x grid
 
      allocate(zx(0:inznbmr))
      zx(0)=ZERO
      zx(inznbmri)=ONE
      zxmin=EPS7

      do ir=1,inznbmri-1
         zx(ir)=((inznbmri-ir)*zx(0) + ir*zx(inznbmri))/(inznbmri)
      enddo
      do ir=inznbmri+1,inznbmr
         zx(ir)=zx(ir-1)+(zx(inznbmri)-zx(inznbmri-1))
      enddo

      ! generate integer array for blackbox containing necessary information
      ! to characterize the grid + some addl info for convenience

      isize = 10 + 2*inznbmr
      allocate(idata(isize))

      !  set grid sizes
 
      idata(1) = inznbmri
      idata(2) = inznbmr
      idata(3) = inphi
      idata(4) = inth0
      if(iudsym.eq.1) then
         idata(5) = 1
         ifac=1
      else
         idata(5) = 0
         ifac=2
      endif
 
      idata(8:10) = 0  ! spares...

      !  get the no. of theta zones per zone row and total no. of zones

      allocate(nthzns(inznbmr),nthzsum(0:inznbmr))
 
      nthzsum(0)=0
 
      do ir=1,inznbmr
 
         nthzns(ir)=ifac*ir*inth0
         nthzsum(ir)=nthzsum(ir-1) + inphi*nthzns(ir)
 
         if(ir.eq.inznbmri) then
            idata(6)=nthzsum(ir)  ! total #zones inside plasma
         endif

         if(ir.eq.inznbmr) then
            idata(7)=nthzsum(ir)  ! total #zones inside & outside
         endif
      enddo
 
      idata(8:10) = 0  ! spares...

      idum=10+inznbmr
      idata(11:idum)=nthzns
      idata(idum+1:idum+inznbmr)=nthzsum(1:inznbmr)

      !  OK: create blackbox

      call xplasma_create_blackbox(s,mcname,xplasma_bbgtype,id,ierr, &
           iarray=idata, r8array=zx, label=label)

      deallocate(zx,idata,nthzns,nthzsum)
      if(ierr.ne.0) return
 
      !  get space for the zone volumes

      call xplasma_mcgrid_info(s,id,ierr, nzons=ifbznsi)
      if(ierr.ne.0) return

      allocate(zvols(ifbznsi))
      call xplasma_mcgrid_volumes(s,id,zvols,ierr)
      deallocate(zvols)

    end subroutine xplasma_mcgrid_define
 
    subroutine xplasma_mcgrid_volumes(s,id_mcgrid,vols,ierr, &
         lstart0_th,ccwflag_th)

      !  return the zone volumes associated with a MCgrid.  Compute them
      !  if necessary.

      type(xplasma), pointer :: s
      integer, intent(in) :: id_mcgrid   ! MC grid id
      real*8, dimension(:) :: vols       ! the array of volume elements
      integer, intent(out) :: ierr       ! status code (0=OK)

      ! theta indexing options
      logical, intent(in), optional :: lstart0_th ! default depends on symmetry
      logical, intent(in), optional :: ccwflag_th ! default: T

      ! lstart0_th=T means user's first theta zone lower bdy is 0; F means it
      ! is -pi.  If defaulted: updown symmetric data assumed to start at 0;
      ! updown asymmetric data assumed to start at -pi.

      ! ccwflag_th=T means zone index increases as one goes along a row of
      ! zones at fixed radial index in a counter-clockwise direction; F means
      ! the opposite.  T is the default.

      !----------------------------------------------
      integer :: ifbznsi  ! number of zones inside plasma
      integer :: inznbmri ! number of zone rows inside plasma
      integer :: inznbmr  ! total number of zone rows
      integer :: idv,idmc,idint,idin2,isizf,inumf,itype,inumchis,iertmp
      integer :: inthz,inths,izp,izn,ir,ith,inzr,id
      character*32 volname,mcname
      character*120 zmsg
      integer, dimension(:), pointer :: idata
      real*8, dimension(:), pointer :: r8data
      real*8, dimension(:), allocatable :: zchis
      real*8, dimension(:,:), allocatable :: voltmp
      real*8 :: x2(2),zfac
      logical :: udsym
      !----------------------------------------------

      id=id_mcgrid

      call xplasma_mcgrid_info(s,id,ierr, &
           nzons=ifbznsi, nzrow=inznbmri, nzrow_ext=inznbmr, name=mcname, &
           udsym=udsym)
      if(ierr.ne.0) return

      if(size(vols).ne.ifbznsi) then
         ierr=9999
         write(zmsg,*) ' ?xplasma_mcgrid_volume: #of MCgrid zones in plasma:',&
              ifbznsi
         call xplasma_errmsg_append(s,trim(zmsg))
         write(zmsg,*) '  size of passed array "vols": ',size(vols), &
              ' does not match.'
         return
      endif

      volname = trim(mcname)//'_DVOL'

      call xplasma_find_item(s,volname,idv,ierr,nf_noerr=.TRUE.)
      if(ierr.ne.0) return

      if(idv.gt.0) then
         call xplasma_get_item_info(s,idv,ierr,itype=itype)
         if(ierr.ne.0) return
         if(itype.ne.xplasma_blkbxType) then
            ierr=9999
            call xplasma_errmsg_append(s,' ?xplasma_mcgrid_volume: internal error looking for volumes: '//trim(volname))
            return
         endif

         call xplasma_blackBox_retrieve(s,idv,ierr, itype=itype, &
              ia_ptr=idata, r8a_ptr=r8data)

         if(itype.ne.xplasma_bbftype) then
            ierr=9999
            call xplasma_errmsg_append(s,' ?xplasma_mcgrid_volume: internal error volumes not a MCgrid profile: '//trim(volname))
            nullify(idata,r8data)
            return
         endif

         idmc=idata(1)
         isizf=idata(2)
         inumf=idata(3)

         if(idmc.ne.id) then
            ierr=9999
            call xplasma_errmsg_append(s,' ?xplasma_mcgrid_volume: internal error: MCgrid id mismatch: '//trim(volname))
            nullify(idata,r8data)
            return
         endif

         if(isizf.ne.ifbznsi) then
            ierr=9999
            call xplasma_errmsg_append(s,' ?xplasma_mcgrid_volume: internal error: profile size mismatch: '//trim(volname))
            nullify(idata,r8data)
            return
         endif

         if(inumf.gt.0) then
            ierr=9999
            call xplasma_errmsg_append(s,' ?xplasma_mcgrid_volume: internal error: profile rank mismatch: '//trim(volname))
            nullify(idata,r8data)
            return
         endif

         !  OK -- actual retrieve below...

         nullify(idata,r8data)

      else

         !  OK, volumes not found: compute them...

         allocate(idata(inznbmri))
         call xplasma_mcgrid_info(s,id,ierr, &
              nths=idata)

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         izp=1
         idint=0
         idin2=0
         do ir=1,inznbmri
 
            x2(1)=(ir-1)*ONE/(inznbmri)
            x2(2)=ir*ONE/(inznbmri)
            call xplasma_create_integ(s,trim(mcname)//'_VOLINT', &
                 x2, idint, ierr, &
                 rhomin=EPS7)

            if(ierr.ne.0) then
               call xplasma_errmsg_append(s, &
                    ' ?xplasma_mcgrid_volumes: integrator setup failed for:'//&
                    ' '//trim(mcname)//' (1d)')
               exit
            endif

            inthz= idata(ir)
            inths= inthz + 1

            allocate(zchis(inths))
            if(udsym) then
               do ith=1,inths
                  zchis(ith)= (ith-1)*CPI/(inthz)
               enddo
               zfac=2.0d0   ! multiply upper half volume by two
            else
               do ith=1,inths
                  zchis(ith)= -CPI + (ith-1)*C2PI/(inthz)
               enddo
               zfac=1.0d0
            endif

            call xplasma_augment_integ(s,trim(mcname)//'_VOLIN2',idint, &
                 zchis, idin2,ierr)

            deallocate(zchis)

            if(ierr.ne.0) then
               call xplasma_errmsg_append(s, &
                    ' ?xplasma_mcgrid_volumes: integrator setup failed for:'//&
                    ' '//trim(mcname)//' (2d)')
               exit
            endif

            izn = izp+inthz-1

            allocate(voltmp(inthz,1))

            call xplasma_2d_zonint(s,idin2,'dVol',voltmp,ierr)

            if(ierr.ne.0) then
               deallocate(voltmp)
               call xplasma_errmsg_append(s, &
                    ' ?xplasma_mcgrid_volumes: integration failed for:'//&
                    ' '//trim(mcname))
               exit
            else
               vols(izp:izn)=voltmp(1:inthz,1)*zfac
               deallocate(voltmp)
            endif

            izp=izn+1
         enddo

         if(idin2.gt.0) call xplasma_remove_item(s,idin2,iertmp)
         if(idint.gt.0) call xplasma_remove_item(s,idint,iertmp)

         if(ierr.eq.0) then

            call xplasma_mcgrid_putobj(s,id,volname,idv,ierr, &
                 data_1d=vols, &
                 label = trim(mcname)//' zone volumes', units = 'm**3')

         endif

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)

         deallocate(idata)

      endif

      ! retrieve here with ordering as specified by user optional
      ! arguments...

      if(ierr.eq.0) then
         call xplasma_mcgrid_getobj(s,idv,ierr, &
              data_1d=vols,lstart0_th=lstart0_th,ccwflag_th=ccwflag_th)
      endif

    end subroutine xplasma_mcgrid_volumes

    subroutine xplasma_mcgrid_info(s,id_mcgrid,ierr, &
         nzrow,nzrow_ext, nzons,nzons_ext, nth0,nphi, nths, udsym, &
         name,label,author)

      ! fetch information on xplasma_mcgrid grid

      type(xplasma), pointer :: s
      integer, intent(in) :: id_mcgrid  ! id of MCgrid object to be queried
      integer, intent(out):: ierr       ! status code 0=OK

      integer, intent(out), optional :: nzrow  ! # of zone rows in plasma
      integer, intent(out), optional :: nzrow_ext  ! # inside and outside

      integer, intent(out), optional :: nzons  ! # of zones in plasma
      integer, intent(out), optional :: nzons_ext  ! # inside and outside

      integer, intent(out), optional :: nth0   ! # of zones @axis cover [0:pi]
      integer, intent(out), optional :: nphi   ! # of phi zones (1:axisymmetry)

      integer, dimension(:), intent(out), optional :: nths
      ! number of theta zones in each zone row

      logical, intent(out), optional :: udsym

      character*(*), intent(out), optional :: name
      character*(*), intent(out), optional :: label
      character*(*), intent(out), optional :: author

      !-------------------------------------------------------
      integer :: itype,inznbmr,inznbmri,isize,id

      integer, dimension(:), pointer :: ia_ptr
      real*8, dimension(:), pointer :: r8a_ptr
      !-------------------------------------------------------

      id = id_mcgrid

      if(present(nzrow)) nzrow = 0
      if(present(nzrow_ext)) nzrow_ext = 0

      if(present(nzons)) nzons = 0
      if(present(nzons_ext)) nzons_ext = 0

      if(present(nth0)) nth0 = 0
      if(present(nphi)) nphi = 0

      if(present(nths)) nths = 0

      if(present(udsym)) udsym=.FALSE.

      if(present(name)) name=' '
      if(present(label)) label=' '
      if(present(author)) author=' '

      ierr = 0

      if(present(name)) call xplasma_get_item_info(s,id,ierr,name=name)
      if(ierr.ne.0) return

      if(present(label)) call xplasma_get_item_info(s,id,ierr,label=label)
      if(ierr.ne.0) return

      if(present(author)) call xplasma_get_item_info(s,id,ierr,author=author)
      if(ierr.ne.0) return

      call xplasma_blackBox_retrieve(s,id,ierr, itype=itype, &
           ia_ptr=ia_ptr, r8a_ptr=r8a_ptr)

      if(ierr.ne.0) return
      if(itype.ne.xplasma_bbgtype) then
         ierr=9999
         call xplasma_errmsg_append(s, &
              ' ?xplasma_mcgrid_info: not a MCgrid object: '//trim(name))
         nullify(ia_ptr,r8a_ptr)
         return
      endif

      inznbmri= ia_ptr(1)
      inznbmr = ia_ptr(2)

      if(present(nzrow)) nzrow=inznbmri
      if(present(nzrow_ext)) nzrow_ext=inznbmr

      if(present(nth0)) nth0=ia_ptr(4)
      if(present(nphi)) nphi=ia_ptr(3)

      if(present(udsym)) udsym = (ia_ptr(5).eq.1)  ! updown symmetry flag

      if(present(nzons)) nzons= ia_ptr(10+inznbmr+inznbmri)
      if(present(nzons_ext)) nzons_ext= ia_ptr(10+inznbmr+inznbmr)

      if(present(nths)) then
         isize=min(size(nths),inznbmr)
         nths(1:isize)=ia_ptr(10+1:10+isize)
      endif

      nullify(ia_ptr,r8a_ptr)

    end subroutine xplasma_mcgrid_info

    subroutine xplasma_mcgrid_profInfo(s,id_mcprof,ierr, &
         id_mcgrid,nzons,irank,idim1,idim2, &
         gridId1,gridId2,normcode, &
         name,label,units,author)

      ! fetch information on a single profile defined over an MCgrid...

      type (xplasma), pointer :: s
      integer, intent(in) :: id_mcprof   ! profile ID
      integer, intent(out) :: ierr       ! status code returned: 0=OK

      integer, intent(out), optional :: id_mcgrid   ! MCgrid ID
      integer, intent(out), optional :: nzons       ! #spatial zones
      integer, intent(out), optional :: irank       ! #non-spatial dimensions
      integer, intent(out), optional :: idim1,idim2 ! non-spatial dim. sizes

      !  grid IDs for non-spatial dimensions (if available)

      integer, intent(out), optional :: gridId1
      integer, intent(out), optional :: gridId2

      ! normalization code

      integer, intent(out), optional :: normcode

      !-------------

      character*(*), intent(out), optional :: name
      character*(*), intent(out), optional :: label
      character*(*), intent(out), optional :: units
      character*(*), intent(out), optional :: author

      !---------------------------------------
      integer, dimension(:), pointer :: iadata
      integer :: itype,iertmp
      character*32 :: zname
      !---------------------------------------

      if(present(id_mcgrid)) id_mcgrid=0
      if(present(nzons)) nzons=0
      if(present(irank)) irank=0
      if(present(idim1)) idim1=0
      if(present(idim2)) idim2=0

      if(present(gridId1)) gridId1=0
      if(present(gridId2)) gridId2=0

      if(present(normcode)) normcode=0

      if(present(name)) name=' '
      if(present(label)) label=' '
      if(present(units)) units=' '
      if(present(author)) author=' '

      call xplasma_blackbox_info(s,id_mcprof,ierr, type=itype, name=zname)
      if((ierr.ne.0).or.(itype.ne.xplasma_bbftype)) then
         if(ierr.ne.0) then
            zname='(no name)'
         else
            ierr=9999
         endif
         call xplasma_errmsg_append(s, &
              ' %xplasma_mcgrid_profInfo: not a MCgrid profile: '//trim(zname))
         return
      endif

      if(present(name)) name=zname

      call xplasma_blackbox_retrieve(s,id_mcprof,iertmp, ia_ptr=iadata)

      if(present(id_mcgrid)) id_mcgrid=iadata(1)
      if(present(nzons)) nzons=iadata(2)
      if(present(irank)) irank=iadata(3)
      if(present(idim1)) idim1=iadata(4)
      if(present(idim2)) idim2=iadata(5)

      call xplasma_get_item_info(s,id_mcprof,iertmp, &
           label=label, units=units, author=author)

      if(size(iadata).gt.6) then
         if(present(gridId1)) gridId1=iadata(10)
         if(present(gridId2)) gridId2=iadata(11)
         if(present(normcode)) normcode=iadata(15)
      endif

    end subroutine xplasma_mcgrid_profInfo

    subroutine xplasma_mcgrid_putobj(s,id_mcgrid,pname,id_p,ierr, &
         lstart0_th,ccwflag_th, &
         data_1d,data_2d,data_3d, &
         gridId1,gridId2,normcode, &
         label,units)

      ! define a profile "object" over a MCgrid.  This can define
      ! a single data value per grid zone (use data_1d), or a vector of
      ! of data values per grid zone (use data_2d) or a 2d array of values
      ! per grid zone (use data_3d).

      type(xplasma), pointer :: s
      integer, intent(in) :: id_mcgrid  ! MCgrid over which profile is defined.
      character*(*), intent(in) :: pname ! name of profile
      integer, intent(out) :: id_p      ! xplasma ID of profile returned
      integer, intent(out) :: ierr      ! status code (0=OK)

      ! theta indexing options
      logical, intent(in), optional :: lstart0_th ! default depends on symmetry
      logical, intent(in), optional :: ccwflag_th ! default: T

      ! lstart0_th=T means user's first theta zone lower bdy is 0; F means it
      ! is -pi.  If defaulted: updown symmetric data assumed to start at 0;
      ! updown asymmetric data assumed to start at -pi.

      ! ccwflag_th=T means zone index increases as one goes along a row of
      ! zones at fixed radial index in a counter-clockwise direction; F means
      ! the opposite.  T is the default.

      ! NOTE: similar options for phi mapping could be defined but that has
      ! not yet been done; may be an issue if/when we get beyond axisymmetry
      ! in XPLASMA...

      ! Although following arguments are optional, exactly one must be
      ! provided.  The last dimension of the array must match the size
      ! of the entire grid, or of the subset of the grid covering the
      ! plasma (either or)

      real*8, intent(in), dimension(:), optional :: data_1d
      real*8, intent(in), dimension(:,:), optional :: data_2d
      real*8, intent(in), dimension(:,:,:), optional :: data_3d

      !  if data_2d or data_3d is specified-- it is desirable to specify
      !  a grid which corresponds to the 1st dimension of the data (default
      !  0, meaning no grid was specified).

      integer, intent(in), optional :: gridId1

      !  if data_3d is specified-- it is desirable to specify a grid which
      !  corresponds to the 2nd dimension of the data 

      integer, intent(in), optional :: gridId2
      
      !  (the actual data_2d or data_3d array dimensions will be one less 
      !  than the corresponding grid sizes, for zone oriented "binned" data).

      !-------------------
      !  normalization code
      !     0 - none -- default -- no information
      !    -1 - divide by vol(...) to get a density
      !    +1 - multiply by vol(...) to integrate a density
      !    +2 -- multiply by vol(...)*d[grid1]*d[grid2]/2 (FBM normalization)
      
      integer, intent(in), optional :: normcode

      !-----------------
      character*(*), intent(in), optional :: label
      character*(*), intent(in), optional :: units

      !---------------------------------------------------
      integer :: icount,irank,isizp,idims(3),isize_all,iertmp
      integer :: nzons,nzons_ext,nzrow_ext,ii,jj,kk,iadr,itot,ir,kth
      real*8, dimension(:), pointer :: r8buf
      integer :: intbuf(15)
      integer, dimension(:), allocatable :: nths
      character*120 zmsg
      logical :: udsym,istart0,ccwflag,pishift,reverse
      !---------------------------------------------------

      ierr=0

      call xplasma_mcgrid_info(s,id_mcgrid,iertmp, &
           udsym=udsym,nzons=nzons,nzons_ext=nzons_ext,nzrow_ext=nzrow_ext)
      if(iertmp.ne.0) then
         call xplasma_errmsg_append(s,' ?xplasma_mcgrid: valid MCgrid ID must be provided.')
         ierr=iertmp
      endif

      icount=0

      if(present(data_1d)) then
         irank=0
         isizp=size(data_1d)
         icount=icount + 1
         isize_all=isizp
      endif

      if(present(data_2d)) then
         irank=1
         isizp=size(data_2d,2)
         idims(1)=size(data_2d,1)
         icount = icount + 1
         isize_all = isizp*idims(1)
      endif

      if(present(data_3d)) then
         irank=2
         isizp=size(data_3d,3)
         idims(1)=size(data_3d,1)
         idims(2)=size(data_3d,2)
         icount = icount + 1
         isize_all = isizp*idims(1)*idims(2)
      endif

      if(icount.eq.0) then
         ierr = 9999
         call xplasma_errmsg_append(s, &
              ' ?xplasma_mcgrid_putobj: no input data provided.')
         call xplasma_errmsg_append(s, &
              '  one of the arguments "data_1d", "data_2d", or "data_3d" must be provided.')

      else if(icount.gt.1) then
         ierr = 9999
         call xplasma_errmsg_append(s, &
              ' ?xplasma_mcgrid_putobj: multiple input data objects provided.')
         call xplasma_errmsg_append(s, &
              '  only one of the arguments "data_1d", "data_2d", "data_3d" can be provided.')

      endif

      if(ierr.ne.0) return

      if((isizp.ne.nzons).and.(isizp.ne.nzons_ext)) then
         ierr=9999
         zmsg=' '
         write(zmsg,*) &
              ' ?xplasma_mcgrid_putobj: MCgrid sizes (in plasma, total): ', &
              nzons,nzons_ext
         call xplasma_errmsg_append(s,zmsg)

         zmsg=' '
         write(zmsg,*) &
              '  data object grid dimension: ',isizp,' does not match.'
         call xplasma_errmsg_append(s,zmsg)
         return
      endif

      intbuf = 0
      intbuf(1)=id_mcgrid
      intbuf(2)=isizp
      intbuf(3)=irank
      if(irank.gt.0) then
         intbuf(4:4+irank-1)=idims(1:irank)
      endif

      if(irank.ge.1) then
         if(present(gridId1)) intbuf(10)=gridId1

         if(irank.ge.2) then
            if(present(gridId2)) intbuf(11)=gridId2
         endif
      endif

      if(present(normcode)) intbuf(15)=normcode

      !  create storage for this data...

      call xplasma_create_blackbox(s,pname,xplasma_bbftype,id_p,ierr, &
           iarray=intbuf, r8asize=isize_all, &
           label=label, units=units, r8a_ptr=r8buf)
      if(ierr.ne.0) return

      allocate(nths(nzrow_ext))
      call xplasma_mcgrid_info(s,id_mcgrid,iertmp,nths=nths)

      ccwflag=.TRUE.
      if(present(ccwflag_th)) ccwflag=ccwflag_th

      if(udsym) then
         istart0=.TRUE.
         if(present(lstart0_th)) istart0=lstart0_th
         pishift=.FALSE.
         reverse=.not.istart0
      else
         istart0=.FALSE.
         if(present(lstart0_th)) istart0=lstart0_th
         pishift=istart0
         reverse=.not.ccwflag
      endif

      itot=0
      do ir=1,nzrow_ext
         if(itot.eq.isizp) exit

         do kk=1,nths(ir)
            kth=kk
            if(pishift) then
               kth=kth+nths(ir)/2
               if(kth.gt.nths(ir)) kth=kth-nths(ir)
            endif
            if(reverse) then
               kth=nths(ir)+1-kth
            endif

            if(present(data_1d)) then
               r8buf(itot+kk)=data_1d(itot+kth)

            else if(present(data_2d)) then
               iadr=itot*idims(1)
               do jj=1,idims(1)
                  r8buf(iadr+idims(1)*(kk-1)+jj)=data_2d(jj,itot+kth)
               enddo

            else if(present(data_3d)) then
               iadr=itot*idims(1)*idims(2)
               do jj=1,idims(2)
                  do ii=1,idims(1)
                     r8buf(iadr+idims(1)*idims(2)*(kk-1)+idims(1)*(jj-1)+ii)=&
                          data_3d(ii,jj,itot+kth)
                  enddo
               enddo
               
            endif
         enddo
         itot=itot+nths(ir)

      enddo

      deallocate(nths)
      nullify(r8buf)

    end subroutine xplasma_mcgrid_putobj

    subroutine xplasma_mcgrid_getobj(s,id_p,ierr, &
         lstart0_th,ccwflag_th, &
         data_1d,data_2d,data_3d,name,label,units,author)

      !  retrieve previously stored MCgrid profile
      !  (see xplasma_mcgrid_putobj)

      type (xplasma), pointer :: s
      integer, intent(in) :: id_p  ! MCgrid profile ID
      integer, intent(out) :: ierr ! exit status code (0=OK)

      ! theta indexing options
      logical, intent(in), optional :: lstart0_th ! default depends on symmetry
      logical, intent(in), optional :: ccwflag_th ! default: T

      ! lstart0_th=T means user's first theta zone lower bdy is 0; F means it
      ! is -pi.  If defaulted: updown symmetric data assumed to start at 0;
      ! updown asymmetric data assumed to start at -pi.

      ! ccwflag_th=T means zone index increases as one goes along a row of
      ! zones at fixed radial index in a counter-clockwise direction; F means
      ! the opposite.

      ! NOTE: similar options for phi mapping could be defined but that has
      ! not yet been done; may be an issue if/when we get beyond axisymmetry
      ! in XPLASMA...

      !-----------------
      !  one of the optional arguments data_1d, data_2d, data_3d must
      !  be provided to receive the data.

      !  the provided data array sizes must match the stored data.

      real*8, dimension(:), intent(out), optional :: data_1d
      real*8, dimension(:,:), intent(out), optional :: data_2d
      real*8, dimension(:,:,:), intent(out), optional :: data_3d

      !-----------------
      character*(*), intent(out), optional :: name
      character*(*), intent(out), optional :: label
      character*(*), intent(out), optional :: units
      character*(*), intent(out), optional :: author

      !---------------------------------------------------
      integer :: icount,irank,isizp,idims(3),isize_all,ipages,ipaged,iertmp
      integer :: nzons,nzons_ext,ii,jj,kk,kth,ir,iadr,itype,isize_r8,isize_int
      integer :: isizd,ialloc,ip

      real*8, dimension(:), pointer :: r8a_ptr,r8b_ptr

      integer :: id_mcgrid,nzrow_ext,itot
      integer, dimension(:), allocatable :: nths
      integer, dimension(:), pointer :: ia_ptr
      character*120 zmsg
      character*32 zname
      logical :: udsym,ccwflag,istart0,pishift,reverse
      !---------------------------------------------------

      ialloc=0
      ierr=0

      call xplasma_blackBox_info(s,id_p,iertmp, type=itype, name=zname, &
           isize=isize_int, r8size=isize_r8)

      if(iertmp.ne.0) then
         ierr=iertmp
      else
         if(itype.ne.xplasma_bbftype) then
            ierr=9999
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_mcgrid_getobj: id not for a MCgrid profile: '// &
                 trim(zname))
         endif
      endif

      if(present(name)) name=' '
      if(present(label)) label=' '
      if(present(units)) units=' '
      if(present(author)) author=' '

      icount=0

      if(present(data_1d)) then
         irank=0
         isizp=size(data_1d)
         icount=icount + 1
         isize_all=isizp
         ipages=1
      endif

      if(present(data_2d)) then
         irank=1
         isizp=size(data_2d,2)
         idims(1)=size(data_2d,1)
         icount = icount + 1
         ipages=idims(1)
         isize_all = isizp*ipages
      endif

      if(present(data_3d)) then
         irank=2
         isizp=size(data_3d,3)
         idims(1)=size(data_3d,1)
         idims(2)=size(data_3d,2)
         icount = icount + 1
         ipages=idims(1)*idims(2)
         isize_all = isizp*ipages
      endif

      if(icount.eq.0) then
         ierr = 9999
         call xplasma_errmsg_append(s, &
              ' ?xplasma_mcgrid_getobj: no output data array provided.')
         call xplasma_errmsg_append(s, &
              '  one of the arguments "data_1d", "data_2d", or "data_3d" must be provided.')

      else if(icount.gt.1) then
         ierr = 9999
         call xplasma_errmsg_append(s, &
              ' ?xplasma_mcgrid_getobj: multiple output data arrays provided.')
         call xplasma_errmsg_append(s, &
              '  only one of the arguments "data_1d", "data_2d", "data_3d" can be provided.')

      endif

      if(ierr.ne.0) return

      call xplasma_blackBox_retrieve(s,id_p,ierr, &
           ia_ptr=ia_ptr,r8a_ptr=r8a_ptr)

      if(ierr.ne.0) return

      id_mcgrid=ia_ptr(1)
      isizd=ia_ptr(2)

      if(irank.ne.ia_ptr(3)) then
         ierr=9999
         call xplasma_errmsg_append(s, &
              ' ?xplasma_mcgrid_getobj: rank mismatch in data.')

         zmsg=' '
         write(zmsg,*) '  item "'//trim(zname)//'" rank is: ',ia_ptr(3), &
              ' (dimensionality per spatial grid pt.)'
         call xplasma_errmsg_append(s,zmsg)

         zmsg=' '
         write(zmsg,*) '  rank of passed array is: ',irank
         call xplasma_errmsg_append(s,zmsg)
      endif

      if(isize_all.eq.isize_r8) then

         r8b_ptr => r8a_ptr   ! size is correct

      else

         ipaged=1
         if(irank.eq.1) ipaged=ia_ptr(4)
         if(irank.eq.2) ipaged=ia_ptr(4)*ia_ptr(5)

         call xplasma_mcgrid_info(s,id_mcgrid,ierr, &
              nzons=nzons,nzons_ext=nzons_ext)

         if((isizp.ne.nzons).and.(isizp.ne.nzons_ext)) then

            ierr=9999

            call xplasma_errmsg_append(s, &
                 ' ?xplasma_mcgrid_getobj: size mismatch in data.')

            zmsg=' '
            write(zmsg,*) '  item "'//trim(zname)//'" size is: ',isize_r8
            call xplasma_errmsg_append(s,zmsg)

            zmsg=' '
            write(zmsg,*) '  size of passed array to receive data: ',isize_all
            call xplasma_errmsg_append(s,zmsg)

            zmsg=' '
            write(zmsg,*) '  rank of passed array to receive data: ',irank+1
            call xplasma_errmsg_append(s,zmsg)

            zmsg=' '
            write(zmsg,*) '  grid size implied by passed array: ',isizp
            call xplasma_errmsg_append(s,zmsg)

            zmsg=' '
            write(zmsg,*) '  available mcgrid sizes: ',nzons,nzons_ext
            call xplasma_errmsg_append(s,zmsg)

            if(present(data_1d)) data_1d=ZERO
            if(present(data_2d)) data_2d=ZERO
            if(present(data_3d)) data_3d=ZERO

         else if(ipages.ne.ipaged) then

            ierr=9999

            call xplasma_errmsg_append(s, &
                 ' ?xplasma_mcgrid_getobj: #data pts per grid pt not correct.')

            zmsg=' '
            write(zmsg,*) '  item "'//trim(zname)//'" has: ',ipaged, &
                 ' data pts per grid pt.'
            call xplasma_errmsg_append(s,zmsg)

            zmsg=' '
            write(zmsg,*) '  passed array has: ',ipages, &
                 ' data pts per grid pt.'
            call xplasma_errmsg_append(s,zmsg)

         else

            ialloc=1
            allocate(r8b_ptr(isizp*ipages))

            ii=0
            jj=0
            do ip=1,ipages
               r8b_ptr(ii+1:ii+isizp)=r8a_ptr(jj+1:jj+isizp)
               ii=ii+isizp
               jj=jj+isizd
            enddo
         endif
      endif

      if(present(name)) call xplasma_get_item_info(s,id_p,iertmp,name=name)
      if(present(label)) call xplasma_get_item_info(s,id_p,iertmp,label=label)
      if(present(units)) call xplasma_get_item_info(s,id_p,iertmp,units=units)
      if(present(author)) call xplasma_get_item_info(s,id_p,iertmp,author=author)

      if(ierr.ne.0) then
         if(present(data_1d)) data_1d=ZERO
         if(present(data_2d)) data_2d=ZERO
         if(present(data_3d)) data_3d=ZERO

      else
         call xplasma_mcgrid_info(s,id_mcgrid,iertmp, udsym=udsym, &
              nzrow_ext=nzrow_ext)
         allocate(nths(nzrow_ext))
         call xplasma_mcgrid_info(s,id_mcgrid,iertmp, nths=nths)

         ccwflag=.TRUE.
         if(present(ccwflag_th)) ccwflag=ccwflag_th

         if(udsym) then
            istart0=.TRUE.
            if(present(lstart0_th)) istart0=lstart0_th
            pishift=.FALSE.
            reverse=.not.istart0
         else
            istart0=.FALSE.
            if(present(lstart0_th)) istart0=lstart0_th
            pishift=istart0
            reverse=.not.ccwflag
         endif

         itot=0
         do ir=1,nzrow_ext
            if(itot.eq.isizp) exit

            do kk=1,nths(ir)
               kth=kk
               if(pishift) then
                  kth=kth+nths(ir)/2
                  if(kth.gt.nths(ir)) kth=kth-nths(ir)
               endif
               if(reverse) then
                  kth=nths(ir)+1-kth
               endif
               
               if(present(data_1d)) then
                  data_1d(itot+kth)=r8b_ptr(itot+kk)
               else if(present(data_2d)) then
                  iadr=itot*idims(1)
                  do jj=1,idims(1)
                     data_2d(jj,itot+kth)=r8b_ptr(iadr+idims(1)*(kk-1)+jj)
                  enddo
               else if(present(data_3d)) then
                  iadr=itot*idims(1)*idims(2)
                  do jj=1,idims(2)
                     do ii=1,idims(1)
                        data_3d(ii,jj,itot+kth)= &
                             r8b_ptr(iadr+idims(1)*idims(2)*(kk-1) + &
                             idims(1)*(jj-1)+ii)
                     enddo
                  enddo
               endif

            enddo
            itot=itot+nths(ir)
         enddo
         deallocate(nths)
      endif

      if(ialloc.eq.1) deallocate(r8b_ptr)
      nullify(r8a_ptr,r8b_ptr)

    end subroutine xplasma_mcgrid_getobj

end module xplasma_mcgrid
