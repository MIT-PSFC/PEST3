subroutine eqm_fRZ(iorder,id_R,id_Z,zlbl,zdata,id1, &
     ibcR0,zbcR0,ibcR1,zbcR1, &
     ibcZ0,zbcZ0,ibcZ1,zbcZ1, &
     id_fun,ierr)

  !  establish a function of "R & Z"

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  !  input:
  integer iorder                    ! desired interpolation order

  !  iorder=0:  bilinear interpolation on data on (nR)x(nZ) grid
  !  iorder=1:  Akima Hermite interpolation on data on(nR)x(nZ) grid
  !  iorder=2 or 3:  Spline interpolation on data on (nR)x(nZ) grid

  !  iorder=99 -- this is a an update of the data in an existing
  !             f(R,Z) object.  The fit order and axes used
  !             cannot change.

  !  iorder=100 -- bilinear data (may be replacement of same sized object)
  !  iorder=101 -- Akima Hermite (may be replacement of same sized object)
  !  iorder=102 or 103 -- spline (may be replacement of same sized object)

  integer id_R                      ! axis id -- 1st dimension, R
  integer id_Z                      ! axis id -- 2nd dimension, Z

  character*(*) zlbl                ! function name (label)
  integer id1                       ! size of 1st dimension of zdata

  REAL*8 zdata(id1,*)               ! data (at least nR*nZ words)

  !  nR and nZ refer to the axes (id_R,id_Z)
  !  id1 must be .ge. size of axis (id_R)

  !  bdy conditions are enforced only for iorder.ge.1
  integer ibcR0,ibcR1               ! BC's at extrema

  !  BC controls are as in "pspline":  0 -- default, 1 -- fix df/dR
  !                                    2 -- fix d2f/dR2, etc...

  REAL*8 zbcR0(*)                   ! if ibc=1:  df/dR @ R(1), nZ pts
                                    ! if ibc=2:  d2f/dR2 @ R(1), nZ pts
  REAL*8 zbcR1(*)                   ! ditto @R(nR)

  !  zbcR0 and zbcR1 ignored if ibcR0=ibcR1=0

  !  bdy conditions are enforced only for iorder.ge.1

  integer ibcZ0,ibcZ1               ! BC's at extrema

  !  BC controls are as in "pspline":  0 -- default, 1 -- fix df/dZ
  !                                    2 -- fix d2f/dZ2, etc...

  REAL*8 zbcZ0(*)                   ! if ibc=1:  df/dZ @ Z(1), nR pts
                                    ! if ibc=2:  d2f/dZ2 @ Z(1), nR pts
  REAL*8 zbcZ1(*)                   ! ditto @Z(nZ)

  !  zbcZ0 and zbcZ1 ignored if ibcZ0=ibcZ1=0

  !  output:

  integer id_fun                    ! function id code
  integer ierr                      ! completion code, 0=OK

  !-----------------------------
  !  local copy

  integer iorderi,iertmp,icoord1,icoord2,ierdum
  integer :: nR,nZ
  character*32 zname

  !-----------------------------
  ierr=0

  iorderi=iorder
  if(iorderi.lt.99) iorderi=max(0,min(2,iorderi))
  if(iorderi.ge.100) iorderi=max(0,min(2,(iorderi-100)))

  !  iorderi=99 gets more special treatment below...
  !--------------------------
  !  check axes

  call xplasma_grid_info(s,id_R,ierr, coord=icoord1)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     write(lunerr,*) &
          ' ?eqm_frz ("',trim(zlbl),'"): id_R=',id_R,' invalid grid.'
     ierr=1
  else if(icoord1.ne.xplasma_R_coord) then
     write(lunerr,*) & 
          ' ?eqm_frz ("',trim(zlbl),'"): id_R=',id_R,' *NOT* major radius.'
     call xplasma_grid_info(s,id_R,ierdum, name=zname)
     write(lunerr,*) &
          '  name of id_R grid: ',trim(zname)
     ierr=2
  endif
  iertmp=ierr

  call xplasma_grid_info(s,id_Z,ierr, coord=icoord1)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     write(lunerr,*) &
          ' ?eqm_frz ("',trim(zlbl),'"): id_Z=',id_Z,' invalid grid.'
     ierr=1
  else if(icoord1.ne.xplasma_Z_coord) then
     write(lunerr,*) & 
          ' ?eqm_frz ("',trim(zlbl),'"): id_Z=',id_Z,' *NOT* Z.'
     call xplasma_grid_info(s,id_Z,ierdum, name=zname)
     write(lunerr,*) &
          '  name of id_Z grid: ',trim(zname)
     ierr=2
  endif
  ierr=max(ierr,iertmp)
  if(ierr.gt.0) return

  !  OK: grid sizes are...

  call xplasma_grid_size(s,id_R,nR,ierr)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     write(lunerr,*) ' ?error in eqm_frz'
  endif

  call xplasma_grid_size(s,id_Z,nZ,ierr)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     write(lunerr,*) ' ?error in eqm_frz'
  endif

  !--------------------------
  !  iorderi=99 means request to use prior splining method for this call;
  !  this implies profile must already exist:  check...

  if(iorderi.eq.99) then
     call xplasma_find_item(s,zlbl,id_fun,iertmp)
     if(iertmp.eq.0) then
        call xplasma_prof_info(s,id_fun,ierr, splineType=iorderi)
     else
        ierr=iertmp
     endif
     if(ierr.ne.0) then
        call xplasma_error(s,ierr,lunerr)
        call eq_errmsg(' ?eqm_frz: iorder=99 call failed: '//trim(zlbl)// &
             ' does not exist.')
        ierr=1
     endif
  endif
  if(ierr.ne.0) return

  !--------------------------
  !  iorderi now contains the desired spline type in all cases, so...

  call xoi_author_set(ierr)

  call xplasma_create_2dprof(s,zlbl,id_R,id_Z,zdata(1:nR,1:nZ), &
       id_fun,ierr, &
       ispline=iorderi, &
       ibcx1a=ibcR0,zbcx1a=zbcR0(1:nR),ibcx1b=ibcR1,zbcx1b=zbcR1(1:nR), &
       ibcx2a=ibcZ0,zbcx2a=zbcZ0(1:nZ),ibcx2b=ibcZ1,zbcx2b=zbcZ1(1:nZ))

  call xoi_author_clear(iertmp)

end subroutine eqm_fRZ
