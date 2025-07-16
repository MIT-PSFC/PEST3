subroutine eqm_frhochi(iorder,id_ax1,id_ax2,zlbl,zdata,id1, &
     ibcrho0,zbcrho0,ibcrho1,zbcrho1,id_fun,ierr)

  !  set up a function f(rho,chi)
  !    rho is the radial coordinate; chi the poloidal angle coordinate
  !    e.g. for a 2d function over an axisymmetric tokamak plasma cross
  !    section.

  !    (lots of people use "theta" instead of "chi" for this angle coordinate)
  !    (but the name is not changed in order to preserve back compatibility
  !    with legacy code)

  !    this routine assumes chi describes the contour of a boundary surface
  !    going in a counterclockwise direction.  For the reverse, set iccw=0
  !    and call eqm_frhochi_ccw(...) -- see below

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  !  input:

  integer iorder                    ! desired interpolation order

  !  [in these comments nrho & nchi refer to the axes of the data passed]

  !  iorder=0:  bilinear interpolation on data on (nrho)x(nchi) grid
  !  iorder=1:  Akima Hermite interpolation on data on(nrho)x(nchi) grid
  !  iorder=2 or 3:  Spline interpolation on data on (nrho)x(nchi) grid

  !  iorder=99 -- this is a an update of the data in an existing
  !             f(rho,chi) object.  The fit order and axes used
  !             cannot change.

  !  iorder=100 -- bilinear data (may be replacement of same sized object)
  !  iorder=101 -- Akima Hermite (may be replacement of same sized object)
  !  iorder=102 or 103 -- spline (may be replacement of same sized object)

  integer id_ax1                    ! axis id -- 1st dimension of zdata
  integer id_ax2                    ! axis id -- 2nd dimension of zdata

  !  id_ax1 or id_ax2 must give a rho coordinate; the other coordinate must
  !  be a chi (poloidal angle) coordinate.

  character*(*) zlbl                ! function name (label)
  integer id1                       ! size of 1st dimension of zdata

  REAL*8 zdata(id1,*)               ! data (at least nrho*nchi words)

  !  id1 must be .ge. size of axis (id_ax1)

  !  bdy conditions are enforced only for iorder.ge.1

  integer ibcrho0,ibcrho1           ! BC's at extrema

  !  BC controls are as in "pspline":  0 -- default, 1 -- fix df/drho

  REAL*8 zbcrho0(*)               ! if ibc=1:  df/drho @ rho(1), nchi pts
                                  ! if ibc=2:  d2f/drho2 @ rho(1), nchi pts
  REAL*8 zbcrho1(*)               ! ditto @rho(nrho)

  !     the BC data zbcrho0,zbcrho1 are used iff ibcrho0,ibcrho1 are 1 or 2
  !        otherwise the BC is derived from the main data array
  !     for Hermite fits, ibc values of 0 or 1 (only) are allowed)
  !  output:

  integer id_fun                    ! function id code
  integer ierr                      ! completion code, 0=OK

  !-----------------------------
  integer :: iccw
  !-----------------------------

  iccw = 1
  call eqm_frhochi_ccw(iccw, &
       iorder,id_ax1,id_ax2,zlbl,zdata,id1, &
       ibcrho0,zbcrho0,ibcrho1,zbcrho1,id_fun,ierr)

end subroutine eqm_frhochi

subroutine eqm_frhochi_ccw(iccw,iorder,id_ax1,id_ax2,zlbl,zdata,id1, &
     ibcrho0,zbcrho0,ibcrho1,zbcrho1,id_fun,ierr)

  !  set up a function f(rho,chi)
  !    rho is the radial coordinate; chi the poloidal angle coordinate
  !    e.g. for a 2d function over an axisymmetric tokamak plasma cross
  !    section.

  !    (lots of people use "theta" instead of "chi" for this angle coordinate)
  !    (but the name is not changed in order to preserve back compatibility
  !    with legacy code)

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  !  input:

  integer :: iccw                   ! orientation of chi: 1 means CCW
  ! (CCW, counter-clockwise); 0 means CW, clockwise.
  ! Xplasma internals ALWAYS use the CCW orientation so an ordinate
  ! reversal is needed if iccw=0

  integer iorder                    ! desired interpolation order

  !  [in these comments nrho & nchi refer to the axes of the data passed]

  !  iorder=0:  bilinear interpolation on data on (nrho)x(nchi) grid
  !  iorder=1:  Akima Hermite interpolation on data on(nrho)x(nchi) grid
  !  iorder=2 or 3:  Spline interpolation on data on (nrho)x(nchi) grid

  !  iorder=99 -- this is a an update of the data in an existing
  !             f(rho,chi) object.

  !  iorder=100 -- bilinear data (may be replacement of same sized object)
  !  iorder=101 -- Akima Hermite (may be replacement of same sized object)
  !  iorder=102 or 103 -- spline (may be replacement of same sized object)

  !  Note: replacement is now *always* allowed; the 100...103 argument
  !  values are maintained for backwards compatibility with the old f77
  !  interface (original xplasma interface).

  !  SPECIAL NAMES:  "R" & "Z" vs. (rho,chi) -- axisymmetric equilibrium
  !  flux surfaces.  For these, the equilibrium order if set to 99 or 0
  !  is taken from a scalar stored in s%rzordr -

  integer id_ax1                    ! axis id -- 1st dimension of zdata
  integer id_ax2                    ! axis id -- 2nd dimension of zdata

  !  id_ax1 or id_ax2 must give a rho coordinate; the other coordinate must
  !  be a chi (poloidal angle) coordinate.

  character*(*) zlbl                ! function name (label)
  integer id1                       ! size of 1st dimension of zdata

  REAL*8 zdata(id1,*)               ! data (at least nrho*nchi words)

  !  id1 must be .ge. size of axis (id_ax1)

  !  bdy conditions are enforced only for iorder.ge.1

  integer ibcrho0,ibcrho1           ! BC's at extrema

  !  BC controls are as in "pspline":  0 -- default, 1 -- fix df/drho

  REAL*8 zbcrho0(*)               ! if ibc=1:  df/drho @ rho(1), nchi pts
                                  ! if ibc=2:  d2f/drho2 @ rho(1), nchi pts
  REAL*8 zbcrho1(*)               ! ditto @rho(nrho)

  !     the BC data zbcrho0,zbcrho1 are used iff ibcrho0,ibcrho1 are 1 or 2
  !        otherwise the BC is derived from the main data array
  !     for Hermite fits, ibc values of 0 or 1 (only) are allowed)
  !  output:

  integer id_fun                    ! function id code
  integer ierr                      ! completion code, 0=OK

  !-----------------------------
  ! local
  integer :: icoord1,icoord2,iorderi,iertmp,inx1,inx2,irho12,irzo
  character*32 zname
  real*8, dimension(:), allocatable :: zbc0,zbc1
  logical :: ccwflag
  !-----------------------------

  ierr=0

  iorderi=iorder
  if(iorderi.lt.99) iorderi=max(0,min(2,iorderi))
  if(iorderi.ge.100) iorderi=max(0,min(2,(iorderi-100)))

  if((zlbl.eq.'r').or.(zlbl.eq.'R').or.(zlbl.eq.'z').or.(zlbl.eq.'Z')) then

     !  R & Z names -- special treatment -- R(rho,theta) & Z(rho,theta)
     !  the equilibrium flux surfaces

     call xplasma_global_info(s,ierr, rzOrder=irzo)
     if(ierr.ne.0) then
        call xplasma_error(s,ierr,lunerr)
        return
     endif
     if(iorderi.eq.99) iorderi=irzo
     iorderi=max(irzo,iorderi)

  endif

  !  iorderi=99 gets more special treatment below...

  !--------------------------
  !  check axes

  call xplasma_grid_info(s,id_ax1,ierr, coord=icoord1)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     write(lunerr,*) &
          ' ?eqm_frhochi ("',trim(zlbl),'"): id_ax1=',id_ax1,' invalid.'
     ierr=1
  endif
  iertmp=ierr

  call xplasma_grid_info(s,id_ax2,ierr, coord=icoord2)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     write(lunerr,*) &
          ' ?eqm_frhochi ("',trim(zlbl),'"): id_ax2=',id_ax2,' invalid.'
     ierr=1
  endif
  ierr=max(ierr,iertmp)
  if(ierr.gt.0) return

  ! OK -- both are x axes at least.  Now, one needs to match the xplasma
  ! module rho and one the xplasma module theta coordinate, order not
  ! important...

  if((icoord1.ne.xplasma_rho_coord).and.(icoord2.ne.xplasma_rho_coord)) then
     ierr=1
     call eq_errmsg(' ?eqm_frhochi: id_ax1 & id_ax2 error: "rho" not found:')
     call xplasma_get_item_info(s,id_ax1,iertmp, name=zname)
     write(lunerr,*) '  id_ax1 = ',id_ax1,' name = ',trim(zname)
     call xplasma_get_item_info(s,id_ax2,iertmp, name=zname)
     write(lunerr,*) '  id_ax2 = ',id_ax2,' name = ',trim(zname)
     return
  endif

  if((icoord1.ne.xplasma_theta_coord).and. &
       (icoord2.ne.xplasma_theta_coord)) then
     ierr=1
     call eq_errmsg(' ?eqm_frhochi: id_ax1 & id_ax2 error: "theta" not found:')
     call xplasma_get_item_info(s,id_ax1,iertmp, name=zname)
     write(lunerr,*) '  id_ax1 = ',id_ax1,' name = ',trim(zname)
     call xplasma_get_item_info(s,id_ax2,iertmp, name=zname)
     write(lunerr,*) '  id_ax2 = ',id_ax2,' name = ',trim(zname)
     return
  endif

  if(icoord1.eq.xplasma_rho_coord) then
     irho12=1
  else
     irho12=2
  endif

  call xplasma_grid_size(s,id_ax1,inx1,iertmp)
  call xplasma_grid_size(s,id_ax2,inx2,iertmp)

  !-----------------
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
        call eq_errmsg(' ?eqm_rhofun: iorder=99 call failed: '//trim(zlbl)// &
             ' does not exist.')
        ierr=1
     endif
  endif

  if(ierr.ne.0) return

  !-----------------
  !  iorderi now contains the desired spline type in all cases, so...

  call xoi_author_set(ierr)

  ccwflag = iccw.eq.1

  if(irho12.eq.1) then
     allocate(zbc0(inx2)); zbc0=0
     allocate(zbc1(inx2)); zbc1=0
     if(ibcrho0.gt.0) zbc0=zbcrho0(1:inx2)
     if(ibcrho1.gt.0) zbc1=zbcrho1(1:inx2)
     call xplasma_create_2dprof(s,zlbl,id_ax1,id_ax2,zdata(1:inx1,1:inx2),&
          id_fun,ierr, &
          iorderi,ibcx1a=ibcrho0,zbcx1a=zbc0,ibcx1b=ibcrho1,zbcx1b=zbc1, &
          ccwflag2=ccwflag)
  else
     allocate(zbc0(inx1)); zbc0=0
     allocate(zbc1(inx1)); zbc1=0
     if(ibcrho0.gt.0) zbc0=zbcrho0(1:inx1)
     if(ibcrho1.gt.0) zbc1=zbcrho1(1:inx1)
     call xplasma_create_2dprof(s,zlbl,id_ax1,id_ax2,zdata(1:inx1,1:inx2),&
          id_fun,ierr, &
          iorderi,ibcx2a=ibcrho0,zbcx2a=zbc0,ibcx2b=ibcrho1,zbcx2b=zbc1, &
          ccwflag1=ccwflag)
  endif

  call xoi_author_clear(iertmp)

  return
end subroutine eqm_frhochi_ccw
