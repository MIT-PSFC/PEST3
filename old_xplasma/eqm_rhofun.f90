subroutine eqm_rhofun(iorder,id_axis,zlbl,zdata, &
     ibc1,zbc1,ibc2,zbc2, &
     id_rhofun,ierr)
  !
  ! set up f(rho); rho = radial flux coordinate, sqrt(Phi/Philim)
  !

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  !  establish a function of "rho" -- axis to plasma bdy
  !  input:

  integer iorder                    ! interpolation order

  !     for all of the following, (nrho) values are to be provided.

  !  iorder = 0 -- piecewise linear
  !  iorder = 1 -- use Hermite piecewise cubics with Akima method for
  !     setting the derivative at the grid points (modified at the bdys)
  !  iorder = 2 or 3 -- use a cubic spline

  !  iorder = 1 -- once differentiable
  !  iorder = 2 or 3 -- twice differentiable (spline)

  !  iorder = 99 -- replacement data for existing f(rho) object
  !    (interpolation order will match that of the original object)
  !    (id_axis must match that of the original object)

  !  iorder=100 -- piecewise linear data, may replace  same sized object
  !  iorder=101 -- Akima Hermite (may be replacement of same sized object)
  !  iorder=102 or 103 -- spline (may be replacement of same sized object)

  integer id_axis                   ! axis id -- must be rho or akin to rho

  !  if id_axis is not rho, then the size of the object is not nrho but the
  !  size of the given axis.  The given axis must be "kin" to rho, i.e. it
  !  must map the same space, but perhaps with a larger or smaller number of
  !  points differently distributed.

  character*(*) zlbl                ! function name (label)
  REAL*8 zdata(*)                   ! data (at least `nrho' words)

  integer ibc1                      ! bc @ rho_axis
  real*8 zbc1                       ! bc parameter
  integer ibc2                      ! bc @ rho_bdy
  real*8 zbc2                       ! bc parameter

  !                =0:  "standard" boundary condition:  "not a knot" for
  !                     splines, Akima for piecewise Hermite.
  !                =1:  assigned df/drho bc, --> df/drho=zbc1|2
  !                =2:  assigned d2f/drho2 bc (splines only), d2f/drho2=zbc1|2

  !  for iorder.eq.0, bc controls are ignored.

  !  output:

  integer id_rhofun                 ! function id, returned
  integer ierr                      ! completion code, 0=OK
  
  !-----------------------------
  integer :: icoord,iorderi,iertmp,inx,icount
  character*32 zname
  logical :: iforce = .FALSE.
  !-----------------------------

  ierr=0

  iorderi=iorder
  if(iorderi.lt.99) iorderi=max(0,min(2,iorderi))
  if(iorderi.ge.100) iorderi=max(0,min(2,(iorderi-100)))

  !  iorderi=99 gets special treatment below...

  !-----------------
  !  check axis id

  call xplasma_grid_info(s,id_axis,ierr, coord=icoord)

  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     write(lunerr,*) &
          ' ?eqm_rhofun ("',trim(zlbl),'"): id_axis=',id_axis,' invalid.'
     ierr=1
     return
  endif

  !  verify that the x axis coordinate is OK ("rho")

  if(icoord.ne.xplasma_rho_coord) then
     ierr=1
     call eq_errmsg(' ?eqm_rhofun:  id_axis invalid:')
     call xplasma_get_item_info(s,id_axis,iertmp, name=zname)
     write(lunerr,*) '  id_axis = ',id_axis,' name = ',trim(zname)
     write(lunerr,*) '  ...not a "rho" coordinate.'
     return
  endif

  call xplasma_grid_size(s,id_axis,inx,iertmp)

  !-----------------
  !  iorderi=99 means request to use prior splining method for this call;
  !  this implies profile must already exist:  check...

  if(iorderi.eq.99) then
     call xplasma_find_item(s,zlbl,id_rhofun,iertmp)
     if(iertmp.eq.0) then
        call xplasma_prof_info(s,id_rhofun,ierr, splineType=iorderi)
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

  call xplasma_create_1dprof(s,zlbl,id_axis,zdata(1:inx),id_rhofun,ierr, &
       iorderi,ibc1,zbc1,ibc2,zbc2)

  call xoi_author_clear(iertmp)

  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     ierr=1
     return
  else
     call xplasma_get_item_info(s,id_rhofun,iertmp, name=zname)
     if((zname.eq.'G').or.(zname.eq.'PSI')) then
        call xplasma_eqcheck(s,iforce,icount,ierr)
        if(icount.gt.0) then
           write(lunerr,*) ' %eqm_rhofun("'//trim(zname)// &
                '",...): equilibrium specification completed.'
        endif
     endif
  endif

end subroutine eqm_rhofun
