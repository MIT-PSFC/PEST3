subroutine eq_fun(id_fun,magdim,magfit,cyldim,cylfit,iwarn,ierr)

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  retrieve information on the specified function

  integer, intent(in) :: id_fun     ! function id code

  !  variation vs. magnetic coordinates...
  !  magdim=0 if no such variation is defined

  integer, intent(out) :: magdim    ! 1=> f(rho), 2=> f(rho,chi), ...
  integer, intent(out) :: magfit    ! 1=> Hermite, 2=> Spline, ...

  !  variation vs. cylindric coordinates
  !  cyldim=0 if no such variation is defined

  integer, intent(out) :: cyldim    ! 2=> f(R,Z), 3=> f(R,Z,phi)
  integer, intent(out) :: cylfit    ! 1=> Hermite, 2=> Spline, ...

  integer, intent(out) :: iwarn     ! 0=OK, 1=[non-standard axes]
  integer, intent(out) :: ierr      ! 0=OK, 1=[bad id_fun code]

  !  iwarn is also set if the name starts with "__" -- reserved to xplasma.
  !  if iwarn is set, then, magdim=cyldim=0 on exit.
  !---------------------------
  integer :: id_funs(2),idi,id,iertmp,irank,igtype,iorder
  character*32 :: fname
  !---------------------------
  
  magdim=0
  magfit=0

  cyldim=0
  cylfit=0

  ierr=0
  iwarn=0

  id_funs(1)=id_fun
  call xplasma_prof_info(s,id_funs(1),ierr, profId1=id_funs(2), name=fname)

  if(ierr.ne.0) then
     iwarn=1
     write(lunerr,*) '?eq_fun:  invalid id code:  ',id_fun
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  if(fname(1:2).eq.'__') then
     write(lunerr,*) ' %eq_fun: reserved name: ',trim(fname)
     iwarn=1
  else
     do idi=1,2
        id=id_funs(idi)
        if(id.eq.0) cycle

        call xplasma_prof_info(s,id,iertmp, &
             rank=irank, gridType=igtype, splineType=iorder)

        if(igtype.eq.1) then
           magdim=irank
           magfit=iorder

        else if(igtype.eq.2) then
           cyldim=irank
           cylfit=iorder
        endif
     enddo
  endif

end subroutine eq_fun

!------------------------------------------------
subroutine eq_fun_rhogrid(id_fun,id_grid)

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id_fun     ! function id (in)
  integer, intent(out) :: id_grid   ! indep. coord. grid id (out)

  !  return rho grid id for function (if these exist), or zero.

  !-------------------
  integer :: ierr
  !-------------------

  call xplasma_prof_gridInfo(s,id_fun,xplasma_rho_coord,id_grid,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_fun_rhogrid: error detected:'
     call xplasma_error(s,ierr,lunerr)
  endif
 
end subroutine eq_fun_rhogrid

!------------------------------------------------
subroutine eq_fun_chigrid(id_fun,id_grid)

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id_fun     ! function id (in)
  integer, intent(out) :: id_grid   ! indep. coord. grid id (out)

  !  return chi (pol. angle) grid id for function (if these exist), or zero.

  !-------------------
  integer :: ierr
  !-------------------

  call xplasma_prof_gridInfo(s,id_fun,xplasma_theta_coord,id_grid,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_fun_rhogrid: error detected:'
     call xplasma_error(s,ierr,lunerr)
  endif
 
end subroutine eq_fun_chigrid

!------------------------------------------------
subroutine eq_fun_phigrid(id_fun,id_grid)

  use xplasma_obj_instance
  use eq_module 
  implicit NONE

  integer, intent(in) :: id_fun     ! function id (in)
  integer, intent(out) :: id_grid   ! indep. coord. grid id (out)

  !  return phi grid id for function (if these exist), or zero.

  !-------------------
  integer :: ierr
  !-------------------

  call xplasma_prof_gridInfo(s,id_fun,xplasma_phi_coord,id_grid,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_fun_rhogrid: error detected:'
     call xplasma_error(s,ierr,lunerr)
  endif
 
end subroutine eq_fun_phigrid

!------------------------------------------------
subroutine eq_fun_Rgrid(id_fun,id_grid)

  use xplasma_obj_instance
  use eq_module 
  implicit NONE

  integer, intent(in) :: id_fun     ! function id (in)
  integer, intent(out) :: id_grid   ! indep. coord. grid id (out)

  !  return R grid id for function (if these exist), or zero.

  !-------------------
  integer :: ierr
  !-------------------

  call xplasma_prof_gridInfo(s,id_fun,xplasma_R_coord,id_grid,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_fun_rhogrid: error detected:'
     call xplasma_error(s,ierr,lunerr)
  endif
 
end subroutine eq_fun_Rgrid

!------------------------------------------------
subroutine eq_fun_Zgrid(id_fun,id_grid)

  use xplasma_obj_instance
  use eq_module 
  implicit NONE

  integer, intent(in) :: id_fun     ! function id (in)
  integer, intent(out) :: id_grid   ! indep. coord. grid id (out)

  !  return Z grid id for function (if these exist), or zero.

  !-------------------
  integer :: ierr
  !-------------------

  call xplasma_prof_gridInfo(s,id_fun,xplasma_Z_coord,id_grid,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_fun_rhogrid: error detected:'
     call xplasma_error(s,ierr,lunerr)
  endif
 
end subroutine eq_fun_Zgrid
