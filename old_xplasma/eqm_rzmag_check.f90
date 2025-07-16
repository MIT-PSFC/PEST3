subroutine eqm_rzmag_check(id1,id2,idsw,id_rho,id_chi,inrho,inchi,ierr)

  ! error check routine shared by eqm_rzmag*.f90
  ! find the ids of __RHO and __CHI or report an error

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  integer, intent(in) :: id1,id2         ! R & Z array dimensions
  integer, intent(in) :: idsw            ! 1: id1->rho dimension;
  ! 2: id2 ->rho dimension

  integer, intent(out) :: id_rho,id_chi  ! the grid ids
  integer, intent(out) :: inrho,inchi    ! the grid sizes
  integer, intent(out) :: ierr           ! completion code, 0=OK

  !------------------------------------
  integer :: iertmp
  !------------------------------------

  inrho=0
  inchi=0

  iertmp=0

  call xplasma_find_item(s,'__RHO',id_rho,ierr)
  if(ierr.ne.0) then
     iertmp=max(iertmp,ierr)
     write(lunerr,*) ' ?eqm_rzmag: "__RHO" grid not found.'
  endif

  call xplasma_find_item(s,'__CHI',id_chi,ierr)
  if(ierr.ne.0) then
     iertmp=max(iertmp,ierr)
     write(lunerr,*) ' ?eqm_rzmag: "__CHI" grid not found.'
  endif

  ierr=iertmp
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  call xplasma_grid_size(s,id_rho,inrho,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_rzmag: "__RHO" grid size fetch failed.'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  call xplasma_grid_size(s,id_chi,inchi,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_rzmag: "__RHO" grid size fetch failed.'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  if(idsw.eq.1) then
     if(id1.lt.inrho) then
        ierr=ierr+1
        write(lunerr,*) ' ?eqm_rzmag: "__RHO" dimension has ',inrho,' pts.'
        write(lunerr,*) '  but data 1st dimension only has: ',id1
     endif
     if(id2.lt.inchi) then
        ierr=ierr+1
        write(lunerr,*) ' ?eqm_rzmag: "__CHI" dimension has ',inchi,' pts.'
        write(lunerr,*) '  but data 2nd dimension only has: ',id2
     endif
  else
     if(id1.lt.inchi) then
        ierr=ierr+1
        write(lunerr,*) ' ?eqm_rzmag: "__CHI" dimension has ',inchi,' pts.'
        write(lunerr,*) '  but data 1st dimension only has: ',id1
     endif
     if(id2.lt.inrho) then
        ierr=ierr+1
        write(lunerr,*) ' ?eqm_rzmag: "__RHO" dimension has ',inrho,' pts.'
        write(lunerr,*) '  but data 2nd dimension only has: ',id2
     endif
  endif

end subroutine eqm_rzmag_check
