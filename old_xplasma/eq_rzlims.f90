subroutine eq_rzlims(Rmin,Rmax,Zmin,Zmax,ierr)

  !  return min/max R & Z of (R,Z) grids
  !  set ierr if grids are not available.

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  !  all output:
  real*8 Rmin,Rmax                  ! min/max R of grids
  real*8 Zmin,Zmax                  ! min/max Z of grids
  integer ierr                      ! completion code:  0 = OK

  !--------------------------------------------------
  integer :: id_Rc,id_Zc,ingrids
  !--------------------------------------------------

  Rmin = 0; Rmax = 0; Zmin = 0; Zmax = 0

  do
     call xplasma_find_item(s,'__R_coord',id_Rc,ierr); if(ierr.ne.0) exit
     call xplasma_find_item(s,'__Z_coord',id_Zc,ierr); if(ierr.ne.0) exit

     call xplasma_coord_info(s,id_Rc,ierr, ngrids=ingrids); if(ierr.ne.0) exit
     if(ingrids.eq.0) then
        write(lunerr,*) ' ?eq_rzlims: R grid has not been defined.'
        exit
     endif

     call xplasma_coord_info(s,id_Zc,ierr, ngrids=ingrids); if(ierr.ne.0) exit
     if(ingrids.eq.0) then
        write(lunerr,*) ' ?eq_rzlims: Z grid has not been defined.'
        exit
     endif

     call xplasma_coord_info(s,id_Rc,ierr, &
          xmin=Rmin,xmax=Rmax); if(ierr.ne.0) exit

     call xplasma_coord_info(s,id_Zc,ierr, &
          xmin=Zmin,xmax=Zmax); if(ierr.ne.0) exit

     exit
  enddo

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_rzlims: unexpected error:'
     call xplasma_error(s,ierr,lunerr)
  else if(ingrids.eq.0) then
     ierr=1
  endif

end subroutine eq_rzlims
