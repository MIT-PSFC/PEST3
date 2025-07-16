subroutine eqm_bset(isnccwb,isnccwi)

  use xplasma_obj_instance
  use eq_module

  !  setup B field signs

  IMPLICIT NONE

  !  input:
  integer isnccwb                   ! toroidal field sign factor...
  integer isnccwi                   ! plasma current sign factor...

  !  CCW = counter-clockwise
  !  isnccwb=+1 -- TF is CCW looking down from top of machine
  !         =-1 -- TF is clockwise

  !  isnccwi=+1 -- Ip is CCW looking down from top of machine
  !         =-1 -- Ip is clockwise

  !  if isnccwi=+1 Bpol points down on the large major radius side
  !  of the machine...

  !-------------------------------------

  INTEGER ierr
  integer :: nsnccwb,nsnccwi,icount
  logical :: iforce = .FALSE.
  !-------------------------------------

  ierr=0
  nsnccwb=0
  if((isnccwb.ne.1).and.(isnccwb.ne.-1)) then
     ierr=1
     call eq_errmsg(' ?eqm_bset:  isnccwb must be +/- 1')
     write(lunerr,*) '  instead, isnccwb = ',isnccwb
  else
     nsnccwb=isnccwb
  endif

  nsnccwi=0
  if((isnccwi.ne.1).and.(isnccwi.ne.-1)) then
     ierr=1
     call eq_errmsg(' ?eqm_bset:  isnccwi must be +/- 1')
     write(lunerr,*) '  instead, isnccwi = ',isnccwi
  else
     nsnccwi=isnccwi
  endif

  call xplasma_field_ccw_set(s,ierr, bphi_ccw=nsnccwb, jphi_ccw=nsnccwi)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_bset: xplasma_field_ccw_set returned status: ',ierr
     call xplasma_error(s,ierr,lunerr)
  else
     call xplasma_eqcheck(s,iforce,icount,ierr)
     if(icount.gt.0) then
        write(lunerr,*) ' %eqm_bset: equilibrium specification completed.'
     endif
  endif

end subroutine eqm_bset

subroutine eq_bchk_sign(isnccwb,isnccwi)

  ! retrieve field and current direction signs e.g. as set by eqm_bset...

  use xplasma_obj_instance
  use eq_module

  !--------------
  integer :: ierr
  !--------------

  call xplasma_global_info(s,ierr, bphi_ccw=isnccwb, jphi_ccw=isnccwi)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_bchk_sign: xplasma_global_info status: ',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_bchk_sign
