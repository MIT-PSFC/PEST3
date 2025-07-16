subroutine eqm_chi_cwdir_set(iswitch,idstring)

  implicit NONE
  integer,intent(in) :: iswitch         ! =0:  CCW (normal), =1:  CW (inverted)
  character*(*), intent(in) :: idstring ! caller's id string

  write(6,*) ' Detected call to defunct xplasma routine "eqm_chi_cwdir_set"'
  write(6,*) '    with idstring= "',trim(idstring),'" and'
  write(6,*) '    with iswitch = ',iswitch,'; NO ACTION TAKEN.'

end subroutine eqm_chi_cwdir_set


subroutine eqm_chi_cwdir_restore(idstring,ierr)
 
  implicit NONE
  character*(*), intent(in) :: idstring ! caller's id string
  integer,intent(out) :: ierr           ! =0:  normal

  write(6,*) ' Detected call to defunct xplasma routine "eqm_chi_cwdir_restore"'
  write(6,*) '    with idstring= "',trim(idstring),'"; NO ACTION TAKEN.'
  ierr=9999

end subroutine eqm_chi_cwdir_restore


subroutine eqm_chi_cwdir(iswitch)
 
  implicit NONE
  integer iswitch         ! =0:  CCW (normal), =1:  CW (inverted)

  write(6,*) ' Detected call to defunct xplasma routine "eqm_chi_cwdir"'
  write(6,*) '   with iswitch = ',iswitch,'; NO ACTION TAKEN.'

end subroutine eqm_chi_cwdir


subroutine eq_get_chi_cwdir(iswitch_value)

  implicit NONE
  integer, intent(out) :: iswitch_value

  iswitch_value = 0
  write(6,*) ' Detected call to defunct xplasma routine "eq_get_chi_cwdir"!'

end subroutine eq_get_chi_cwdir

!-------------------------------------------------

subroutine eqm_phi_cwdir_set(iswitch,idstring)

  implicit NONE
  integer,intent(in) :: iswitch         ! =0:  CCW (normal), =1:  CW (inverted)
  character*(*), intent(in) :: idstring ! caller's id string

  write(6,*) ' Detected call to defunct xplasma routine "eqm_phi_cwdir_set"'
  write(6,*) '    with idstring= "',trim(idstring),'" and'
  write(6,*) '    with iswitch = ',iswitch,'; NO ACTION TAKEN.'

end subroutine eqm_phi_cwdir_set


subroutine eqm_phi_cwdir_restore(idstring,ierr)
 
  implicit NONE
  character*(*), intent(in) :: idstring ! caller's id string
  integer,intent(out) :: ierr           ! =0:  normal

  write(6,*) ' Detected call to defunct xplasma routine "eqm_phi_cwdir_restore"'
  write(6,*) '    with idstring= "',trim(idstring),'"; NO ACTION TAKEN.'
  ierr=9999

end subroutine eqm_phi_cwdir_restore


subroutine eqm_phi_cwdir(iswitch)
 
  implicit NONE
  integer iswitch         ! =0:  CCW (normal), =1:  CW (inverted)

  write(6,*) ' Detected call to defunct xplasma routine "eqm_phi_cwdir"'
  write(6,*) '   with iswitch = ',iswitch,'; NO ACTION TAKEN.'

end subroutine eqm_phi_cwdir


subroutine eq_get_phi_cwdir(iswitch_value)

  implicit NONE
  integer, intent(out) :: iswitch_value

  iswitch_value = 0
  write(6,*) ' Detected call to defunct xplasma routine "eq_get_phi_cwdir"!'

end subroutine eq_get_phi_cwdir

!-------------------------------------------------

subroutine eqm_flxbdy(iwarn)

  implicit NONE
  integer, intent(out) :: iwarn

  iwarn=9999

  write(6,*) ' Detected call to defunct xplasma routine "eqm_flxbdy".'

end subroutine eqm_flxbdy


subroutine eqm_rbdy(iwarn)

  implicit NONE
  integer, intent(out) :: iwarn

  iwarn=9999

  write(6,*) ' Detected call to defunct xplasma routine "eqm_rbdy".'

end subroutine eqm_rbdy

!-------------------------------------------------

subroutine eqm_spordr_set(ival)

  implicit NONE

  integer, intent(in) :: ival

  write(6,*) ' Detected call to defunct xplasma routine "eqm_spordr_set"'
  write(6,*) '    with ival = ',ival,'; NO ACTION TAKEN.'

end subroutine eqm_spordr_set


subroutine eqm_spordr_get(ival)

  implicit NONE

  integer, intent(out) :: ival

  write(6,*) ' Detected call to defunct xplasma routine "eqm_spordr_get"'
  write(6,*) ' (All cubic splines are in compact form now).'

  ival=0

end subroutine eqm_spordr_get

!-------------------------------------------------

subroutine eqm_nfast(inR,inZ,ierr)

  implicit NONE

  integer inR,inZ                   ! desired grid size in R and Z
  integer ierr                      ! completion code, 0=OK

  !  the "nfast" map no longer exists.

  ierr=9999
  write(6,*) ' Detected call to defunct xplasma routine "eqm_nfast"'
  write(6,*) '    with inR,inZ = ',inR,inZ,'; NO ACTION TAKEN.'

end subroutine eqm_nfast
