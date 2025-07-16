!-----------------------------------------------------------------------
subroutine eq_flxint(name,noption,result,inumchi,inumrho,ierr)

  !  return results of numerical integrations
  !  fortran-77 interface

  !  ported to fortran-90 rewritten XPLASMA dmc March 2006

  !  the integration option "noption" no longer has any effect but is
  !  left in place for backwards compatibility reasons.

  use xplasma_obj_instance
  use eq_module
  use eqi_intg_module

  implicit NONE

  !input:

  character*(*) name           ! name of item wanted
  integer noption              ! integration method option (*ignored*)
  integer inumchi              ! 1st dimension, chi dimension of output array
  integer inumrho              ! 2nd dimension, rho dimension of output array

  !output:

  real*8 result(inumchi,inumrho) ! result of integrations

  integer ierr                 ! completion code (0=OK)

  !-------------------------------------------------
  integer :: idi,inum,inrhoi,inthi,inthzons,inzons,inrho_need
  logical :: ilsurf,dvol_need,dvdrho_need,bmaxa_need
  character*32 zname
  !-------------------------------------------------
  result = 0

  call eq_find_intg(inumchi,idi,inrhoi,inthi,ierr)
  if(ierr.ne.0) return

  inthzons = inthi-1
  inzons = inrhoi-1

  zname=name
  call uupper(zname)

  sp => s

  call xplasma_integ_num(zname,inum,ilsurf, &
       dvol_need,dvdrho_need,bmaxa_need)

  if(inum.eq.0) then
     write(lunerr,*) &
          '  unrecognized integration name:  "'//trim(zname)//'".'
     ierr=1
     return
  endif

  if((.not.ilsurf).or.(inum.eq.eqi_intg_dvol).or.(inum.eq.eqi_intg_darea)) then
     inrho_need = inzons
  else
     inrho_need = inzons + 1
  endif

  if(inumchi.lt.inthzons) then
     write(lunerr,*) ' ?eq_flxint: too few poloidal zones in result array:'
     write(lunerr,*) '  received: ',inumchi,'; need at least: ',inthzons
     ierr = 1
  endif

  if(inumrho.lt.inrho_need) then
     write(lunerr,*) ' ?eq_flxint: too few radial zones in result array:'
     write(lunerr,*) '  received: ',inumrho,'; need at least: ',inrho_need
     ierr = 1
  endif

  !---------------------
  !  done with error checking; do integration

  if(inumchi.eq.1) then
     call xplasma_rho_zonint(s,idi,zname,result(1,1:inrho_need),ierr)

  else
     call xplasma_2d_zonint(s,idi,zname,result(1:inthzons,1:inrho_need),ierr)

  endif

end subroutine eq_flxint
