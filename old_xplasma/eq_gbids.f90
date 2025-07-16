subroutine eq_gbids(id_g,id_R,id_BR,id_BZ,id_Bmod)

  !  get the field related Xplasma profile IDs
  !   g(rho)  BR(theta,rho)  BZ(theta,rho)  Bmod(theta,rho)

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer, intent(out) :: id_g,id_R,id_BR,id_BZ,id_Bmod

  !-----------------
  integer :: iertmp
  !-----------------

  call xplasma_common_ids(s,iertmp, &
       id_g=id_g, id_R=id_R, id_BR=id_BR, id_BZ=id_BZ, id_Bmod=id_Bmod)

  if(iertmp.ne.0) then
     write(lunerr,*) ' ?eq_gbids -- failed to get B field IDs'
     call xplasma_error(s,iertmp,lunerr)

     id_g = 0
     id_R = 0
     id_BR = 0
     id_BZ = 0
     id_Bmod = 0
  endif

end subroutine eq_gbids
