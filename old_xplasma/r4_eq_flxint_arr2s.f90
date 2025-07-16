!-------------------------------------------------------------
 
subroutine r4_eq_flxint_arr2(zdata,iwant,result,inumchi,inumrho,ierr)

  ! integration of user defined data -- legacy f77 interface REAL*4

  use xplasma_obj_instance
  IMPLICIT NONE

  !  REAL (*4) precision vsn of eq_flxint_arr2

  !  **input**

  real zdata(*)

  integer iwant                ! type of integrations wanted

  integer inumchi              ! 1st (chi) dim of result
  integer inumrho              ! 2nd (rho) dim of result

  !  if inumchi=1, do integrations over the entire flux zone / flux surface
  !  if inumchi>1, do integrations over chi segments.

  !  **output**

  real result(inumchi,inumrho) ! results of integrations
  integer ierr                 ! completion code, 0=OK
  !                                  ! see description, above...

  !  this routine has a simplified argument list; eq_flxint_arr2s,
  !  eq_flxint_arr2d are more complex, with more control options.

  !  **note** the number of values returned by result is set by a prior
  !  call to eq_flxint_init -- #integration zones if iwant.gt.1;
  !  #integration zones + 1 if iwant.eq.1
  !-------------------------------------
  integer :: ilun
  real :: zdum = 0.0
  !-------------------------------------

  call r4_eq_flxint_arr2s(zdata, &
       0,zdum,0,zdum, &
       iwant,result,inumchi,inumrho,ierr)

  if(ierr.ne.0) then
     call eq_get_lunerr(ilun)
     write(ilun,*) ' r4_eq_flxint_arr2: r4_eq_flxint_arr2s returned ierr=',ierr
  endif

end subroutine r4_eq_flxint_arr2

 
subroutine r4_eq_flxint_arr2s(zdata, &
     ibcrho0,zbcrho0,ibcrho1,zbcrho1, &
     iwant,result,inumchi,inumrho,ierr)

  !  integration of user defined data with user specified boundary conditions
  !  legacy f77 xplasma interface REAL*4

  use xplasma_obj_instance
  IMPLICIT NONE

  !  given zdata(1:nchi,1:nrho) represent f(chi,rho), use eq_flxint_arr2d
  !  to (a) setup a spline representation, and (b) do the chosen set of
  !  numeric integrations using this representation.

  !  *** iwant *** specifies the type of integration desired:
  !    here "f" denotes the interpolating function generated from the
  !    passed array data.

  !    iwant=1:  surface integral:
  !       int[phi=0 to 2pi] dphi * {
  !            int[chi=0 to 2pi] dchi * f }

  !    iwant=5:  volume weighted surface integral:
  !        int[phi=0 to 2pi] dphi * {
  !            int[chi=0 to 2pi] dchi * f * det(Jacobian) }

  !      **caution** if rho(1) corresponds to the mag. axis, a
  !        det(Jacobian)=0 and a value of 0.0 is always returned!

  !    iwant=6:  volume weighted surface average:

  !        int[phi=0 to 2pi] dphi * {
  !            int[chi=0 to 2pi] dchi * f * det(Jacobian)}
  !        -----------------------------------------------
  !        int[phi=0 to 2pi] dphi * {
  !            int[chi=0 to 2pi] dchi * det(Jacobian)}

  !      **caution** if rho(1) corresponds to the mag. axis,
  !        the value "f" evaluated at the axis is returned there.

  !  zonal integrals:  result(1:nrho-1) are computed; result(j)
  !                    corresponds to the integral from rho(j) to rho(j+1)

  !    iwant=2:  integrate f*drho*dtheta*dphi:

  !       int[rho=rho(j) to rho(j+1)] drho * {
  !          int[phi=0 to 2pi] dphi * {
  !               int[chi=0 to 2pi] dchi * f } }

  !    iwant=3:  integrate f*dV (use det(Jacobian))

  !       int[rho=rho(j) to rho(j+1)] drho * {
  !          int[phi=0 to 2pi] dphi * {
  !               int[chi=0 to 2pi] dchi * f * det(Jacobian) } }

  !    iwant=4:  compute <f> = (integrated f*dV)/dVol
  !              where dVol is the volume btw rho(j) and rho(j+1)

  !  **caution** in axisymmetric case there is no variation with phi,
  !  but a factor of 2pi still arises from the phi integration...

  !  **input**

  real zdata(*)

  integer ibcrho0,ibcrho1      ! BC controls
  real zbcrho0(*),zbcrho1(*)   ! BC data values

  !  BC type values:  for ibcrho0:
  !     =0 -- use "not a knot", zbcrho0(...) ignored
  !     =1 -- match slope, given at x(1),chi(ichi) by zbcrho0(ichi)
  !     =2 -- match 2nd deriv., given at x(1),chi(ichi) by zbcrho0(ichi)
  !     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all chi(j)
  !           (zbcrho0(...) ignored).
  !  and similarly for (ibcrho1,zbcrho1(...)).

  integer iwant                ! type of integrations wanted

  integer inumchi              ! 1st (chi) dim of result
  integer inumrho              ! 2nd (rho) dim of result

  !  if inumchi=1, do integrations over the entire flux zone / flux surface
  !  if inumchi>1, do integrations over chi segments.

  !  **output**

  real result(inumchi,inumrho) ! results of integrations
  integer ierr                 ! completion code, 0=OK

  !  this routine has a simplified argument list; eq_flxint_arr2d
  !  is more complex, with more control options.

  !  **note** the number of values returned by result is set by a prior
  !  call to eq_flxint_init -- #integration zones if iwant.gt.1;
  !  #integration zones + 1 if iwant.eq.1
  !-------------------------------------
  integer idr,idbc0,idbc1
  integer :: nchi,nrho,iertmp
  integer :: id_R,id_rho,id_chi
  integer :: ilun
  !-------------------------------------

  call eq_get_lunerr(ilun)

  call xplasma_common_ids(s,ierr, id_R=id_R)
  if(ierr.ne.0) then
     write(ilun,*) ' ** error detected in eq_flxint_arr2s:'
     call xplasma_error(s,ierr,ilun)
     return
  else if(id_R.eq.0) then
     write(ilun,*) ' ** equilibrium R(x,theta) not found in ex_flxint_arr2s.'
     ierr=1
     return
  endif

  call xplasma_prof_info(s,id_R,iertmp, gridId1=id_chi, gridId2=id_rho)
  call xplasma_grid_size(s,id_chi,nchi,iertmp)
  call xplasma_grid_size(s,id_rho,nrho,iertmp)

  if((ibcrho0.eq.1).or.(ibcrho0.eq.2)) then
     idbc0=nchi
  else
     idbc0=1
  endif

  if((ibcrho1.eq.1).or.(ibcrho1.eq.2)) then
     idbc1=nchi
  else
     idbc1=1
  endif

  call r4_eq_flxint_arr2d(2,2,zdata,nchi,nrho, &
       ibcrho0,idbc0,zbcrho0,ibcrho1,idbc1,zbcrho1, &
       iwant,result,inumchi,inumrho,ierr)

  if(ierr.ne.0) then
     write(ilun,*) ' ** eq_flxint_arr2s: error in eq_flxint_arr2d: '
     call xplasma_error(s,ierr,ilun)
  endif

end subroutine r4_eq_flxint_arr2s
