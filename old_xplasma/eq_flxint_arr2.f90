subroutine eq_flxint_arr2(zdata,iwant,result,inumchi,inumrho,ierr)

  ! integration of user defined data -- legacy f77 interface

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE
  
  !  given zdata(1:nchi,1:nrho) represent f(chi,rho), use eq_flxint_arr2s
  !  to (a) setup a spline representation, and (b) do the chosen set of
  !  numeric integrations using this representation.

  !  "not a knot" boundary condition will be applied to the spline.
  !  for more control of the boundary condition see eq_flxint_arr2s.

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

  real*8 zdata(*)

  integer iwant                ! type of integrations wanted

  integer inumchi              ! 1st (chi) dim of result
  integer inumrho              ! 2nd (rho) dim of result

  !  if inumchi=1, do integrations over the entire flux zone / flux surface
  !  if inumchi>1, do integrations over chi segments.

  !  **output**

  real*8 result(inumchi,inumrho) ! results of integrations
  integer ierr                   ! completion code, 0=OK
  !                                  ! see description, above...

  !  this routine has a simplified argument list; eq_flxint_arr2s
  !  and eq_flxint_arr2d are more complex, with more control options.

  !  **note** the number of values returned by result is set by a prior
  !  call to eq_flxint_init -- #integration zones if iwant.gt.1;
  !  #integration zones + 1 if iwant.eq.1
  !-------------------------------------
  real*8 :: zdum(1)
  !-------------------------------------

  zdum(1)=0.0d0

  call eq_flxint_arr2s(zdata, &
       0,zdum,0,zdum, &
       iwant,result,inumchi,inumrho,ierr)

  if(ierr.ne.0) then
     write(lunerr,*) ' eq_flxint_arr2: eq_flxint_arr2s returned ierr=',ierr
  endif

end subroutine eq_flxint_arr2
