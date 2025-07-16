subroutine eq_gaxis(ivec,phi,b_axis,R_axis,Z_axis)

  !  (vectorized)
  !  given a toroidal angle coordinate (phi), return B and (R,Z) at the
  !  magnetic axis.

  !  if the equilibrium is AXISYMMETRIC then the value of the input variable
  !  phi is immaterial

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer ivec
  REAL*8 phi(ivec)                  ! (input) toroidal angle coordinate

  REAL*8 b_axis(ivec)               ! (output) B @ magnetic axis
  REAL*8 R_axis(ivec)               ! (output) R @ magnetic axis
  REAL*8 Z_axis(ivec)               ! (output) Z @ magnetic axis

  !-------------------------
  integer :: ier
  logical :: axisymm
  integer :: iv
  !-------------------------

  if(ivec.le.0) return

  call xplasma_global_info(s,ier,axisymm=axisymm)
  if(ier.ne.0) then
     write(lunerr,*) ' ?eq_gaxis: xplasma error:'
     call xplasma_error(s,ier,lunerr)
     call all_zero

  else if(.not.axisymm) then
     write(lunerr,*) ' ?eq_gaxis: non-axisymmetry not yet supported.'
     call all_zero

  else

     call eq_sgaxis(b_axis(1),r_axis(1),z_axis(1))
     do iv=2,ivec
        b_axis(iv)=b_axis(1)
        r_axis(iv)=r_axis(1)
        z_axis(iv)=z_axis(1)
     enddo
  endif

  contains

    subroutine all_zero

      b_axis = 0
      r_axis = 0
      z_axis = 0

    end subroutine all_zero
end subroutine eq_gaxis

subroutine eq_sgaxis(zB,zR,zZ)

  !  get mod(B), R, and Z on the magnetic axis -- axisymmetry expected.

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  real*8, intent(out) :: zB,zR,zZ  ! axial values of mod(B), R, Z

  !-------------------------
  integer :: ier
  logical :: axisymm
  !-------------------------

  call xplasma_global_info(s,ier,axisymm=axisymm)
  if(.not.axisymm) then
     write(lunerr,*) ' ?eq_sgaxis: non-axisymmetry, axial values of mod(B)', &
          ' R and Z not unique.'
     zB=0
     zR=0
     zZ=0
  else

     call xplasma_mag_axis(s,ier, raxis=zR,zaxis=zZ,baxis=zB)
     if(ier.ne.0) then
        write(lunerr,*) ' ?eq_sgaxis: xplasma error:'
        call xplasma_error(s,ier,lunerr)
        zB=0
        zR=0
        zZ=0
     endif

  endif

end subroutine eq_sgaxis
