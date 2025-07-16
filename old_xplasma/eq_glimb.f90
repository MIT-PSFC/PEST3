subroutine eq_glimb(rho,b_min,b_max,ierr)

  ! find min & max mod(B) on specified flux surface...

  IMPLICIT NONE
  REAL*8 rho                        ! rho of surface (in)
  REAL*8 b_min                      ! minimum mod(B) found (out)
  REAL*8 b_max                      ! maximum mod(B) found (out)
  integer ierr                      ! completion code, 0=OK

  !---------------------------
  real*8 zdum1,zdum2
  real*8, parameter :: ZERO = 0.0d0
  !---------------------------

  call eq_glimb1(rho,ZERO,b_min,zdum1,b_max,zdum2,ierr)

end subroutine eq_glimb
