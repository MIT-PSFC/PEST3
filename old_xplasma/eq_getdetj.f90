subroutine eq_getdetj(ivec,zrho,zchi,zphi,zscale,detj,ierr)

  !  get determinant of Jacobian (metric) without
  !  returning the tensor itself
  !    Note: if calling code uses a clockwise oriented theta, pass ivec < 0
  !          if calling code uses a counter-clockwise oriented theta, 
  !              pass ivec > 0; abs(ivec) = vector size.

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  input:
  integer ivec                      ! vector dimensioning
  real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi)
  real*8 zscale                     ! scale factor to apply to det(J)

  !  output:
  real*8 detj(abs(ivec))                 ! determinant of jtens

  integer ierr                      ! error code, 0=OK
 
  !--------------------------------
 
  logical :: ccwflag,axisymm
  integer :: jvec
  real*8 :: zr(abs(ivec)),zdetj(abs(ivec))
 
  !--------------------------------
 
  jvec=abs(ivec)
  ccwflag = (ivec.gt.0)

  call xplasma_global_info(s,ierr, axisymm=axisymm)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?error detected in eq_getdetj:'
     call xplasma_error(s,ierr,lunerr)
  else if(.not.axisymm) then
     ierr=1
     write(lunerr,*) ' ?eq_getdetj: only available for axisymmetric geometry!'
  endif
  if(ierr.eq.0) then

     call xplasma_rzjac(s,zrho,zchi,ierr, &
          ccwflag=ccwflag, &
          r=zr,rzdetj=zdetj)

     if(ierr.eq.0) then

        detj=zr*zdetj*zscale

     else
        write(lunerr,*) ' ?error detected in eq_getdetj xplasm_rzjac call:'
        call xplasma_error(s,ierr,lunerr)
     endif
  endif

  if(ierr.ne.0) detj = 0.0d0

end subroutine eq_getdetj
