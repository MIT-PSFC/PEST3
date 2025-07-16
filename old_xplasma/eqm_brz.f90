subroutine eqm_brz(userbvec,zsm,ierr)

  !  Legacy f77 Xplasma interface...

  !  set up Bphi(R,Z), BR(R,Z), and BZ(R,Z) -- Akima-Hermite interpolation
  !  for axisymmetric case, setup psi(pol) also.

  !    data values from user supplied function  of form:

  !        subroutine userbvec(iv,zR,zZ,zphi,init,BR,BZ,BPHI,aux,kaux,ierr)
  !          iv -- vector dimension
  !          zR(iv),zZ(iv),zphi(iv)   ! coordinates at which to evaluate
  !          init -- action code
  !          BR(iv),BZ(iv),BPHI(iv)   ! field components returned
  !          aux(iv)                  ! aux. values returned
  !          kaux                     ! =1 if aux. values were set
  !          ierr                     ! completion code, 0=OK

  !  unlike the eqm_rzfunf and eqm_rzfun2 "userfcn" routines, this routine
  !  is called and must return good values for all (Ri,Zj) on the (R,Z) grid
  !  (established by the user with and eqm_rzgrid call, earlier)--where (Ri,Zj)
  !  are beyond the plasma boundary; it may return Bphi=BR=BZ=Psi=0 for (Ri,Zj)
  !  which are inside the plasma.

  !     normal calls to userbvec use init=0
  !         returns BR,BZ,BPHI, and (possibly) psi(pol) in aux.
  !     an initialization call uses init=1, returns BR=BZ=BPHI=czero
  !         also returns naux=1 if a poloidal flux function is to be
  !         returned.
  !         **note** for axisymmetric equilibria, naux=1 & psi(pol.)
  !         is expected.

  !     a cleanup call uses init=2, returns BR=BZ=BPHI=czero

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  !  input:
  external userbvec                 ! user supplied Bvec subroutine
  real*8 zsm                        ! smoothing parameter for edge (m)

  !  output:
  integer ierr                      ! completion code, 0=OK

  !---------------------------------------------
  integer :: iertmp
  !---------------------------------------------

  call xplasma_brz(s,userbvec,ierr, sm_edge=zsm)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_brz: error in xplasma_brz:'
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eqm_brz
