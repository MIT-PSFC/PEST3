subroutine eqm_rzmagb(Rarr,Zarr,id1,id2,idrho,  &
     ibcR0,zbcR0,ibcR1,zbcR1,  &
     ibcZ0,zbcZ0,ibcZ1,zbcZ1,  &
     id_R,id_Z,ierr)
 
  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  establish R(chi,rho) and Z(chi,rho) bicubic splines.
  !  specify spline boundary conditions @ rho(axis) & @ rho(bdy)
  !  independently for R, Z

  !  input arguments:

  integer id1,id2            ! R,Z array dimensions
  real*8 Rarr(id1,id2)       ! R array
  real*8 Zarr(id1,id2)       ! Z array

  integer idrho              ! =1:  id1 is "rho" dimension; =2: id2 is.

  !  ** boundary conditions **

  !  BC controls are as in "pspline":  0 -- default, 1 -- fix df/drho
  !                                    2 -- fix d2f/drho2  ... f is R or Z

  !  for BC types 1 & 2, data must be supplied:  nchi 1st or 2nd derivative
  !  values, repsectively;  for other types the "zbc" arrays are ignored.

  integer ibcR0              ! R BC type @ rho(1)
  real*8 zbcR0(*)            ! R BC data @ rho(1) (if needed)
  integer ibcR1              ! R BC type @ rho(nrho)
  real*8 zbcR1(*)            ! R BC data @ rho(nrho) (if needed)

  integer ibcZ0              ! Z BC type @ rho(1)
  real*8 zbcZ0(*)            ! Z BC data @ rho(1) (if needed)
  integer ibcZ1              ! Z BC type @ rho(nrho)
  real*8 zbcZ1(*)            ! Z BC data @ rho(nrho) (if needed)
 
  !  output arguments:

  integer id_R               ! xplasma id for "R" bicubic spline (returned)
  integer id_Z               ! xplasma id for "Z" bicubic spline (returned)

  integer ierr               ! completion code (0 = OK)

  !------------------------------------------------
  !  local:

  integer inx1,inx2,iertmp,icount
  integer inrho,inchi,id_rho,id_chi
  logical :: iforce=.FALSE.

  !------------------------------------------------

  id_R=0
  id_Z=0
  ierr=0

  call eqm_rzmag_check(id1,id2,idrho,id_rho,id_chi,inrho,inchi,ierr)
  if(ierr.ne.0) return

  if(ibcR0.ne.ibcZ0) then
     write(lunerr,*) &
          ' ?eqm_rzmagb: different rho=0 boundary condition for R & Z:'
     write(lunerr,*) &
          '  ibcR0=',ibcR0,' ibcZ0=',ibcZ0
     write(lunerr,*) &
          '  These must match, i.e. have the same value.'
     ierr=1
  endif

  if(ibcR1.ne.ibcZ1) then
     write(lunerr,*) &
          ' ?eqm_rzmagb: different rho=1 boundary condition for R & Z:'
     write(lunerr,*) &
          '  ibcR0=',ibcR0,' ibcZ0=',ibcZ0
     write(lunerr,*) &
          '  These must match, i.e. have the same value.'
     ierr=1
  endif

  if(ierr.ne.0) return

  call xoi_author_set(iertmp)
  
  if(id_rho.eq.1) then
     call xplasma_rzmag(s,id_rho,id_chi, &
          Rarr(1:inrho,1:inchi),Zarr(1:inrho,1:inchi), &
          id_R,id_Z,ierr, &
          ibc0=ibcR0,Rbc0=zbcR0(1:inchi),Zbc0=zbcZ0(1:inchi), &
          ibc1=ibcR1,Rbc1=zbcR1(1:inchi),Zbc1=zbcZ1(1:inchi))
  else
     call xplasma_rzmag(s,id_chi,id_rho, &
          Rarr(1:inchi,1:inrho),Zarr(1:inchi,1:inrho), &
          id_R,id_Z,ierr, &
          ibc0=ibcR0,Rbc0=zbcR0(1:inchi),Zbc0=zbcZ0(1:inchi), &
          ibc1=ibcR1,Rbc1=zbcR1(1:inchi),Zbc1=zbcZ1(1:inchi))
  endif

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_rzmagb: xplasma_rzmag returned ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  else
     call xplasma_eqcheck(s,iforce,icount,ierr)
     if(icount.gt.0) then
        write(lunerr,*) ' %eqm_rzmagb: equilibrium specification completed.'
     endif
  endif

  call xoi_author_clear(iertmp)

end subroutine eqm_rzmagb
