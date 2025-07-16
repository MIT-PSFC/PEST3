subroutine eq_flxint_init(iauto,inum,zrho,zrhomin,ierr)

  !  fortran-77 interface -- initialize flux surface numerical integration
 
  use xplasma_obj_instance
  use eq_module

  implicit NONE
 
  !   initialize computation of flux surface / flux zone averages / integrals
  !   using xplasma f90 routines...

  ! input:

  integer iauto            ! =1 for automatic generation of rho partitioning
  integer inum             ! number of rho surfaces partitioning region

  ! input/output:
  !   (output if iauto=1, input otherwise)

  real*8 zrho(inum)        ! sequence of rho surfaces (strict ascending)
  real*8 zrhomin           ! minimum safe rho (approach to axial singularity)

  ! output:

  integer ierr             ! completion code, 0=OK

  !-----------------------------------
  !  note (dmc) -- may need to upgrade to support some integrals in 
  !  extrapolated region beyond plasma boundary...
  !-----------------------------------
  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE  = 1.0d0
  real*8, parameter :: CEPS7 = 1.0d-7
  real*8 :: zrhominu

  integer :: i,idi,iertmp
  character*32 :: iiname,authname
  !-----------------------------------
  !  error check...

  if(inum.lt.2) then
     write(lunerr,*) ' ?eq_flxint_init:  #surfaces, inum.ge.2 required!'
     ierr=1
     return
  endif

  !  iauto option...

  if(iauto.eq.1) then
     do i=1,inum
        zrho(i)=(i-1)*ONE/(inum-1)
     enddo
     zrhominu=min(CEPS4*(zrho(2)),CEPS7)

  else
     zrhominu=max(CEPS4*CEPS7,min(10*CEPS7,zrhomin))  ! btw 1.e-6 and 1.e-11
  endif

  !  create integrator information...

  iiname = "__F77_INTEGRATOR_1D"
  call xplasma_blkbxId(s,iiname,idi)
  if(idi.ne.0) then
     call xplasma_get_item_info(s,idi,ierr, author=authname)
     if(authname.ne.xplasma_xmhd) then
        write(lunerr,*) ' %eq_flxint_init: old integrator owned by "', &
             trim(authname),'" removed (OK).'
        call xplasma_author_set(s,authname,ierr)
        call xplasma_remove_item(s,idi,ierr)
        call xplasma_author_clear(s,authname,ierr)
     endif
  endif

  call xplasma_author_set(s,xplasma_xmhd,ierr)
  if(ierr.ne.0) return

  call xplasma_create_integ(s, iiname, zrho, idi, ierr, &
       rhomin=zrhominu, cache_enable = .TRUE.)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_flxint_init-- integration setup failed:'
     call xplasma_error(s,ierr,lunerr)
  endif
 
  call xplasma_author_clear(s,xplasma_xmhd,iertmp)

end subroutine eq_flxint_init

subroutine eq_flxint_chinit(iauto,inumchis,zchi,ierr)

  !  fortran-77 xplasma interface:
  !  supplemental routine:  initialize chi break points for integrations
  !  over chi intervals

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  integer iauto                 ! =1 for automatic generation of zchi values
  integer inumchis              ! no. of chi interval boundaries (at least 3).
  real*8 zchi(inumchis)         ! chi intervals (output if iauto=1).

  integer ierr                  ! completion code:  0 = OK

  !---------------------------------------------------------------------
  integer :: i,idi,inrhoi,inthi,idi_2d,iertmp
  !---------------------------------------------------------------------
  !  error check...

  if(inumchis.lt.2) then
     write(lunerr,*) ' ?eq_flxint_chinit:  #surfaces, inumchis.ge.2 required!'
     ierr=1
     return
  endif

  ! find 1d integrator dataset (which will be extended to form 2d dataset).

  call eq_find_intg(1,idi,inrhoi,inthi,ierr)
  if((ierr.ne.0).or.(idi.eq.0)) then
     write(lunerr,*) ' ?eq_flxint_chinit: no 1d integrator dataset found.'
     write(lunerr,*) '  Perhaps there was no prior eq_flxint_init call.'
     write(lunerr,*) '  More information:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  !  iauto option...

  if(iauto.eq.1) then
     do i=1,inumchis
        zchi(i)=(i-1)*C2PI/(inumchis-1)
     enddo
  endif
 
  ! OK...
  !  create integrator information...

  call xoi_author_set(ierr)
  if(ierr.ne.0) return
 
  call xplasma_augment_integ(s,"__F77_INTEGRATOR_2D",idi,zchi,idi_2d,ierr)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_flxint_chinit-- integration setup failed:'
     call xplasma_error(s,ierr,lunerr)
  endif
 
  call xoi_author_clear(iertmp)

end subroutine eq_flxint_chinit
