subroutine eqm_dbdy_grid(zdist,zRmini,zRmaxi,zZmini,zZmaxi, &
     inumRi,inumZi,id_R,id_Z,ierr)
  !
  !  f77 xplasma legacy interface:
  !  set up limiter which is a combination of a box and a fixed distance
  !  set back from the plasma boundary (eqm_dbdy); then set up evenly spaced
  !  grid to cover vacuum region (eqm_rzgrid).

  !  mod DMC Jan 2007 -- force Zmin = -Zmax
  !    some physics modules require this Z grid symmetry

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  !  set up evenly spaced (R,Z) grid of specified size, and set up
  !  a limiter which a fixed distance (zdist) from the plasma, except
  !  where a box [zRmini,zRmaxi]x[zZmini,zZmaxi] cuts in

  !  input:

  real*8 zdist           ! fixed distance of limiter from plasma edge
  real*8 zRmini          ! Rmin of enclosing box
  real*8 zRmaxi          ! Rmax of enclosing box
  real*8 zZmini          ! Zmin of enclosing box
  real*8 zZmaxi          ! Zmax of enclosing box

  integer inumRi         ! size of R grid to generate (minimum: 32)
  integer inumZi         ! size of Z grid to generate (minimum: 32)

  !  output:

  integer id_R                      ! R grid id
  integer id_Z                      ! Z grid id

  integer ierr           ! completion code, 0 = OK

  !--------------------------------------
  !  local:
  real*8, dimension(:), allocatable :: zRgrid   ! R grid (generated)
  real*8, dimension(:), allocatable :: zZgrid   ! Z grid (generated)

  real*8 zRmin,zRmax,zZmin,zZmax,zRlim
  real*8 zmidR,zmidZ,zdelR,zdelZ
  integer i,inumR,inumZ

  real*8, parameter :: zrho_bdy = 1.0d0

  !--------------------------------------

  ierr=0

  inumR=inumRi
  inumZ=inumZi

  if(inumR.lt.32) then
     write(lunerr,*) ' %eqm_dbdy_grid:  inumR = ',inumR,' too small.'
     write(lunerr,*) '  minimum grid size allowed:  32; fixup applied.'
     inumR=32
  endif
  if(inumZ.lt.32) then
     write(lunerr,*) ' %eqm_dbdy_grid:  inumZ = ',inumZ,' too small.'
     write(lunerr,*) '  minimum grid size allowed:  32; fixup applied.'
     inumZ=32
  endif

  if(ierr.ne.0) return

  call eq_glimrz(zrho_bdy,zRmin,zRmax,zZmin,zZmax,ierr)
  if(ierr.ne.0) then
     call eq_errmsg(' ?eqm_dbdy_grid:  eq_glimRZ error occurred.')
     return
  endif

  zRlim=0.5*zRmin
  zRmin=max(zRlim,min(zRmin,zRmini)) ! btw .5 & 1.0*Rmin(plasma bdy)
  zRmax=max(zRmax,zRmaxi)
  zZmin=min(zZmin,zZmini)
  zZmax=max(zZmax,zZmaxi)

  ! force Z grid up-down symmetry
  zZmax = max(abs(zZmin),abs(zZmax))
  zZmin = -zZmax

  ! expand mapped region...

  zmidR=chalf*(zRmin+zRmax)
  zdelR=chalf*(inumR-1)*(zRmax-zRmin)/(inumR-4)
  zmidZ=chalf*(zZmin+zZmax)
  zdelZ=chalf*(inumZ-1)*(zZmax-zZmin)/(inumZ-4)

  allocate(zRgrid(inumR),zZgrid(inumZ))

  do i=1,inumR
     zRgrid(i)=zmidR-zdelR + (i-1)*ctwo*zdelR/(inumR-1)
  enddo

  do i=1,inumZ
     zZgrid(i)=zmidZ-zdelZ + (i-1)*ctwo*zdelZ/(inumZ-1)
  enddo

  call eqm_dbdy(zdist,zRmini,zRmaxi,zZmini,zZmaxi,ierr)

  if(ierr.eq.0) then
     call eqm_rzgrid(zRgrid,zZgrid,0,0,inumR,inumZ,ceps4,id_R,id_Z,ierr)
  endif

  deallocate(zRgrid,zZgrid)

end subroutine eqm_dbdy_grid
