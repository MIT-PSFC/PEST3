subroutine eqm_cbdy_grid(ipts,rlim,zlim,inumRi,inumZi,id_R,id_Z,ierr)

  !  (a)
  !  set up an axisymmetric piecewise linear boundary / limiter contour.
  !  the contour must close on itself (rlim(1)=rlim(ipts) &
  !   zlim(1)=zlim(ipts)) and not cross itself.

  !  (b)
  !  set up a rectangular (R,Z) grid, evenly spaced and with the dimensioning
  !  given, which covers the plasma core and "scrapeoff layer" i.e. the space
  !  between the core plasma and the limiters.  The grid is "padded" with
  !  1.5 zones around each edge.  The eqm_rzgrid routine is used to actually
  !  set the grid.

  !  the padding serves two functions:  (1) provides some mapped space
  !  in all directions beyond the limiters, and (2) ameliorates possible
  !  boundary effects in Akima Hermite bicubic fitting over the (R,Z) grid.

  !  ...for more details on the contours, see subroutine eqm_cbdy

  !  mod DMC Jan 2007 -- force Zmin = -Zmax
  !    some physics modules require this Z grid symmetry
  !-----------------------

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  !  input:

  integer ipts                      ! # of (R,Z) points defining contour
  real*8 rlim(ipts),zlim(ipts)      ! (R,Z) points

  integer inumRi                    ! desired no. of points in R grid
  integer inumZi                    ! desired no. of points in Z grid

  !  output:

  integer id_R                      ! R grid xplasma id code
  integer id_Z                      ! Z grid xplasma id code

  integer ierr                      ! completion code, 0=OK

  !-----------------------
  !  local...

  real*8, dimension(:), allocatable :: zRgrid,zZgrid  ! generated R & Z grids

  real*8 zRmin,zRmax,zZmin,zZmax,zRlim
  integer itype

  real*8 zRmid,zZmid
  real*8 zdelR,zdelZ

  integer i,inumR,inumZ

  !-----------------------

  inumR = inumRi
  inumZ = inumZi

  ierr=0
  if(inumR.lt.32) then
     write(lunerr,*) ' ?eqm_tbdy_grid: inumR=',inumR,' .lt. 32: fixup applied.'
     inumR=32
  endif
  if(inumZ.lt.32) then
     write(lunerr,*) ' ?eqm_tbdy_grid: inumZ=',inumZ,' .lt. 32: fixup applied.'
     inumZ=32
  endif
  if(ierr.ne.0) return

  !  1.  set up limiters

  call eqm_cbdy(ipts,rlim,zlim,ierr)
  if(ierr.ne.0) return

  !  2.  fetch the limiter extrema

  call eq_bdlims(itype,zRmin,zRmax,zZmin,zZmax,ierr)

  ! force Z grid up-down symmetry
  zZmax = max(abs(zZmin),abs(zZmax))
  zZmin = -zZmax

  zRmid=chalf*(zRmin+zRmax)
  zZmid=chalf*(zZmin+zZmax)

  !  2a: expand slightly, but check...

  zdelR=(inumR-1)*chalf*(zRmax-zRmin)/(inumR-4)
  zdelZ=(inumZ-1)*chalf*(zZmax-zZmin)/(inumZ-4)

  zRmin=max((zRmin/2),zRmid-zdelR)
  zRmax=zRmid+zdelR

  zZmin=zZmid-zdelZ
  zZmax=zZmid+zdelZ

  !  3.  generate the R,Z grids

  allocate(zRgrid(inumR),zZgrid(inumZ))

  do i=1,inumR
     zRgrid(i)=(zRmin*(inumR-i)+zRmax*(i-1))/(inumR-1)
  enddo

  do i=1,inumZ
     zZgrid(i)=(zZmin*(inumZ-i)+zZmax*(i-1))/(inumZ-1)
  enddo

  call eqm_rzgrid(zRgrid,zZgrid,0,0,inumR,inumZ,ceps4,id_R,id_Z,ierr)

  deallocate(zRgrid,zZgrid)

end subroutine eqm_cbdy_grid
