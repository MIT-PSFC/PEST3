subroutine eqm_tbdy_grid(ilines,rl,zl,thl,icircs,rc,zc,rad, &
     inumRi,inumZi,id_R,id_Z,ierr)

  !  a)
  !  set up a TRANSP-style machine boundary specification.  This consists
  !  two lists:  a list of circles (R,Z) & radius, and a list of infinite
  !  lines (R,Z) and inclination angle, thl (degrees).

  !  ...see eqm_tbdy subroutine for details on the limiter definition

  !  b)
  !  set up a rectangular (R,Z) grid, evenly spaced and with the dimensioning
  !  given, which covers the plasma core and "scrapeoff layer" i.e. the space
  !  between the core plasma and the limiters.  The grid is "padded" with
  !  1.5 zones around each edge.  The eqm_rzgrid routine is used to actually
  !  set the grid.

  !  the padding serves two functions:  (1) provides some mapped space
  !  in all directions beyond the limiters, and (2) ameliorates possible
  !  boundary effects in Akima Hermite bicubic fitting over the (R,Z) grid.

  !  mod DMC Jan 2007 -- force Zmin = -Zmax
  !    some physics modules require this Z grid symmetry

  !----------------

  !  details of the circle and lines limiters:  see eqm_tbdy subroutine.
  !  the limiters define an excluded space -- regions behind limiters
  !  with respect to the plasma magnetic axis.  The plasma boundary must
  !  not intersect the excluded space.

  !-----------------------

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  !  input:

  integer ilines                    ! # of line limiters (can be zero)
  real*8 rl(ilines),zl(ilines)      ! (R,Z) starting point of line(s)
  real*8 thl(ilines)                ! in DEGREES, line orientation

  integer icircs                    ! # of circle limiters (can be zero)
  real*8 rc(icircs),zc(icircs)      ! (R,Z) center of circle(s)
  real*8 rad(icircs)                ! radius of circle(s)

  integer inumRi                    ! size of R grid -- at least 32 pts
  integer inumZi                    ! size of Z grid -- at least 32 pts

  !  output:

  integer id_R                      ! R grid xplasma id code
  integer id_Z                      ! Z grid xplasma id code

  integer ierr                      ! completion code, 0=OK

  !-----------------------
  !  local...

  real*8, dimension(:), allocatable :: zRgrid,zZgrid ! evenly spaced R&Z grids

  real*8 zRmin,zRmax,zZmin,zZmax
  integer itype

  real*8 zRmid,zZmid,zRlim
  real*8 zdelR,zdelZ

  integer i,inumR,inumZ

  !-----------------------

  ierr=0
  inumR=inumRi
  inumZ=inumZi

  if(inumR.lt.32) then
     write(lunerr,*) ' ?eqm_tbdy_grid:  inumR=',inumR,' .lt. 32: fixed.'
     inumR=32
  endif
  if(inumZ.lt.32) then
     write(lunerr,*) ' ?eqm_tbdy_grid:  inumZ=',inumZ,' .lt. 32: fixed.'
     inumZ=32
  endif
  if(ierr.ne.0) return

  !  1.  set up limiters

  call eqm_tbdy(ilines,rl,zl,thl,icircs,rc,zc,rad,ierr)
  if(ierr.ne.0) return

  !  2.  fetch the limiter extrema

  call eq_bdlims(itype,zRmin,zRmax,zZmin,zZmax,ierr)

  ! force Z grid up-down symmetry
  zZmax = max(abs(zZmin),abs(zZmax))
  zZmin = -zZmax

  zRmid=chalf*(zRmin+zRmax)
  zZmid=chalf*(zZmin+zZmax)

  zdelR=(inumR-1)*chalf*(zRmax-zRmin)/(inumR-4)
  zdelZ=(inumZ-1)*chalf*(zZmax-zZmin)/(inumZ-4)
  
  !  3.  generate the R,Z grids

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

end subroutine eqm_tbdy_grid
