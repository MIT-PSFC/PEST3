subroutine eqi_rzsmedg(R,inumR,Z,inumZ,fdata,zsm,ierr)

  !  given f(R,Z) data,
  !  smooth in a band of width (zsm) along the plasma boundary

  !  the smoothing of f(Rj,Zk) is a bi-triangular weighted convolution
  !  over Rj +/- zsmloc, Zk +/- zsmloc, where zsmloc=zsm if distance from
  !  bdy is .le. 0.5*zsm, and a value linearly declining to zero if distance
  !  from bdy is between 0.5*zsm and 1.5*zsm.

  implicit NONE

  integer inumR,inumZ               ! R,Z grid sizes
  real*8 R(inumR),Z(inumZ)          ! R,Z grids
  real*8 fdata(inumR,inumZ)         ! data to be smoothed
  real*8 zsm                        ! smoothing bandwidth

  integer ierr                      ! completion code returned, 0=OK

  !---------------------------------
  real*8, dimension(:,:), allocatable :: wk1,wk2
  real*8, parameter :: HALF = 0.5d0
  !---------------------------------

  ierr=0

  allocate(wk1(inumR,inumZ))
  allocate(wk2(inumR,inumZ))

  !  the bidirectional smooth with varying convolution width is not quite
  !  linear, so compute it both ways and average the result

  wk1=fdata
  wk2=fdata

  call eqi_rzsmedR(R,inumR,Z,inumZ,wk1,zsm)
  call eqi_rzsmedZ(R,inumR,Z,inumZ,wk1,zsm)

  call eqi_rzsmedZ(R,inumR,Z,inumZ,wk2,zsm)
  call eqi_rzsmedR(R,inumR,Z,inumZ,wk2,zsm)

  fdata=HALF*(wk1+wk2)

  deallocate(wk1,wk2)

end subroutine eqi_rzsmedg
