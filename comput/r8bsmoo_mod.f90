module r8bsmoo_mod
 
  SAVE
 
  !  mock up the part of TRCOM used by bsmoo.for,bsmooh.for,bsmoox.for
  !  which are being made independent of TRCOM
 
  integer mj                  ! grid array size
  integer lcentr,ledge        ! active grid limits
  integer lcp1,lep1           ! as above, +1
  integer nzones              ! no. of active zones
 
  !  the normalized grid
 
  real*8, dimension(:,:), allocatable :: xi
 
  !  the smoothing parameter
 
  real*8 dxbsmoo
 
end module r8bsmoo_mod
