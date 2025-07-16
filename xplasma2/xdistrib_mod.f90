module xdistrib_mod
 
  !  module for xdistrib.f90 -- share information with function passed
  !  to root finder
 
  integer, dimension(-1:3), parameter :: ncdim = (/ 1, 1, 2, 2, 4 /)
 
  !  to aid interpolation evaluations...
  real*8, dimension(:,:), allocatable :: xpkg
 
  !  interpolating function coefficients array...
  real*8, dimension(:,:), allocatable :: splcoef
 
  !  table of integral results at node points...
  real*8, dimension(:), allocatable :: ftabl
 
  integer nx_spl,iorder_spl   ! no. of pts in spline, & order (type of spline)
 
end module xdistrib_mod
