subroutine eq_bdlims(itype,zRminLim,zRmaxLim,zZminLim,zZmaxLim,ierr)

  !  fetch limiter type and extrema of vacuum region enclosed by
  !  limiter

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  integer, intent(out) :: itype  ! type: 100=contour
  !  101 = lines & circles with possible plasma bdy distance constraint

  real*8, intent(out) :: zRminLim,zRmaxLim  ! R range of vacuum region
  real*8, intent(out) :: zZminLim,zZmaxLim  ! Z range of vacuum region

  integer, intent(out) :: ierr ! status code (0=OK)

  !-----------------------------------------

  call xplasma_lim_RZminmax(s,zRminLim,zRmaxLim,zZminLim,zZmaxLim,ierr, &
       itype=itype)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_bdlims: xplasma_lim_RZminmax failed, ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_bdlims
