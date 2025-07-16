subroutine eqi_evenx(x,nx,xtol)

  ! if spacing of x grid varies by less than xtol from the averaged,
  ! rebuild the grid to be precisely evenly spaced to real*8 precision

  implicit NONE

  integer, intent(in) :: nx
  real*8 :: x(nx)  ! can be modified on exit
  real*8, intent(in) :: xtol

  !-----------------------------------
  ! local:

  real*8 :: dx_avg,zdx,ztol,ztest
  integer :: i,iredo

  !-----------------------------------
  !  no error messages, but do sanity checks...

  if(nx.le.1) return
  if(x(nx).le.x(1)) return

  ztol=min(1.d-4,max(1.d-12,xtol))

  dx_avg = (x(nx)-x(1))/(nx-1)

  iredo=1
  do i=1,nx-1
     zdx=x(i+1)-x(i)
     ztest = abs(zdx-dx_avg)/dx_avg
     if(ztest.gt.ztol) then
        iredo=0
        exit
     endif
  enddo

  if(iredo.eq.0) return

  do i=2,nx-1
     x(i)=x(1)+(i-1)*dx_avg
  enddo

end subroutine eqi_evenx
