subroutine pstCheckMetric(i, j, nt1, ns, t, p)
  USE pstcom
  USE newmet
  USE mtrik1
  
  implicit none
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
  integer, intent(in) :: i ! poloidal ray index
  integer, intent(in) :: j ! radial index
  integer, intent(in) :: nt1 ! max poloidal index
  integer, intent(in) :: ns ! max radial index
  real(r8), intent(in) :: t(nt1), p(ns) ! the grids

  real(r8)  :: dp, dt

  if(i<=1 .or. i>=nt1 .or. j<=1 .or. j>=ns) then
     print *,'error: indices must lie between [2,nt1-1] and [2, ns-1]'
     return
  end if

  dp = p(j+1) - p(j-1)
  dt = t(i+1) - t(i-1)

  write(*,'(a,2i3   )') 'INDICES   ', i, j
  write(*,'(a,2f10.6)') 'POSITIONS ', t(i), p(j)
  write(*,*) 'Check PP'
  write(*,'(2e17.9)') ppa(j), (pa(j+1)-pa(j-1))/dp
  write(*,*) 'Check QP'
  write(*,'(2e17.9)') qpa(j), (qa(j+1)-qa(j-1))/dp
  write(*,*) 'Check GP'
  write(*,'(2e17.9)') gpa(j), (ga(j+1)-ga(j-1))/dp

  write(*,*) 'Check xdth'
  write(*,'(2e17.9)') xdth(i,j), (xa(i+1,j  )-xa(i-1,j  ))/dt
  write(*,*) 'Check xdps'
  write(*,'(2e17.9)') xdps(i,j), (xa(i  ,j+1)-xa(i  ,j-1))/dp
  write(*,*) 'Check zdth'
  write(*,'(2e17.9)') zdth(i,j), (za(i+1,j  )-za(i-1,j  ))/dt
  write(*,*) 'Check zdps'
  write(*,'(2e17.9)') zdps(i,j), (za(i  ,j+1)-za(i  ,j-1))/dp
  
  write(*,*) 'Check xsq'
  write(*,'(2e17.9)') xsq(i,j), xa(i,j)**2
  write(*,*) 'Check xjacob'
  write(*,'(2e17.9)') xjacob(i,j), xa(i,j)*( &
       & (xa(i+1,j  )-xa(i-1,j  ))*(za(i  ,j+1)-za(i  ,j-1)) - &
       & (za(i+1,j  )-za(i-1,j  ))*(xa(i  ,j+1)-xa(i  ,j-1)) &
       & )/(dt*dp)
  write(*,*) 'Check xjprym'
  write(*,'(2e17.9)') xjprym(i,j), (xjacob(i  ,j+1)-xjacob(i  ,j-1))/dp
  write(*,*) 'grpssq'
  write(*,'(e17.9)') grpssq(i,j)
  write(*,*) 'grpsth'
  write(*,'(e17.9)') grpsth(i,j)
  write(*,*) 'delta'
  write(*,'(e17.9)') delta(i,j)
  write(*,*) 'qdelp'
  write(*,'(e17.9)') qdelp(i,j)
  
  
end subroutine pstCheckMetric

  
