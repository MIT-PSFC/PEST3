subroutine pstIntegral(nt1, zx, zy, zyint)
  !
  ! given the set of independent zx and dependent zy points, compute
  ! the integrals of zy from zx(1) to zx(i). zy must be periodic
  ! ie zy(1) = zy(nt1)
  !
  ! A. Pletzer April 2000
  !
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  integer, intent(in) :: nt1
  real(r8), intent(in) :: zx(*), zy(*)
  real(r8), intent(out) :: zyint(*) ! the result

  integer i, im1, j, ier, nwk, ibcxmin, ibcxmax, ilinx 
  real(r8) :: bcxmin, bcxmax, dx
  real(r8), dimension(:), allocatable :: wk
  real(r8), dimension(:,:), allocatable :: fspl


  ! periodic boundary conditions
  ibcxmin=-1
  ibcxmax=-1
  nwk = 5*nt1

  allocate( wk(nwk), fspl(4,nt1))

  fspl(1,1:nt1) = zy(1:nt1)
  call r8cspline(zx, nt1, fspl, ibcxmin,bcxmin,ibcxmax,bcxmax, &
       &   wk, nwk, ilinx, ier)

  zyint(1) = 0._r8
  do i=2, nt1
     im1 = i - 1
     dx = zx(i) - zx(im1)
     zyint(i) = zyint(im1) + &
          & dx*( fspl(1,im1) + &
          & dx*( fspl(2,im1)/2.0_r8 + &
          & dx*( fspl(3,im1)/3.0_r8 + &
          & dx*  fspl(4,im1)/4.0_r8   &
          &    ) ) )
  enddo
  deallocate(wk, fspl)   

end subroutine pstIntegral

