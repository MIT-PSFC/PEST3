subroutine r8_qksmooth(n,x,f,delta,bc0,bc1,ier)
 
  ! use R8FILTR6 for quick smooth
  ! on output f is smoothed using moving triangular weighting function
  ! convolution; half width of weighting function is delta.
  ! x(1:n) is expected to be strict ascending.
 
  implicit NONE
 
  integer, intent(in) :: n     ! size of x & f
  real*8, intent(in) :: x(n)     ! independent coordinate of f(x)
  real*8, intent(inout) :: f(n)  ! dependent coordinate of f(x)
                               ! smoothed on output
 
  real*8, intent(in) :: delta    ! smoothing parameter, in units of x
                               ! x-half-width of triangular weighting fcn
                               ! delta <= 0 ==> no smoothing.
                               ! delta <= (x(n)-x(1))/4 is enforced.
 
  real*8, intent(in) :: bc0      ! bdy cond at x(1), btw 0 and 1
  real*8, intent(in) :: bc1      ! bdy cond at x(n), btw 0 and 1
 
  integer, intent(out) :: ier  ! =0: OK; =1: ERROR
 
  ! possible errors: x(1:n) not strict ascending; BC not btw 0 and 1.
 
  ! boundary conditions: =0 means endpoint is fixed during smoothing
  !                         (minimal smoothing at endpoint)
  !                      =1 means df/dx -> 0 at endpoint
  !                         (maximal smoothing at endpoint)
  ! separate boundary condition control for each boundary or endpoint.
 
  !----------------------------------
  !  local
 
  real*8 deltax,epsx,zdum1
  integer i,idum1,idum2
 
  real*8, parameter :: ZERO = 0
  real*8, parameter :: ONE = 1
 
  real*8, dimension(:), allocatable :: del,eps,fwk
  !----------------------------------
  !  error checking...
 
  ier=0
  if(delta.le.ZERO) then
     return                ! no smoothing
  else if((bc0.lt.ZERO).or.(bc0.gt.ONE)) then
     ier=1
  else if((bc1.lt.ZERO).or.(bc1.gt.ONE)) then
     ier=1
  else
     do i=1,n-1
        if(x(i+1).le.x(i)) ier=1
     enddo
  endif
 
  if(ier.ne.0) return
 
  deltax=max(ZERO,min( (x(n)-x(1))/4, delta))
  epsx=max(ONE, 10*maxval(abs(f)))
 
  allocate(del(n),eps(n),fwk(n))
 
  fwk=f
  del=deltax
  eps=epsx
 
  call r8filtr6(x,fwk,f,n,eps,n,eps,0,del, &
       0,ZERO,0,ZERO,bc0,bc1,zdum1,idum1,idum2)
 
  deallocate(del,eps,fwk)
 
end subroutine r8_qksmooth
