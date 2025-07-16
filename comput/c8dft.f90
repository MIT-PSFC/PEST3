subroutine c8dft(x,y,n)
 
  !  compute the Discrete Fourier Transform of a complex sequence.
  !  Replacement for NAG c06fcf; may be less efficient.  DMC 15 Oct 2002
 
  implicit NONE
 
  integer, intent(in) :: n
  real*8, intent(inout) :: x(0:n-1)  ! real parts
  real*8, intent(inout) :: y(0:n-1)  ! imaginary parts
 
  integer i,j
  real*8 x0(0:n-1),y0(0:n-1)
  real*8 zsin(n-1),zcos(n-1)
 
  real*8, parameter :: C2PI = 6.2831853071795862D+00
  real*8, parameter :: ZERO = 0.0D+00
  real*8, parameter :: ONE = 1.0D+00
 
  real*8 zang0,znorm
 
  !---------------------------
 
  znorm=ONE/sqrt(ONE*n)
 
  x0(0)=x(0)
  y0(0)=y(0)
  do i=1,n-1
     x0(i)=x(i)
     y0(i)=y(i)
     x(0)=x(0)+x(i)    ! form a0
     y(0)=y(0)+y(i)    ! form a0
  enddo
 
  do i=1,n-1
     x(i)=x0(0)
     y(i)=y0(0)
 
     zang0=i*C2PI/n
     call r8sincos(zang0,n-1,zsin,zcos)
     do j=1,n-1
        x(i)=x(i)+x0(j)*zcos(j)+y0(j)*zsin(j)
        y(i)=y(i)+y0(j)*zcos(j)-x0(j)*zsin(j)
     enddo
  enddo
 
  x=znorm*x
  y=znorm*y
 
end subroutine c8dft
