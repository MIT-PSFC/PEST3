subroutine r4dft(f,n)
 
  !  compute the Discrete Fourier Transform of a real sequence; return
  !  the transform in Hermetian form.  Replacement for NAG c06eae; may
  !  be less efficient.  DMC 15 Oct 2002
 
  implicit NONE
 
  integer, intent(in) :: n
  real, intent(inout) :: f(0:n-1)
 
  integer i,ib,j
  real f0(0:n-1)
  real zsin(n-1),zcos(n-1)
 
  real, parameter :: C2PI = 6.2831853071795862D+00
  real, parameter :: ZERO = 0.0D+00
  real, parameter :: ONE = 1.0D+00
 
  real zang0,znorm
 
  !---------------------------
 
  znorm=ONE/sqrt(ONE*n)
 
  f0(0)=f(0)
  do i=1,n-1
     f0(i)=f(i)
     f(0)=f(0)+f(i)  ! form a0
  enddo
 
  do i=1,n/2
     ib=n-i
     f(i)=f0(0)
     if(ib.gt.(n/2)) f(ib)=ZERO
 
     zang0=i*C2PI/n
     call sincos(zang0,n-1,zsin,zcos)
     do j=1,n-1
        f(i)=f(i)+f0(j)*zcos(j)
        if(ib.gt.(n/2)) f(ib)=f(ib)-f0(j)*zsin(j)
     enddo
  enddo
 
  f=znorm*f
 
end subroutine r4dft
