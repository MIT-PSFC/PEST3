subroutine integ_wts(x,w)

  !  return x evaluation points and weights for 
  !  a definite integral over the region [0,1].
  !  Note the x evaluation points returned are not in ascending order.

  real*8, intent(out) :: x(10),w(10)

  !----------------
  real*8 :: x1(5),w10(5)

  ! gauss-kronrod-patterson quadrature coefficients for use in
  ! quadpack routine qng.  these coefficients were calculated with
  ! 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.

  data x1    (  1) / 0.973906528517171720077964012084452D0/
  data x1    (  2) / 0.865063366688984510732096688423493D0/
  data x1    (  3) / 0.679409568299024406234327365114874D0/
  data x1    (  4) / 0.433395394129247190799265943165784D0/
  data x1    (  5) / 0.148874338981631210884826001129720D0/

  data w10   (  1) / 0.066671344308688137593568809893332D0/
  data w10   (  2) / 0.149451349150580593145776339657697D0/
  data w10   (  3) / 0.219086362515982043995534934228163D0/
  data w10   (  4) / 0.269266719309996355091226921569469D0/
  data w10   (  5) / 0.295524224714752870173892994651338D0/
  !----------------

  integer :: ii
  real*8, parameter :: HALF = 0.5d0
  real*8 :: xinc,ww

  do ii=1,5
     xinc = x1(ii)*HALF
     x(ii) = HALF + xinc
     x(5+ii) = HALF - xinc
     ww = w10(ii)*HALF
     w(ii) = ww
     w(5+ii) = ww
  enddo

end subroutine integ_wts
