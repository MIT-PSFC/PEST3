subroutine eq_arc_scrunch(r,z,nth,ns, rmc1d,rms1d,zmc1d,zms1d, imom,ier)

  ! remap {R(theta,rho),Z(theta,rho)} to an evenly spaced theta grid [0:2pi]
  ! defined by equal arclengh on each surface, i.e.
  !
  !     L(theta1,theta2)/L(0,2pi) = (theta2-theta)/2pi

  ! on each surface contour, theta=0 is defined as the point where the
  ! contour passes Z=0 on the large major radius side.

  ! algorithm:  the contour {R(:,j),Z(:,j)} is remapped surface-by-surface;
  !             moments are only calculated at the end.  The FFT is not used
  !             in any intermediate calculation, so there is no moments 
  !             truncation error in the {R,Z} contour regeneration.

  ! minimal error checking here -- this is a "private" routine of the
  ! nscrunch library -- this should be called from scrunch_rz.f90 only

  use eq_arc_scrunch_mod
  implicit NONE

  !-------------------------------------------------------
  !  arguments

  integer, intent(in) :: nth    ! no. of theta pts
  integer, intent(in) :: ns     ! no. of surfaces
 
  ! numerical surfaces, modified on output:
  real*8, intent(inout) :: r(nth,ns),z(nth,ns)

  integer, intent(in) :: imom   ! no. of moments not counting 0'th moment
 
  real*8, intent(out) :: rmc1d(0:imom,ns)
  real*8, intent(out) :: rms1d(0:imom,ns)
  real*8, intent(out) :: zmc1d(0:imom,ns)
  real*8, intent(out) :: zms1d(0:imom,ns)
 
  integer, intent(out) :: ier   ! completion code: 0 = OK
 
  !-------------------------------------------------------
  !  local

  real*8 :: Raxis, Zaxis
  integer :: is
  !-------------------------------------------------------

  Raxis = R(1,1)
  Zaxis = Z(1,1)
  R(2:nth,1)=Raxis
  Z(2:nth,1)=Zaxis

  rmc1d(0:imom,1)=0
  rms1d(0:imom,1)=0
  zmc1d(0:imom,1)=0
  zms1d(0:imom,1)=0

  rmc1d(0,1)=Raxis
  zmc1d(0,1)=Zaxis

  do is = 2,ns

     call eq_arc_scrunch1(Raxis, Zaxis, &
          R(1:nth,is), Z(1:nth,is), &
          rmc1d(0:imom,is), rms1d(0:imom,is), &
          zmc1d(0:imom,is), zms1d(0:imom,is), &
          imom, ier)

     if(ier.ne.0) then
        R(1:nth,is:ns)=0
        Z(1:nth,is:ns)=0
        rmc1d(0:imom,is:ns)=0
        rms1d(0:imom,is:ns)=0
        zmc1d(0:imom,is:ns)=0
        zms1d(0:imom,is:ns)=0
        exit
     endif
  enddo

end subroutine eq_arc_scrunch
