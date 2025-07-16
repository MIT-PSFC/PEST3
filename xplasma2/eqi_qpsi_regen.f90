subroutine eqi_qpsi_regen(gval,inum,rs,zs,qval)

  !  estimate q(psi) from loop integral formula  
  !     q(psi) = 2pi*(R*B_phi)*[integral dl * (R/grad(psi)) * (1/R**2)]
  !  where the integral goes around the flux surface at psi.

  use eqi_geq_mod   ! contains Psi(R,Z) bicubic spline
  implicit NONE
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  !-----------------------------------------------------------
  !  arguments:

  real*8, intent(in) :: gval    ! "g", (R*B_phi), Tesla*m
  integer, intent(in) :: inum   ! #pts on surface
  real*8, intent(in) :: rs(inum),zs(inum)  ! closed flux surface contour
  !                     with rs(1)=rs(inum) and zs(1)=zs(inum)

  real*8, intent(out) :: qval   ! "q" value estimated

  !-----------------------------------------------------------
  !  local:

  real*8, dimension(:,:), allocatable :: grad_psi
  real*8 :: rz1(2),rz2(2),rz3(2),rcen(2),rad,rscale,arc,dl,zrgpsii
  real*8 :: zdot,zcross

  integer :: ii,iok

  real*8, parameter :: zlarge = 1.0e5_R8
  real*8, parameter :: z2pi = 6.2831853071795862_R8

  !-----------------------------------------------------------

  allocate(grad_psi(inum,2))

  call ezspline_gradient(pspl,inum,rs,zs,grad_psi,iok)
  call ezspline_error(iok)

  qval = 0

  rz2(1) = rs(inum-1)   ! (R,Z) point #2
  rz2(2) = zs(inum-1)

  rz3(1) = rs(1) ! = rs(inum)  ! (R,Z) point #3
  rz3(2) = zs(1) ! = zs(inum)

  !  loop around the contour

  do ii=1,inum-1

     rz1 = rz2          ! #2 -> #1
     rz2 = rz3          ! #3 -> #2
     rz3(1) = rs(ii+1)  ! new #3
     rz3(2) = zs(ii+1)

     zdot = (rz3(1)-rz2(1))*(rz2(1)-rz1(1)) + (rz3(2)-rz2(2))*(rz2(2)-rz1(2))
     if(zdot.le.0.0_R8) then
        ! negative dot product; two consective sections at acute angle
        ! for dL use 1/2 sum of lengths

        dl = 0.5_R8*( sqrt((rz3(1)-rz2(1))**2 + (rz3(2)-rz2(2))**2) + &
             sqrt((rz2(1)-rz1(1))**2 + (rz2(2)-rz1(2))**2) )
     else
        rscale = max(abs(rz3(1)-rz1(1)),abs(rz3(2)-rz1(2)))

        !  find radius of curvature; treat points as colinear if iok=1
        !  or if rad > zlarge*rscale

        call r8_trirad(rz1,rz2,rz3, rcen, rad, iok)
        if(iok.eq.0) then
           if(rad.gt.zlarge*rscale) iok=1
        endif

        if(iok.eq.1) then
           ! colinear formula -- half the distance rz1 to rz3
           dl = 0.5_R8*sqrt((rz3(1)-rz1(1))**2 + (rz3(2)-rz1(2))**2)

        else
           ! arc formula: atan2(sin(alpha),cos(alpha));
           !   v1.v3=|v1|*|v3|*cos(alpha); v1^v3=|v1|*|v3|*sin(alpha)

           zcross = (rz1(1)-rcen(1))*(rz3(2)-rcen(2)) - &
                (rz3(1)-rcen(1))*(rz1(2)-rcen(2))
           zdot = (rz1(1)-rcen(1))*(rz3(1)-rcen(1)) + &
                (rz1(2)-rcen(2))*(rz3(2)-rcen(2))

           arc = abs(atan2(zcross,zdot))

           dl = 0.5_R8*arc*rad   ! use half the arc rz1 to rz3
        endif
     endif

     !  R/(R**2*|grad(psi)|) = 1/(R*|grad(psi)|)

     zrgpsii = 1.0_R8/(rs(ii)* &
          max(1.0e-20_R8,(sqrt(grad_psi(ii,1)**2 + grad_psi(ii,2)**2))) )

     qval = qval + dl*zrgpsii

  enddo

  qval = qval * abs(gval)/z2pi

end subroutine eqi_qpsi_regen
