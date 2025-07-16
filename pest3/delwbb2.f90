!------------------------------------------------------
      SUBROUTINE pstdelwbb2
!------------------------------------------------------
! computes the big-big contribution to the energy...
!
! a. pletzer october  91
! modified for up-down asymmetric plasmas on 19/1/94
! perform cubic spline quadrature ap 3/11/99
!
 USE pstcom

 USE r33com

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

  integer jms, j1
  COMPLEX*16, dimension(nusurf,ms) :: zwoo, zwoe, zweo, zwee
  COMPLEX*16, PARAMETER :: IMAG = (0.0_r8,1.0_r8)
  REAL(r8) zintegral(nusurf)
  REAL(r8), dimension(nusurf) :: zy

  zwoo = (0.0_r8, 0.0_r8)
  zwoe = (0.0_r8, 0.0_r8)
  zweo = (0.0_r8, 0.0_r8)
  zwee = (0.0_r8, 0.0_r8)
  DO jms=1, ms
     DO surf=1, nusurf
        zwoo(surf,jms) = SUM( &
            &   conjg( x1frbo(surf,1:JMAX1,jms) )*alxi1o(surf,1:JMAX1) &
            & - conjg( x1frbo(surf,1:JMAX1,jms) )*qdchio(surf,1:JMAX1) &
            & - conjg( x1fpro(surf,1:JMAX1,jms) )*ddchio(surf,1:JMAX1) &
            &                )
        zwoe(surf,jms) = SUM( &
            &   conjg( x1frbo(surf,1:JMAX1,jms) )*alxi1e(surf,1:JMAX1) &
            & - conjg( x1frbo(surf,1:JMAX1,jms) )*qdchie(surf,1:JMAX1) &
            & - conjg( x1fpro(surf,1:JMAX1,jms) )*ddchie(surf,1:JMAX1) &
            &                )
        zweo(surf,jms) = SUM( &
            &   conjg( x1frbe(surf,1:JMAX1,jms) )*alxi1o(surf,1:JMAX1) &
            & - conjg( x1frbe(surf,1:JMAX1,jms) )*qdchio(surf,1:JMAX1) &
            & - conjg( x1fpre(surf,1:JMAX1,jms) )*ddchio(surf,1:JMAX1) &
            &                )
        zwee(surf,jms) = SUM( &
            &   conjg( x1frbe(surf,1:JMAX1,jms) )*alxi1e(surf,1:JMAX1) &
            & - conjg( x1frbe(surf,1:JMAX1,jms) )*qdchie(surf,1:JMAX1) &
            & - conjg( x1fpre(surf,1:JMAX1,jms) )*ddchie(surf,1:JMAX1) &
            &                )

        ENDDO

        zy(1:nusurf) = REAL(      zwoo(1:nusurf,jms), r8)
        CALL i2mex_integrate1d(nusurf,psinew,zy,nusurf,psinew,zintegral)
        w1l1oo(jms,ms) = zintegral(nusurf)
        zy(1:nusurf) = REAL(-IMAG*zwoo(1:nusurf,jms), r8)
        CALL i2mex_integrate1d(nusurf,psinew,zy,nusurf,psinew,zintegral)
        w1l1oo(jms,ms) = w1l1oo(jms,ms) + (0.0_r8,1.0_r8)*zintegral(nusurf)

        zy(1:nusurf) = REAL(      zwoe(1:nusurf,jms), r8)
        CALL i2mex_integrate1d(nusurf,psinew,zy,nusurf,psinew,zintegral)
        w1l1oe(jms,ms) = zintegral(nusurf)
        zy(1:nusurf) = REAL(-IMAG*zwoe(1:nusurf,jms), r8)
        CALL i2mex_integrate1d(nusurf,psinew,zy,nusurf,psinew,zintegral)
        w1l1oe(jms,ms) = w1l1oe(jms,ms) + (0.0_r8,1.0_r8)*zintegral(nusurf)

        zy(1:nusurf) = REAL(      zweo(1:nusurf,jms), r8)
        CALL i2mex_integrate1d(nusurf,psinew,zy,nusurf,psinew,zintegral)
        w1l1eo(jms,ms) = zintegral(nusurf)
        zy(1:nusurf) = REAL(-IMAG*zweo(1:nusurf,jms), r8)
        CALL i2mex_integrate1d(nusurf,psinew,zy,nusurf,psinew,zintegral)
        w1l1eo(jms,ms) = w1l1eo(jms,ms) + (0.0_r8,1.0_r8)*zintegral(nusurf)

        zy(1:nusurf) = REAL(      zwee(1:nusurf,jms), r8)
        CALL i2mex_integrate1d(nusurf,psinew,zy,nusurf,psinew,zintegral)
        w1l1ee(jms,ms) = zintegral(nusurf)
        zy(1:nusurf) = REAL(-IMAG*zwee(1:nusurf,jms), r8)
        CALL i2mex_integrate1d(nusurf,psinew,zy,nusurf,psinew,zintegral)
        w1l1ee(jms,ms) = w1l1ee(jms,ms) + (0.0_r8,1.0_r8)*zintegral(nusurf)
     ENDDO

!
      w1l1oo(ms,1:ms) = conjg( w1l1oo(1:ms,ms) )
      w1l1oe(ms,1:ms) = conjg( w1l1eo(1:ms,ms) )
      w1l1eo(ms,1:ms) = conjg( w1l1oe(1:ms,ms) )
      w1l1ee(ms,1:ms) = conjg( w1l1ee(1:ms,ms) )
!
    end SUBROUTINE pstdelwbb2

