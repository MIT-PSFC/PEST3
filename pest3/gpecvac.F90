!...................................................................
      SUBROUTINE pstgpecvac
!...................................................................
!
!     Vacuum matrix from GPEC's VACUUM (mscvac), as used by DCON/STRIDE.
!     Wall comes from vac.in in the working directory. Replaces pstvacuum
!     when lgpecvac is set (-V<mthvac>).
!
 USE pstcom
 USE newmet
 USE l22com
#ifdef PEST3_GPEC_VACUUM
 USE vacuum_mod, ONLY: mscvac
#endif
 IMPLICIT NONE
 integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
#ifdef PEST3_GPEC_VACUUM
 INTEGER :: j1, j2
 REAL(r8) :: q, nq, alnq1, alnq2, clockwise
 REAL(r8), DIMENSION(mth+1) :: xb, zb, db
 COMPLEX(r8), DIMENSION(jmax1,jmax1) :: wv
 REAL(r8), DIMENSION(:,:), ALLOCATABLE :: grri, xzpts
 REAL(r8), EXTERNAL :: pstclockwise

 q = qa(nusurf)
 nq = n*q
! boundary on PEST3's theta grid, theta_j = (j-1)*dth, closed.
 xb = xa(1:mth+1,nusurf)
 zb = za(1:mth+1,nusurf)
 db = delta(1:mth+1,nusurf)
 clockwise = pstclockwise(xb, zb, mth+1)
 ! VACUUM needs clockwise theta (as DCON passes it): map theta -> -theta, delta -> -delta,
 ! which turns each basis function into its conjugate; conjugate wv back below.
 if (clockwise < 0.0_r8) then
    xb = xb(mth+1:1:-1)
    zb = zb(mth+1:1:-1)
    db = -db(mth+1:1:-1)
 endif

 write(outmod,'(/,1x,a,i6)') 'GPEC vacuum (mscvac), mthvac =', mthvac
 write(itty  ,'(/,1x,a,i6)') 'GPEC vacuum (mscvac), mthvac =', mthvac

 CALL set_dcon_params(mth, lmin(1), lmax(1), NINT(n), q, xb, zb, db)
 ALLOCATE(grri(2*(mthvac+5),2*jmax1), xzpts(mthvac+5,4))
 CALL mscvac(wv, jmax1, mth, mthvac, .TRUE., 1.0_r8, .FALSE., .FALSE., &
      grri, xzpts)
 ! wall used (none if farwal): mthvac unique points, Julia GPEC wall-file format
 if (maxval(abs(xzpts(1:mthvac,3))) > 0.0_r8) then
    open(unit=97, file='gpecvac_wall.dat', status='replace')
    write(97,'(i6)') mthvac
    write(97,'(es24.15)') 0.5_r8*(maxval(xzpts(1:mthvac,3)) + minval(xzpts(1:mthvac,3)))
    write(97,'(a)') '  i  R  Z (VACUUM wall, from PEST3 -V)'
    do j1 = 1, mthvac
       write(97,'(i6,2es24.15)') j1, xzpts(j1,3), xzpts(j1,4)
    enddo
    close(97)
 endif
 DEALLOCATE(grri, xzpts)
 CALL unset_dcon_params
 if (clockwise < 0.0_r8) wv = CONJG(wv)

 do j1 = 1, jmax1
    alnq1 = j1 + lmin(1) - 1 - nq
    do j2 = 1, jmax1
       alnq2 = j2 + lmin(1) - 1 - nq
       vacmat(j1,j2) = alnq1*alnq2*REAL(wv(j1,j2),r8)
       vacmti(j1,j2) = alnq1*alnq2*AIMAG(wv(j1,j2))
    enddo
 enddo
#else
 write(itty,*) 'ERROR: pest3 built without GPEC vacuum (PEST3_GPEC_VACUUM)'
 stop 'pest3 built without GPEC vacuum'
#endif
      END SUBROUTINE pstgpecvac

!...................................................................
      FUNCTION pstclockwise(x, z, np1)
!...................................................................
!     +1 if the closed curve (x,z)(1:np1) runs clockwise in (R,Z), else -1.
!     Signed area, so it does not depend on where theta = 0 lies.
 IMPLICIT NONE
 integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
 INTEGER, INTENT(IN) :: np1
 REAL(r8), INTENT(IN) :: x(np1), z(np1)
 REAL(r8) :: pstclockwise
 pstclockwise = -SIGN(1.0_r8, SUM(x(1:np1-1)*z(2:np1) - x(2:np1)*z(1:np1-1)))
      END FUNCTION pstclockwise
