!
       SUBROUTINE pstexsolx(kms,kparit)
!
! extract solution from right-hand side vector ut and store into
! xisolo and xisole...
!
 USE pstcom

 USE combla

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER KMS
      INTEGER KPARIT
      INTEGER JBEG
      INTEGER JMAX3
      INTEGER J1
      INTEGER J
      INTEGER J0
      INTEGER J0M


!
ibas = 1
jmax3 = ibas*jmax1
!
       if( kparit == 1 ) then
!
! odd responses...
!
       xisolo(:,:,kms) =  0.0_r8 
!
          do m1 = 1,mp
          jbeg = jtot(m1)
          jmax1 = lmax(1) - lmin(1) + 1
          do j1 = 1,jmax3
             j = jbeg + j1 - 1
             xisolo(m1,j1,kms) = ut(j)
          end do
       end do
!
 13   continue
!
      else
!
! even responses...
!
       xisole(:,:,kms) =  0.0_r8 
!
       do m1 = 1,mp
          jbeg = jtot(m1)
          jmax1 = lmax(1) - lmin(1) + 1
          do j1 = 1,jmax3
             j = jbeg + j1 - 1
             xisole(m1,j1,kms) = ut(j)
          end do
       end do
!
      j0 = jmax3/2 + 1
      j0m = j0 + lmin(1) - 1
!
       end if
!
       return
!
 100  format(1x,'rational surface: ',i3,' odd parity' &
/ 1x,'psinod ',6(5x,'re',i1,7x,'im',i1,2x))
!
!
 101  format(1x,'rational surface: ',i3,' even parity' &
/ 1x,'psinod',6(5x,'re',i1,7x,'im',i1,2x))
 102  format(1x, f9.7,6(1x,f9.3,1x,f9.3))
 105  format(1x,f11.8,9(1x,e12.5))
       end


