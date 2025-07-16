!........................................................................
       SUBROUTINE pstwallck
!........................................................................
!      fix the off-diagonal m-1 x m block. diagonal m x m will be fixed
!       in addvac .
!       special treatment if plasma terminates at a conducting wall
!       this version added feb 17,77 jm.
!
!      off-diagonal block...
!
 USE pstcom

 USE l21com

 USE l22com

!     use l33com                  (pstcom)
 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER JMAX
      INTEGER J1
      INTEGER JM1
      INTEGER J2
      INTEGER JM2


!
      jmax = lmax(1) - lmin(1) + 1

!      modified for hermite cubic elements 4/10/79 rg
!
      if ( lcub )  go to 420
!
!      for linear elements...
!
      jtot1 = 1
      jtot2 = ipol * jsub(m) + 1
!
      do 401 j1 = 1, jmax
      jm1 = jtot1 + j1 - 1
      do 401 j2 = 1, jmax
      jm2 = jtot2 + j2 - 1
      amat(jm1,jm2) =  0.0_r8 
  401 continue
!
      return
!
  420 continue
!
!      for hermite cubic elements...
!
      jmax = 2 * jmax
      jtot1 = 1
      jtot2 = ipol * jsub(mp) + 1
      do 422 j1 = 1, jmax
      jm1 = jtot1 + j1 - 1
      do 421 j2 = 1, jmax, 2
      jm2 = jtot2 + j2 - 1
      amat ( jm1,jm2 ) =  0.0_r8 
  421 continue
  422 continue
!
      return
      end

