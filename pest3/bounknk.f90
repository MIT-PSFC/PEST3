!**********************************************************************
       SUBROUTINE pstbounknk
!**********************************************************************
!     selectively apply ideal wall boundary conditions to some fourier
!     components to simulate selected helical feed-back coils.
!     feb. 1992 jm
 USE pstcom

 USE l21com

 USE l22com

!     use l33com                  (pstcom)
 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER JMAX
      INTEGER JTOT3
      INTEGER I
      INTEGER JM2
      INTEGER J1
      INTEGER JM1
      INTEGER J2
      INTEGER JD


!     
      if ( .not. kinkstab) return
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
      jtot3 = jtot2 + jmax - 1
!     
      do 405 i = 1, lstab(1)
      jm2 = jtot2 + lstab(i+1)-lmin(1) 
      do 401 j1 = 1, jtot3
      amat(j1,jm2) =  0.0_r8 
  401 continue
      jm1 = jtot2 + lstab(i+1) -lmin(1)
      do 402 j2 = 1, jtot3
      amat(jm1,j2) =  0.0_r8 
 402  continue
!     
!     fill the diagonal with unity
!     
      jd = jtot2 + lstab(i+1) - lmin(1)
      amat(jd,jd) =  1.0_r8 
 405  continue
!
      return
!
 420  continue
!     
!     error single mode stabilization not implemented here as yet.
!     
      CALL psterrmes(outpst,'bounknk')
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

