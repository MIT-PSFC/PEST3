!
!      4.2    addvac - add contributions from vacuum and surface terms
!           to total potential energy associated with perturbations.
!.......................................................................
      SUBROUTINE pstaddvac(is)
!.......................................................................
 USE pstcom
 USE mtrik1
 USE l21com
 USE l22com
 USE r33com

 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IS
      INTEGER IOK
      INTEGER JMAX
      INTEGER LM1
      real*8 QN
      INTEGER J1
      INTEGER JM1
      INTEGER L1
      real*8 ALNQ1
      INTEGER J2
      INTEGER JM2
      INTEGER L2
      real*8 ALNQ2
      INTEGER LGIVUP


 REAL*8, DIMENSION(nths) 	 :: sfft
 REAL*8, DIMENSION(nths) 	 :: omfft
 REAL*8, DIMENSION(nths) 	 :: onfft
 INTEGER, DIMENSION(3) 	 :: lfft
 INTEGER, DIMENSION(nths) 	 :: invfft
!
      if ( wall ) go to 11
      if(is /= 1) return


      jmax = lmax(1) - lmin(1) + 1
!
!      adjust vacmat for general toroidal coordinates expansion...
!
      lm1 = lmin(1) - 1
      qn = n * qa(nusurf)
!
      jtot1 = ipol * jsub(m) + 1
      jtot2 = jtot1
      do 102 j1 = 1, jmax
      jm1 = jtot1 + ibas * ( j1-1 )
      l1 = j1 + lm1
      alnq1 = l1 - qn
      do 102 j2 = 1, jmax
      jm2 = jtot2 + ibas * ( j2-1 )
      l2 = j2 + lm1
      alnq2 = l2 - qn
!
      amat(jm1,jm2) = amat(jm1,jm2) &
  +  0.5_r8       *( vacmat(j1,j2) + vacmat(j2,j1) ) &
  + ( 0._r8 , 0.5_r8 )  *( vacmti(j1,j2) - vacmti(j2,j1) )

!
 102  continue

      return
   11 continue
!      special treatment if plasma terminates at a conducting wall
!
!      modified for hermite cubics 4/10/79 rg
!
      jmax = lmax(1) - lmin(1) + 1
      if ( lcub )  go to 120
!
!      for linear piecewise elements...
!
      jtot1 = ipol * jsub(m) + 1
      jtot2 = jtot1
      do 112 j1 = 1, jmax
      jm1 = jtot1 + j1 - 1
      do 112 j2 = 1, jmax
      jm2 = jtot2 + j2 - 1
      amat(jm1,jm2) =  0.0_r8 

      if ( jm1  ==  jm2 )then
         amat(jm1,jm2) =  1.0_r8 
      endif
  112 continue
!
      return
  120 continue
!
      write(outmod,*)'vacuum calculation not implemented' &
               ,' for lcub = ',lcub
      write(itty  ,*)'vacuum calculation not implemented' &
               ,' for lcub = ',lcub
      stop 'ERROR in addvac'
!
 7000 CALL psterrmes ( outpst,'addvac',lgivup )
      end

