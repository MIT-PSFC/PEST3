!**********************************************************************
!      entry point for segment governing vacuum and boundary conditions
!
!**********************************************************************
!
      SUBROUTINE pstent34
!
!      in this version of 34, external boundary conditions are imposed
!      in modular form after delta-w has been computed and stored on disk.
!
!      author: j manickam pppl june 1978
!.......................................................................
!      last revision date:
!      mar 27, 1979      modified bestpest version for pestr code rg.
!      feb,1980          modified for pest 2.3_r8    rg.
!......................................................................
!
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IPK
      INTEGER I
      INTEGER IS

!
!      calculate vacuum contribution if required
!
!aplet  19.5_r8 .94      if ( wall ) go to 10
      if ( wall  .OR.  imode >  1 ) go to 10
      CALL pstvacuum
   10 continue
!
!      fetch the last two sub-blocks for manipulation
!
!      repeat for potential and kinetic energy
!
!aplet  17.5_r8 .94
      ipk = 2
      if(lmarg) ipk = 1
      do 20 i = 1, ipk
!      do 20 i = 1, 2
      is = i
      CALL pstfetchm(is)
!
!      now for appropriate boundary conditions
!
      CALL pstboundy(is)
!
!      invoke wall boundary conditions for selected fourier comps.
!
! aplet  18.5_r8 .94 lstab not defined
!      if(is  ==  1) CALL pstbounknk
!
!      now write out the corrected delta-w block
!
      CALL pstwritem(is)
!
   20 continue
!
      return
      end

