!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_BRZ
      SUBROUTINE R4_EQM_BRZ(
     > USERBVEC,R4_ZSM,IERR)
 
      external EQM_BRZ
 
! argument declarations
      EXTERNAL USERBVEC
 ! floating type, input only:
      REAL R4_ZSM
 
      INTEGER IERR
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: ZSM_act
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      ZSM_act=R4_ZSM
! call to original routine:  EQM_BRZ
 
      CALL EQM_BRZ(
     > USERBVEC,ZSM_act,IERR)
 
! copy back outputs if modified.
 
! exit
      return
      end
