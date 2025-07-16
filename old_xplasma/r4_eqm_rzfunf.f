!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_RZFUNF
      SUBROUTINE R4_EQM_RZFUNF(
     > IFUN,USERFCN,IARG,IORDER,R4_ZSM,IERR)
 
      external EQM_RZFUNF
 
! argument declarations
      INTEGER IFUN
      REAL*8 USERFCN
      EXTERNAL USERFCN
      INTEGER IARG
      INTEGER IORDER
 ! floating type, input only:
      REAL R4_ZSM
 
      INTEGER IERR
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: ZSM_act
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      ZSM_act=R4_ZSM
! call to original routine:  EQM_RZFUNF
 
      CALL EQM_RZFUNF(
     > IFUN,USERFCN,IARG,IORDER,ZSM_act,IERR)
 
! copy back outputs if modified.
 
! exit
      return
      end
