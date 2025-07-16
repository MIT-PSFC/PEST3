!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_RZFUN
      SUBROUTINE R4_EQM_RZFUN(
     > IFUN,R4_LAMDA,IORDER,R4_ZSM,IERR)
 
      external EQM_RZFUN
 
! argument declarations
      INTEGER IFUN
 ! floating type, input only:
      REAL R4_LAMDA
 
      INTEGER IORDER
 ! floating type, input only:
      REAL R4_ZSM
 
      INTEGER IERR
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: LAMDA_act
      REAL*8 :: ZSM_act
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      LAMDA_act=R4_LAMDA
      ZSM_act=R4_ZSM
! call to original routine:  EQM_RZFUN
 
      CALL EQM_RZFUN(
     > IFUN,LAMDA_act,IORDER,ZSM_act,IERR)
 
! copy back outputs if modified.
 
! exit
      return
      end
