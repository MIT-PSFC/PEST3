!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_DBDY
      SUBROUTINE R4_EQM_DBDY(
     > R4_ZDIST,R4_ZRMINI,R4_ZRMAXI,R4_ZZMINI,R4_ZZMAXI,IERR)
 
      external EQM_DBDY
 
! argument declarations
 ! floating type, input only:
      REAL R4_ZDIST
 
 ! floating type, input only:
      REAL R4_ZRMINI
 
 ! floating type, input only:
      REAL R4_ZRMAXI
 
 ! floating type, input only:
      REAL R4_ZZMINI
 
 ! floating type, input only:
      REAL R4_ZZMAXI
 
      INTEGER IERR
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: ZDIST_act
      REAL*8 :: ZRMINI_act
      REAL*8 :: ZRMAXI_act
      REAL*8 :: ZZMINI_act
      REAL*8 :: ZZMAXI_act
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      ZDIST_act=R4_ZDIST
      ZRMINI_act=R4_ZRMINI
      ZRMAXI_act=R4_ZRMAXI
      ZZMINI_act=R4_ZZMINI
      ZZMAXI_act=R4_ZZMAXI
! call to original routine:  EQM_DBDY
 
      CALL EQM_DBDY(
     > ZDIST_act,ZRMINI_act,ZRMAXI_act,ZZMINI_act,ZZMAXI_act,IERR)
 
! copy back outputs if modified.
 
! exit
      return
      end
