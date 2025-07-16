 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQM_RZFUN2
      SUBROUTINE  C_EQM_RZFUN2(
     > C_ZNAME_LEN,C_ZNAME,IFUN,USERFCN,IARG,IORDER,ZSM,IERR)
 
      external EQM_RZFUN2
 
! argument declarations
 ! scalar char*(*) type, input only:
      INTEGER C_ZNAME_LEN
      BYTE_DECLARE C_ZNAME(C_ZNAME_LEN)
 
      INTEGER IFUN
      REAL*8 USERFCN
      EXTERNAL USERFCN
      INTEGER IARG
      INTEGER IORDER
      REAL*8 ZSM
      INTEGER IERR
 
! local (automatic array) declarations
      CHARACTER(C_ZNAME_LEN) :: ZNAME_act
 
! executable code:  copy for input
 
      call cstring(ZNAME_act, C_ZNAME,'2F')
! call to original routine:  EQM_RZFUN2
 
      CALL EQM_RZFUN2(
     > ZNAME_act,IFUN,USERFCN,IARG,IORDER,ZSM,IERR)
 
! copy back outputs if modified.
 
! exit
      return
      end
