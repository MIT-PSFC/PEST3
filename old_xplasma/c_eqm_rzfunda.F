 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQM_RZFUNDA
      SUBROUTINE  C_EQM_RZFUNDA(
     > C_ZNAME_LEN,C_ZNAME,IFUN,ZDATA,INR,INZ,IORDER,ZSM,IERR)
 
      external EQM_RZFUNDA
 
! argument declarations
 ! scalar char*(*) type, input only:
      INTEGER C_ZNAME_LEN
      BYTE_DECLARE C_ZNAME(C_ZNAME_LEN)
 
      INTEGER IFUN
      INTEGER INR
      INTEGER INZ
      INTEGER IORDER
      REAL*8 ZSM
      INTEGER IERR
      REAL*8 ZDATA(INR,INZ)
 
! local (automatic array) declarations
      CHARACTER(C_ZNAME_LEN) :: ZNAME_act
 
! executable code:  copy for input
 
      call cstring(ZNAME_act, C_ZNAME,'2F')
! call to original routine:  EQM_RZFUNDA
 
      CALL EQM_RZFUNDA(
     > ZNAME_act,IFUN,ZDATA,INR,INZ,IORDER,ZSM,IERR)
 
! copy back outputs if modified.
 
! exit
      return
      end
