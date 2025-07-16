 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQ_FLXINT
      SUBROUTINE  C_EQ_FLXINT(
     > C_NAME_LEN,C_NAME,NOPTION,RESULT,INUMCHI,INUMRHO,IERR)
 
      external EQ_FLXINT
 
! argument declarations
 ! scalar char*(*) type, input/output:
      INTEGER C_NAME_LEN
      BYTE_DECLARE C_NAME(C_NAME_LEN)
 
      INTEGER NOPTION
      INTEGER INUMCHI
      INTEGER INUMRHO
      INTEGER IERR
      REAL*8 RESULT(INUMCHI,INUMRHO)
 
! local (automatic array) declarations
      CHARACTER(C_NAME_LEN) :: NAME_act
 
! executable code:  copy for input
 
      call cstring(NAME_act, C_NAME,'2F')
! call to original routine:  EQ_FLXINT
 
      CALL EQ_FLXINT(
     > NAME_act,NOPTION,RESULT,INUMCHI,INUMRHO,IERR)
 
! copy back outputs if modified.
      call cstring(NAME_act, C_NAME,'2C')
 
! exit
      return
      end
