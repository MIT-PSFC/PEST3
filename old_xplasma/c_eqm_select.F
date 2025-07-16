 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQM_SELECT
      SUBROUTINE  C_EQM_SELECT(
     > C_MODEL_LEN,C_MODEL,ILAXI)
 
      external EQM_SELECT
 
! argument declarations
 ! scalar char*(*) type, input only:
      INTEGER C_MODEL_LEN
      BYTE_DECLARE C_MODEL(C_MODEL_LEN)
 
      INTEGER ILAXI
 
! local (automatic array) declarations
      CHARACTER(C_MODEL_LEN) :: MODEL_act
 
! executable code:  copy for input
 
      call cstring(MODEL_act, C_MODEL,'2F')
! call to original routine:  EQM_SELECT
 
      CALL EQM_SELECT(
     > MODEL_act,ILAXI)
 
! copy back outputs if modified.
 
! exit
      return
      end
