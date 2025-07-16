 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQM_WR_GEQDSK
      SUBROUTINE  C_EQM_WR_GEQDSK(
     > C_FILENAME_LEN,C_FILENAME,IER)
 
      external EQM_WR_GEQDSK
 
! argument declarations
 ! scalar char*(*) type, input/output:
      INTEGER C_FILENAME_LEN
      BYTE_DECLARE C_FILENAME(C_FILENAME_LEN)
 
      INTEGER IER
 
! local (automatic array) declarations
      CHARACTER(C_FILENAME_LEN) :: FILENAME_act
 
! executable code:  copy for input
 
      call cstring(FILENAME_act, C_FILENAME,'2F')
! call to original routine:  EQM_WR_GEQDSK
 
      CALL EQM_WR_GEQDSK(
     > FILENAME_act,IER)
 
! copy back outputs if modified.
      call cstring(FILENAME_act, C_FILENAME,'2C')
 
! exit
      return
      end
