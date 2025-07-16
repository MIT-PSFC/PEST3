 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQ_RESTORE
      SUBROUTINE  C_EQ_RESTORE(
     > C_FILENAME_LEN,C_FILENAME,IER)
 
      external EQ_RESTORE
 
! argument declarations
 ! scalar char*(*) type, input only:
      INTEGER C_FILENAME_LEN
      BYTE_DECLARE C_FILENAME(C_FILENAME_LEN)
 
      INTEGER IER
 
! local (automatic array) declarations
      CHARACTER(C_FILENAME_LEN) :: FILENAME_act
 
! executable code:  copy for input
 
      call cstring(FILENAME_act, C_FILENAME,'2F')
! call to original routine:  EQ_RESTORE
 
      CALL EQ_RESTORE(
     > FILENAME_act,IER)
 
! copy back outputs if modified.
 
! exit
      return
      end
