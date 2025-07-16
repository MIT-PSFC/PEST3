 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQM_FROMGEQDSK
      SUBROUTINE  C_EQM_FROMGEQDSK(
     > C_LABEL_LEN,C_LABEL,C_ZFILE_LEN,C_ZFILE,NS,NT1,FBDY,CBDY,IRZ,IER)
 
      external EQM_FROMGEQDSK
 
! argument declarations
 ! scalar char*(*) type, input/output:
      INTEGER C_LABEL_LEN
      BYTE_DECLARE C_LABEL(C_LABEL_LEN)
 
 ! scalar char*(*) type, input/output:
      INTEGER C_ZFILE_LEN
      BYTE_DECLARE C_ZFILE(C_ZFILE_LEN)
 
      INTEGER NS
      INTEGER NT1
      REAL*8 FBDY
      REAL*8 CBDY
      INTEGER IRZ
      INTEGER IER
 
! local (automatic array) declarations
      CHARACTER(C_LABEL_LEN) :: LABEL_act
      CHARACTER(C_ZFILE_LEN) :: ZFILE_act
 
! executable code:  copy for input
 
      call cstring(LABEL_act, C_LABEL,'2F')
      call cstring(ZFILE_act, C_ZFILE,'2F')
! call to original routine:  EQM_FROMGEQDSK
 
      CALL EQM_FROMGEQDSK(
     > LABEL_act,ZFILE_act,NS,NT1,FBDY,CBDY,IRZ,IER)
 
! copy back outputs if modified.
      call cstring(LABEL_act, C_LABEL,'2C')
      call cstring(ZFILE_act, C_ZFILE,'2C')
 
! exit
      return
      end
