 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQM_IRHOFUN
      SUBROUTINE  C_EQM_IRHOFUN(
     > ID_AXIS,C_ZLBL_LEN,C_ZLBL,INPROF,ZPROF,IFLAG,ID,IERR)
 
      external EQM_IRHOFUN
 
! argument declarations
      INTEGER ID_AXIS
 ! scalar char*(*) type, input/output:
      INTEGER C_ZLBL_LEN
      BYTE_DECLARE C_ZLBL(C_ZLBL_LEN)
 
      INTEGER INPROF
      INTEGER IFLAG
      INTEGER ID
      INTEGER IERR
      REAL*8 ZPROF(INPROF)
 
! local (automatic array) declarations
      CHARACTER(C_ZLBL_LEN) :: ZLBL_act
 
! executable code:  copy for input
 
      call cstring(ZLBL_act, C_ZLBL,'2F')
! call to original routine:  EQM_IRHOFUN
 
      CALL EQM_IRHOFUN(
     > ID_AXIS,ZLBL_act,INPROF,ZPROF,IFLAG,ID,IERR)
 
! copy back outputs if modified.
      call cstring(ZLBL_act, C_ZLBL,'2C')
 
! exit
      return
      end
