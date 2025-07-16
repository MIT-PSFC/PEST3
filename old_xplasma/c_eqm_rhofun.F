 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQM_RHOFUN
      SUBROUTINE  C_EQM_RHOFUN(
     > IORDER,ID_AXIS,C_ZLBL_LEN,C_ZLBL,ZDATA,IBC1,ZBC1,IBC2,ZBC2,
     > ID_RHOFUN,IERR)
 
      external EQM_RHOFUN
 
! argument declarations
      INTEGER IORDER
      INTEGER ID_AXIS
 ! scalar char*(*) type, input only:
      INTEGER C_ZLBL_LEN
      BYTE_DECLARE C_ZLBL(C_ZLBL_LEN)
 
      INTEGER IBC1
      REAL*8 ZBC1
      INTEGER IBC2
      REAL*8 ZBC2
      INTEGER ID_RHOFUN
      INTEGER IERR
      REAL*8 ZDATA(*)
 
! local (automatic array) declarations
      CHARACTER(C_ZLBL_LEN) :: ZLBL_act
 
! executable code:  copy for input
 
      call cstring(ZLBL_act, C_ZLBL,'2F')
! call to original routine:  EQM_RHOFUN
 
      CALL EQM_RHOFUN(
     > IORDER,ID_AXIS,ZLBL_act,ZDATA,IBC1,ZBC1,IBC2,ZBC2,ID_RHOFUN,IERR)
 
! copy back outputs if modified.
 
! exit
      return
      end
