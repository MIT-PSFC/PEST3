 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQM_FRHOCHI
      SUBROUTINE  C_EQM_FRHOCHI(
     > IORDER,ID_AX1,ID_AX2,C_ZLBL_LEN,C_ZLBL,ZDATA,ID1,IBCRHO0,ZBCRHO0,
     > IBCRHO1,ZBCRHO1,ID_FUN,IERR)
 
      external EQM_FRHOCHI
 
! argument declarations
      INTEGER IORDER
      INTEGER ID_AX1
      INTEGER ID_AX2
 ! scalar char*(*) type, input only:
      INTEGER C_ZLBL_LEN
      BYTE_DECLARE C_ZLBL(C_ZLBL_LEN)
 
      INTEGER ID1
      INTEGER IBCRHO0
      INTEGER IBCRHO1
      INTEGER ID_FUN
      INTEGER IERR
      REAL*8 ZDATA(ID1,*)
      REAL*8 ZBCRHO0(*)
      REAL*8 ZBCRHO1(*)
 
! local (automatic array) declarations
      CHARACTER(C_ZLBL_LEN) :: ZLBL_act
 
! executable code:  copy for input
 
      call cstring(ZLBL_act, C_ZLBL,'2F')
! call to original routine:  EQM_FRHOCHI
 
      CALL EQM_FRHOCHI(
     > IORDER,ID_AX1,ID_AX2,ZLBL_act,ZDATA,ID1,IBCRHO0,ZBCRHO0,IBCRHO1,
     > ZBCRHO1,ID_FUN,IERR)
 
! copy back outputs if modified.
 
! exit
      return
      end
