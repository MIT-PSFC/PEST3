 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQM_UAXIS
      SUBROUTINE  C_EQM_UAXIS(
     > C_ANAME_LEN,C_ANAME,IKIN,IPER,ZAXIS,INUM,ZTOL,ID,IERR)
 
      external EQM_UAXIS
 
! argument declarations
 ! scalar char*(*) type, input only:
      INTEGER C_ANAME_LEN
      BYTE_DECLARE C_ANAME(C_ANAME_LEN)
 
      INTEGER IKIN
      INTEGER IPER
      INTEGER INUM
      REAL*8 ZTOL
      INTEGER ID
      INTEGER IERR
      REAL*8 ZAXIS(INUM)
 
! local (automatic array) declarations
      CHARACTER(C_ANAME_LEN) :: ANAME_act
 
! executable code:  copy for input
 
      call cstring(ANAME_act, C_ANAME,'2F')
! call to original routine:  EQM_UAXIS
 
      CALL EQM_UAXIS(
     > ANAME_act,IKIN,IPER,ZAXIS,INUM,ZTOL,ID,IERR)
 
! copy back outputs if modified.
 
! exit
      return
      end
