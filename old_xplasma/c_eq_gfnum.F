 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQ_GFNUM
      SUBROUTINE  C_EQ_GFNUM(
     > C_ZNAME_LEN,C_ZNAME,IFNUM)
 
      external EQ_GFNUM
 
! argument declarations
 ! scalar char*(*) type, input only:
      INTEGER C_ZNAME_LEN
      BYTE_DECLARE C_ZNAME(C_ZNAME_LEN)
 
      INTEGER IFNUM
 
! local (automatic array) declarations
      CHARACTER(C_ZNAME_LEN) :: ZNAME_act
 
! executable code:  copy for input
 
      call cstring(ZNAME_act, C_ZNAME,'2F')
! call to original routine:  EQ_GFNUM
 
      CALL EQ_GFNUM(
     > ZNAME_act,IFNUM)
 
! copy back outputs if modified.
 
! exit
      return
      end
 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQ_GET_FNAME
      SUBROUTINE  C_EQ_GET_FNAME(
     > IFNUM,C_ZNAME_LEN,C_ZNAME)
 
      external EQ_GET_FNAME
 
! argument declarations
      INTEGER IFNUM
 ! scalar char*(*) type, output only:
      INTEGER C_ZNAME_LEN
      BYTE_DECLARE C_ZNAME(C_ZNAME_LEN)
 
 
! local (automatic array) declarations
      CHARACTER(C_ZNAME_LEN) :: ZNAME_act
 
! executable code:  copy for input
 
      ZNAME_act=' '
! call to original routine:  EQ_GET_FNAME
 
      CALL EQ_GET_FNAME(
     > IFNUM,ZNAME_act)
 
! copy back outputs if modified.
      call cstring(ZNAME_act, C_ZNAME,'2C')
 
! exit
      return
      end
