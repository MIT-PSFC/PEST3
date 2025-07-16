 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQ_GANUM
      SUBROUTINE  C_EQ_GANUM(
     > C_ZNAME_LEN,C_ZNAME,IANUM)
 
      external EQ_GANUM
 
! argument declarations
 ! scalar char*(*) type, input only:
      INTEGER C_ZNAME_LEN
      BYTE_DECLARE C_ZNAME(C_ZNAME_LEN)
 
      INTEGER IANUM
 
! local (automatic array) declarations
      CHARACTER(C_ZNAME_LEN) :: ZNAME_act
 
! executable code:  copy for input
 
      call cstring(ZNAME_act, C_ZNAME,'2F')
! call to original routine:  EQ_GANUM
 
      CALL EQ_GANUM(
     > ZNAME_act,IANUM)
 
! copy back outputs if modified.
 
! exit
      return
      end
 
#include "fpreproc/byte_declare.h"
 
!----------------------------------------------------------
! cface  -- generated (f90 fixed form) C interface to:  EQ_GET_ANAME
      SUBROUTINE  C_EQ_GET_ANAME(
     > IANUM,C_ZNAME_LEN,C_ZNAME)
 
      external EQ_GET_ANAME
 
! argument declarations
      INTEGER IANUM
 ! scalar char*(*) type, output only:
      INTEGER C_ZNAME_LEN
      BYTE_DECLARE C_ZNAME(C_ZNAME_LEN)
 
 
! local (automatic array) declarations
      CHARACTER(C_ZNAME_LEN) :: ZNAME_act
 
! executable code:  copy for input
 
      ZNAME_act=' '
! call to original routine:  EQ_GET_ANAME
 
      CALL EQ_GET_ANAME(
     > IANUM,ZNAME_act)
 
! copy back outputs if modified.
      call cstring(ZNAME_act, C_ZNAME,'2C')
 
! exit
      return
      end
