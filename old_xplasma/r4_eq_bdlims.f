!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_BDLIMS
      SUBROUTINE R4_EQ_BDLIMS(
     > ITYPE,R4_ZRMIN,R4_ZRMAX,R4_ZZMIN,R4_ZZMAX,IERR)
 
      external EQ_BDLIMS
 
! argument declarations
      INTEGER ITYPE
 ! floating type, output only:
      REAL R4_ZRMIN
 
 ! floating type, output only:
      REAL R4_ZRMAX
 
 ! floating type, output only:
      REAL R4_ZZMIN
 
 ! floating type, output only:
      REAL R4_ZZMAX
 
      INTEGER IERR
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: ZRMIN_act
      REAL*8 :: ZRMAX_act
      REAL*8 :: ZZMIN_act
      REAL*8 :: ZZMAX_act
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      ZRMIN_act=0
      ZRMAX_act=0
      ZZMIN_act=0
      ZZMAX_act=0
! call to original routine:  EQ_BDLIMS
 
      CALL EQ_BDLIMS(
     > ITYPE,ZRMIN_act,ZRMAX_act,ZZMIN_act,ZZMAX_act,IERR)
 
! copy back outputs if modified.
      R4_ZRMIN = ZRMIN_act
      R4_ZRMAX = ZRMAX_act
      R4_ZZMIN = ZZMIN_act
      R4_ZZMAX = ZZMAX_act
 
! exit
      return
      end
