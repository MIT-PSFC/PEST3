!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_FROMGEQDSK
      SUBROUTINE R4_EQM_FROMGEQDSK(
     > LABEL,ZFILE,NS,NT1,R4_FBDY,R4_CBDY,IRZ,IER)
 
      external EQM_FROMGEQDSK
 
! argument declarations
      CHARACTER*(*) LABEL
      CHARACTER*(*) ZFILE
      INTEGER NS
      INTEGER NT1
 ! floating type, input/output:
      REAL R4_FBDY
 
 ! floating type, input/output:
      REAL R4_CBDY
 
      INTEGER IRZ
      INTEGER IER
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: FBDY_act
      REAL*8 :: FBDY_ref
      REAL*8 :: CBDY_act
      REAL*8 :: CBDY_ref
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      FBDY_act=R4_FBDY
      FBDY_ref=R4_FBDY
 
      CBDY_act=R4_CBDY
      CBDY_ref=R4_CBDY
 
! call to original routine:  EQM_FROMGEQDSK
 
      CALL EQM_FROMGEQDSK(
     > LABEL,ZFILE,NS,NT1,FBDY_act,CBDY_act,IRZ,IER)
 
! copy back outputs if modified.
 
      if(FBDY_act.ne.FBDY_ref) then
         R4_FBDY = FBDY_act
      endif
 
      if(CBDY_act.ne.CBDY_ref) then
         R4_CBDY = CBDY_act
      endif
 
! exit
      return
      end

!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_GEQ_MAXGSERR
      SUBROUTINE R4_EQM_GEQ_MAXGSERR(
     > R4_ZMAXERR)
 
      external EQM_GEQ_MAXGSERR
 
! argument declarations
 ! floating type, input/output:
      REAL R4_ZMAXERR
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8 :: ZMAXERR_act
      REAL*8 :: ZMAXERR_ref
 
! allocation of working arrays...
 
 
! executable code:  copy for input
 
      ZMAXERR_act=R4_ZMAXERR
      ZMAXERR_ref=R4_ZMAXERR
 
! call to original routine:  EQM_GEQ_MAXGSERR
 
      CALL EQM_GEQ_MAXGSERR(
     > ZMAXERR_act)
 
! copy back outputs if modified.
 
      if(ZMAXERR_act.ne.ZMAXERR_ref) then
         R4_ZMAXERR = ZMAXERR_act
      endif
 
! exit
      return
      end
