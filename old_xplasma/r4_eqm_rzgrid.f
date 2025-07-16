!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_RZGRID
      SUBROUTINE R4_EQM_RZGRID(
     > R4_ZR,R4_ZZ,IAUTOR,IAUTOZ,INUMR,INUMZ,R4_ZTOL,ID_R,ID_Z,IERR)
 
      external EQM_RZGRID
 
! argument declarations
      INTEGER IAUTOR
      INTEGER IAUTOZ
      INTEGER INUMR
      INTEGER INUMZ
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER ID_R
      INTEGER ID_Z
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_ZR(INUMR)
 
 ! floating type, input/output:
      REAL R4_ZZ(INUMZ)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZR_ref(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZZ_ref(:)
      REAL*8 :: ZTOL_act
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(INUMR)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZGRID -- ZR_act ALLOCATE error!')
      ALLOCATE(ZR_ref(INUMR)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZGRID -- ZR_ref ALLOCATE error!')
      ALLOCATE(ZZ_act(INUMZ)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZGRID -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZZ_ref(INUMZ)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZGRID -- ZZ_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZR_ref=R4_ZR
 
      ZZ_act=R4_ZZ
      ZZ_ref=R4_ZZ
 
      ZTOL_act=R4_ZTOL
! call to original routine:  EQM_RZGRID
 
      CALL EQM_RZGRID(
     > ZR_act(1),ZZ_act(1),IAUTOR,IAUTOZ,INUMR,INUMZ,ZTOL_act,ID_R,ID_Z,
     > IERR)
 
! copy back outputs if modified.
 
      where(ZR_act.ne.ZR_ref)
         R4_ZR = ZR_act
      end where
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZR_ref)
 
      where(ZZ_act.ne.ZZ_ref)
         R4_ZZ = ZZ_act
      end where
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZZ_ref)
 
! exit
      return
      end
