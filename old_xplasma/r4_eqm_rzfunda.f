!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_RZFUNDA
      SUBROUTINE R4_EQM_RZFUNDA(
     > ZNAME,IFUN,R4_ZDATA,INR,INZ,IORDER,R4_ZSM,IERR)
 
      external EQM_RZFUNDA
 
! argument declarations
      CHARACTER*(*) ZNAME
      INTEGER IFUN
      INTEGER INR
      INTEGER INZ
      INTEGER IORDER
 ! floating type, input only:
      REAL R4_ZSM
 
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_ZDATA(INR,INZ)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZDATA_act(:,:)
      REAL*8, allocatable :: ZDATA_ref(:,:)
      REAL*8 :: ZSM_act
 
! allocation of working arrays...
 
      ALLOCATE(ZDATA_act(INR,INZ)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZFUNDA -- ZDATA_act ALLOCATE error!')
      ALLOCATE(ZDATA_ref(INR,INZ)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZFUNDA -- ZDATA_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZDATA_act=R4_ZDATA
      ZDATA_ref=R4_ZDATA
 
      ZSM_act=R4_ZSM
! call to original routine:  EQM_RZFUNDA
 
      CALL EQM_RZFUNDA(
     > ZNAME,IFUN,ZDATA_act(1,1),INR,INZ,IORDER,ZSM_act,IERR)
 
! copy back outputs if modified.
 
      where(ZDATA_act.ne.ZDATA_ref)
         R4_ZDATA = ZDATA_act
      end where
      DEALLOCATE(ZDATA_act)
      DEALLOCATE(ZDATA_ref)
 
! exit
      return
      end
