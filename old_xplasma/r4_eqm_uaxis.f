!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_UAXIS
      SUBROUTINE R4_EQM_UAXIS(
     > ANAME,IKIN,IPER,R4_ZAXIS,INUM,R4_ZTOL,ID,IERR)
 
      external EQM_UAXIS
 
! argument declarations
      CHARACTER*(*) ANAME
      INTEGER IKIN
      INTEGER IPER
      INTEGER INUM
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER ID
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZAXIS(INUM)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZAXIS_act(:)
      REAL*8 :: ZTOL_act
 
! allocation of working arrays...
 
      ALLOCATE(ZAXIS_act(INUM)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_UAXIS -- ZAXIS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZAXIS_act=R4_ZAXIS
      ZTOL_act=R4_ZTOL
! call to original routine:  EQM_UAXIS
 
      CALL EQM_UAXIS(
     > ANAME,IKIN,IPER,ZAXIS_act(1),INUM,ZTOL_act,ID,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZAXIS_act)
 
! exit
      return
      end
