!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GRID
      SUBROUTINE R4_EQ_GRID(
     > ID_AXIS,R4_ZDATA,NMAX,NGOT,IERR)
 
      external EQ_GRID
 
! argument declarations
      INTEGER ID_AXIS
      INTEGER NMAX
      INTEGER NGOT
      INTEGER IERR
 ! floating type, output only:
      REAL R4_ZDATA(NMAX)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZDATA_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZDATA_act(NMAX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GRID -- ZDATA_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZDATA_act=0
! call to original routine:  EQ_GRID
 
      CALL EQ_GRID(
     > ID_AXIS,ZDATA_act(1),NMAX,NGOT,IERR)
 
! copy back outputs if modified.
      R4_ZDATA = ZDATA_act
      DEALLOCATE(ZDATA_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GRID_ZC
      SUBROUTINE R4_EQ_GRID_ZC(
     > ID_AXIS,R4_ZDATA,NMAX,NGOT,IERR)
 
      external EQ_GRID_ZC
 
! argument declarations
      INTEGER ID_AXIS
      INTEGER NMAX
      INTEGER NGOT
      INTEGER IERR
 ! floating type, output only:
      REAL R4_ZDATA(NMAX)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZDATA_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZDATA_act(NMAX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GRID_ZC -- ZDATA_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZDATA_act=0
! call to original routine:  EQ_GRID_ZC
 
      CALL EQ_GRID_ZC(
     > ID_AXIS,ZDATA_act(1),NMAX,NGOT,IERR)
 
! copy back outputs if modified.
      R4_ZDATA = ZDATA_act
      DEALLOCATE(ZDATA_act)
 
! exit
      return
      end
