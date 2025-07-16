!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_RHO
      SUBROUTINE R4_EQM_RHO(
     > R4_ZRHO,IBDY,R4_ZTOL,ID,IERR)
 
      external EQM_RHO
 
! argument declarations
      INTEGER IBDY
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER ID
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(IBDY)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8 :: ZTOL_act
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(IBDY)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RHO -- ZRHO_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZTOL_act=R4_ZTOL
! call to original routine:  EQM_RHO
 
      CALL EQM_RHO(
     > ZRHO_act(1),IBDY,ZTOL_act,ID,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
 
! exit
      return
      end
