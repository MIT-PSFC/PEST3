!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQX_DFAST
      SUBROUTINE R4_EQX_DFAST(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZDIST,IERR)
 
      external EQX_DFAST
 
! argument declarations
      INTEGER IVEC
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_ZDIST(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZDIST_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DFAST -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DFAST -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DFAST -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZDIST_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DFAST -- ZDIST_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZDIST_act=0
! call to original routine:  EQX_DFAST
 
      CALL EQX_DFAST(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZDIST_act(1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_ZDIST = ZDIST_act
      DEALLOCATE(ZDIST_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_DFAST
      SUBROUTINE R4_EQ_DFAST(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZDIST,IERR)
 
      external EQ_DFAST
 
! argument declarations
      INTEGER IVEC
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZDIST(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZDIST_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_DFAST -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_DFAST -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_DFAST -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZDIST_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_DFAST -- ZDIST_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZDIST_act=0
! call to original routine:  EQ_DFAST
 
      CALL EQ_DFAST(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZDIST_act(1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZDIST = ZDIST_act
      DEALLOCATE(ZDIST_act)
 
! exit
      return
      end
