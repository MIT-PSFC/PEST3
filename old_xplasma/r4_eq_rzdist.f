!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQX_DIST
      SUBROUTINE R4_EQX_DIST(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZDIST,ILPT,R4_RLIM,R4_ZLIM,R4_PHILIM,
     > IERR)
 
      external EQX_DIST
 
! argument declarations
      INTEGER IVEC
      INTEGER ILPT
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_ZDIST(IVEC)
 
 ! floating type, output only:
      REAL R4_RLIM(IVEC)
 
 ! floating type, output only:
      REAL R4_ZLIM(IVEC)
 
 ! floating type, output only:
      REAL R4_PHILIM(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZDIST_act(:)
      REAL*8, allocatable :: RLIM_act(:)
      REAL*8, allocatable :: ZLIM_act(:)
      REAL*8, allocatable :: PHILIM_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DIST -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DIST -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DIST -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZDIST_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DIST -- ZDIST_act ALLOCATE error!')
      ALLOCATE(RLIM_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DIST -- RLIM_act ALLOCATE error!')
      ALLOCATE(ZLIM_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DIST -- ZLIM_act ALLOCATE error!')
      ALLOCATE(PHILIM_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_DIST -- PHILIM_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZDIST_act=0
      RLIM_act=0
      ZLIM_act=0
      PHILIM_act=0
! call to original routine:  EQX_DIST
 
      CALL EQX_DIST(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZDIST_act(1),ILPT,RLIM_act(1),
     > ZLIM_act(1),PHILIM_act(1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_ZDIST = ZDIST_act
      DEALLOCATE(ZDIST_act)
      R4_RLIM = RLIM_act
      DEALLOCATE(RLIM_act)
      R4_ZLIM = ZLIM_act
      DEALLOCATE(ZLIM_act)
      R4_PHILIM = PHILIM_act
      DEALLOCATE(PHILIM_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RZDIST
      SUBROUTINE R4_EQ_RZDIST(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZDIST,ILPT,R4_RLIM,R4_ZLIM,R4_PHILIM,
     > IERR)
 
      external EQ_RZDIST
 
! argument declarations
      INTEGER IVEC
      INTEGER ILPT
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZDIST(IVEC)
 
 ! floating type, output only:
      REAL R4_RLIM(IVEC)
 
 ! floating type, output only:
      REAL R4_ZLIM(IVEC)
 
 ! floating type, output only:
      REAL R4_PHILIM(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZDIST_act(:)
      REAL*8, allocatable :: RLIM_act(:)
      REAL*8, allocatable :: ZLIM_act(:)
      REAL*8, allocatable :: PHILIM_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZDIST -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZDIST -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZDIST -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZDIST_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZDIST -- ZDIST_act ALLOCATE error!')
      ALLOCATE(RLIM_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZDIST -- RLIM_act ALLOCATE error!')
      ALLOCATE(ZLIM_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZDIST -- ZLIM_act ALLOCATE error!')
      ALLOCATE(PHILIM_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZDIST -- PHILIM_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZDIST_act=0
      RLIM_act=0
      ZLIM_act=0
      PHILIM_act=0
! call to original routine:  EQ_RZDIST
 
      CALL EQ_RZDIST(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZDIST_act(1),ILPT,
     > RLIM_act(1),ZLIM_act(1),PHILIM_act(1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZDIST = ZDIST_act
      DEALLOCATE(ZDIST_act)
      R4_RLIM = RLIM_act
      DEALLOCATE(RLIM_act)
      R4_ZLIM = ZLIM_act
      DEALLOCATE(ZLIM_act)
      R4_PHILIM = PHILIM_act
      DEALLOCATE(PHILIM_act)
 
! exit
      return
      end
