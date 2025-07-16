!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQX_BDFIND
      SUBROUTINE R4_EQX_BDFIND(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZRHO_IN,R4_ZCHI_OUT,R4_ZPHI_OUT,
     > R4_ZDIST,IERR)
 
      external EQX_BDFIND
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZRHO_IN
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI_OUT(IVEC)
 
 ! floating type, output only:
      REAL R4_ZPHI_OUT(IVEC)
 
 ! floating type, output only:
      REAL R4_ZDIST(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8 :: ZRHO_IN_act
      REAL*8, allocatable :: ZCHI_OUT_act(:)
      REAL*8, allocatable :: ZPHI_OUT_act(:)
      REAL*8, allocatable :: ZDIST_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_BDFIND -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_BDFIND -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_BDFIND -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZCHI_OUT_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_BDFIND -- ZCHI_OUT_act ALLOCATE error!')
      ALLOCATE(ZPHI_OUT_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_BDFIND -- ZPHI_OUT_act ALLOCATE error!')
      ALLOCATE(ZDIST_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQX_BDFIND -- ZDIST_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZRHO_IN_act=R4_ZRHO_IN
      ZCHI_OUT_act=0
      ZPHI_OUT_act=0
      ZDIST_act=0
! call to original routine:  EQX_BDFIND
 
      CALL EQX_BDFIND(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZRHO_IN_act,ZCHI_OUT_act(1),
     > ZPHI_OUT_act(1),ZDIST_act(1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_ZCHI_OUT = ZCHI_OUT_act
      DEALLOCATE(ZCHI_OUT_act)
      R4_ZPHI_OUT = ZPHI_OUT_act
      DEALLOCATE(ZPHI_OUT_act)
      R4_ZDIST = ZDIST_act
      DEALLOCATE(ZDIST_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_BDFIND
      SUBROUTINE R4_EQ_BDFIND(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZRHO_IN,R4_ZCHI_OUT,R4_ZPHI_OUT,
     > R4_ZDIST,IERR)
 
      external EQ_BDFIND
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZRHO_IN
 
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI_OUT(IVEC)
 
 ! floating type, output only:
      REAL R4_ZPHI_OUT(IVEC)
 
 ! floating type, output only:
      REAL R4_ZDIST(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZRHO_IN_act
      REAL*8, allocatable :: ZCHI_OUT_act(:)
      REAL*8, allocatable :: ZPHI_OUT_act(:)
      REAL*8, allocatable :: ZDIST_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BDFIND -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BDFIND -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BDFIND -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZCHI_OUT_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BDFIND -- ZCHI_OUT_act ALLOCATE error!')
      ALLOCATE(ZPHI_OUT_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BDFIND -- ZPHI_OUT_act ALLOCATE error!')
      ALLOCATE(ZDIST_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BDFIND -- ZDIST_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZRHO_IN_act=R4_ZRHO_IN
      ZCHI_OUT_act=0
      ZPHI_OUT_act=0
      ZDIST_act=0
! call to original routine:  EQ_BDFIND
 
      CALL EQ_BDFIND(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZRHO_IN_act,ZCHI_OUT_act(1),
     > ZPHI_OUT_act(1),ZDIST_act(1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZCHI_OUT = ZCHI_OUT_act
      DEALLOCATE(ZCHI_OUT_act)
      R4_ZPHI_OUT = ZPHI_OUT_act
      DEALLOCATE(ZPHI_OUT_act)
      R4_ZDIST = ZDIST_act
      DEALLOCATE(ZDIST_act)
 
! exit
      return
      end
