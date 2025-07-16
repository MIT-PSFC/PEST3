!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RCYL
      SUBROUTINE R4_EQ_RCYL(
     > IVEC,R4_XX,R4_YY,R4_ZR,R4_ZPHI)
 
      external EQ_RCYL
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, output only:
      REAL R4_ZR(IVEC)
 
 ! floating type, output only:
      REAL R4_ZPHI(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYL -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYL -- YY_act ALLOCATE error!')
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYL -- ZR_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYL -- ZPHI_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZR_act=0
      ZPHI_act=0
! call to original routine:  EQ_RCYL
 
      CALL EQ_RCYL(
     > IVEC,XX_act(1),YY_act(1),ZR_act(1),ZPHI_act(1))
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZPHI = ZPHI_act
      DEALLOCATE(ZPHI_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RCYLCS
      SUBROUTINE R4_EQ_RCYLCS(
     > IVEC,R4_XX,R4_YY,R4_ZR,R4_ZPHI,R4_ZCOSPHI,R4_ZSINPHI)
 
      external EQ_RCYLCS
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, output only:
      REAL R4_ZR(IVEC)
 
 ! floating type, output only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCOSPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZSINPHI(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZCOSPHI_act(:)
      REAL*8, allocatable :: ZSINPHI_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYLCS -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYLCS -- YY_act ALLOCATE error!')
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYLCS -- ZR_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYLCS -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZCOSPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYLCS -- ZCOSPHI_act ALLOCATE error!')
      ALLOCATE(ZSINPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RCYLCS -- ZSINPHI_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZR_act=0
      ZPHI_act=0
      ZCOSPHI_act=0
      ZSINPHI_act=0
! call to original routine:  EQ_RCYLCS
 
      CALL EQ_RCYLCS(
     > IVEC,XX_act(1),YY_act(1),ZR_act(1),ZPHI_act(1),ZCOSPHI_act(1),
     > ZSINPHI_act(1))
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZPHI = ZPHI_act
      DEALLOCATE(ZPHI_act)
      R4_ZCOSPHI = ZCOSPHI_act
      DEALLOCATE(ZCOSPHI_act)
      R4_ZSINPHI = ZSINPHI_act
      DEALLOCATE(ZSINPHI_act)
 
! exit
      return
      end
