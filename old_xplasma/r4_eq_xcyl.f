!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_XCYL
      SUBROUTINE R4_EQ_XCYL(
     > IVEC,R4_ZR,R4_ZPHI,R4_XX,R4_YY)
 
      external EQ_XCYL
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_XX(IVEC)
 
 ! floating type, output only:
      REAL R4_YY(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYL -- ZR_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYL -- ZPHI_act ALLOCATE error!')
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYL -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYL -- YY_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZPHI_act=R4_ZPHI
      XX_act=0
      YY_act=0
! call to original routine:  EQ_XCYL
 
      CALL EQ_XCYL(
     > IVEC,ZR_act(1),ZPHI_act(1),XX_act(1),YY_act(1))
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZPHI_act)
      R4_XX = XX_act
      DEALLOCATE(XX_act)
      R4_YY = YY_act
      DEALLOCATE(YY_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_XCYLCS
      SUBROUTINE R4_EQ_XCYLCS(
     > IVEC,R4_ZR,R4_ZPHI,R4_XX,R4_YY,R4_ZCOSPHI,R4_ZSINPHI)
 
      external EQ_XCYLCS
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_XX(IVEC)
 
 ! floating type, output only:
      REAL R4_YY(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCOSPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZSINPHI(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZCOSPHI_act(:)
      REAL*8, allocatable :: ZSINPHI_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYLCS -- ZR_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYLCS -- ZPHI_act ALLOCATE error!')
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYLCS -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYLCS -- YY_act ALLOCATE error!')
      ALLOCATE(ZCOSPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYLCS -- ZCOSPHI_act ALLOCATE error!')
      ALLOCATE(ZSINPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XCYLCS -- ZSINPHI_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZPHI_act=R4_ZPHI
      XX_act=0
      YY_act=0
      ZCOSPHI_act=0
      ZSINPHI_act=0
! call to original routine:  EQ_XCYLCS
 
      CALL EQ_XCYLCS(
     > IVEC,ZR_act(1),ZPHI_act(1),XX_act(1),YY_act(1),ZCOSPHI_act(1),
     > ZSINPHI_act(1))
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZPHI_act)
      R4_XX = XX_act
      DEALLOCATE(XX_act)
      R4_YY = YY_act
      DEALLOCATE(YY_act)
      R4_ZCOSPHI = ZCOSPHI_act
      DEALLOCATE(ZCOSPHI_act)
      R4_ZSINPHI = ZSINPHI_act
      DEALLOCATE(ZSINPHI_act)
 
! exit
      return
      end
