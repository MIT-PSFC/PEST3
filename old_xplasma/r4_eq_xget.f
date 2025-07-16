!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_XGET
      SUBROUTINE R4_EQ_XGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_XX,R4_YY,R4_ZZ,NREGION,IERR)
 
      external EQ_XGET
 
! argument declarations
      INTEGER IVEC
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, input only:
      REAL R4_ZCHI(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_XX(IVEC)
 
 ! floating type, output only:
      REAL R4_YY(IVEC)
 
 ! floating type, output only:
      REAL R4_ZZ(IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XGET -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XGET -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_XGET -- ZZ_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      XX_act=0
      YY_act=0
      ZZ_act=0
! call to original routine:  EQ_XGET
 
      CALL EQ_XGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),XX_act(1),YY_act(1),
     > ZZ_act(1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_XX = XX_act
      DEALLOCATE(XX_act)
      R4_YY = YY_act
      DEALLOCATE(YY_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RZGET
      SUBROUTINE R4_EQ_RZGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZR,R4_ZZ,NREGION,IERR)
 
      external EQ_RZGET
 
! argument declarations
      INTEGER IVEC
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, input only:
      REAL R4_ZCHI(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZR(IVEC)
 
 ! floating type, output only:
      REAL R4_ZZ(IVEC)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RZGET -- ZZ_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      ZR_act=0
      ZZ_act=0
! call to original routine:  EQ_RZGET
 
      CALL EQ_RZGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZR_act(1),ZZ_act(1),
     > NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
 
! exit
      return
      end
