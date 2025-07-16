!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FRZ
      SUBROUTINE R4_EQ_FRZ(
     > IVEC,R4_ZR,R4_ZZ,NLIST,IFLIST,IVECD,R4_ZANS,IERR)
 
      external EQ_FRZ
 
! argument declarations
      INTEGER IVEC
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
      INTEGER IFLIST(NLIST)
 ! floating type, output only:
      REAL R4_ZANS(IVECD,NLIST)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZANS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FRZ -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FRZ -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZANS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FRZ -- ZANS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZANS_act=0
! call to original routine:  EQ_FRZ
 
      CALL EQ_FRZ(
     > IVEC,ZR_act(1),ZZ_act(1),NLIST,IFLIST,IVECD,ZANS_act(1,1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      R4_ZANS = ZANS_act
      DEALLOCATE(ZANS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GRZ
      SUBROUTINE R4_EQ_GRZ(
     > IVEC,R4_ZR,R4_ZZ,NLIST,IFLIST,IVECD,R4_ZANS,IERR)
 
      external EQ_GRZ
 
! argument declarations
      INTEGER IVEC
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
      INTEGER IFLIST(NLIST)
 ! floating type, output only:
      REAL R4_ZANS(IVECD,2,NLIST)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZANS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GRZ -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GRZ -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZANS_act(IVECD,2,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GRZ -- ZANS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZANS_act=0
! call to original routine:  EQ_GRZ
 
      CALL EQ_GRZ(
     > IVEC,ZR_act(1),ZZ_act(1),NLIST,IFLIST,IVECD,ZANS_act(1,1,1),IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      R4_ZANS = ZANS_act
      DEALLOCATE(ZANS_act)
 
! exit
      return
      end
