!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FRHOCHI
      SUBROUTINE R4_EQ_FRHOCHI(
     > IVEC,R4_ZRHO,R4_ZCHI,NLIST,IFLIST,IVECD,R4_ZANS,IERR)
 
      external EQ_FRHOCHI
 
! argument declarations
      INTEGER IVEC
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
      INTEGER IFLIST(NLIST)
 ! floating type, output only:
      REAL R4_ZANS(IVECD,NLIST)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZANS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FRHOCHI -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FRHOCHI -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZANS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FRHOCHI -- ZANS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZANS_act=0
! call to original routine:  EQ_FRHOCHI
 
      CALL EQ_FRHOCHI(
     > IVEC,ZRHO_act(1),ZCHI_act(1),NLIST,IFLIST,IVECD,ZANS_act(1,1),
     > IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      R4_ZANS = ZANS_act
      DEALLOCATE(ZANS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GRHOCHI
      SUBROUTINE R4_EQ_GRHOCHI(
     > IVEC,R4_ZRHO,R4_ZCHI,NLIST,IFLIST,IVECD,R4_ZANS,IERR)
 
      external EQ_GRHOCHI
 
! argument declarations
      INTEGER IVEC
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
      INTEGER IFLIST(NLIST)
 ! floating type, output only:
      REAL R4_ZANS(IVECD,2,NLIST)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZANS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GRHOCHI -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GRHOCHI -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZANS_act(IVECD,2,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GRHOCHI -- ZANS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZANS_act=0
! call to original routine:  EQ_GRHOCHI
 
      CALL EQ_GRHOCHI(
     > IVEC,ZRHO_act(1),ZCHI_act(1),NLIST,IFLIST,IVECD,ZANS_act(1,1,1),
     > IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      R4_ZANS = ZANS_act
      DEALLOCATE(ZANS_act)
 
! exit
      return
      end
