!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FFGET
      SUBROUTINE R4_EQ_FFGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,NLIST,IFCNS,IVECD,R4_FVALS,NREGION,
     > IERR)
 
      external EQ_FFGET
 
! argument declarations
      INTEGER IVEC
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZPHI(abs(IVEC))
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FFGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FFGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FFGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      FVALS_act=0
! call to original routine:  EQ_FFGET
 
      CALL EQ_FFGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),NLIST,IFCNS,IVECD,
     > FVALS_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_FFGRAD
      SUBROUTINE R4_EQ_FFGRAD(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,IRZMODE,NLIST,IFCNS,IVECD,R4_FVALS,
     > R4_FGRADS,NREGION,IERR)
 
      external EQ_FFGRAD
 
! argument declarations
      INTEGER IVEC
      INTEGER IRZMODE
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZPHI(abs(IVEC))
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
 ! floating type, output only:
      REAL R4_FGRADS(3,IVECD,NLIST)
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FFGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FFGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FFGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_FFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQ_FFGRAD
 
      CALL EQ_FFGRAD(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),IRZMODE,NLIST,IFCNS,
     > IVECD,FVALS_act(1,1),FGRADS_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
      R4_FGRADS = FGRADS_act
      DEALLOCATE(FGRADS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMAP_FFGET
      SUBROUTINE R4_EQMAP_FFGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZR,R4_ZZ,NLIST,IFCNS,IVECD,
     > R4_FVALS,NREGION,IERR)
 
      external EQMAP_FFGET
 
! argument declarations
      INTEGER IVEC
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZPHI(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZR(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZZ(abs(IVEC))
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZR_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      ZR_act=0
      ZZ_act=0
      FVALS_act=0
! call to original routine:  EQMAP_FFGET
 
      CALL EQMAP_FFGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZR_act(1),ZZ_act(1),
     > NLIST,IFCNS,IVECD,FVALS_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMAP_FFGRAD
      SUBROUTINE R4_EQMAP_FFGRAD(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZR,R4_ZZ,IRZMODE,NLIST,IFCNS,
     > IVECD,R4_FVALS,R4_FGRADS,NREGION,IERR)
 
      external EQMAP_FFGRAD
 
! argument declarations
      INTEGER IVEC
      INTEGER IRZMODE
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input only:
      REAL R4_ZPHI(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZR(abs(IVEC))
 
 ! floating type, output only:
      REAL R4_ZZ(abs(IVEC))
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
 ! floating type, output only:
      REAL R4_FGRADS(3,IVECD,NLIST)
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZR_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_FFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      ZR_act=0
      ZZ_act=0
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQMAP_FFGRAD
 
      CALL EQMAP_FFGRAD(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZR_act(1),ZZ_act(1),
     > IRZMODE,NLIST,IFCNS,IVECD,FVALS_act(1,1),FGRADS_act(1,1,1),
     > NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
      R4_FGRADS = FGRADS_act
      DEALLOCATE(FGRADS_act)
 
! exit
      return
      end
