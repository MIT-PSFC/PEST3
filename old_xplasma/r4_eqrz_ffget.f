!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZ_FFGET
      SUBROUTINE R4_EQRZ_FFGET(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZTOL,NLIST,IFCNS,IVECD,R4_FVALS,
     > NREGION,IERR)
 
      external EQRZ_FFGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FFGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FFGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FFGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZTOL_act=R4_ZTOL
      FVALS_act=0
! call to original routine:  EQRZ_FFGET
 
      CALL EQRZ_FFGET(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZTOL_act,NLIST,IFCNS,IVECD,
     > FVALS_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZ_FFGRAD
      SUBROUTINE R4_EQRZ_FFGRAD(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZTOL,NLIST,IFCNS,IVECD,R4_FVALS,
     > R4_FGRADS,NREGION,IERR)
 
      external EQRZ_FFGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
 ! floating type, output only:
      REAL R4_FGRADS(3,IVECD,NLIST)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FFGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FFGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FFGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_FFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZTOL_act=R4_ZTOL
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQRZ_FFGRAD
 
      CALL EQRZ_FFGRAD(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZTOL_act,NLIST,IFCNS,IVECD,
     > FVALS_act(1,1),FGRADS_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
      R4_FGRADS = FGRADS_act
      DEALLOCATE(FGRADS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZMAP_FFGET
      SUBROUTINE R4_EQRZMAP_FFGET(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZRHO,R4_ZCHI,R4_ZTOL,NLIST,IFCNS,
     > IVECD,R4_FVALS,NREGION,IERR)
 
      external EQRZMAP_FFGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI(IVEC)
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZRHO_act=0
      ZCHI_act=0
      ZTOL_act=R4_ZTOL
      FVALS_act=0
! call to original routine:  EQRZMAP_FFGET
 
      CALL EQRZMAP_FFGET(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZTOL_act,NLIST,IFCNS,IVECD,FVALS_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZMAP_FFGRAD
      SUBROUTINE R4_EQRZMAP_FFGRAD(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZRHO,R4_ZCHI,R4_ZTOL,NLIST,IFCNS,
     > IVECD,R4_FVALS,R4_FGRADS,NREGION,IERR)
 
      external EQRZMAP_FFGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZR(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, input only:
      REAL R4_ZPHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI(IVEC)
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
 ! floating type, output only:
      REAL R4_FGRADS(3,IVECD,NLIST)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZR_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_FFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZRHO_act=0
      ZCHI_act=0
      ZTOL_act=R4_ZTOL
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQRZMAP_FFGRAD
 
      CALL EQRZMAP_FFGRAD(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZTOL_act,NLIST,IFCNS,IVECD,FVALS_act(1,1),FGRADS_act(1,1,1),
     > NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
      R4_FGRADS = FGRADS_act
      DEALLOCATE(FGRADS_act)
 
! exit
      return
      end
