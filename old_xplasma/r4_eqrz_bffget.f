!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZ_BFFGET
      SUBROUTINE R4_EQRZ_BFFGET(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZTOL,R4_BVEC,NLIST,IFCNS,IVECD,
     > R4_FVALS,NREGION,IERR)
 
      external EQRZ_BFFGET
 
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
      REAL R4_BVEC(3,IVEC)
 
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
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGET -- BVEC_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZTOL_act=R4_ZTOL
      BVEC_act=0
      FVALS_act=0
! call to original routine:  EQRZ_BFFGET
 
      CALL EQRZ_BFFGET(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZTOL_act,BVEC_act(1,1),
     > NLIST,IFCNS,IVECD,FVALS_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZ_BFFGRAD
      SUBROUTINE R4_EQRZ_BFFGRAD(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZTOL,R4_BVEC,R4_GBTENSR,NLIST,IFCNS,
     > IVECD,R4_FVALS,R4_FGRADS,NREGION,IERR)
 
      external EQRZ_BFFGRAD
 
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
      REAL R4_BVEC(3,IVEC)
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,IVEC)
 
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
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGRAD -- GBTENSR_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZ_BFFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZTOL_act=R4_ZTOL
      BVEC_act=0
      GBTENSR_act=0
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQRZ_BFFGRAD
 
      CALL EQRZ_BFFGRAD(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZTOL_act,BVEC_act(1,1),
     > GBTENSR_act(1,1,1),NLIST,IFCNS,IVECD,FVALS_act(1,1),
     > FGRADS_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_GBTENSR = GBTENSR_act
      DEALLOCATE(GBTENSR_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
      R4_FGRADS = FGRADS_act
      DEALLOCATE(FGRADS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZMAP_BFFGET
      SUBROUTINE R4_EQRZMAP_BFFGET(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZRHO,R4_ZCHI,R4_ZTOL,R4_BVEC,NLIST,
     > IFCNS,IVECD,R4_FVALS,NREGION,IERR)
 
      external EQRZMAP_BFFGET
 
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
 
 ! floating type, output only:
      REAL R4_BVEC(3,IVEC)
 
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
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGET -- BVEC_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZRHO_act=0
      ZCHI_act=0
      ZTOL_act=R4_ZTOL
      BVEC_act=0
      FVALS_act=0
! call to original routine:  EQRZMAP_BFFGET
 
      CALL EQRZMAP_BFFGET(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZTOL_act,BVEC_act(1,1),NLIST,IFCNS,IVECD,FVALS_act(1,1),NREGION,
     > IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQRZMAP_BFFGRAD
      SUBROUTINE R4_EQRZMAP_BFFGRAD(
     > IVEC,R4_ZR,R4_ZZ,R4_ZPHI,R4_ZRHO,R4_ZCHI,R4_ZTOL,R4_BVEC,
     > R4_GBTENSR,NLIST,IFCNS,IVECD,R4_FVALS,R4_FGRADS,NREGION,IERR)
 
      external EQRZMAP_BFFGRAD
 
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
 
 ! floating type, output only:
      REAL R4_BVEC(3,IVEC)
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,IVEC)
 
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
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGRAD -- GBTENSR_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQRZMAP_BFFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZR_act=R4_ZR
      ZZ_act=R4_ZZ
      ZPHI_act=R4_ZPHI
      ZRHO_act=0
      ZCHI_act=0
      ZTOL_act=R4_ZTOL
      BVEC_act=0
      GBTENSR_act=0
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQRZMAP_BFFGRAD
 
      CALL EQRZMAP_BFFGRAD(
     > IVEC,ZR_act(1),ZZ_act(1),ZPHI_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZTOL_act,BVEC_act(1,1),GBTENSR_act(1,1,1),NLIST,IFCNS,IVECD,
     > FVALS_act(1,1),FGRADS_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZR_act)
      DEALLOCATE(ZZ_act)
      DEALLOCATE(ZPHI_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_GBTENSR = GBTENSR_act
      DEALLOCATE(GBTENSR_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
      R4_FGRADS = FGRADS_act
      DEALLOCATE(FGRADS_act)
 
! exit
      return
      end
