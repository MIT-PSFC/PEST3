!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_BFFGET
      SUBROUTINE R4_EQ_BFFGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,IRZMODE,R4_BVEC,NLIST,IFCNS,IVECD,
     > R4_FVALS,NREGION,IERR)
 
      external EQ_BFFGET
 
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
      REAL R4_BVEC(3,abs(IVEC))
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
      INTEGER NREGION(abs(IVEC))
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGET -- BVEC_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      BVEC_act=0
      FVALS_act=0
! call to original routine:  EQ_BFFGET
 
      CALL EQ_BFFGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),IRZMODE,BVEC_act(1,1),
     > NLIST,IFCNS,IVECD,FVALS_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_BFFGRAD
      SUBROUTINE R4_EQ_BFFGRAD(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,IRZMODE,R4_BVEC,R4_GBTENSR,NLIST,
     > IFCNS,IVECD,R4_FVALS,R4_FGRADS,NREGION,IERR)
 
      external EQ_BFFGRAD
 
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
      REAL R4_BVEC(3,abs(IVEC))
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,abs(IVEC))
 
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
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGRAD -- GBTENSR_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_BFFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      BVEC_act=0
      GBTENSR_act=0
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQ_BFFGRAD
 
      CALL EQ_BFFGRAD(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),IRZMODE,BVEC_act(1,1),
     > GBTENSR_act(1,1,1),NLIST,IFCNS,IVECD,FVALS_act(1,1),
     > FGRADS_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
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
! r8real -- generated (f90 fixed form) REAL interface to:  EQMAP_BFFGET
      SUBROUTINE R4_EQMAP_BFFGET(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZR,R4_ZZ,IRZMODE,R4_BVEC,NLIST,
     > IFCNS,IVECD,R4_FVALS,NREGION,IERR)
 
      external EQMAP_BFFGET
 
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
 
 ! floating type, output only:
      REAL R4_BVEC(3,abs(IVEC))
 
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
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZR_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGET -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGET -- BVEC_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      ZR_act=0
      ZZ_act=0
      BVEC_act=0
      FVALS_act=0
! call to original routine:  EQMAP_BFFGET
 
      CALL EQMAP_BFFGET(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZR_act(1),ZZ_act(1),
     > IRZMODE,BVEC_act(1,1),NLIST,IFCNS,IVECD,FVALS_act(1,1),NREGION,
     > IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
      R4_BVEC = BVEC_act
      DEALLOCATE(BVEC_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMAP_BFFGRAD
      SUBROUTINE R4_EQMAP_BFFGRAD(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZR,R4_ZZ,IRZMODE,R4_BVEC,
     > R4_GBTENSR,NLIST,IFCNS,IVECD,R4_FVALS,R4_FGRADS,NREGION,IERR)
 
      external EQMAP_BFFGRAD
 
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
 
 ! floating type, output only:
      REAL R4_BVEC(3,abs(IVEC))
 
 ! floating type, output only:
      REAL R4_GBTENSR(3,3,abs(IVEC))
 
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
      REAL*8, allocatable :: BVEC_act(:,:)
      REAL*8, allocatable :: GBTENSR_act(:,:,:)
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZR_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGRAD -- ZR_act ALLOCATE error!')
      ALLOCATE(ZZ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(BVEC_act(3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGRAD -- BVEC_act ALLOCATE error!')
      ALLOCATE(GBTENSR_act(3,3,abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGRAD -- GBTENSR_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMAP_BFFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZCHI_act=R4_ZCHI
      ZPHI_act=R4_ZPHI
      ZR_act=0
      ZZ_act=0
      BVEC_act=0
      GBTENSR_act=0
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQMAP_BFFGRAD
 
      CALL EQMAP_BFFGRAD(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZR_act(1),ZZ_act(1),
     > IRZMODE,BVEC_act(1,1),GBTENSR_act(1,1,1),NLIST,IFCNS,IVECD,
     > FVALS_act(1,1),FGRADS_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZPHI_act)
      R4_ZR = ZR_act
      DEALLOCATE(ZR_act)
      R4_ZZ = ZZ_act
      DEALLOCATE(ZZ_act)
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
