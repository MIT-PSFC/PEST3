!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZ_FFGET
      SUBROUTINE R4_EQXYZ_FFGET(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZTOL,NLIST,IFCNS,IVECD,R4_FVALS,
     > NREGION,IERR)
 
      external EQXYZ_FFGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FFGET -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FFGET -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FFGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZTOL_act=R4_ZTOL
      FVALS_act=0
! call to original routine:  EQXYZ_FFGET
 
      CALL EQXYZ_FFGET(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZTOL_act,NLIST,IFCNS,IVECD,
     > FVALS_act(1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZ_FFGRAD
      SUBROUTINE R4_EQXYZ_FFGRAD(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZTOL,NLIST,IFCNS,IVECD,R4_FVALS,
     > R4_FGRADS,NREGION,IERR)
 
      external EQXYZ_FFGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
 ! floating type, output only:
      REAL R4_FGRADS(3,IVECD,NLIST)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FFGRAD -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FFGRAD -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FFGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZ_FFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZTOL_act=R4_ZTOL
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQXYZ_FFGRAD
 
      CALL EQXYZ_FFGRAD(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZTOL_act,NLIST,IFCNS,IVECD,
     > FVALS_act(1,1),FGRADS_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
      R4_FGRADS = FGRADS_act
      DEALLOCATE(FGRADS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZMAP_FFGET
      SUBROUTINE R4_EQXYZMAP_FFGET(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZTOL,NLIST,
     > IFCNS,IVECD,R4_FVALS,NREGION,IERR)
 
      external EQXYZMAP_FFGET
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZPHI(IVEC)
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVALS_act(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGET -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGET -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGET -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGET -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGET -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGET -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGET -- FVALS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZRHO_act=0
      ZCHI_act=0
      ZPHI_act=0
      ZTOL_act=R4_ZTOL
      FVALS_act=0
! call to original routine:  EQXYZMAP_FFGET
 
      CALL EQXYZMAP_FFGET(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZPHI_act(1),ZTOL_act,NLIST,IFCNS,IVECD,FVALS_act(1,1),NREGION,
     > IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_ZPHI = ZPHI_act
      DEALLOCATE(ZPHI_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQXYZMAP_FFGRAD
      SUBROUTINE R4_EQXYZMAP_FFGRAD(
     > IVEC,R4_XX,R4_YY,R4_ZZ,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZTOL,NLIST,
     > IFCNS,IVECD,R4_FVALS,R4_FGRADS,NREGION,IERR)
 
      external EQXYZMAP_FFGRAD
 
! argument declarations
      INTEGER IVEC
 ! floating type, input only:
      REAL R4_ZTOL
 
      INTEGER NLIST
      INTEGER IVECD
      INTEGER IERR
 ! floating type, input only:
      REAL R4_XX(IVEC)
 
 ! floating type, input only:
      REAL R4_YY(IVEC)
 
 ! floating type, input only:
      REAL R4_ZZ(IVEC)
 
 ! floating type, output only:
      REAL R4_ZRHO(IVEC)
 
 ! floating type, output only:
      REAL R4_ZCHI(IVEC)
 
 ! floating type, output only:
      REAL R4_ZPHI(IVEC)
 
      INTEGER IFCNS(NLIST)
 ! floating type, output only:
      REAL R4_FVALS(IVECD,NLIST)
 
 ! floating type, output only:
      REAL R4_FGRADS(3,IVECD,NLIST)
 
      INTEGER NREGION(IVEC)
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: XX_act(:)
      REAL*8, allocatable :: YY_act(:)
      REAL*8, allocatable :: ZZ_act(:)
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8 :: ZTOL_act
      REAL*8, allocatable :: FVALS_act(:,:)
      REAL*8, allocatable :: FGRADS_act(:,:,:)
 
! allocation of working arrays...
 
      ALLOCATE(XX_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGRAD -- XX_act ALLOCATE error!')
      ALLOCATE(YY_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGRAD -- YY_act ALLOCATE error!')
      ALLOCATE(ZZ_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGRAD -- ZZ_act ALLOCATE error!')
      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGRAD -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZCHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGRAD -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGRAD -- ZPHI_act ALLOCATE error!')
      ALLOCATE(FVALS_act(IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGRAD -- FVALS_act ALLOCATE error!')
      ALLOCATE(FGRADS_act(3,IVECD,NLIST)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQXYZMAP_FFGRAD -- FGRADS_act ALLOCATE error!')
 
! executable code:  copy for input
 
      XX_act=R4_XX
      YY_act=R4_YY
      ZZ_act=R4_ZZ
      ZRHO_act=0
      ZCHI_act=0
      ZPHI_act=0
      ZTOL_act=R4_ZTOL
      FVALS_act=0
      FGRADS_act=0
! call to original routine:  EQXYZMAP_FFGRAD
 
      CALL EQXYZMAP_FFGRAD(
     > IVEC,XX_act(1),YY_act(1),ZZ_act(1),ZRHO_act(1),ZCHI_act(1),
     > ZPHI_act(1),ZTOL_act,NLIST,IFCNS,IVECD,FVALS_act(1,1),
     > FGRADS_act(1,1,1),NREGION,IERR)
 
! copy back outputs if modified.
      DEALLOCATE(XX_act)
      DEALLOCATE(YY_act)
      DEALLOCATE(ZZ_act)
      R4_ZRHO = ZRHO_act
      DEALLOCATE(ZRHO_act)
      R4_ZCHI = ZCHI_act
      DEALLOCATE(ZCHI_act)
      R4_ZPHI = ZPHI_act
      DEALLOCATE(ZPHI_act)
      R4_FVALS = FVALS_act
      DEALLOCATE(FVALS_act)
      R4_FGRADS = FGRADS_act
      DEALLOCATE(FGRADS_act)
 
! exit
      return
      end
