!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_GETDETJ
      SUBROUTINE R4_EQ_GETDETJ(
     > IVEC,R4_ZRHO,R4_ZCHI,R4_ZPHI,R4_ZSCALE,R4_DETJ,IERR)
 
      external EQ_GETDETJ
 
! argument declarations
      INTEGER IVEC
 ! floating type, input/output:
      REAL R4_ZSCALE
 
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_ZRHO(abs(IVEC))
 
 ! floating type, input/output:
      REAL R4_ZCHI(abs(IVEC))
 
 ! floating type, input/output:
      REAL R4_ZPHI(abs(IVEC))
 
 ! floating type, input/output:
      REAL R4_DETJ(abs(IVEC))
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZRHO_ref(:)
      REAL*8, allocatable :: ZCHI_act(:)
      REAL*8, allocatable :: ZCHI_ref(:)
      REAL*8, allocatable :: ZPHI_act(:)
      REAL*8, allocatable :: ZPHI_ref(:)
      REAL*8 :: ZSCALE_act
      REAL*8 :: ZSCALE_ref
      REAL*8, allocatable :: DETJ_act(:)
      REAL*8, allocatable :: DETJ_ref(:)
 
! allocation of working arrays...
 
      ALLOCATE(ZRHO_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETDETJ -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZRHO_ref(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETDETJ -- ZRHO_ref ALLOCATE error!')
      ALLOCATE(ZCHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETDETJ -- ZCHI_act ALLOCATE error!')
      ALLOCATE(ZCHI_ref(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETDETJ -- ZCHI_ref ALLOCATE error!')
      ALLOCATE(ZPHI_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETDETJ -- ZPHI_act ALLOCATE error!')
      ALLOCATE(ZPHI_ref(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETDETJ -- ZPHI_ref ALLOCATE error!')
      ALLOCATE(DETJ_act(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETDETJ -- DETJ_act ALLOCATE error!')
      ALLOCATE(DETJ_ref(abs(IVEC))
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_GETDETJ -- DETJ_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      ZRHO_act=R4_ZRHO
      ZRHO_ref=R4_ZRHO
 
      ZCHI_act=R4_ZCHI
      ZCHI_ref=R4_ZCHI
 
      ZPHI_act=R4_ZPHI
      ZPHI_ref=R4_ZPHI
 
      ZSCALE_act=R4_ZSCALE
      ZSCALE_ref=R4_ZSCALE
 
      DETJ_act=R4_DETJ
      DETJ_ref=R4_DETJ
 
! call to original routine:  EQ_GETDETJ
 
      CALL EQ_GETDETJ(
     > IVEC,ZRHO_act(1),ZCHI_act(1),ZPHI_act(1),ZSCALE_act,DETJ_act(1),
     > IERR)
 
! copy back outputs if modified.
 
      where(ZRHO_act.ne.ZRHO_ref)
         R4_ZRHO = ZRHO_act
      end where
      DEALLOCATE(ZRHO_act)
      DEALLOCATE(ZRHO_ref)
 
      where(ZCHI_act.ne.ZCHI_ref)
         R4_ZCHI = ZCHI_act
      end where
      DEALLOCATE(ZCHI_act)
      DEALLOCATE(ZCHI_ref)
 
      where(ZPHI_act.ne.ZPHI_ref)
         R4_ZPHI = ZPHI_act
      end where
      DEALLOCATE(ZPHI_act)
      DEALLOCATE(ZPHI_ref)
 
      if(ZSCALE_act.ne.ZSCALE_ref) then
         R4_ZSCALE = ZSCALE_act
      endif
 
      where(DETJ_act.ne.DETJ_ref)
         R4_DETJ = DETJ_act
      end where
      DEALLOCATE(DETJ_act)
      DEALLOCATE(DETJ_ref)
 
! exit
      return
      end
