!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_SURFGEO
      SUBROUTINE R4_EQ_SURFGEO(
     > IVEC,R4_RHO,R4_PHI,NTHETA,R4_ELONG,R4_TRIANG,R4_INDENT,R4_ZMIDP,
     > R4_RMNMIDP,R4_RMJMIDP,R4_LIMITS,IER)
 
      external EQ_SURFGEO
 
! argument declarations
      INTEGER IVEC
      INTEGER NTHETA
      INTEGER IER
 ! floating type, input/output:
      REAL R4_RHO(IVEC)
 
 ! floating type, input/output:
      REAL R4_PHI(IVEC)
 
 ! floating type, input/output:
      REAL R4_ELONG(IVEC)
 
 ! floating type, input/output:
      REAL R4_TRIANG(IVEC)
 
 ! floating type, input/output:
      REAL R4_INDENT(IVEC)
 
 ! floating type, input/output:
      REAL R4_ZMIDP(IVEC)
 
 ! floating type, input/output:
      REAL R4_RMNMIDP(IVEC)
 
 ! floating type, input/output:
      REAL R4_RMJMIDP(IVEC)
 
 ! floating type, input/output:
      REAL R4_LIMITS(20,IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: RHO_act(:)
      REAL*8, allocatable :: RHO_ref(:)
      REAL*8, allocatable :: PHI_act(:)
      REAL*8, allocatable :: PHI_ref(:)
      REAL*8, allocatable :: ELONG_act(:)
      REAL*8, allocatable :: ELONG_ref(:)
      REAL*8, allocatable :: TRIANG_act(:)
      REAL*8, allocatable :: TRIANG_ref(:)
      REAL*8, allocatable :: INDENT_act(:)
      REAL*8, allocatable :: INDENT_ref(:)
      REAL*8, allocatable :: ZMIDP_act(:)
      REAL*8, allocatable :: ZMIDP_ref(:)
      REAL*8, allocatable :: RMNMIDP_act(:)
      REAL*8, allocatable :: RMNMIDP_ref(:)
      REAL*8, allocatable :: RMJMIDP_act(:)
      REAL*8, allocatable :: RMJMIDP_ref(:)
      REAL*8, allocatable :: LIMITS_act(:,:)
      REAL*8, allocatable :: LIMITS_ref(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(RHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- RHO_act ALLOCATE error!')
      ALLOCATE(RHO_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- RHO_ref ALLOCATE error!')
      ALLOCATE(PHI_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- PHI_act ALLOCATE error!')
      ALLOCATE(PHI_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- PHI_ref ALLOCATE error!')
      ALLOCATE(ELONG_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- ELONG_act ALLOCATE error!')
      ALLOCATE(ELONG_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- ELONG_ref ALLOCATE error!')
      ALLOCATE(TRIANG_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- TRIANG_act ALLOCATE error!')
      ALLOCATE(TRIANG_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- TRIANG_ref ALLOCATE error!')
      ALLOCATE(INDENT_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- INDENT_act ALLOCATE error!')
      ALLOCATE(INDENT_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- INDENT_ref ALLOCATE error!')
      ALLOCATE(ZMIDP_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- ZMIDP_act ALLOCATE error!')
      ALLOCATE(ZMIDP_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- ZMIDP_ref ALLOCATE error!')
      ALLOCATE(RMNMIDP_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- RMNMIDP_act ALLOCATE error!')
      ALLOCATE(RMNMIDP_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- RMNMIDP_ref ALLOCATE error!')
      ALLOCATE(RMJMIDP_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- RMJMIDP_act ALLOCATE error!')
      ALLOCATE(RMJMIDP_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- RMJMIDP_ref ALLOCATE error!')
      ALLOCATE(LIMITS_act(20,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- LIMITS_act ALLOCATE error!')
      ALLOCATE(LIMITS_ref(20,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_SURFGEO -- LIMITS_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      RHO_act=R4_RHO
      RHO_ref=R4_RHO
 
      PHI_act=R4_PHI
      PHI_ref=R4_PHI
 
      ELONG_act=R4_ELONG
      ELONG_ref=R4_ELONG
 
      TRIANG_act=R4_TRIANG
      TRIANG_ref=R4_TRIANG
 
      INDENT_act=R4_INDENT
      INDENT_ref=R4_INDENT
 
      ZMIDP_act=R4_ZMIDP
      ZMIDP_ref=R4_ZMIDP
 
      RMNMIDP_act=R4_RMNMIDP
      RMNMIDP_ref=R4_RMNMIDP
 
      RMJMIDP_act=R4_RMJMIDP
      RMJMIDP_ref=R4_RMJMIDP
 
      LIMITS_act=R4_LIMITS
      LIMITS_ref=R4_LIMITS
 
! call to original routine:  EQ_SURFGEO
 
      CALL EQ_SURFGEO(
     > IVEC,RHO_act(1),PHI_act(1),NTHETA,ELONG_act(1),TRIANG_act(1),
     > INDENT_act(1),ZMIDP_act(1),RMNMIDP_act(1),RMJMIDP_act(1),
     > LIMITS_act(1,1),IER)
 
! copy back outputs if modified.
 
      where(RHO_act.ne.RHO_ref)
         R4_RHO = RHO_act
      end where
      DEALLOCATE(RHO_act)
      DEALLOCATE(RHO_ref)
 
      where(PHI_act.ne.PHI_ref)
         R4_PHI = PHI_act
      end where
      DEALLOCATE(PHI_act)
      DEALLOCATE(PHI_ref)
 
      where(ELONG_act.ne.ELONG_ref)
         R4_ELONG = ELONG_act
      end where
      DEALLOCATE(ELONG_act)
      DEALLOCATE(ELONG_ref)
 
      where(TRIANG_act.ne.TRIANG_ref)
         R4_TRIANG = TRIANG_act
      end where
      DEALLOCATE(TRIANG_act)
      DEALLOCATE(TRIANG_ref)
 
      where(INDENT_act.ne.INDENT_ref)
         R4_INDENT = INDENT_act
      end where
      DEALLOCATE(INDENT_act)
      DEALLOCATE(INDENT_ref)
 
      where(ZMIDP_act.ne.ZMIDP_ref)
         R4_ZMIDP = ZMIDP_act
      end where
      DEALLOCATE(ZMIDP_act)
      DEALLOCATE(ZMIDP_ref)
 
      where(RMNMIDP_act.ne.RMNMIDP_ref)
         R4_RMNMIDP = RMNMIDP_act
      end where
      DEALLOCATE(RMNMIDP_act)
      DEALLOCATE(RMNMIDP_ref)
 
      where(RMJMIDP_act.ne.RMJMIDP_ref)
         R4_RMJMIDP = RMJMIDP_act
      end where
      DEALLOCATE(RMJMIDP_act)
      DEALLOCATE(RMJMIDP_ref)
 
      where(LIMITS_act.ne.LIMITS_ref)
         R4_LIMITS = LIMITS_act
      end where
      DEALLOCATE(LIMITS_act)
      DEALLOCATE(LIMITS_ref)
 
! exit
      return
      end
