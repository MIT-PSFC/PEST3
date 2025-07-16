!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_PSMOM
      SUBROUTINE R4_EQ_PSMOM(
     > IVEC,R4_RHO,R4_RHOMIN,NMOM,R4_PHI1,IDMOM,R4_PSMOM,R4_GAMMA,
     > R4_NGRDB2,R4_B2,R4_ZERROR,IER)
 
      external EQ_PSMOM
 
! argument declarations
      INTEGER IVEC
 ! floating type, input/output:
      REAL R4_RHOMIN
 
      INTEGER NMOM
 ! floating type, input/output:
      REAL R4_PHI1
 
      INTEGER IDMOM
      INTEGER IER
 ! floating type, input/output:
      REAL R4_RHO(IVEC)
 
 ! floating type, input/output:
      REAL R4_PSMOM(IDMOM,IVEC)
 
 ! floating type, input/output:
      REAL R4_GAMMA(IVEC)
 
 ! floating type, input/output:
      REAL R4_NGRDB2(IVEC)
 
 ! floating type, input/output:
      REAL R4_B2(IVEC)
 
 ! floating type, input/output:
      REAL R4_ZERROR(IVEC)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: RHO_act(:)
      REAL*8, allocatable :: RHO_ref(:)
      REAL*8 :: RHOMIN_act
      REAL*8 :: RHOMIN_ref
      REAL*8 :: PHI1_act
      REAL*8 :: PHI1_ref
      REAL*8, allocatable :: PSMOM_act(:,:)
      REAL*8, allocatable :: PSMOM_ref(:,:)
      REAL*8, allocatable :: GAMMA_act(:)
      REAL*8, allocatable :: GAMMA_ref(:)
      REAL*8, allocatable :: NGRDB2_act(:)
      REAL*8, allocatable :: NGRDB2_ref(:)
      REAL*8, allocatable :: B2_act(:)
      REAL*8, allocatable :: B2_ref(:)
      REAL*8, allocatable :: ZERROR_act(:)
      REAL*8, allocatable :: ZERROR_ref(:)
 
! allocation of working arrays...
 
      ALLOCATE(RHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- RHO_act ALLOCATE error!')
      ALLOCATE(RHO_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- RHO_ref ALLOCATE error!')
      ALLOCATE(PSMOM_act(IDMOM,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- PSMOM_act ALLOCATE error!')
      ALLOCATE(PSMOM_ref(IDMOM,IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- PSMOM_ref ALLOCATE error!')
      ALLOCATE(GAMMA_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- GAMMA_act ALLOCATE error!')
      ALLOCATE(GAMMA_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- GAMMA_ref ALLOCATE error!')
      ALLOCATE(NGRDB2_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- NGRDB2_act ALLOCATE error!')
      ALLOCATE(NGRDB2_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- NGRDB2_ref ALLOCATE error!')
      ALLOCATE(B2_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- B2_act ALLOCATE error!')
      ALLOCATE(B2_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- B2_ref ALLOCATE error!')
      ALLOCATE(ZERROR_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- ZERROR_act ALLOCATE error!')
      ALLOCATE(ZERROR_ref(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_PSMOM -- ZERROR_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      RHO_act=R4_RHO
      RHO_ref=R4_RHO
 
      RHOMIN_act=R4_RHOMIN
      RHOMIN_ref=R4_RHOMIN
 
      PHI1_act=R4_PHI1
      PHI1_ref=R4_PHI1
 
      PSMOM_act=R4_PSMOM
      PSMOM_ref=R4_PSMOM
 
      GAMMA_act=R4_GAMMA
      GAMMA_ref=R4_GAMMA
 
      NGRDB2_act=R4_NGRDB2
      NGRDB2_ref=R4_NGRDB2
 
      B2_act=R4_B2
      B2_ref=R4_B2
 
      ZERROR_act=R4_ZERROR
      ZERROR_ref=R4_ZERROR
 
! call to original routine:  EQ_PSMOM
 
      CALL EQ_PSMOM(
     > IVEC,RHO_act(1),RHOMIN_act,NMOM,PHI1_act,IDMOM,PSMOM_act(1,1),
     > GAMMA_act(1),NGRDB2_act(1),B2_act(1),ZERROR_act(1),IER)
 
! copy back outputs if modified.
 
      where(RHO_act.ne.RHO_ref)
         R4_RHO = RHO_act
      end where
      DEALLOCATE(RHO_act)
      DEALLOCATE(RHO_ref)
 
      if(RHOMIN_act.ne.RHOMIN_ref) then
         R4_RHOMIN = RHOMIN_act
      endif
 
      if(PHI1_act.ne.PHI1_ref) then
         R4_PHI1 = PHI1_act
      endif
 
      where(PSMOM_act.ne.PSMOM_ref)
         R4_PSMOM = PSMOM_act
      end where
      DEALLOCATE(PSMOM_act)
      DEALLOCATE(PSMOM_ref)
 
      where(GAMMA_act.ne.GAMMA_ref)
         R4_GAMMA = GAMMA_act
      end where
      DEALLOCATE(GAMMA_act)
      DEALLOCATE(GAMMA_ref)
 
      where(NGRDB2_act.ne.NGRDB2_ref)
         R4_NGRDB2 = NGRDB2_act
      end where
      DEALLOCATE(NGRDB2_act)
      DEALLOCATE(NGRDB2_ref)
 
      where(B2_act.ne.B2_ref)
         R4_B2 = B2_act
      end where
      DEALLOCATE(B2_act)
      DEALLOCATE(B2_ref)
 
      where(ZERROR_act.ne.ZERROR_ref)
         R4_ZERROR = ZERROR_act
      end where
      DEALLOCATE(ZERROR_act)
      DEALLOCATE(ZERROR_ref)
 
! exit
      return
      end
