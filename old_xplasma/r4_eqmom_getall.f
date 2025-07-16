!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMOM_RCOS
      SUBROUTINE R4_EQMOM_RCOS(
     > R4_X,NX,INORM,IM,R4_MOMARR,IWARN)
 
      external EQMOM_RCOS
 
! argument declarations
      INTEGER NX
      INTEGER INORM
      INTEGER IM
      INTEGER IWARN
 ! floating type, input/output:
      REAL R4_X(NX)
 
 ! floating type, output only:
      REAL R4_MOMARR(NX)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: X_act(:)
      REAL*8, allocatable :: X_ref(:)
      REAL*8, allocatable :: MOMARR_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(X_act(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_RCOS -- X_act ALLOCATE error!')
      ALLOCATE(X_ref(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_RCOS -- X_ref ALLOCATE error!')
      ALLOCATE(MOMARR_act(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_RCOS -- MOMARR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      X_act=R4_X
      X_ref=R4_X
 
      MOMARR_act=0
! call to original routine:  EQMOM_RCOS
 
      CALL EQMOM_RCOS(
     > X_act(1),NX,INORM,IM,MOMARR_act(1),IWARN)
 
! copy back outputs if modified.
 
      where(X_act.ne.X_ref)
         R4_X = X_act
      end where
      DEALLOCATE(X_act)
      DEALLOCATE(X_ref)
      R4_MOMARR = MOMARR_act
      DEALLOCATE(MOMARR_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMOM_RSIN
      SUBROUTINE R4_EQMOM_RSIN(
     > R4_X,NX,INORM,IM,R4_MOMARR,IWARN)
 
      external EQMOM_RSIN
 
! argument declarations
      INTEGER NX
      INTEGER INORM
      INTEGER IM
      INTEGER IWARN
 ! floating type, input/output:
      REAL R4_X(NX)
 
 ! floating type, output only:
      REAL R4_MOMARR(NX)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: X_act(:)
      REAL*8, allocatable :: X_ref(:)
      REAL*8, allocatable :: MOMARR_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(X_act(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_RSIN -- X_act ALLOCATE error!')
      ALLOCATE(X_ref(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_RSIN -- X_ref ALLOCATE error!')
      ALLOCATE(MOMARR_act(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_RSIN -- MOMARR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      X_act=R4_X
      X_ref=R4_X
 
      MOMARR_act=0
! call to original routine:  EQMOM_RSIN
 
      CALL EQMOM_RSIN(
     > X_act(1),NX,INORM,IM,MOMARR_act(1),IWARN)
 
! copy back outputs if modified.
 
      where(X_act.ne.X_ref)
         R4_X = X_act
      end where
      DEALLOCATE(X_act)
      DEALLOCATE(X_ref)
      R4_MOMARR = MOMARR_act
      DEALLOCATE(MOMARR_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMOM_ZCOS
      SUBROUTINE R4_EQMOM_ZCOS(
     > R4_X,NX,INORM,IM,R4_MOMARR,IWARN)
 
      external EQMOM_ZCOS
 
! argument declarations
      INTEGER NX
      INTEGER INORM
      INTEGER IM
      INTEGER IWARN
 ! floating type, input/output:
      REAL R4_X(NX)
 
 ! floating type, output only:
      REAL R4_MOMARR(NX)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: X_act(:)
      REAL*8, allocatable :: X_ref(:)
      REAL*8, allocatable :: MOMARR_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(X_act(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_ZCOS -- X_act ALLOCATE error!')
      ALLOCATE(X_ref(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_ZCOS -- X_ref ALLOCATE error!')
      ALLOCATE(MOMARR_act(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_ZCOS -- MOMARR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      X_act=R4_X
      X_ref=R4_X
 
      MOMARR_act=0
! call to original routine:  EQMOM_ZCOS
 
      CALL EQMOM_ZCOS(
     > X_act(1),NX,INORM,IM,MOMARR_act(1),IWARN)
 
! copy back outputs if modified.
 
      where(X_act.ne.X_ref)
         R4_X = X_act
      end where
      DEALLOCATE(X_act)
      DEALLOCATE(X_ref)
      R4_MOMARR = MOMARR_act
      DEALLOCATE(MOMARR_act)
 
! exit
      return
      end
!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQMOM_ZSIN
      SUBROUTINE R4_EQMOM_ZSIN(
     > R4_X,NX,INORM,IM,R4_MOMARR,IWARN)
 
      external EQMOM_ZSIN
 
! argument declarations
      INTEGER NX
      INTEGER INORM
      INTEGER IM
      INTEGER IWARN
 ! floating type, input/output:
      REAL R4_X(NX)
 
 ! floating type, output only:
      REAL R4_MOMARR(NX)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: X_act(:)
      REAL*8, allocatable :: X_ref(:)
      REAL*8, allocatable :: MOMARR_act(:)
 
! allocation of working arrays...
 
      ALLOCATE(X_act(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_ZSIN -- X_act ALLOCATE error!')
      ALLOCATE(X_ref(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_ZSIN -- X_ref ALLOCATE error!')
      ALLOCATE(MOMARR_act(NX)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQMOM_ZSIN -- MOMARR_act ALLOCATE error!')
 
! executable code:  copy for input
 
      X_act=R4_X
      X_ref=R4_X
 
      MOMARR_act=0
! call to original routine:  EQMOM_ZSIN
 
      CALL EQMOM_ZSIN(
     > X_act(1),NX,INORM,IM,MOMARR_act(1),IWARN)
 
! copy back outputs if modified.
 
      where(X_act.ne.X_ref)
         R4_X = X_act
      end where
      DEALLOCATE(X_act)
      DEALLOCATE(X_ref)
      R4_MOMARR = MOMARR_act
      DEALLOCATE(MOMARR_act)
 
! exit
      return
      end
