!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQM_RZMAG
      SUBROUTINE R4_EQM_RZMAG(
     > R4_RARR,R4_ZARR,ID1,ID2,IDRHO,ID_R,ID_Z,IERR)
 
      external EQM_RZMAG
 
! argument declarations
      INTEGER ID1
      INTEGER ID2
      INTEGER IDRHO
      INTEGER ID_R
      INTEGER ID_Z
      INTEGER IERR
 ! floating type, input/output:
      REAL R4_RARR(ID1,ID2)
 
 ! floating type, input/output:
      REAL R4_ZARR(ID1,ID2)
 
 
! local (automatic array) declarations
 
      integer alloc_stat
      REAL*8, allocatable :: RARR_act(:,:)
      REAL*8, allocatable :: RARR_ref(:,:)
      REAL*8, allocatable :: ZARR_act(:,:)
      REAL*8, allocatable :: ZARR_ref(:,:)
 
! allocation of working arrays...
 
      ALLOCATE(RARR_act(ID1,ID2)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZMAG -- RARR_act ALLOCATE error!')
      ALLOCATE(RARR_ref(ID1,ID2)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZMAG -- RARR_ref ALLOCATE error!')
      ALLOCATE(ZARR_act(ID1,ID2)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZMAG -- ZARR_act ALLOCATE error!')
      ALLOCATE(ZARR_ref(ID1,ID2)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQM_RZMAG -- ZARR_ref ALLOCATE error!')
 
! executable code:  copy for input
 
      RARR_act=R4_RARR
      RARR_ref=R4_RARR
 
      ZARR_act=R4_ZARR
      ZARR_ref=R4_ZARR
 
! call to original routine:  EQM_RZMAG
 
      CALL EQM_RZMAG(
     > RARR_act(1,1),ZARR_act(1,1),ID1,ID2,IDRHO,ID_R,ID_Z,IERR)
 
! copy back outputs if modified.
 
      where(RARR_act.ne.RARR_ref)
         R4_RARR = RARR_act
      end where
      DEALLOCATE(RARR_act)
      DEALLOCATE(RARR_ref)
 
      where(ZARR_act.ne.ZARR_ref)
         R4_ZARR = ZARR_act
      end where
      DEALLOCATE(ZARR_act)
      DEALLOCATE(ZARR_ref)
 
! exit
      return
      end
