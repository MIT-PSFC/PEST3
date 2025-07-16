subroutine eq_fRZ(ivec,zR,zZ,nlist,iflist,ivecd,zans,ierr)

  !  evaluate 1 or more functions f(R,Z)
  !  see also eq_gRZ (below) to evaluate 1st derivatives

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  input:

  integer ivec                      ! input vector dimension
  REAL*8 zR(ivec)                   ! argument -- where to evaluate
  REAL*8 zZ(ivec)                   ! argument -- where to evaluate
  integer nlist                     ! number of functions
  integer iflist(nlist)             ! functions to evaluate

  integer ivecd                     ! output vector dimension

  !  output:

  REAL*8 zans(ivecd,nlist)          ! evaluation results
  integer ierr                      ! completion code (0=OK)

  integer :: iRcoord,iZcoord
  !--------------------------------------------

  call eqi_fRZ_errck(nlist,iflist,iRcoord,iZcoord,ierr)
  if(ierr.ne.0) return

  call xplasma_eval_prof(s,iflist(1:nlist), &
       iRcoord,zR, iZcoord,zZ, zans(1:ivec,1:nlist), &
       ierr)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?xplasma_eval_prof error in eq_frz.f90:'
     call xplasma_error(s,ierr,lunerr)
     zans = 0
  endif

end subroutine eq_fRZ

!===================================================================
subroutine eq_gRZ(ivec,zR,zZ,nlist,iflist,ivecd,zans,ierr)

  !  evaluate 1 or more functions f(R,Z) derivatives df/dR, df/dZ

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  input:

  integer ivec                      ! input vector dimension
  REAL*8 zR(ivec)                   ! argument -- where to evaluate
  REAL*8 zZ(ivec)                   ! argument -- where to evaluate
  integer nlist                     ! number of functions
  integer iflist(nlist)             ! functions to evaluate
  
  integer ivecd                     ! output vector dimension

  !  output:

  REAL*8 zans(ivecd,2,nlist)        ! evaluation results
  integer ierr                      ! completion code (0=OK)

  !  zans(1:ivec,1,j) = df/dR for function id iflist(j)
  !  zans(1:ivec,2,j) = df/dZ for function id iflist(j)
  !--------------------------------------------

  integer :: iRcoord,iZcoord
  integer :: ii,jj
  integer :: iflist2(2*nlist),iRderiv(2*nlist),iZderiv(2*nlist)
  real*8, dimension(:,:), allocatable :: zansb

  !--------------------------------------------

  call eqi_fRZ_errck(nlist,iflist,iRcoord,iZcoord,ierr)
  if(ierr.ne.0) return

  allocate(zansb(ivec,2*nlist))

  jj=0
  do ii=1,nlist
     jj=jj+1
     iflist2(jj)=iflist(ii)
     iRderiv(jj)=1
     iZderiv(jj)=0
     jj=jj+1
     iflist2(jj)=iflist(ii)
     iRderiv(jj)=0
     iZderiv(jj)=1
  enddo

  call xplasma_eval_prof(s,iflist2, &
       iRcoord,zR, iZcoord,zZ, zansb, &
       ierr, &
       ideriv1s=iRderiv, ideriv2s=iZderiv)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?xplasma_eval_prof error in eq_frz.f90:'
     call xplasma_error(s,ierr,lunerr)
     zans = 0

  else
     ierr = 0

     jj=0
     do ii=1,nlist
        jj=jj+1
        zans(1:ivec,1,ii)=zansb(1:ivec,jj)
        jj=jj+1
        zans(1:ivec,2,ii)=zansb(1:ivec,jj)
     enddo
  endif

end subroutine eq_gRZ

subroutine eqi_fRZ_errck(nlist,iflist,iRcoord,iZcoord,ierr)

  !  check that function list arguments are all ok

  !  the first function on the list is probed to determine the [R,Z]
  !  coordinate pair that is in use for this evaluation; the pair may
  !  be reached indirectly through an "associated" profile ID

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  input:
  integer nlist                     ! no. of fcns
  integer iflist(nlist)             ! fcn id list

  !  output:
  integer iRcoord                   ! R coordinate ID
  integer iZcoord                   ! Z coordinate ID
  integer ierr                      ! completion code: 0=OK

  !------------------------------------------------------
  integer ii,id,idg1,idg2,irank,icoord1,icoord2,idgR,idgZ,iertmp
  integer ida,iranka,idg1a,idg2a
  character*32 fname,cname1,cname2
  !------------------------------------------------------

  iRcoord=0
  iZcoord=0

  id = iflist(1)

  ! look for grids & associated profile ID

  call xplasma_prof_info(s,id,ierr,rank=irank,gridId1=idg1,gridId2=idg2, &
       profId1=ida)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqi_frz_errck (eq_frz.f90): unexpected info error 1:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  ! if primary profile is not rank 2, try associated profile

  if(irank.ne.2) then
     if(ida.ne.0) then
        id=ida
        ida=0
        call xplasma_prof_info(s,id,ierr,rank=irank,gridId1=idg1,gridId2=idg2)
     endif
  endif

  ! if associated profile ID is non-zero, then primary profile ID is of rank 2
  ! but still might want associated profile ID data

  if(ida.ne.0) then
     ! rank=2 f(rho,theta) could have associated f(R,Z) @ida
     call xplasma_prof_info(s,ida,ierr,rank=iranka,gridId1=idg1a,gridId2=idg2a)
     if(iranka.ne.2) ida=0  ! ignore assoc. ID if not rank 2
  endif

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqi_frz_errck (eq_frz.f90): unexpected info error 2:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  if(irank.ne.2) then
     write(lunerr,*) &
          ' ?eqi_frz_errck (eq_frz.f90): not a function of (R,Z): rank=', &
          irank
     ierr=1

  else
     call xplasma_grid_info(s,idg1,iertmp, coord=icoord1)
     ierr=max(ierr,iertmp)
     if((ierr.eq.0).and.(ida.ne.0)) then
        
        ! decide now whether to use ID or IDA data

        call xplasma_coord_info(s,icoord1,iertmp, name=cname1)
        ierr=max(ierr,iertmp)

        if((cname1.ne.'__R').and.(cname1.ne.'__Z').and. &
             (cname1(1:4).ne.'__R_').and.(cname1(1:4).ne.'__Z_')) then

           ! switch to assoc ID grids

           idg1=idg1a
           idg2=idg2a

           call xplasma_grid_info(s,idg1,iertmp, coord=icoord1)
           ierr=max(ierr,iertmp)

        endif
     endif

     call xplasma_grid_info(s,idg2,iertmp, coord=icoord2)
     ierr=max(ierr,iertmp)

     if(ierr.eq.0) then
        call xplasma_coord_info(s,icoord1,iertmp, name=cname1)
        ierr=max(ierr,iertmp)
        call xplasma_coord_info(s,icoord2,iertmp, name=cname2)
        ierr=max(ierr,iertmp)
     endif

     if(ierr.ne.0) then
        write(lunerr,*) ' ?eqi_frz_errck (eq_frz.f90): unexpected coord error:'
        call xplasma_error(s,ierr,lunerr)
     else
        if((cname1.eq.'__R').or.(cname1(1:4).eq.'__R_')) then
           iRcoord=icoord1
        else if((cname1.eq.'__Z').or.(cname1(1:4).eq.'__Z_')) then
           iZcoord=icoord1
        else
           ierr=1
           write(lunerr,*) ' ?eqi_frz_errck: could not decipher as R or Z coordinate: ',trim(cname1)
        endif

        if((cname2.eq.'__R').or.(cname2(1:4).eq.'__R_')) then
           iRcoord=icoord2
        else if((cname2.eq.'__Z').or.(cname2(1:4).eq.'__Z_')) then
           iZcoord=icoord2
        else
           ierr=1
           write(lunerr,*) ' ?eqi_frz_errck: could not decipher as R or Z coordinate: ',trim(cname2)
        endif

        if(iRcoord.eq.0) then
           ierr=1
           write(lunerr,*) ' ?eqi_frz_errck: failed to find R dimension.'
        endif

        if(iZcoord.eq.0) then
           ierr=1
           write(lunerr,*) ' ?eqi_frz_errck: failed to find Z dimension.'
        endif
     endif
  endif

  if(ierr.ne.0) then
     call xplasma_get_item_info(s,id,iertmp, name=fname)
     write(lunerr,*) &
          ' ?eqi_frz_errck (eq_frz.f90): not a function of (R,Z): ', &
          trim(fname)
     RETURN
  endif

  !  OK -- R & Z coords of 1st function found.  All functions must use
  !        these specific R & Z coords (there can be multiple pairs of
  !        R & Z coords, functions cannot be mixed across them).

  ierr=0
  do ii=1,nlist
     id=iflist(ii)

     call xplasma_prof_gridInfo(s,id,iRcoord,idgR,ierr)
     if(ierr.ne.0) exit

     call xplasma_prof_gridInfo(s,id,iZcoord,idgZ,ierr)
     if(ierr.ne.0) exit

     if(min(idgR,idgZ).eq.0) exit
  enddo

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqi_frz_errck (eq_frz.f90): unexpected error:'
     call xplasma_error(s,ierr,lunerr)

  else if(min(idgR,idgZ).eq.0) then
     call xplasma_get_item_info(s,id,iertmp, name=fname)
     write(lunerr,*) &
          ' ?eqi_frz_errck (eq_frz.f90): not a function of (R,Z): ', &
          trim(fname)
     ierr=1

  else
     ierr=0
  endif

end subroutine eqi_fRZ_errck
