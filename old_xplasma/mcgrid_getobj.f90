subroutine mcgrid_getobj(id_gridi,zname,ist_th,ist_ph,zconv,zdata,idim,ierr)
 
  !  get (fetch) mcgrid scalar field -- one number per grid zone
 
  !  the data can cover the plasma out to the boundary, or, beyond the
  !  boundary into an extrapolated space.  If the stored data goes beyond
  !  the boundary AND the user data has space for it, the beyond-the-boundary
  !  data will be copied, otherwise zeroes will be written.
 
  !  the output data is on the user's (theta,phi) grids; stored data is on
  !  the standardixed xplasma mcgrid:
  !    ist_th(ph) = -1 -- user grid starts at -pi (theta, phi)
  !    ist_th(ph) =  0 -- user grid starts at 0   (theta, phi)
 
  !  if the mcgrid is updown symmetric, only data covering the range
  !  [0,pi] or [-pi,0] are copied, depending on the setting of ist_th.
 
  !  the mcgrid# id_gridi refers to a grid set up by a prior mcgrid_define
  !  call.
 
  !  if id_gridi is negative:  if the named object does not exist,
  !  set the output array to zero without setting an error flag.
 
  !  if id_gridi is positive:  non-existance of named object is an error.
 
  use xplasma_obj_instance
  use eq_module

  implicit NONE
 
  !  NOTE: for compatibility with f77 xplasma: if zname="ZONE_VOLUME" then
  !  the name [gridname]//'_DVOL' is substituted...

  integer, intent(in) :: id_gridi    ! id from prior mcgrid_define call
  character*(*), intent(in) :: zname ! name of object to fetch
  integer, intent(in) :: ist_th,ist_ph  ! angle grid start flags
  integer, intent(in) :: idim        ! dimension of zdata array
  real*8, intent(in) :: zconv        ! units conversion factor
  !  set zconv.lt.0.0 to specify clockwise orientation of poloidal zones
  !  data scaled by a factor of abs(zconv)

  real*8 :: zdata(idim)              ! data to fetch
 
  integer, intent(out) :: ierr       ! completion code (0=OK)
 
  !----------------------------------
  integer :: idtmp,id_grid,inzons,inzonu,nzons,nzons_ext,itype_bb,iertmp
  integer, dimension(:), pointer :: ia_ptr
  character*32 :: mcname,fname
  !----------------------------------
 
  id_grid=abs(id_gridi)

  call mcgrid_argcheck('mcgrid_getobj',ist_th,ist_ph,0,ierr)
  if(ierr.ne.0) return
 
  call xplasma_mcgrid_info(s,id_grid,ierr, &
       nzons=nzons,nzons_ext=nzons_ext,name=mcname)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getobj: xplasma_mcgrid_info error.'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  fname=zname
  call uupper(fname)
  if((fname.eq.'ZONE_VOLUME').or.(fname.eq.'ZONE_VOLUMES')) then 
     fname=trim(mcname)//'_DVOL'
  endif

  call xplasma_find_item(s,fname,idtmp,ierr,nf_noerr=.TRUE.)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getobj: xplasma_find_item error.'
     call xplasma_error(s,ierr,lunerr)
     return
  endif
  if(idtmp.eq.0) then
     zdata=0.0d0
     if(id_gridi.gt.0) then
        write(lunerr,*) ' ?mcgrid_getobj: "',trim(fname),'" not found.'
        ierr=1
     endif
     return
  endif
        
  call xplasma_blackBox_retrieve(s,idtmp,ierr,ia_ptr=ia_ptr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getobj: xplasma_blackBox_retrieve error.'
     call xplasma_error(s,ierr,lunerr)
     nullify(ia_ptr)
     return
  endif

  inzons=ia_ptr(2)

  !  nzons = #interior grid pts
  !  nzons_ext = #total grid pts
  !  inzons = #in data (could be nzons or nzons_ext).
  !  idim = #in passed array

  if(idim.lt.nzons) then
     write(lunerr,*) ' ?mcgrid_getobj: array size too small: idim=',idim, &
          '; need=',inzons
     ierr=1
     nullify(ia_ptr)
     return
  else if(idim.lt.inzons) then
     inzonu=nzons
     if(idim.gt.inzonu) then
        zdata(inzonu+1:idim)=0.0d0
     endif
  else
     inzonu=inzons
     if(idim.gt.inzonu) then
        zdata(inzonu+1:idim)=0.0d0
     endif
  endif

  call xplasma_mcgrid_getobj(s,idtmp,ierr, &
       lstart0_th=(ist_th.eq.0), ccwflag_th=(zconv.ge.0.0d0), &
       data_1d=zdata(1:inzonu))

  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getobj: xplasma_mcgrid_getobj error:'
     call xplasma_error(s,ierr,lunerr)
     zdata(1:inzonu)=0.0d0
  else
     zdata(1:inzonu)=zdata(1:inzonu)*abs(zconv)
  endif

  nullify(ia_ptr)

end subroutine mcgrid_getobj

subroutine mcgrid_getobj_r4(id_gridi,zname,ist_th,ist_ph,zconv,zdata,idim,ierr)
 
  !  get (fetch) mcgrid scalar field -- one number per grid zone
  !    _r4 (single precision) version -- user's array is declared REAL.
 
  !  the data can cover the plasma out to the boundary, or, beyond the
  !  boundary into an extrapolated space.  If the stored data goes beyond
  !  the boundary AND the user data has space for it, the beyond-the-boundary
  !  data will be copied, otherwise zeroes will be written.
 
  !  the output data is on the user's (theta,phi) grids; stored data is on
  !  the standardixed xplasma mcgrid:
  !    ist_th(ph) = -1 -- user grid starts at -pi (theta, phi)
  !    ist_th(ph) =  0 -- user grid starts at 0   (theta, phi)
 
  !  if the mcgrid is updown symmetric, only data covering the range
  !  [0,pi] or [-pi,0] are copied, depending on the setting of ist_th.
 
  !  the mcgrid# id_gridi refers to a grid set up by a prior mcgrid_define
  !  call.
 
  !  if id_gridi is negative:  if the named object does not exist,
  !  set the output array to zero without setting an error flag.
 
  !  if id_gridi is positive:  non-existance of named object is an error.
 
  implicit NONE
 
  integer, intent(in) :: id_gridi     ! id from prior mcgrid_define call
  character*(*), intent(in) :: zname ! name of object to fetch
  integer, intent(in) :: ist_th,ist_ph  ! angle grid start flags
  integer, intent(in) :: idim        ! dimension of zdata array
  real, intent(in) :: zconv          ! units conversion factor
  !  set zconv.lt.0.0 to specify clockwise orientation of poloidal zones
  !  data scaled by a factor of abs(zconv)

  real :: zdata(idim)                ! data to fetch
 
  integer, intent(out) :: ierr       ! completion code (0=OK)
 
  !----------------------------------
 
  real*8 :: zconv8
  real*8, dimension(:), allocatable :: zdata8
  !----------------------------------
 
  zconv8=zconv

  allocate(zdata8(idim)); zdata=0.0d0

  call mcgrid_getobj(id_gridi,zname,ist_th,ist_ph,zconv8,zdata8,idim,ierr)

  zdata = zdata8
  deallocate(zdata8)
 
end subroutine mcgrid_getobj_r4
 
subroutine mcgrid_getarr2(id_gridi,zname,ist_th,ist_ph, &
     zconv,zdata,idim1,idim2,idim3,in1,in2,ierr)
 
  !  get (fetch) mcgrid 2d array field -- array(in1,in2) per grid zone
 
  !  the data can cover the plasma out to the boundary, or, beyond the
  !  boundary into an extrapolated space.  If the stored data goes beyond
  !  the boundary AND the user data has space for it, the beyond-the-boundary
  !  data will be copied, otherwise zeroes will be written.
 
  !  the output data is on the user's (theta,phi) grids; stored data is on
  !  the standardixed xplasma mcgrid:
  !    ist_th(ph) = -1 -- user grid starts at -pi (theta, phi)
  !    ist_th(ph) =  0 -- user grid starts at 0   (theta, phi)
 
  !  if the mcgrid is updown symmetric, only data covering the range
  !  [0,pi] or [-pi,0] are copied, depending on the setting of ist_th.
 
  !  the mcgrid# id_grid refers to a grid set up by a prior mcgrid_define
  !  call.
 
  !  if id_gridi is negative:  if the named object does not exist,
  !  set the output array to zero without setting an error flag.
 
  !  if id_gridi is positive:  non-existance of named object is an error.
 
  use xplasma_obj_instance
  use eq_module

  implicit NONE
 
  integer, intent(in) :: id_gridi     ! id from prior mcgrid_define call
  character*(*), intent(in) :: zname ! name of object to fetch
  integer, intent(in) :: ist_th,ist_ph  ! angle grid start flags
  integer, intent(in) :: idim1,idim2,idim3        ! dimension of zdata array
  real*8, intent(in) :: zconv        ! units conversion factor
  !  set zconv.lt.0.0 to specify clockwise orientation of poloidal zones
  !  data scaled by a factor of abs(zconv)

  real*8 :: zdata(idim1,idim2,idim3) ! data to fetch
  integer, intent(inout) :: in1,in2  ! activ dimensions zdata(1:in1,1:in2,j)
 
  integer, intent(out) :: ierr       ! completion code (0=OK)
 
  !----------------------------------
  integer :: idtmp,id_grid,inzons,nzons_ext,itype_bb,iertmp,ii,jj,kk
  integer :: inzonu,nzons
  integer :: irank
  integer, dimension(:), pointer :: ia_ptr
  !----------------------------------
  
  id_grid=abs(id_gridi)

  call mcgrid_argcheck('mcgrid_getarr2',ist_th,ist_ph,0,ierr)
  if(ierr.ne.0) return
  
  call xplasma_mcgrid_info(s,id_grid,ierr, &
       nzons=nzons,nzons_ext=nzons_ext)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getarr2: xplasma_mcgrid_info error.'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  call xplasma_find_item(s,zname,idtmp,ierr,nf_noerr=.TRUE.)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getarr2: xplasma_find_item error.'
     call xplasma_error(s,ierr,lunerr)
     return
  endif
  if(idtmp.eq.0) then
     zdata=0.0d0
     if(id_gridi.gt.0) then
        write(lunerr,*) ' ?mcgrid_getarr2: "',trim(zname),'" not found.'
        ierr=1
     endif
     return
  endif
        
  call xplasma_blackBox_retrieve(s,idtmp,ierr,ia_ptr=ia_ptr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getarr2: xplasma_blackBox_retrieve error.'
     call xplasma_error(s,ierr,lunerr)
     ierr=1
     nullify(ia_ptr)
     return
  endif

  inzons=ia_ptr(2)

  !  nzons = #interior grid pts
  !  nzons_ext = #total grid pts
  !  inzons = #in data (could be nzons or nzons_ext).
  !  idim3 = #in passed array

  if(idim3.lt.nzons) then
     write(lunerr,*) ' ?mcgrid_getarr2: array idim3 too small: idim3=',idim3, &
          '; need=',inzons
     ierr=1
     nullify(ia_ptr)
     return
  else if(idim3.lt.inzons) then
     inzonu=nzons
  else
     inzonu=inzons
  endif

  irank=ia_ptr(3)
  if(irank.ne.2) then
     write(lunerr,*) ' ?mcgrid_getarr2: expected rank2 MCgrid array; irank=',&
          irank
     ierr=1
     nullify(ia_ptr)
     return
  endif

  if(in1.le.0) in1=ia_ptr(4)
  if(in2.le.0) in2=ia_ptr(5)

  if(idim1.lt.in1) then
     write(lunerr,*) ' ?mcgrid_getarr2: array idim1 too small: idim1=',idim1, &
          '; need=',in1
     ierr=1
  else if(in1.ne.ia_ptr(4)) then
     write(lunerr,*) ' ?mcgrid_getarr2: 1st dimension of data is: ',ia_ptr(4)
     write(lunerr,*) '  subroutine argument (in1) says it should be: ',in1
     ierr=1
  endif

  if(idim2.lt.in2) then
     write(lunerr,*) ' ?mcgrid_getarr2: array idim2 too small: idim2=',idim2, &
          '; need=',in2
     ierr=1
  else if(in2.ne.ia_ptr(5)) then
     write(lunerr,*) ' ?mcgrid_getarr2: 1st dimension of data is: ',ia_ptr(5)
     write(lunerr,*) '  subroutine argument (in2) says it should be: ',in2
     ierr=1
  endif

  if(ierr.ne.0) then
     nullify(ia_ptr)
     return
  endif

  call xplasma_mcgrid_getobj(s,idtmp,ierr, &
       lstart0_th=(ist_th.eq.0), ccwflag_th=(zconv.ge.0.0d0), &
       data_3d=zdata(1:in1,1:in2,1:inzonu))

  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getarr2: xplasma_mcgrid_getobj error:'
     call xplasma_error(s,ierr,lunerr)
     zdata=0.0d0
  else
     do kk=1,idim3
        do jj=1,idim2
           do ii=1,idim1
              if((ii.le.in1).and.(jj.le.in2).and.(kk.le.inzons)) then
                 zdata(ii,jj,kk)=zdata(ii,jj,kk)*abs(zconv)
              else
                 zdata(ii,jj,ii)=0.0d0
              endif
           enddo
        enddo
     enddo
  endif

  nullify(ia_ptr)

end subroutine mcgrid_getarr2
 
subroutine mcgrid_getarr2_r4(id_gridi,zname,ist_th,ist_ph, &
     zconv,zdata,idim1,idim2,idim3,in1,in2,ierr)
 
  !  get (fetch) mcgrid 2d array field -- array(in1,in2) per grid zone
  !    _r4 (single precision) version -- user's array is declared REAL.
 
  !  the data can cover the plasma out to the boundary, or, beyond the
  !  boundary into an extrapolated space.  If the stored data goes beyond
  !  the boundary AND the user data has space for it, the beyond-the-boundary
  !  data will be copied, otherwise zeroes will be written.
 
  !  the output data is on the user's (theta,phi) grids; stored data is on
  !  the standardixed xplasma mcgrid:
  !    ist_th(ph) = -1 -- user grid starts at -pi (theta, phi)
  !    ist_th(ph) =  0 -- user grid starts at 0   (theta, phi)
 
  !  if the mcgrid is updown symmetric, only data covering the range
  !  [0,pi] or [-pi,0] are copied, depending on the setting of ist_th.
 
  !  the mcgrid# id_gridi refers to a grid set up by a prior mcgrid_define
  !  call.
 
  !  if id_gridi is negative:  if the named object does not exist,
  !  set the output array to zero without setting an error flag.
 
  !  if id_gridi is positive:  non-existance of named object is an error.
 
  implicit NONE
 
  integer, intent(in) :: id_gridi     ! id from prior mcgrid_define call
  character*(*), intent(in) :: zname ! name of object to fetch
  integer, intent(in) :: ist_th,ist_ph  ! angle grid start flags
  integer, intent(in) :: idim1,idim2,idim3        ! dimension of zdata array
  real, intent(in) :: zconv          ! units conversion factor
  !  set zconv.lt.0.0 to specify clockwise orientation of poloidal zones
  !  data scaled by a factor of abs(zconv)

  real :: zdata(idim1,idim2,idim3)   ! data to fetch
  integer, intent(inout) :: in1,in2  ! activ dimensions zdata(1:in1,1:in2,j)
 
  integer, intent(out) :: ierr       ! completion code (0=OK)
 
  !----------------------------------
 
  real*8 :: zconv8
  real*8, dimension(:,:,:), allocatable :: zdata8
 
  !----------------------------------
 
  zconv8 = zconv

  allocate(zdata8(idim1,idim2,idim3)); zdata8=0.0d0

  call mcgrid_getarr2(id_gridi,zname,ist_th,ist_ph, &
       zconv8,zdata8,idim1,idim2,idim3,in1,in2,ierr)

  zdata = zdata8
  deallocate(zdata8)
 
end subroutine mcgrid_getarr2_r4
