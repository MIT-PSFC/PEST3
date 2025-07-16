subroutine mcgrid_putobj(id_grid,zname, &
     ist_th,ist_ph,iextend,zmks,zdata,idimi,ierr)
 
  !  (dmc: intent(in) attribute removed from zdata array arguments)
  !  (to make debugger happy: mcgrid_put uses intent(inout)).

  !  put (store) mcgrid scalar field -- one number per grid zone
 
  !  this can cover the plasma out to the boundary (iextend=0) or
  !  beyond the boundary into the extrapolated region (iextend=1)
 
  !  the input data's theta and phi grid are as per user's definitions;
  !  ist_th(ph) = -1 -- user grid starts at -pi (theta, phi)
  !  ist_th(ph) =  0 -- user grid starts at 0   (theta, phi)
  !  theta, phi grids are evenly spaced; no. of phi points and no. of
  !  theta pts per zone row are as specified in mcgrid #id_grid
  !  set up by prior call to mcgrid_define.
 
  !  specify a negative number for "idimi" if object replacement
  !  is to generate a warning.

  use xplasma_obj_instance
  use eq_module

  implicit NONE
 
  integer, intent(in) :: id_grid     ! id from prior mcgrid_define call
  character*(*), intent(in) :: zname ! (unique) name of object to store
  integer, intent(in) :: ist_th,ist_ph  ! angle grid start flags
  integer, intent(in) :: iextend     ! =0: nznbmri rows; =1: nznbmr rows
  integer, intent(in) :: idimi       ! dimension of zdata array

  real*8, intent(in) :: zmks         ! MKS conversion factor to apply to zdata
  !  set zmks.lt.0.0 to specify clockwise orientation of poloidal zones
  !  data scaled by a factor of abs(zmks)

  real*8 :: zdata(abs(idimi))             ! data to store
 
  integer, intent(out) :: ierr       ! completion code (0=OK)
 
  ! -------------------------------------
  integer :: idtmp,idima,nzons,nzons_ext,inzons,iertmp
  ! -------------------------------------
 
  call mcgrid_argcheck('mcgrid_putobj',ist_th,ist_ph,iextend,ierr)
  if(ierr.ne.0) return

  call xplasma_mcgrid_info(s,id_grid,ierr, nzons=nzons,nzons_ext=nzons_ext)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_putobj: xplasma_mcgrid_info error.'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  if(idimi.lt.0) then
     call xplasma_find_item(s,zname,idtmp,ierr,nf_noerr=.TRUE.)
     if(ierr.ne.0) then
        write(lunerr,*) ' ?mcgrid_putobj: xplasma_find_item error.'
        call xplasma_error(s,ierr,lunerr)
        return
     endif
     if(idtmp.gt.0) then
        write(lunerr,*) &
             ' %mcgrid_putobj: warning: "',trim(zname),'" redefined.'
     endif
  endif

  idima=abs(idimi)
  if(iextend.eq.1) then
     inzons=nzons_ext
  else
     inzons=nzons
  endif

  if(inzons.gt.idima) then
     write(lunerr,*) &
          ' ?mcgrid_putobj: array size too small for iextend=',iextend
     write(lunerr,*) '  have: ',idima,'; need: ',inzons
     ierr=1
     return
  endif

  call xoi_author_set(iertmp)

  call xplasma_mcgrid_putobj(s,id_grid,zname,idtmp,ierr, &
       lstart0_th = (ist_th.eq.0), ccwflag_th = (zmks.ge.0.0d0), &
       data_1d = zdata(1:inzons)*abs(zmks))

  call xoi_author_clear(iertmp)
 
end subroutine mcgrid_putobj
 
subroutine mcgrid_putobj_r4(id_grid,zname, &
     ist_th,ist_ph,iextend,zmks,zdata,idimi,ierr)
 
  !  put (store) mcgrid scalar field -- one number per grid zone
 
  !  this can cover the plasma out to the boundary (iextend=0) or
  !  beyond the boundary into the extrapolated region (iextend=1)
 
  !  the input data's theta and phi grid are as per user's definitions;
  !  ist_th(ph) = -1 -- user grid starts at -pi (theta, phi)
  !  ist_th(ph) =  0 -- user grid starts at 0   (theta, phi)
  !  theta, phi grids are evenly spaced; no. of phi points and no. of
  !  theta pts per zone row are as specified in mcgrid #id_grid
  !  set up by prior call to mcgrid_define.
 
  implicit NONE
 
  integer, intent(in) :: id_grid     ! id from prior mcgrid_define call
  character*(*), intent(in) :: zname ! (unique) name of object to store
  integer, intent(in) :: ist_th,ist_ph  ! angle grid start flags
  integer, intent(in) :: iextend     ! =0: nznbmri rows; =1: nznbmr rows
  integer, intent(in) :: idimi       ! dimension of zdata array
  real, intent(in) :: zmks           ! MKS conversion factor to apply to zdata
  !  set zmks.lt.0.0 to specify clockwise orientation of poloidal zones
  !  data scaled by a factor of abs(zmks)

  real :: zdata(abs(idimi))                ! data to store
 
  integer, intent(out) :: ierr       ! completion code (0=OK)
  ! -------------------------------------
  ! local
 
  real*8 zmks8
  real*8, dimension(:), allocatable :: zdata8
  ! -------------------------------------

  allocate(zdata8(abs(idimi)))

  zmks8 = zmks
  zdata8 = zdata

  call mcgrid_putobj(id_grid,zname, &
     ist_th,ist_ph,iextend,zmks8,zdata8,idimi,ierr)

  deallocate(zdata8)

end subroutine mcgrid_putobj_r4
 
subroutine mcgrid_putarr2(id_grid,zname, &
     ist_th,ist_ph,iextend,zmks,zdata,in1,in2,idimi1,idim2,idim3,ierr)
 
  !  put (store) mcgrid 2d array field -- array(in1,in2) per grid zone
 
  !  this can cover the plasma out to the boundary (iextend=0) or
  !  beyond the boundary into the extrapolated region (iextend=1)
 
  !  the input data's theta and phi grid are as per user's definitions;
  !  ist_th(ph) = -1 -- user grid starts at -pi (theta, phi)
  !  ist_th(ph) =  0 -- user grid starts at 0   (theta, phi)
  !  theta, phi grids are evenly spaced; no. of phi points and no. of
  !  theta pts per zone row are as specified in mcgrid #id_grid
  !  set up by prior call to mcgrid_define.
 
  use xplasma_obj_instance
  use eq_module

  implicit NONE
 
  integer, intent(in) :: id_grid     ! id from prior mcgrid_define call
  character*(*), intent(in) :: zname ! (unique) name of object to store
  integer, intent(in) :: ist_th,ist_ph  ! angle grid start flags
  integer, intent(in) :: iextend     ! =0: nznbmri rows; =1: nznbmr rows
  integer, intent(in) :: idimi1,idim2,idim3 ! dimensions of zdata array
  real*8, intent(in) :: zmks         ! MKS conversion factor to apply to zdata
  !  set zmks.lt.0.0 to specify clockwise orientation of poloidal zones
  !  data scaled by a factor of abs(zmks)

  real*8 :: zdata(abs(idimi1),idim2,idim3) ! data to store
  integer, intent(in) :: in1,in2     ! store:  zdata(1:in1,1:in2,j)
 
  integer, intent(out) :: ierr       ! completion code (0=OK)
 
  ! -------------------------------------
  integer :: idtmp,idima,nzons,nzons_ext,inzons,ii,jj,kk,iertmp
  real*8, dimension(:,:,:), allocatable :: zdbuf
  ! -------------------------------------
 
  call mcgrid_argcheck('mcgrid_putarr2',ist_th,ist_ph,iextend,ierr)
  if(ierr.ne.0) return

  call xplasma_mcgrid_info(s,id_grid,ierr, nzons=nzons,nzons_ext=nzons_ext)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_putarr2: xplasma_mcgrid_info error.'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  if(idimi1.lt.0) then
     call xplasma_find_item(s,zname,idtmp,ierr,nf_noerr=.TRUE.)
     if(ierr.ne.0) then
        write(lunerr,*) ' ?mcgrid_putarr2: xplasma_find_item error.'
        call xplasma_error(s,ierr,lunerr)
        return
     endif
     if(idtmp.gt.0) then
        write(lunerr,*) &
             ' %mcgrid_putarr2: warning: "',trim(zname),'" redefined.'
     endif
  endif

  idima=idim3
  if(iextend.eq.1) then
     inzons=nzons_ext
  else
     inzons=nzons
  endif

  if(inzons.gt.idima) then
     write(lunerr,*) &
          ' ?mcgrid_putarr2: array size too small for iextend=',iextend
     write(lunerr,*) '  have: ',idima,'; need: ',inzons
     ierr=1
     return
  endif

  allocate(zdbuf(in1,in2,inzons))

  do kk=1,inzons
     do jj=1,in2
        do ii=1,in1
           zdbuf(ii,jj,kk) = abs(zmks)*zdata(ii,jj,kk)
        enddo
     enddo
  enddo

  call xoi_author_set(iertmp)

  call xplasma_mcgrid_putobj(s,id_grid,zname,idtmp,ierr, &
       lstart0_th = (ist_th.eq.0), ccwflag_th = (zmks.ge.0.0d0), &
       data_3d = zdbuf)

  call xoi_author_clear(iertmp)

  deallocate(zdbuf)

end subroutine mcgrid_putarr2
 
subroutine mcgrid_putarr2_r4(id_grid,zname, &
     ist_th,ist_ph,iextend,zmks,zdata,in1,in2,idimi1,idim2,idim3,ierr)
 
  !  put (store) mcgrid 2d array field -- array(in1,in2) per grid zone
 
  !  this can cover the plasma out to the boundary (iextend=0) or
  !  beyond the boundary into the extrapolated region (iextend=1)
 
  !  the input data's theta and phi grid are as per user's definitions;
  !  ist_th(ph) = -1 -- user grid starts at -pi (theta, phi)
  !  ist_th(ph) =  0 -- user grid starts at 0   (theta, phi)
  !  theta, phi grids are evenly spaced; no. of phi points and no. of
  !  theta pts per zone row are as specified in mcgrid #id_grid
  !  set up by prior call to mcgrid_define.
 
  implicit NONE
 
  integer, intent(in) :: id_grid     ! id from prior mcgrid_define call
  character*(*), intent(in) :: zname ! (unique) name of object to store
  integer, intent(in) :: ist_th,ist_ph  ! angle grid start flags
  integer, intent(in) :: iextend     ! =0: nznbmri rows; =1: nznbmr rows
  integer, intent(in) :: idimi1,idim2,idim3 ! dimensions of zdata array
  real, intent(in) :: zmks           ! MKS conversion factor to apply to zdata
  !  set zmks.lt.0.0 to specify clockwise orientation of poloidal zones
  !  data scaled by a factor of abs(zmks)

  real :: zdata(abs(idimi1),idim2,idim3)   ! data to store
  integer, intent(in) :: in1,in2     ! store:  zdata(1:in1,1:in2,j)
 
  integer, intent(out) :: ierr       ! completion code (0=OK)
 
  ! -------------------------------------
  ! local
 
  real*8 :: zmks8
  real*8, dimension(:,:,:), allocatable :: zdata8
  integer :: idima1
  ! -------------------------------------

  idima1=abs(idimi1)
  allocate(zdata8(idima1,idim2,idim3))

  zmks8 = zmks
  zdata8 = zdata

  call mcgrid_putarr2(id_grid,zname, &
     ist_th,ist_ph,iextend,zmks8,zdata8,in1,in2,idimi1,idim2,idim3,ierr)

  deallocate(zdata8)

end subroutine mcgrid_putarr2_r4
