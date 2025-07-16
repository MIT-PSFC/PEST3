subroutine eqi_rzbox_binfo_set(zrho,zrzbuf,id,ierr)

  !  store bdy R & Z data in a labeled xplasma black box

  use xplasma_definitions
  use eqi_rzbox_module
  implicit NONE

  real*8, intent(in) :: zrho
  real*8 :: zrzbuf(nsrch1*4)
  integer, intent(out) :: id
  integer, intent(out) :: ierr

  !-------------------------
  integer :: idata(1),iertmp
  real*8 :: zrho_bdy
  real*8, dimension(:), pointer :: zrz_ptr
  !-------------------------

  idata = nsrch1

  call xplasma_author_set(sp, xplasma_xmhd, ierr)
  if(ierr.ne.0) return

  call eqi_extrap_rho_bdy(zrho_bdy)

  if(abs(zrho-1.0d0).lt.1.0d-14) then

     !  create

     allocate(zrz_ptr(4*nsrch1+1))
     zrz_ptr(1:4*nsrch1)=zrzbuf
     zrz_ptr(4*nsrch1+1)=zrho
     call xplasma_create_blackbox(sp,'__RZ_BDY_INFO',99,id,ierr,&
          iarray=idata, r8array=zrz_ptr, &
          label='(R,Z) plasma bdy search contour', units='m')
     deallocate(zrz_ptr)
  else if(abs(zrho-zrho_bdy).lt.1.0d-14) then

     !  create

     allocate(zrz_ptr(4*nsrch1+1))
     zrz_ptr(1:4*nsrch1)=zrzbuf
     zrz_ptr(4*nsrch1+1)=zrho
     call xplasma_create_blackbox(sp,'__RZ_EXTRAP_BDY_INFO',99,id,ierr, &
          iarray=idata, r8array=zrz_ptr, &
          label='(R,Z) extrapolated bdy search contour', units='m')
     deallocate(zrz_ptr)
  else

     ! update if exists; otherwise: create...

     call xplasma_find_item(sp,'__RZ_SURF_INFO',id,ierr,nf_noerr=.TRUE.)
     if(id.ne.0) then
        call xplasma_blackbox_retrieve(sp,id,ierr, r8a_ptr = zrz_ptr)
        zrz_ptr(1:4*nsrch1) = zrzbuf
        zrz_ptr(4*nsrch1+1) = zrho
     else
        allocate(zrz_ptr(4*nsrch1+1))
        zrz_ptr(1:4*nsrch1)=zrzbuf
        zrz_ptr(4*nsrch1+1)=zrho
        call xplasma_create_blackbox(sp,'__RZ_SURF_INFO',99,id,ierr, &
             iarray=idata, r8array=zrz_ptr, &
             label='(R,Z) plasma surface search contour', units='m')
        deallocate(zrz_ptr)
     endif
  endif

  call xplasma_author_clear(sp, xplasma_xmhd, iertmp)

  nullify(zrz_ptr)

end subroutine eqi_rzbox_binfo_set

subroutine eqi_rzbox_binfo_get(zrho,zrzbuf,ierr)

  !  restore bdy R & Z data from xplasma black box

  use xplasma_definitions
  use eqi_rzbox_module
  implicit NONE

  real*8, intent(in) :: zrho  ! expect 1.00 for bdy;  >1.00 for extrap. bdy.
  real*8 :: zrzbuf(nsrch1*4)
  integer, intent(out) :: ierr

  !-------------------------
  integer :: idata(1),id,icheck
  real*8 :: zrho_bdy
  real*8, dimension(:), pointer :: zrz_ptr
  !-------------------------

  call eqi_extrap_rho_bdy(zrho_bdy)

  if(abs(zrho-1.0d0).lt.1.0d-14) then
     call xplasma_find_item(sp,'__RZ_BDY_INFO',id,ierr)
     icheck=0
  else if(abs(zrho-zrho_bdy).lt.1.0d-14) then
     call xplasma_find_item(sp,'__RZ_EXTRAP_BDY_INFO',id,ierr)
     icheck=0
  else
     call xplasma_find_item(sp,'__RZ_SURF_INFO',id,ierr,nf_noerr=.TRUE.)
     icheck=1
  endif

  if(id.eq.0) ierr=1

  if(ierr.ne.0) then
     zrzbuf = 0
     return
  endif

  call xplasma_blackbox_retrieve(sp,id,ierr, iarray=idata, r8a_ptr=zrz_ptr)
  if(icheck.eq.1) then
     if(zrz_ptr(nsrch1*4+1).ne.zrho) then
        ierr=1
        zrzbuf=0
     endif
  endif
  if(ierr.eq.0) then
     zrzbuf=zrz_ptr(1:nsrch1*4)
  endif

  nullify(zrz_ptr)

end subroutine eqi_rzbox_binfo_get
