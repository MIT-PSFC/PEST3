subroutine eqmom_rcos(x,nx,inorm,im,momarr,iwarn)
 
  !  fetch im'th R cosine moment profile
  !  R(x,th) = R0(x)+sum[im=1 to kmom](Rc(im;x)*cos(i*th)+Rs(im;x)*sin(i*th))
  !  Z(x,th) = Z0(x)+sum[im=1 to kmom](Zc(im;x)*cos(i*th)+Zs(im;x)*sin(i*th))
  !
  !  for inorm=0 this routine returns (in momarr) Rc(im;x) or R0(x) if im=0.
  !  if im is out of range, momarr=0 and iwarn=1 are returned.
  !
  !  for inorm=1 return Rc(im;x)/x**im (a function which is extended smoothly
  !  to x=0).
 
  use xplasma_obj_instance
  use eq_module
  implicit NONE
 
  integer, intent(in) :: nx   ! no. of x values for which moment is desired.
  real*8, intent(in) :: x(nx) ! the x values
 
  integer, intent(in) :: inorm   ! =0: return moment; =1: return moment/x**im
  integer, intent(in) :: im      ! moment desired this call.
 
  real*8, intent(out) :: momarr(nx)  ! the moment coefficients returned.
  integer, intent(out) :: iwarn      ! 0=OK, -1 = im out of range.
                                     ! >0 means an some other error.
 
  !----------------------------------------
  integer :: ier,kmom
  !----------------------------------------
 
  call xmoments_kmom_get(kmom)

  iwarn=0
  if((im.lt.0).or.(im.gt.kmom)) then
     iwarn=-1
     momarr=0.0d0
  else
     call xpmom_get1(s,x,inorm,im,ier, rcmom=momarr)
     if(ier.ne.0) then
        write(lunerr,*) ' %eqmom_rcos: spline evaluation => ier=',ier
        call xplasma_error(s,ier,lunerr)
        iwarn=ier
     endif
  endif

end subroutine eqmom_rcos
 
subroutine eqmom_rsin(x,nx,inorm,im,momarr,iwarn)
 
  !  fetch im'th R sine moment profile
  !  R(x,th) = R0(x)+sum[im=1 to kmom](Rc(im;x)*cos(i*th)+Rs(im;x)*sin(i*th))
  !  Z(x,th) = Z0(x)+sum[im=1 to kmom](Zc(im;x)*cos(i*th)+Zs(im;x)*sin(i*th))
  !
  !  for inorm=0 this routine returns (in momarr) Rs(im;x)
  !  if im is out of range, momarr=0 and iwarn=1 are returned.
  !
  !  for inorm=1 return Rs(im;x)/x**im (a function which is extended smoothly
  !  to x=0).
 
  use xplasma_obj_instance
  use eq_module
  implicit NONE
 
  integer, intent(in) :: nx   ! no. of x values for which moment is desired.
  real*8, intent(in) :: x(nx) ! the x values
 
  integer, intent(in) :: inorm   ! =0: return moment; =1: return moment/x**im
  integer, intent(in) :: im      ! moment desired this call.
 
  real*8, intent(out) :: momarr(nx)  ! the moment coefficients returned.
  integer, intent(out) :: iwarn      ! 0=OK, -1 = im out of range.
                                     ! >0 means an some other error.
 
  !----------------------------------------
  integer :: ier,kmom
  !----------------------------------------
 
  call xmoments_kmom_get(kmom)

  iwarn=0
  if((im.lt.0).or.(im.gt.kmom)) then
     iwarn=-1
     momarr=0.0d0
  else
     call xpmom_get1(s,x,inorm,im,ier, rsmom=momarr)
     if(ier.ne.0) then
        write(lunerr,*) ' %eqmom_rsin: spline evaluation => ier=',ier
        call xplasma_error(s,ier,lunerr)
        iwarn=ier
     endif
  endif

end subroutine eqmom_rsin
 
subroutine eqmom_zcos(x,nx,inorm,im,momarr,iwarn)
 
  !  fetch im'th Z cosine moment profile
  !  R(x,th) = R0(x)+sum[im=1 to kmom](Rc(im;x)*cos(i*th)+Rs(im;x)*sin(i*th))
  !  Z(x,th) = Z0(x)+sum[im=1 to kmom](Zc(im;x)*cos(i*th)+Zs(im;x)*sin(i*th))
  !
  !  for inorm=0 this routine returns (in momarr) Zc(im;x) or Z0(x) if im=0.
  !  if im is out of range, momarr=0 and iwarn=1 are returned.
  !
  !  for inorm=1 return Zc(im;x)/x**im (a function which is extended smoothly
  !  to x=0).
 
  use xplasma_obj_instance
  use eq_module
  implicit NONE
 
  integer, intent(in) :: nx   ! no. of x values for which moment is desired.
  real*8, intent(in) :: x(nx) ! the x values
 
  integer, intent(in) :: inorm   ! =0: return moment; =1: return moment/x**im
  integer, intent(in) :: im      ! moment desired this call.
 
  real*8, intent(out) :: momarr(nx)  ! the moment coefficients returned.
  integer, intent(out) :: iwarn      ! 0=OK, -1 = im out of range.
                                     ! >0 means an some other error.
 
  !----------------------------------------
  integer :: ier,kmom
  !----------------------------------------
 
  call xmoments_kmom_get(kmom)

  iwarn=0
  if((im.lt.0).or.(im.gt.kmom)) then
     iwarn=-1
     momarr=0.0d0
  else
     call xpmom_get1(s,x,inorm,im,ier, zcmom=momarr)
     if(ier.ne.0) then
        write(lunerr,*) ' %eqmom_zcos: spline evaluation => ier=',ier
        call xplasma_error(s,ier,lunerr)
        iwarn=ier
     endif
  endif

end subroutine eqmom_zcos
 
subroutine eqmom_zsin(x,nx,inorm,im,momarr,iwarn)
 
  !  fetch im'th Z sine moment profile
  !  R(x,th) = R0(x)+sum[im=1 to kmom](Rc(im;x)*cos(i*th)+Rs(im;x)*sin(i*th))
  !  Z(x,th) = Z0(x)+sum[im=1 to kmom](Zc(im;x)*cos(i*th)+Zs(im;x)*sin(i*th))
  !
  !  for inorm=0 this routine returns (in momarr) Zs(im;x)
  !  if im is out of range, momarr=0 and iwarn=1 are returned.
  !
  !  for inorm=1 return Zs(im;x)/x**im (a function which is extended smoothly
  !  to x=0).
 
  use xplasma_obj_instance
  use eq_module
  implicit NONE
 
  integer, intent(in) :: nx   ! no. of x values for which moment is desired.
  real*8, intent(in) :: x(nx) ! the x values
 
  integer, intent(in) :: inorm   ! =0: return moment; =1: return moment/x**im
  integer, intent(in) :: im      ! moment desired this call.
 
  real*8, intent(out) :: momarr(nx)  ! the moment coefficients returned.
  integer, intent(out) :: iwarn      ! 0=OK, -1 = im out of range.
                                     ! >0 means an some other error.
 
  !----------------------------------------
  integer :: ier,kmom
  !----------------------------------------
 
  call xmoments_kmom_get(kmom)
 
  iwarn=0
  if((im.lt.0).or.(im.gt.kmom)) then
     iwarn=-1
     momarr=0.0d0
  else
     call xpmom_get1(s,x,inorm,im,ier, zsmom=momarr)
     if(ier.ne.0) then
        write(lunerr,*) ' %eqmom_zsin: spline evaluation => ier=',ier
        call xplasma_error(s,ier,lunerr)
        iwarn=ier
     endif
  endif

end subroutine eqmom_zsin
 
