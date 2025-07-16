subroutine r8bsmoo_init(imj,icentr,iedge,zsmoo,zxi)
 
  use r8bsmoo_mod
 
  implicit NONE
 
  integer, intent(in) :: imj   ! grid array size
  integer, intent(in) :: icentr,iedge   ! active grid limits
  real*8, intent(in) :: zsmoo    ! smoothing parameter
  real*8, intent(in) :: zxi(imj,2)        ! the grid (TRCOM style)
 
  !  ***** initialize the module *****
 
  mj=imj
  lcentr=icentr
  ledge=iedge
 
  if(allocated(xi)) deallocate(xi)
  allocate(xi(mj,2))
 
  xi = zxi
 
  dxbsmoo = zsmoo
 
  lcp1 = lcentr+1
  lep1 = ledge+1
 
  nzones = ledge-lcentr + 1
 
end subroutine r8bsmoo_init
 
subroutine r8bsmoo_get(zsmoo)
 
  use r8bsmoo_mod
 
  implicit NONE
 
  real*8, intent(out) :: zsmoo
 
  !------------------
 
  zsmoo = dxbsmoo
 
end subroutine r8bsmoo_get
