subroutine eqi_geq_axis_zsrch(ivec,iok,Z,dPsidZ,ivecd, &
     zinput,ninput,zoutput,noutput)

  ! see comput/zridderx.for -- root finder functional
  ! return Z s.t. dPsidZ=0, along fixed R (see eqi_geq_axis_Rsrch.f90)

  use eqi_geq_axis_data
  implicit NONE

  !-----------------------------

  integer,intent(in) :: ivec
  logical,intent(in) :: iok(ivec)
  real*8, intent(in) :: Z(ivec)
  real*8,intent(out) :: dPsidZ(ivec)

  integer,intent(in) :: ivecd
  integer,intent(in) :: ninput,noutput
  real*8, intent(in) :: zinput(ivecd,ninput)      ! not used
  real*8, intent(inout) :: zoutput(ivecd,noutput) ! not used

  !-----------------------------
  !  local...

  real*8 :: zspl_data(1,1)
  integer :: ifail,idum

  integer :: ict_Zsrch(6) = (/ 0, 0, 1, 0, 0, 0 /)
  !-----------------------------
  !  expecting ivec=1 !!!

  if(ivec.ne.1) call errmsg_exit(' ?? eqi_geq_axis_zsrch: ivec=1 expected !! ')

  call r8xlookup(1,Z,size(zpkg,1),zpkg,2,jj,zparam,hz,hzi,idum)
  
  ! evaluate spline for dpsi/dR

  call r8fvbicub(ict_Zsrch,1,1, &
       zspl_data,ii,jj,rparam,zparam,hr,hri,hz,hzi, &
       psirz_spl,size(psirz_spl,2),size(psirz_spl,3))

  dPsidZ = zspl_data(1,1)

end subroutine eqi_geq_axis_zsrch
