subroutine eqi_geq_axis_rsrch(ivec,iok,R,dPsidR,ivecd, &
     zinput,ninput,zoutput,noutput)

  ! see comput/zriddery.for -- root finder functional
  ! return dPsidR where dPsidZ=0

  use eqi_geq_axis_data
  implicit NONE

  !-----------------------------

  integer,intent(in) :: ivec
  logical,intent(in) :: iok(ivec)
  real*8, intent(in) :: R(ivec)
  real*8,intent(out) :: dPsidR(ivec)

  integer,intent(in) :: ivecd
  integer,intent(in) :: ninput,noutput
  real*8, intent(in) :: zinput(ivecd,ninput)
  real*8, intent(inout) :: zoutput(ivecd,noutput)

  !-----------------------------
  !  local...

  real*8 :: zmin(ivec),zmax(ivec),Z(ivec),zspl_data(1,1)
  integer :: ifail,idum

  integer :: ict_Rsrch(6) = (/ 0, 1, 0, 0, 0, 0 /)

  EXTERNAL eqi_geq_axis_Zsrch

  !-----------------------------
  !  expecting ivec=1 !!!

  if(ivec.ne.1) call errmsg_exit(' ?? eqi_geq_axis_rsrch: ivec=1 expected !! ')

  zmin(1) = zinput(1,1)
  zmax(1) = zinput(1,2)

  call r8xlookup(1,R,size(rpkg,1),rpkg,2,ii,rparam,hr,hri,idum)
  
  call zridderx(ivec,iok,zmin,zmax,eqi_eps,eqi_eta, &
       eqi_geq_axis_Zsrch,Z,ifail, &
       ivecd,zinput,ninput,zoutput,noutput)

  if(ifail.ne.0) then
#ifdef __DEBUG
     write(6,*) ' ?? eq8_geq_axis_rsrch: zridderx ifail = ',ifail
#endif

     dPsidR = 1000.0d0

  else
     ! search for dPsi/Z = 0 succeeded

     zoutput(1,1) = Z(1)

     ! evaluate spline for dpsi/dR

     call r8fvbicub(ict_Rsrch,1,1, &
          zspl_data,ii,jj,rparam,zparam,hr,hri,hz,hzi, &
          psirz_spl,size(psirz_spl,2),size(psirz_spl,3))

     dPsidR = zspl_data(1,1)

  endif

end subroutine eqi_geq_axis_rsrch
