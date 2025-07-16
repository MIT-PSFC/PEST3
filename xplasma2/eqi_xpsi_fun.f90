subroutine eqi_xpsi_fun(ivec,iok,x,ansr,ivecd,zinput,ninput,zoutput,noutput)

  ! subroutine passed to "zridderx" root finder -- 
  !   looking for rho to satisfy Psi(rho) - Psi_target = 0

  use xplasma_definitions
  use eqi_rzbox_module
  implicit NONE

  integer ivec
  logical iok(ivec)  ! if iok(i) TRUE then SKIP element i of vector
  real*8 x(ivec)     ! input vector {x(i)}
  real*8 ansr(ivec)  ! output vector {f(x(i))}
  integer ivecd      ! auxilliary information vector size (.ge.ivec)
  integer ninput,noutput
  real*8 zinput(ivecd,ninput),zoutput(ivecd,noutput)

  !-----------------
  integer, parameter :: ilun_dbg=6
  integer :: id_psi,i,iertmp,itry
  real*8 :: xtry(ivec),pstry(ivec),pstarg(ivec)
  !-----------------

  call xplasma_common_ids(sp,iertmp, id_psi=id_psi) ! assume this is OK

  itry = 0
  do i=1,ivec
     if(.not.iok(i)) then
        itry=itry+1
        xtry(itry)=x(i)
        pstarg(itry)=zinput(i,1)
     endif
  enddo

  call xplasma_eval_prof(sp,id_psi,xtry,pstry,iertmp, force_bounds=.TRUE.)
  if(iertmp.ne.0) then
     write(ilun_dbg,*) ' ? xplasma2/eqi_xpsi_fun.f90: unexpected error:'
     call xplasma_error(sp,iertmp,6)
  endif

  itry = 0
  do i=1,ivec
     if(.not.iok(i)) then
        itry=itry+1
        ansr(i)=pstry(itry)-pstarg(itry)
     endif
  enddo

end subroutine eqi_xpsi_fun
