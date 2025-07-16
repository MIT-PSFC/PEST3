subroutine eqi_fbmod(ivec,zrho,zchii,ict,ivecd,zans,ier)

  use xplasma_definitions
  use eqi_rzbox_module

  !  evaluate 2d spline (and derivatives) of mod(B) vs. rho,chi

  IMPLICIT NONE

  !  input:

  integer ivec
  REAL*8 zrho(ivec)                 ! argument -- where to evaluate
  REAL*8 zchii(ivec)                ! argument -- where to evaluate
  integer ict                       ! 1: return B; 2: return B & dB/dtheta

  !  chi is periodic and will be adjusted, but, the user should avoid
  !  passing chi arguments that are out of the normal range by more than
  !  a few periods.

  !  output:

  integer ivecd
  REAL*8 zans(ivecd,ict)            ! result of evaluation(s)
  integer ier                       ! exit code, 0= normal

  !----------------------------------
  integer :: ilist(ict)
  integer :: ideriv(ict)
  integer :: id_Bmod
  !----------------------------------

  call xplasma_common_ids(sp,ier,id_Bmod=id_Bmod)
  if(ier.ne.0) return

  if(id_Bmod.eq.0) then
     call xplasma_errmsg_append(sp, &
          ' ?eqi_fBod: no Bmod(rho,theta) spline found.')
     ier=106
     return
  endif

  ilist=id_Bmod
  ideriv(1)=0
  if(ict.eq.2) ideriv(2)=1

  call xplasma_eval_prof(sp,ilist, &
       xplasma_theta_coord,zchii, xplasma_rho_coord,zrho, &
       zans(1:ivec,1:ict),ier, ideriv1s=ideriv)

end subroutine eqi_fbmod
