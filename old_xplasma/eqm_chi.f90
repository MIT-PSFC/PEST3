subroutine eqm_chi(zchi,iauto,inum,ztol,id,ierr)

  use xplasma_obj_instance
  use eq_module

  !  set up chi grid -- must be strict ascending
  !  periodic grid, usually -pi to pi or 0 to 2pi.
  !  zchi(inum) = zchi(1) + 2pi.

  !  if iauto=1, the grid zchi(...) is generated here.

  !  ** THIS ROUTINE ASSUMES the poloidal angle coordinate chi is oriented
  !  so that chi=0 is near Rmax of a flux surface and dZ/dchi @chi=0 > 0, i.e.
  !  polar coordinate goes counter-clockwise around the flux surface for
  !  a plasma cross section drawn right of the machine centerline.  If this
  !  assumption is not desired, use eqm_chi_setccw instead!

  IMPLICIT NONE

  !  input:

  integer, intent(in) :: inum           ! #of grid bdy pts
  integer, intent(in) :: iauto          ! =1 to generate the chi grid -pi to pi
  REAL*8 zchi(inum)                     ! the grid (zone bdys, inum-1 zones).
  REAL*8, intent(in) :: ztol            ! even spacing tolerance

!  zchi is OUTPUT if iauto=1 is set: an evenly spaced grid spanning [-pi,pi]
!  output:

  integer id                        ! axis id code (returned)
  integer ierr                      ! =0: OK

  !--------------------------
  integer :: iccw
  !--------------------------

  iccw = 1

  call eqm_chi_setccw(iccw,zchi,iauto,inum,ztol,id,ierr)

end subroutine eqm_chi

subroutine eqm_chi_setccw(iccw,zchi,iauto,inum,ztol,id,ierr)

  use xplasma_obj_instance
  use eq_module

  !  set up chi grid -- must be strict ascending
  !  periodic grid, usually -pi to pi or 0 to 2pi.
  !  zchi(inum) = zchi(1) + 2pi.

  !  if iauto=1, the grid zchi(...) is generated here.

  !  ** ICCW=1 means ** the poloidal angle coordinate chi is oriented
  !  so that chi=0 is near Rmax of a flux surface and dZ/dchi @chi=0 > 0, i.e.
  !  polar coordinate goes counter-clockwise around the flux surface for
  !  a plasma cross section drawn right of the machine centerline.
  ! 
  !  ** ICCW<>1 the poloidal angle coordinate goes in reverse direction,
  !  clockwise around the flux surface.

  IMPLICIT NONE

  !  input:

  integer, intent(in) :: iccw           ! =1: CCW orientation of CHI

  integer, intent(in) :: inum           ! #of grid bdy pts
  integer, intent(in) :: iauto          ! =1 to generate the chi grid -pi to pi
  REAL*8 zchi(inum)                     ! the grid (zone bdys, inum-1 zones).
  REAL*8, intent(in) :: ztol            ! even spacing tolerance

!  zchi is OUTPUT if iauto=1 is set: an evenly spaced grid spanning [-pi,pi]
!  output:

  integer id                        ! axis id code (returned)
  integer ierr                      ! =0: OK

!-----------------

  real*8 zzchi(inum)
  integer i,inchi,iertmp
  logical :: lccw

!-------------------------------------------------------------

!  error checks:

  ierr=0
  id=0

  if(.not.eq_module_init) then
     write(lunerr,*) ' ?eqm_chi: initialization: call eqm_select first!'
     ierr=1
     return
  else 
     call xplasma_coord_info(s,xplasma_theta_coord,ierr, ngrids=inchi)
     if(ierr.ne.0) then
        call xplasma_error(s,ierr,lunerr)
     else if(inchi.gt.0) then
        write(lunerr,*) ' ?eqm_rho: rho already defined!'
        write(lunerr,*) '  re-initialization: call eqm_select first.'
        ierr=1
        return
     endif
  endif
  if(ierr.ne.0) return

  if(iauto.eq.1) then
     do i=1,inum
        zzchi(i)=-cpi + (i-1)*c2pi/(inum-1)
     enddo
     zchi=zzchi
  else

     if(abs((zchi(inum)-zchi(1))-c2pi).gt.ceps4) then
        call eq_errmsg(' ?eqm_chi: zchi(nchi)-zchi(1) .ne. 2*pi')
        write(lunerr,*) ' zchi(1)=',zchi(1),' zchi(nchi)=',zchi(inum)
        ierr=1
     else
        zzchi=zchi
        zzchi(inum)=zzchi(1)+c2pi
     endif
  endif

  lccw = (iccw.eq.1) 

  if(ierr.gt.0) return

  call eqi_evenx(zzchi,inum,ztol)

  call xoi_author_set(ierr)

  call xplasma_create_grid(s,'__CHI',xplasma_theta_coord,zzchi, &
       id,ierr, lccw)

  call xoi_author_clear(iertmp)

  return
end subroutine eqm_chi_setccw
