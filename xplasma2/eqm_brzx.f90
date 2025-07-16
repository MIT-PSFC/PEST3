subroutine eqm_brzx(ivec,zR,zZ,zphi,init,BR,BZ,BPHI,zpsi,ipsi,ierr)

  !  **vectorized**
  !  compute B vector on (R,Z,phi) grid -- using psi extrapolation technique
  !  also return psi(pol) -- this routine assumes an axisymmetric geometry

  !  **NOTE** algorithm replaced: call eqm_brz_adhoc instead!!

  !  COMMENTS DESCRIBING OLD ALGORITHM:
  !  extrapolation uses del.B=0, delxB=j=0 until current has to be re-intro-
  !  duced to force smoothness.  (Reality is that there are coils out there
  !  which control the plasma position and shape.  Their presence is essential
  !  and delxB=0 extrapolation is exponentially unstable without them.  But,
  !  as their positions are unknown, j is re-introduced in such a way as to
  !  smooth out the field extrapolation and act in their place.  Even so, 
  !  the field extrapolation can still be rough.  Therefore, the alternative
  !    "eqm_brz_adhoc"
  !  is also provided -- even less physical, but numerically better behaved).

  !  see subroutine eqm_brz_adhoc, below...

  !  axisymmetry is assumed... phi argument is ignored, but left in so
  !  that the routine can be used as an argument of the eqm_brz subroutine
  !--------------------------

  implicit NONE

  integer ivec                      ! vector lengths
  real*8 zR(ivec),zZ(ivec),zphi(ivec) ! locations to evaluate B vector

  integer init                      ! init-cleanup-compute control flag

  !  init=1  -- initialize
  !  init=2  -- cleanup
  !  init=0  -- compute

  !  ...a vector of B vectors...
  real*8 BR(ivec),BZ(ivec),BPHI(ivec) ! field vectors components (returned)

  real*8 zpsi(ivec)                 ! psi(poloidal) returned
  integer ipsi                      ! =1 returned, psi(poloidal) flag

  integer ierr                      ! completion code, 0=OK

  !--------------------------------

  call eqm_brz_adhoc(ivec,zR,zZ,zphi,init,BR,BZ,BPHI,zpsi,ipsi,ierr)

end subroutine eqm_brzx

!----------------------------------------------------------------------
subroutine eqm_brz_adhoc(ivec,zR,zZ,zphi,init,BR,BZ,BPHI,zpsi,ipsi,ierr)

  !  **vectorized**
  !  compute B vector on (R,Z,phi) grid -- using ad hoc extrapolation technique
  !  also return psi(pol) -- this routine assumes an axisymmetric geometry

  !  axisymmetry is assumed... phi argument is ignored, but left in so
  !  that the routine can be used as an argument of the eqm_brz subroutine

  !  the flux surface grid is extrapolated using a TRANSP-style "bloat"
  !  algorithm-- creating a fiction of closed "flux surfaces" beyond the
  !  plasma boundary.  These aren't real flux surfaces but serve as a
  !  convenient numerical grid.

  !  for rho > rhobdy, Bpol(rho,theta) = Bpol(rhobdy,theta)*L(rhobdy)/L(rho)
  !  is assumed (this determines the magnitude; the direction is tangent to
  !  the external (pseudo) flux surface with sign consistent with the interior
  !  plasma current.  L(rhobdy) = poloidal path length around plasma boundary;
  !  L(rho) = poloidal path length around rho surface.

  !  also, for rho > rhobdy, psi(rho) is integrated from dR*2pi*R*Bz, selected
  !  at R=Rmax of each extrapolated pseudo-flux surface

  implicit NONE

  integer ivec                      ! vector lengths
  real*8 zR(ivec),zZ(ivec),zphi(ivec) ! locations to evaluate B vector

  integer init                      ! init-cleanup-compute control flag

  !  init=1  -- initialize
  !  init=2  -- cleanup
  !  init=0  -- compute

  !  ...a vector of B vectors...
  real*8 BR(ivec),BZ(ivec),BPHI(ivec) ! field vectors components (returned)

  real*8 zpsi(ivec)                 ! psi(poloidal) returned
  integer ipsi                      ! =1 returned, psi(poloidal) flag

  integer ierr                      ! completion code, 0=OK

  !--------------------------------

  call eqm_brz_exec(ivec,zR,zZ,zphi,init,BR,BZ,BPHI,zpsi,ipsi,ierr)

end subroutine eqm_brz_adhoc

!--------------------------------
subroutine eqm_brz_exec(ivec,zR,zZ,zphi,init,BR,BZ,BPHI,zpsi,ipsi,ierr)

  !  **vectorized**
  !  compute B vector on (R,Z,phi) grid -- using psi extrapolation technique
  !  also return psi(pol) -- this routine assumes an axisymmetric geometry

  !  axisymmetry is assumed... phi argument is ignored, but left in so
  !  that the routine can be used as an argument of the eqm_brz subroutine

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE

  integer ivec                      ! vector lengths
  real*8 zR(ivec),zZ(ivec),zphi(ivec) ! locations to evaluate B vector

  integer init                      ! init-cleanup-compute control flag

  !  init=1  -- initialize
  !  init=2  -- cleanup
  !  init=0  -- compute

  !  ...a vector of B vectors...
  real*8 BR(ivec),BZ(ivec),BPHI(ivec) ! field vectors components (returned)

  real*8 zpsi(ivec)                 ! psi(poloidal) returned
  integer ipsi                      ! =1 returned, psi(poloidal) flag

  integer ierr                      ! completion code, 0=OK

  !--------------------------------

  real*8 zrho,zindx,zdi
  integer :: indxR,indxZ,itype,i,iadr
  real*8 :: zrmin,zrmax,zzmin,zzmax

  integer :: id_Rg,id_Zg,id_g,id_psi,id_bmod
  integer :: bphi_ccw,jphi_ccw,iertmp,id_bb,nR,nZ
  logical :: axisymm,sol

  real*8, pointer, dimension(:) :: eqbuf
  integer, pointer, dimension(:) :: idata
  
  character*32 :: bbxname

  !--------------------------------

  ipsi=1
  bbxname='__BBXTRAP_TMP'

  call xplasma_global_info(sp,ierr, axisymm=axisymm, scrapeoff=sol, &
       bphi_ccw=bphi_ccw, jphi_ccw=jphi_ccw)

  if(ierr.ne.0) return

  call xplasma_find_item(sp,bbxname,id_bb,iertmp,nf_noerr=.TRUE.)
  if(ierr.ne.0) return

  !  cleanup call-- remove temporary black box

  if(init.eq.2) then
     call xplasma_author_set(sp,xplasma_xmhd,iertmp)
     if(id_bb.gt.0) then
        call xplasma_remove_item(sp,id_bb,iertmp)
     endif
     call xplasma_author_clear(sp,xplasma_xmhd,iertmp)

     BR=0
     BZ=0
     Bphi=0
     zpsi=0

     return
  endif

  !  initialization; error checks done once...

  if(init.eq.1) then

     if(.not.axisymm) then
        call xplasma_errmsg_append(sp,'  ...error in eqm_brz_exec')
        ierr=107
        return
     endif

     if(.not.sol) then
        call xplasma_errmsg_append(sp,'  ...error in eqm_brz_exec')
        ierr=109
        return
     endif

     !  B field SIGN factors must be set
     if((jphi_ccw.eq.0).or.(bphi_ccw.eq.0)) then
        call xplasma_errmsg_append(sp, &
             ' ?eqm_brzx/eqm_brz_adhoc:  set B sign factors first!')
        ierr=9999
        return
     endif

     call xplasma_find_item(sp,'BMOD',id_bmod,ierr,nf_noerr=.TRUE.)
     if(id_bmod.eq.0) then
        call xplasma_errmsg_append(sp, &
             ' ?eqm_brzx/eqm_brz_adhoc:  "BMOD" undefined.')
        call xplasma_errmsg_append(sp, &
             '  it is likely that equilibrium is not defined.')
        ierr=9999
        return
     endif

  !  create B(R,Z) & psi temporary space

     nR=0
     nZ=0

     call xplasma_find_item(sp,'__Rgrid',id_Rg,ierr)
     if(ierr.ne.0) return

     call xplasma_find_item(sp,'__Zgrid',id_Zg,ierr)
     if(ierr.ne.0) return

     call xplasma_grid_size(sp,id_Rg,nR,ierr)
     if(ierr.ne.0) return

     call xplasma_grid_size(sp,id_Zg,nZ,ierr)
     if(ierr.ne.0) return

     call xplasma_author_set(sp,xplasma_xmhd,iertmp)

     itype=9999
     call xplasma_create_blackBox(sp,bbxname,itype,id_bb,ierr, &
          iasize=5, r8asize= 4*nR*nZ+nR+nZ+2, ia_ptr=idata, r8a_ptr=eqbuf)

     call xplasma_author_clear(sp,xplasma_xmhd,iertmp)

     if(ierr.ne.0) then
        nullify(eqbuf,idata)
        return
     endif

     idata(1)=nR
     idata(2)=nZ
     idata(3)=4*nR*nZ+1
     idata(4)=idata(3)+nR
     idata(5)=idata(4)+nZ

     call xplasma_grid(sp,id_Rg,eqbuf(idata(3):idata(3)+nR-1),ierr)
     if(ierr.ne.0) then
        nullify(eqbuf,idata)
        return
     endif

     call xplasma_grid(sp,id_Zg,eqbuf(idata(4):idata(4)+nZ-1),ierr)
     if(ierr.ne.0) then
        nullify(eqbuf,idata)
        return
     endif

     !  g(edge)

     call xplasma_common_ids(sp,ierr, id_g=id_g, id_psi=id_psi)

     zrho=1.0d0
     call xplasma_eval_prof(sp,id_g,zrho,eqbuf(idata(5)),ierr)
     if(ierr.ne.0) then
        nullify(eqbuf,idata)
        return
     endif

     zrho=1.0d0
     call xplasma_eval_prof(sp,id_psi,zrho,eqbuf(idata(5)+1),ierr)
     if(ierr.ne.0) then
        nullify(eqbuf,idata)
        return
     endif

     if(ierr.ne.0) return

     call xplasma_blackBox_retrieve(sp,id_bb,ierr, &
          ia_ptr=idata, r8a_ptr=eqbuf)

     !  fill in the fields  BR  BZ  Bphi  Psi -- outside plasma only...

     call eqi_xpsi(nR,nZ,eqbuf(idata(3)),eqbuf(idata(4)), &
          eqbuf(idata(5)),eqbuf(idata(5)+1),bphi_ccw, jphi_ccw, &
          eqbuf(1),eqbuf(1+nR*nZ),eqbuf(1+2*nR*nZ),eqbuf(1+3*nR*nZ),ierr)

     BR=0
     BZ=0
     Bphi=0
     zpsi=0

     nullify(idata)
     nullify(eqbuf)

     return
  endif

  !-----------------------------
  !  normal evaluation, using precomputed data (init=0)

  if(id_bb.eq.0) then
     call xplasma_errmsg_append(sp,' ?eqm_brzx-- never initialized.')
     ierr=9999
     return
  endif

  call xplasma_blackBox_retrieve(sp,id_bb,ierr, &
       ia_ptr=idata, r8a_ptr=eqbuf)

  nR=idata(1)
  nZ=idata(2)

  zrmin=eqbuf(idata(3)); zrmax=eqbuf(idata(3)+nR-1)
  zzmin=eqbuf(idata(4)); zzmax=eqbuf(idata(4)+nZ-1)

  do i=1,ivec
     zindx=1.0d0+(nR-1)*(zR(i)-zRmin)/(zRmax-zRmin)
     indxR=zindx+0.01d0
     zdi=zindx-indxR
     zindx=1.0d0+(nZ-1)*(zZ(i)-zZmin)/(zZmax-zZmin)
     indxZ=zindx+0.01d0
     zdi=max(zdi,(zindx-indxZ))
     if(zdi.gt.1.0d-5) then
        call xplasma_errmsg_append(sp, &
             '  ?eqm_brzx: expected (R,Z) to be at a grid point.')
        ierr=9999
        exit
     endif

     iadr=indxR+nR*(indxZ-1)
     BR(i)=eqbuf(iadr)
     BZ(i)=eqbuf(iadr+nR*nZ)
     BPHI(i)=eqbuf(iadr+2*nR*nZ)
     zpsi(i)=eqbuf(iadr+3*nR*nZ)
  enddo

  nullify(eqbuf)
  nullify(idata)

end subroutine eqm_brz_exec
