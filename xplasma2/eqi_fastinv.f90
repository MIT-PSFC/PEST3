subroutine eqi_fastinv(ivec,zR,zZ,zphi,zrhomin,zrho_out,zchi_out,iregion,ierr)

  !  eval scrape-off-region (R,Z,phi)->(rho,chi,phi) map.
  !  rho,chi form a nested nonsingular coordinate system, but, no
  !  particular relation to field or flux surfaces.

  !  iregion output:
  !  The SOL is defined by a RxZ grid rectangle.  If an input point
  !  (zR(i),zZ(i)) is not inside this rectangle, iregion(i)=3 is returned.
  !  Otherwise iregion(i)=1 or iregion(i)=2 according as zrho_out(i) is or
  !  is not greater than one.

  !  *** if the map is undefined, create it ***

  use xplasma_definitions
  use eqi_rzbox_module
  implicit NONE

  integer ivec                      ! vector dimension
  real*8 zR(ivec),zZ(ivec),zphi(ivec) ! input positions from which to seek
  real*8 zrhomin                    ! min. allowed rho output

  real*8 zrho_out(ivec)             ! output rho of sol pt
  real*8 zchi_out(ivec)             ! output chi of sol pt
  integer :: iregion(ivec)          ! output region code as described

  integer ierr                      ! completion code, 0=OK

  !-----------------------------------------
  real*8 zRtarg(ivec),zZtarg(ivec)
  real*8 zrmid,zzmid,zcorr(ivec)

  integer, dimension(:), pointer :: ii,jj
  real*8, dimension(:), pointer :: zRoff,zZoff

  integer iadr00,iadr01,iadr10,iadr11,i,nR,nZ
  real*8 zrmin,zrmax,zzmin,zzmax
  logical :: sol

  real*8 zchibr00,zchibr01,zchibr10,zchibr11,zchimin,zchitest

  integer :: id_map,id_Rgrid,id_Zgrid,lrhomap,lchimap
  integer, dimension(:), pointer :: idata
  real*8, dimension(:), pointer :: eqbuf

  type (xpeval) :: R_intrp, Z_intrp

  real*8, parameter :: czero=0.0d0
  real*8, parameter :: chalf=0.5d0
  real*8, parameter :: cone =1.0d0
  real*8, parameter :: ceps =1.0d-15
  real*8, parameter :: c2pi =6.2831853071795862D+00
  real*8, parameter :: cpi  =3.1415926535897931D+00
  !-----------------------------------------

  ierr=0

  call xplasma_RZminmax_extended(sp, zrmin,zrmax, zzmin,zzmax, ierr, sol=sol)
  if(ierr.ne.0) return

  if(.not.sol) then
     ierr=109
     call xplasma_errmsg_append(sp, &
         ' ?? eqi_fastinv.f90 was called, but no scrapeoff region is defined.')
     return
  endif

  call xplasma_find_item(sp,'__FASTMAP',id_map,ierr, nf_noerr=.TRUE.)
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'__Rgrid',id_Rgrid,ierr)
  if(ierr.ne.0) return

  call xplasma_grid_size(sp,id_Rgrid,nR,ierr)
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'__Zgrid',id_Zgrid,ierr)
  if(ierr.ne.0) return

  call xplasma_grid_size(sp,id_Zgrid,nZ,ierr)
  if(ierr.ne.0) return

  if(id_map.eq.0) call eqi_fastinv_gen(id_Rgrid,id_Zgrid,nR,nZ,id_map,ierr)
  if(ierr.ne.0) return

  call xplasma_blackbox_retrieve(sp, id_map, ierr, &
       ia_ptr = idata, r8a_ptr = eqbuf)

  lrhomap = idata(1)
  lchimap = idata(2)

  !  SOL (rho,chi) -- use bilinear approximation

  zrmid=chalf*(zrmin+zrmax)

  zzmid=chalf*(zzmin+zzmax)

  !  rough correction for points beyond the (R,Z) grid

  do i=1,ivec
     zRtarg(i)=max(zrmin,min(zrmax,zR(i)))
     zZtarg(i)=max(zzmin,min(zzmax,zZ(i)))
     if((zR(i).ne.zRtarg(i)).or.(zZ(i).ne.zZtarg(i))) then
        zcorr(i)=sqrt((zR(i)-zrmid)**2+(zZ(i)-zzmid)**2)/ &
             sqrt((zrtarg(i)-zrmid)**2+(zztarg(i)-zzmid)**2)
        zcorr(i)=max(cone+ceps,zcorr(i))
        iregion(i)=3
     else
        zcorr(i)=cone
        iregion(i)=2
     endif
  enddo

  !  lookup...

  call xplasma_x_lookup(sp, id_Rgrid, zRtarg, R_intrp, ierr)
  if(ierr.ne.0) then
     call xpeval_free(R_intrp)
     return
  endif

  call xplasma_x_lookup(sp, id_Zgrid, zZtarg, Z_intrp, ierr)
  if(ierr.ne.0) then
     call xpeval_free(R_intrp)
     call xpeval_free(Z_intrp)
     return
  endif

  call xpeval_info(R_intrp, ierr, ix=ii, dxn=zRoff)
  call xpeval_info(Z_intrp, ierr, ix=jj, dxn=zZoff)

  !  evaluation...

  do i=1,ivec

     !  rho

     iadr00=lrhomap +(jj(i)-1)*nR + ii(i)-1
     iadr01=iadr00+1
     iadr10=iadr00+nR
     iadr11=iadr01+nR

     zrho_out(i)= &
          (cone-zZoff(i))* &
          ((cone-zRoff(i))*eqbuf(iadr00)+ &
                           zRoff(i)*eqbuf(iadr01)) + &
          zZoff(i)* &
          ((cone-zRoff(i))*eqbuf(iadr10)+ &
                           zRoff(i)*eqbuf(iadr11))

     zrho_out(i)=max(zrhomin,zrho_out(i)*zcorr(i))
     if(zrho_out(i).le.cone) iregion(i)=1

     !  chi -- watch for branch cut

     iadr00=lchimap +(jj(i)-1)*nR + ii(i)-1
     iadr01=iadr00+1
     iadr10=iadr00+nR
     iadr11=iadr01+nR

     zchimin=min(eqbuf(iadr00),eqbuf(iadr01),eqbuf(iadr10),eqbuf(iadr11))
     zchitest=zchimin+cpi

     if(eqbuf(iadr00).gt.zchitest) then
        zchibr00=-c2pi
     else
        zchibr00=czero
     endif

     if(eqbuf(iadr01).gt.zchitest) then
        zchibr01=-c2pi
     else
        zchibr01=czero
     endif

     if(eqbuf(iadr10).gt.zchitest) then
        zchibr10=-c2pi
     else
        zchibr10=czero
     endif

     if(eqbuf(iadr11).gt.zchitest) then
        zchibr11=-c2pi
     else
        zchibr11=czero
     endif

     zchi_out(i)= &
          (cone-zZoff(i))* &
          ((cone-zRoff(i))*(eqbuf(iadr00)+zchibr00)+ &
                           zRoff(i)*(eqbuf(iadr01)+zchibr01)) + &
          zZoff(i)* &
          ((cone-zRoff(i))*(eqbuf(iadr10)+zchibr10)+ &
                           zRoff(i)*(eqbuf(iadr11)+zchibr11))

  enddo

  nullify(ii,jj,zRoff,zZoff,idata,eqbuf)
  call xpeval_free(R_intrp)
  call xpeval_free(Z_intrp)

end subroutine eqi_fastinv

subroutine eqi_fastinv_gen(id_Rgrid,id_Zgrid,nR,nZ,id_map,ierr)

  ! generate fast inverse map on RxZ extended grids

  use xplasma_definitions
  use eqi_rzbox_module
  implicit NONE

  !-----------

  integer, intent(in) :: id_Rgrid,id_Zgrid  ! grid IDs
  integer, intent(in) :: nR,nZ              ! grid sizes

  integer, intent(out) :: id_map            ! "blackbox" id
  integer, intent(out) :: ierr              ! status code, 0=OK

  !-----------
  real*8 :: Rgrid(nR),Zgrid(nZ),Ztmp(nR)
  real*8, dimension(:), allocatable :: Rtarg,Ztarg,phidum
  real*8, dimension(:), pointer :: r8data
  integer, dimension(:), allocatable :: iregion
  real*8 :: reltol,rtol
  integer :: idata(2),ii,jj,iR,iZ,iertmp,itype
  real*8 :: zrho_bdy
  !-----------

  call xplasma_grid(sp, id_Rgrid, Rgrid, ierr)
  if(ierr.ne.0) return

  call xplasma_grid(sp, id_Zgrid, Zgrid, ierr)
  if(ierr.ne.0) return

  call xplasma_global_info(sp, ierr, bdytol=reltol)
  if(ierr.ne.0) return

  if(reltol.lt.1.0d-13) reltol=1.0d-13

  allocate(Rtarg(nR*nZ),Ztarg(nR*nZ))
  allocate(phidum(nR*nZ),iregion(nR*nZ))

  call eqi_extrap_rho_bdy(zrho_bdy)

  call xplasma_author_set(sp, xplasma_xmhd,iertmp)

  do 
     rtol=Rgrid(nR)*reltol

     ii=0
     jj=nR*nZ

     idata(1)=ii+1
     idata(2)=jj+1

     itype=77 ! arbitrary type code
     call xplasma_create_blackBox(sp,'__FASTMAP',itype,id_map,ierr, &
          iarray=idata, r8asize=2*nR*nZ, label='fast inverse map data', &
          r8a_ptr=r8data)
     if(ierr.ne.0) exit

     do iZ=1,nZ
        Ztmp=Zgrid(iZ)
        call xplasma_rhoth_inv(sp,Rgrid,Ztmp, &
             r8data(ii+1:ii+nR),r8data(jj+1:jj+nR),iregion(1:nR),ierr)
        if(ierr.ne.0) exit
        Rtarg(ii+1:ii+nR)=Rgrid
        Ztarg(ii+1:ii+nR)=Ztmp
        Phidum(ii+1:ii+nR)=0
        ii=ii+nR
        jj=jj+nR
     enddo
     if(ierr.ne.0) exit

     ! refine this...
     call eqi_xinv2d(nR*nZ,r8data(1:nR*nZ),r8data(nR*nZ+1:2*nR*nZ),Phidum, &
          Rtarg,Ztarg,iregion,reltol,rtol,zrho_bdy,ierr)
     if(ierr.ne.0) exit

     exit
  enddo

  call xplasma_author_clear(sp, xplasma_xmhd,iertmp)
  nullify(r8data)

end subroutine eqi_fastinv_gen
