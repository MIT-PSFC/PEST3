subroutine eqi_bdfast(ivec,zR,zZ,iwant,want_th,zdist,zth_out,ierr)

  !  evaluate distance map: (R,Z)->(distance,nearest-theta)
  !    distance from rho=1 surface, i.e. the plasma boundary

  !  SOL = scrape off layer region enclosing a plasma.

  !  The SOL is defined by a RxZ grid rectangle.  If an input point
  !  (zR(i),zZ(i)) is not inside this rectangle, a rough correction is
  !  provided to the distance and nearest-theta

  !  *** if the map is undefined, create it ***

  use xplasma_definitions
  use eqi_rzbox_module
  implicit NONE

  integer ivec                      ! vector dimension
  real*8 zR(ivec),zZ(ivec)          ! input positions for which 

  integer, intent(in) :: iwant      ! =2: want entire fast map defined
  ! =1: only need points with d>0, i.e. outside the plasma.

  logical :: want_th                ! .TRUE. if zth_out(...) output is wanted.

  real*8 zdist(ivec)                ! output distance to nearest pt on bdy
  real*8 zth_out(ivec)              ! output theta of nearest pt on bdy

  integer ierr                      ! completion code, 0=OK

  !-----------------------------------------
  real*8 zRtarg(ivec),zZtarg(ivec)
  real*8 zcorr(ivec)

  integer, dimension(:), pointer :: ii,jj
  real*8, dimension(:), pointer :: zRoff,zZoff

  integer iadr00,iadr01,iadr10,iadr11,i,nR,nZ
  real*8 zrmin,zrmax,zzmin,zzmax
  logical :: sol

  real*8 zthbr00,zthbr01,zthbr10,zthbr11,zthmin,zthtest

  integer :: id_map,id_Rgrid,id_Zgrid,ldmap,lthmap,istat
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
         ' ?? eqi_bdfast.f90 was called, but no scrapeoff region is defined.')
     return
  endif

  call xplasma_find_item(sp,'__DISTMAP',id_map,ierr, nf_noerr=.TRUE.)
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'__Rgrid',id_Rgrid,ierr)
  if(ierr.ne.0) return

  call xplasma_grid_size(sp,id_Rgrid,nR,ierr)
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'__Zgrid',id_Zgrid,ierr)
  if(ierr.ne.0) return

  call xplasma_grid_size(sp,id_Zgrid,nZ,ierr)
  if(ierr.ne.0) return

  ! create full map, if needed...
  if(id_map.eq.0) call eqi_bdfast_gen(id_Rgrid,id_Zgrid,nR,nZ,id_map,iwant, &
       ierr)
  if(ierr.ne.0) return

  call xplasma_blackbox_retrieve(sp, id_map, ierr, &
       ia_ptr = idata, r8a_ptr = eqbuf)

  ldmap = idata(1)
  lthmap = idata(2)
  istat = idata(3)

  ! complete map, if needed...
  if(istat.lt.iwant) then
     call eqi_bdfast_gen(id_Rgrid,id_Zgrid,nR,nZ,id_map,iwant,ierr)
  endif

  !  SOL (rho,th) -- use bilinear approximation
  !  rough correction for points beyond the (R,Z) grid

  do i=1,ivec
     zRtarg(i)=max(zrmin,min(zrmax,zR(i)))
     zZtarg(i)=max(zzmin,min(zzmax,zZ(i)))
     if((zR(i).ne.zRtarg(i)).or.(zZ(i).ne.zZtarg(i))) then
        zcorr(i)=sqrt((zR(i)-zRtarg(i))**2+(zZ(i)-zZtarg(i))**2)
        zcorr(i)=max(ceps,zcorr(i))
     else
        zcorr(i)=czero
     endif
  enddo

  !  lookup...

  call xplasma_x_lookup(sp, id_Rgrid, zRtarg, R_intrp, ierr)
  if(ierr.ne.0) then
     nullify(eqbuf,idata)
     call xpeval_free(R_intrp)
     return
  endif

  call xplasma_x_lookup(sp, id_Zgrid, zZtarg, Z_intrp, ierr)
  if(ierr.ne.0) then
     nullify(eqbuf,idata)
     call xpeval_free(R_intrp)
     call xpeval_free(Z_intrp)
     return
  endif 

  call xpeval_info(R_intrp, ierr, ix=ii, dxn=zRoff)
  call xpeval_info(Z_intrp, ierr, ix=jj, dxn=zZoff)

  !  evaluation...

  do i=1,ivec

     !  rho

     iadr00=ldmap +(jj(i)-1)*nR + ii(i)-1
     iadr01=iadr00+1
     iadr10=iadr00+nR
     iadr11=iadr01+nR

     zdist(i)= &
          (cone-zZoff(i))* &
          ((cone-zRoff(i))*eqbuf(iadr00)+ &
                           zRoff(i)*eqbuf(iadr01)) + &
          zZoff(i)* &
          ((cone-zRoff(i))*eqbuf(iadr10)+ &
                           zRoff(i)*eqbuf(iadr11))

     zdist(i)=zdist(i)+zcorr(i)

     !  th -- watch for branch cut  ** if wanted **

     if(.not.want_th) cycle

     iadr00=lthmap +(jj(i)-1)*nR + ii(i)-1
     iadr01=iadr00+1
     iadr10=iadr00+nR
     iadr11=iadr01+nR

     zthmin=min(eqbuf(iadr00),eqbuf(iadr01),eqbuf(iadr10),eqbuf(iadr11))
     zthtest=zthmin+cpi

     if(eqbuf(iadr00).gt.zthtest) then
        zthbr00=-c2pi
     else
        zthbr00=czero
     endif

     if(eqbuf(iadr01).gt.zthtest) then
        zthbr01=-c2pi
     else
        zthbr01=czero
     endif

     if(eqbuf(iadr10).gt.zthtest) then
        zthbr10=-c2pi
     else
        zthbr10=czero
     endif

     if(eqbuf(iadr11).gt.zthtest) then
        zthbr11=-c2pi
     else
        zthbr11=czero
     endif

     zth_out(i)= &
          (cone-zZoff(i))* &
          ((cone-zRoff(i))*(eqbuf(iadr00)+zthbr00)+ &
                           zRoff(i)*(eqbuf(iadr01)+zthbr01)) + &
          zZoff(i)* &
          ((cone-zRoff(i))*(eqbuf(iadr10)+zthbr10)+ &
                           zRoff(i)*(eqbuf(iadr11)+zthbr11))

  enddo

  nullify(ii,jj,zRoff,zZoff,idata,eqbuf)
  call xpeval_free(R_intrp)
  call xpeval_free(Z_intrp)

end subroutine eqi_bdfast

subroutine eqi_bdfast_gen(id_Rgrid,id_Zgrid,nR,nZ,id_map,iwant,ierr)

  ! generate or complete fast distance map on RxZ extended grids

  !   id_map=0 on input means map does not exist at all;
  !            (it is created and id_map receives the id as output).
  !   id_map=1 on input means map exists but may not be complete (check).

  !   iwant=2 on input means entire map is needed.
  !   iwant=1 on input means only points beyond the plasma boundary are
  !   needed.

  use xplasma_definitions
  use eqi_rzbox_module
  implicit NONE

  !-----------

  integer, intent(in) :: id_Rgrid,id_Zgrid  ! grid IDs
  integer, intent(in) :: nR,nZ              ! grid sizes

  integer, intent(inout) :: id_map          ! "blackbox" id
  integer, intent(in) :: iwant              ! map section desired: 2 for all

  integer, intent(out) :: ierr              ! status code, 0=OK

  !-----------
  real*8 :: Rgrid(nR),Zgrid(nZ),Ztmp(nR),zrho_bdy,zrho0
  real*8, dimension(:), allocatable :: Rvec,Zvec,phidum,phitmp
  real*8, dimension(:), allocatable :: dmap,thmap,rhotmp,chtmp  ! 1 row only
  real*8, dimension(:), pointer :: r8data
  logical, dimension(:), allocatable :: iuse
  integer, dimension(:), pointer :: idata
  integer :: ii,jj,iR,iZ,iRegion(nR),iertmp,itype
  integer :: id_dmap,id_thmap,ihave,jvec
  logical :: iskipi,iskipo,iexact
  !-----------

  iskipi=.FALSE.
  iskipo=.FALSE.
  iexact=.FALSE.
  zrho_bdy = 1.0d0
  zrho0 = 0.0d0

  call xplasma_grid(sp, id_Rgrid, Rgrid, ierr)
  if(ierr.ne.0) return

  call xplasma_grid(sp, id_Zgrid, Zgrid, ierr)
  if(ierr.ne.0) return

  if(id_map.eq.0) then
     allocate(idata(3))
     idata(1)=1
     idata(2)=nR*nZ+1
     idata(3)=0

     itype=77 ! arbitrary type code
     call xplasma_author_set(sp, xplasma_xmhd,iertmp)

     call xplasma_create_blackBox(sp,'__DISTMAP',itype,id_map,ierr, &
          iarray=idata,r8asize=2*nR*nZ, &
          label='fast distance map data')

     call xplasma_author_clear(sp, xplasma_xmhd,iertmp)

     deallocate(idata)
  endif
  if(ierr.ne.0) return

  call xplasma_blackBox_retrieve(sp, id_map, ierr, &
       ia_ptr=idata, r8a_ptr=r8data)
  if(ierr.ne.0) then
     nullify(r8data,idata)
     return
  endif

  ihave = idata(3)

  !  ihave = 0 means: have no data;
  !  ihave = 1 means: have data for pts outside plasma only;
  !  ihave = 2 means: have all data.

  if(ihave.ge.iwant) return

  !  ihave < iwant:
  !    ihave=0, iwant=2: compute entire map
  !    ihave=0, iwant=1: compute map outside plasma only
  !    ihave=1, iwant=2: compute map inside plasma only

  allocate(dmap(nR*nZ),thmap(nR*nZ),rhotmp(nR),chtmp(nR))
  allocate(Phidum(nR*nZ),Phitmp(nR*nZ),Rvec(nR*nZ),Zvec(nR*nZ),iuse(nR*nZ))

  do 

     Phidum = 0

     if((ihave.eq.0).and.(iwant.eq.2)) then
        !  compute entire map...
        ii=0
        do iZ=1,nZ
           Rvec(ii+1:ii+nR)=Rgrid
           Zvec(ii+1:ii+nR)=Zgrid(iZ)
           ii=ii+nR
        enddo
        call eqi_bdfind(nR*nZ,Rvec,Zvec,Phidum,zrho_bdy,iskipi,iskipo,iexact,&
             r8data(nR*nZ+1:2*nR*nZ),Phitmp,r8data(1:nR*nZ),ierr)
        if(ierr.ne.0) exit

     else
        !  compute necessary subset of map
        jvec=0
        ii=0
        do iZ=1,nZ
           Ztmp=Zgrid(iZ)
           call eqi_fastinv(nR,Rgrid,Ztmp,Phidum,zrho0,rhotmp,chtmp,iregion, &
                ierr)
           if(ierr.ne.0) exit

           do iR=1,nR
              ii=ii+1
              if(rhotmp(iR).gt.1.0d0) then
                 if(ihave.eq.0) then
                    !  compute outside points only
                    jvec=jvec+1
                    Rvec(jvec)=Rgrid(iR)
                    Zvec(jvec)=Zgrid(iZ)
                    iuse(ii)=.TRUE.
                 else
                    iuse(ii)=.FALSE.
                 endif
              else
                 if(ihave.eq.1) then
                    !  compute inside points only
                    jvec=jvec+1
                    Rvec(jvec)=Rgrid(iR)
                    Zvec(jvec)=Zgrid(iZ)
                    iuse(ii)=.TRUE.
                 else
                    iuse(ii)=.FALSE.
                 endif
              endif
           enddo
        enddo

        call eqi_bdfind(jvec,Rvec,Zvec,Phidum,zrho_bdy, &
             iskipi,iskipo,iexact,thmap,Phitmp,dmap,ierr)
        if(ierr.ne.0) exit

        jvec=0
        ii=0
        jj=nR*nZ

        do iZ=1,nZ
           do iR=1,nR
              ii=ii+1
              jj=jj+1

              if(iuse(ii)) then
                 jvec=jvec+1
                 r8data(ii)=dmap(jvec)
                 r8data(jj)=thmap(jvec)
              endif
           enddo
        enddo
     endif

     exit
  enddo

  deallocate(dmap,thmap,rhotmp,chtmp)
  deallocate(Phidum,Phitmp,Rvec,Zvec,iuse)

  if(ierr.eq.0) then
     ihave = iwant
     idata(3) = ihave
  endif

  nullify(r8data,idata)

end subroutine eqi_bdfast_gen
