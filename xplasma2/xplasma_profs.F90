module xplasma_profs

  ! module with routines that help to create various types of xplasma profiles

  use xplasma_obj
  use xplasma_ctran
  use xplasma_sol
  use xplasma_rzgeo
  use xplasma_flxint
  use eqi_rzbox_module

  implicit NONE

  private

  public :: xplasma_rzprof,xplasma_rzprof_fun
  public :: xplasma_brz,xplasma_brz_extrap
  public :: xplasma_irhofun
  public :: xplasma_geqdsk_rewrite
  public :: xplasma_wr_geqdsk,xplasma_rhopsi_gen,xplasma_rhopsi_find

  contains

    subroutine xplasma_brz_extrap(s,ier, ispline,sm_edge)

      ! create extrapolated B(R,Z) and Psi(R,Z) by standard method, when
      ! starting only with a prescribed boundary core (inverse representation)
      ! equilibrium.

      ! do not use this if a free boundary equilibrium is available...

      type(xplasma), pointer :: s
      integer, intent(out) :: ier ! exit status code 0=OK

      integer, intent(in), optional :: ispline   ! spline fit control
      !  DEFAULT =1, C1 Akima Hermite
      !   ...use 0 for piecewise linear; use 2 for C2 bicubic splines

      real*8, intent(in),optional :: sm_edge   ! edge smoothing
      !  smoothin in vicinity +/- sm_edge (meters) from plasma boundary
      !  i.e. option to smooth transition btw field internal to plasma and
      !  the field outside; default = no smoothing.

      !------------------------
      external eqm_brz_adhoc
      !------------------------

      call xplasma_brz(s,eqm_brz_adhoc,ier, ispline,sm_edge)

    end subroutine xplasma_brz_extrap

    subroutine xplasma_brz(s,userbvec,ier, ispline,sm_edge)

      ! create B(R,Z) and Psi(R,Z), extending fields already defined
      ! inside the plasma to the region outside-- using user provided
      ! subroutine...

      type(xplasma), pointer :: s

      external userbvec
      integer, intent(out) :: ier ! exit status code 0=OK

      integer, intent(in), optional :: ispline   ! spline fit control
      !  DEFAULT =1, C1 Akima Hermite
      !   ...use 0 for piecewise linear; use 2 for C2 bicubic splines

      real*8, intent(in),optional :: sm_edge   ! edge smoothing
      !  smoothin in vicinity +/- sm_edge (meters) from plasma boundary
      !  i.e. option to smooth transition btw field internal to plasma and
      !  the field outside; default = no smoothing.

      !-------------------------
      !  set up Bphi(R,Z), BR(R,Z), and BZ(R,Z) -- Akima-Hermite interpolation
      !  for axisymmetric case, setup psi(pol) also.

      !    data values from user supplied function  of form:

      !        subroutine userbvec(iv,zR,zZ,zphi,init,BR,BZ,BPHI,Psi,kpsi,ierr)
      !          iv -- vector dimension
      !          zR(iv),zZ(iv),zphi(iv)   ! coordinates at which to evaluate
      !          init -- action code
      !          BR(iv),BZ(iv),BPHI(iv)   ! field components returned
      !          Psi(iv)                  ! Psi values returned
      !          kpsi                     ! =1 if Psi values were set
      !          ierr                     ! completion code, 0=OK

      !  the subroutine must be capable of providing BR,BZ,BPHI, and Psi on
      !  points (R(i),Z(j)) which are on the extrapolated (R,Z) grids
      !  i.e. __RGRID and __ZGRID.  Only the points which are beyond the 
      !  plasma boundary need be filled; others can be zero.
      !

      !     normal calls to userbvec use init=0
      !         returns BR,BZ,BPHI, and Psi
      !     an initialization call uses init=1, returns BR=BZ=BPHI=czero
      !         also returns kpsi=1 if a poloidal flux function is to be
      !         returned. 

      !     a cleanup call uses init=2, returns BR=BZ=BPHI=czero

      !------------------------------
      !  This routine is not required for setting up B and Psi on an extended
      !  (R,Z) grid; the alternative is to call xplasma_rzprof separately for
      !  Psi and for each component of B-- as this routine itself does at the
      !  end...
      !------------------------------

      integer :: id_Rg,id_Zg
      integer :: nR,nZ,iR,iZ
      real*8, dimension(:), allocatable :: Rg,Zg,Rtmp,Ztmp,Phidum
      real*8, dimension(:,:), allocatable :: BRx,BZx,Bphix,Bmodx,Psix
      real*8 :: zdum1(1),zdum2(1),zdum3(1),zdum4(1)
      real*8 :: zsm_edge
      integer :: idum,jvec,ii,iersum,kpsi

      integer :: id_psi,id_BR,id_BZ,id_Bmod,id_Bphi,id_out,jspline,iertmp
      integer :: id_g,bphi_ccw,id_map,lrhomap
      logical :: axisymm,scrapeoff,outside_only

      real*8, dimension(:), pointer :: eqbuf
      integer, dimension(:), pointer :: idata

      real*8, parameter :: ZERO = 0.0d0
      !------------------------------

      sp => s
               
      call xplasma_global_info(s,ier, axisymm=axisymm,scrapeoff=scrapeoff, &
           bphi_ccw=bphi_ccw)
      if(ier.ne.0) return

      if(.not.axisymm) then
         ier=107
         call xplasma_errmsg_append(s,'  xplasma_brz still requires axisymmetry!')
         return
      endif

      if(.not.scrapeoff) then
         ier=109
         call xplasma_errmsg_append(s,'  xplasma_brz needs scrapeoff region defined.')
         return
      endif

      !  find needed resources...

      call xplasma_common_ids(s,ier, &
           id_psi=id_psi,id_BR=id_BR,id_BZ=id_BZ,id_Bmod=id_Bmod,id_g=id_g)
      if(ier.ne.0) return

      call xplasma_find_item(s,'__RGRID',id_Rg,ier)
      if(ier.ne.0) return

      call xplasma_find_item(s,'__ZGRID',id_Zg,ier)
      if(ier.ne.0) return

      call xplasma_find_item(s,'Bphi',id_Bphi,ier)
      if(ier.ne.0) return

      !  fetch grids

      call xplasma_grid_size(s,id_Rg,nR,ier)
      if(ier.ne.0) return

      call xplasma_grid_size(s,id_Zg,nZ,ier)
      if(ier.ne.0) return

      !  initialize userbvec

      zdum1=ZERO; zdum2=ZERO; zdum3=ZERO
      call userbvec(1,zdum1,zdum2,zdum3,1,zdum1,zdum2,zdum3,zdum4,kpsi,ier)

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_brz: external call to userbvec failed.')

         return
      endif

      if(kpsi.eq.1) then
         zsm_edge=ZERO  ! if (smooth Psi) is constructed, no further smoothing
      else
         ! If Psi(ext) imposed, smoothing is an option.
         zsm_edge=ZERO
         if(present(sm_edge)) then
            zsm_edge=sm_edge
         endif
      endif

      call xplasma_find_item(s,'__FASTMAP',id_map,iertmp,nf_noerr=.TRUE.)
      if(id_map.eq.0) then
         call eqi_fastinv_gen(id_Rg,id_Zg,nR,nZ,id_map,ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_brz: fast map setup failed.')
            return
         endif
      endif

      call xplasma_blackbox_retrieve(s, id_map, iertmp, &
           ia_ptr=idata, r8a_ptr=eqbuf)
      lrhomap=idata(1)

      jspline=1
      if(present(ispline)) jspline=ispline

      !  OK...

      iersum=0

      allocate(Rg(nR),Zg(nZ))
      allocate(Rtmp(nR*nZ),Ztmp(nR*nZ),Phidum(nR*nZ))
      allocate(BRx(nR,nZ),BZx(nR,nZ),Bphix(nR,nZ),Bmodx(nR,nZ),Psix(nR,nZ))

      call xplasma_grid(s,id_Rg,Rg,iertmp)
      call xplasma_grid(s,id_Zg,Zg,iertmp)

      Phidum=ZERO

      ii=0
      do iZ=1,nZ
         Rtmp(ii+1:ii+nR)=Rg     ! (1:nR)
         Ztmp(ii+1:ii+nR)=Zg(iZ)
         ii=ii+nR
      enddo

      call userbvec(nR*nZ,Rtmp,Ztmp,Phidum,0, &
              BRx,BZx,Bphix,Psix,idum,ier)
      iersum=iersum+ier

      !  generally profiles may have been evaluated at points beyond plasma
      !  boundary only.  If so Bphix(iR,iZ) will have zeroes; it is 
      !  advantageous to fill in Bphi now using g(rho)/R formula...

      jvec=0
      do iZ=1,nZ
         do iR=1,nR
            if(Bphix(iR,iZ).eq.ZERO) then
               jvec=jvec+1
               exit
            endif
         enddo
         if(jvec.gt.0) exit
      enddo

      outside_only = (jvec.gt.0)

      !  cleanup call for userbvec...

      zdum1=ZERO; zdum2=ZERO; zdum3=ZERO
      call userbvec(1,zdum1,zdum2,zdum3,2,zdum1,zdum2,zdum3,zdum4,idum,ier)
      iersum=iersum+ier

      if(iersum.ne.0) then
         call xplasma_errmsg_append(s,' ?xplasma_brz: extrapolated field setup failed.')
      else

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         if(outside_only) then
            ! code "id_fun_in=-1" indicates g(rho)/R
            call chk_auth('Bphi_RZ')
            if(id_Bphi.gt.0) then
               call xplasma_rzprof(s,'Bphi_RZ',id_out,ier, &
                    ispline=jspline, sm_edge=zsm_edge, id_fun_in=id_Bphi, &
                    data=Bphix,label='B_phi on (R,Z) grid',units='T')
            else
               call xplasma_rzprof(s,'Bphi_RZ',id_out,ier, &
                    ispline=jspline, sm_edge=zsm_edge, id_fun_in=-1, &
                    data=Bphix,label='B_phi on (R,Z) grid',units='T')
            endif
            iersum=iersum+ier

            call chk_auth('BR_RZ')
            call xplasma_rzprof(s,'BR_RZ',id_out,ier, &
                 ispline=jspline, sm_edge=zsm_edge, id_fun_in=id_BR, &
                 data=BRx,label='B_R on (R,Z) grid',units='T')
            iersum=iersum+ier

            call chk_auth('BZ_RZ')
            call xplasma_rzprof(s,'BZ_RZ',id_out,ier, &
                 ispline=jspline, sm_edge=zsm_edge, id_fun_in=id_BZ, &
                 data=BZx,label='B_Z on (R,Z) grid',units='T')
            iersum=iersum+ier

            if(kpsi.eq.1) then
               call chk_auth('Psi_RZ')
               call xplasma_rzprof(s,'Psi_RZ',id_out,ier, &
                    ispline=jspline, sm_edge=zsm_edge, id_fun_in=id_Psi, &
                    data=Psix, &
                    label='Poloidal flux on (R,Z) grid',units='Wb/rad')
               iersum=iersum+ier
            endif

         endif

         if(zsm_edge.gt.ZERO) then
            if(.not.outside_only) then
               call eqi_rzsmedg(Rg,nR,Zg,nZ,Bphix,zsm_edge,iertmp)
               call eqi_rzsmedg(Rg,nR,Zg,nZ,BRx,zsm_edge,iertmp)
               call eqi_rzsmedg(Rg,nR,Zg,nZ,BZx,zsm_edge,iertmp)
            endif
         endif

         Bmodx=sqrt(BRx**2+BZx**2+Bphix**2)

         if(.not.outside_only) then
            call chk_auth('BR_RZ')
            call xplasma_create_prof(s,'BR_RZ',id_Rg,id_Zg,BRx,id_out,ier, &
                 ispline=jspline,assoc_id=id_BR, &
                 label='B_R on (R,Z) grid',units='T')
            iersum=iersum+ier

            call chk_auth('BZ_RZ')
            call xplasma_create_prof(s,'BZ_RZ',id_Rg,id_Zg,BZx,id_out,ier, &
                 ispline=jspline,assoc_id=id_BZ, &
                 label='B_Z on (R,Z) grid',units='T')
            iersum=iersum+ier

            call chk_auth('BPhi_RZ')
            call xplasma_create_prof(s,'Bphi_RZ',id_Rg,id_Zg,Bphix,id_out,ier,&
                 ispline=jspline,assoc_id=id_Bphi, &
                 label='Btoroidal on (R,Z) grid',units='T')
            iersum=iersum+ier
         endif

         call chk_auth('BMod_RZ')
         call xplasma_create_prof(s,'Bmod_RZ',id_Rg,id_Zg,Bmodx,id_out,ier, &
              ispline=jspline,assoc_id=id_Bmod, &
              label='mod(B) on (R,Z) grid',units='T')
         iersum=iersum+ier

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)
      endif

      deallocate(Rg,Zg,Rtmp,Ztmp,Phidum)
      deallocate(BRx,BZx,Bphix,Bmodx,Psix)

      if(iersum.gt.0) ier=9999

      nullify(eqbuf,idata)

    CONTAINS
      subroutine chk_auth(zname)
        character*(*), intent(in) :: zname

        !  acquire ownership of existing profile, if necessary...

        !----------------------
        integer :: idp,iertmp
        character*32 :: zauth
        integer :: ildbg = 6
        !----------------------

        call xplasma_profId(s,zname,idp)
        if(idp.eq.0) return

        call xplasma_prof_info(s,idp,iertmp, author=zauth)
        if(zauth.ne.xplasma_xmhd) then
#ifdef __DEBUG
           write(ildbg,*) &
                ' %xplasma_brz(chk_auth): resetting author/owner of '// &
                trim(zname)
           write(ildbg,*) '  from "'//trim(zauth)//'" to "'// &
                trim(xplasma_xmhd)//'".'
#endif
           call xplasma_reset_author(s,idp,zauth,xplasma_xmhd,iertmp)
        endif
      end subroutine chk_auth

    end subroutine xplasma_brz

    !-------------------------------------------------------
    subroutine xplasma_rzprof(s,fname,id_out,ier, &
         id_Rgrid,id_Zgrid, &
         ispline,sm_edge, &
         id_fun_in,lamda, &
         label,units,data)

      ! create a profile f(R,Z) from existing data with extrapolation,
      ! or with array data provided.

      ! either an existing function (id_fun_in) or input data (data) or
      ! both must be provided.

      ! modes of use:
      !   (id_fun_in omitted, data(:,:) provided) -- just use the data given
      !   (id_fun_in provided, data(:,:) omitted) -- use existing profile
      !       at id_fun_in to give variation inside plasma; use extrapolation
      !       based on distance map: f_outside = f_bdy*exp(-d/lamda) where
      !       d is the distance from the plasma
      !   (both id_fun_in and data(:,:) provided) -- use existing profile
      !       at id_fun_in to give variation inside plasma: use data to 
      !       specify the variation beyond the plasma.

      !   ispline -- 0 = Bilinear, 1 = C1 Hermite, 2 = C2 Spline

      !   sm_edge -- Meters -- gives smoothing in vicinity of boundary
      !       a hat function convolution of half width sm_edge is applied
      !       in the vicinity of the boundary

      type(xplasma), pointer :: s

      character*(*), intent(in) :: fname  ! name of profile to create...

      integer, intent(out) :: id_out   ! ID of function just created
      integer, intent(out) :: ier      ! completion code 0=OK

      !------

      integer, intent(in), optional :: id_Rgrid,id_Zgrid  ! R & Z grids to use
      !  (default: __RGRID & __ZGRID)

      integer, intent(in), optional :: ispline  ! interpolation order
      ! 0 (default): piecewise linear; 1: C1 Hermite; 2: C2 Spline

      real*8, intent(in), optional :: sm_edge   ! edge smoothing control
      ! default: no smoothing

      integer,intent(in), optional :: id_fun_in ! f(rho) from which to 
      ! form f(R,Z) (default: 0, none)
      ! note: if id_fun_in = -1, use the formula bphi_ccw*g(rho)/R(rho,theta)

      real*8, intent(in), optional :: lamda     ! scrape off distance
      ! default: huge, making for flat extrapolation

      character*(*), intent(in), optional :: label,units  ! labeling info

      real*8, intent(inout), dimension(:,:), optional :: data  ! f(R,Z) data
      ! default: NONE
      ! data must be defined in extrapolated region-- internal region
      ! will be filled in...

      !--------------------------------
      integer :: id_Rg,id_Zg,nR,nZ,i,iZ,jvec
      integer :: jspline,isource,idf,id_dmap,maptype,iertmp,bphi_ccw
      real*8 :: zsm_edge,zlamda
      logical :: standard_RZ,outside_only,need_dmap,axisymm,scrapeoff,iRflag

      real*8, dimension(:,:), allocatable :: zdata
      real*8, dimension(:), allocatable :: zR,zZ,zZtmp,zwk1,zwk2,zwk3,zwk4
      real*8, dimension(:), allocatable :: zdist,zrho,zth
      logical, dimension(:), allocatable :: inside

      integer :: lrhomap,lchimap,ii,jj,kk,id_map
      real*8, dimension(:), pointer :: eqbuf
      integer, dimension(:), pointer :: idata

      real*8, parameter :: lamda_min=0.0001d0
      real*8, parameter :: zlarge = 1.0d35
      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter :: ONE = 1.0d0
      !--------------------------------

      id_out=0

      id_Rg=0; id_Zg=0
      if(present(id_Rgrid)) id_Rg=id_Rgrid
      if(present(id_Zgrid)) id_Zg=id_Zgrid

      call xplasma_ck_rzgrid(s,id_Rg,nR,id_Zg,nZ,standard_RZ,ier)
      if(ier.ne.0) return
               
      call xplasma_global_info(s,ier, axisymm=axisymm,scrapeoff=scrapeoff, &
           bphi_ccw=bphi_ccw)
      if(ier.ne.0) return

      jspline=0
      if(present(ispline)) jspline=ispline

      zsm_edge=ZERO
      if(present(sm_edge)) zsm_edge=sm_edge

      idf=0
      iRflag=.FALSE.
      if(present(id_fun_in)) then
         if(id_fun_in.gt.0) then
            call xplasma_ck_fun(s,id_fun_in,ier)
            if(ier.ne.0) return
            idf=id_fun_in
         else if(id_fun_in.eq.-1) then
            iRflag=.TRUE.  ! sign*g/R
            call xplasma_common_ids(s,ier, id_g=idf)
         endif
      else
         if(present(lamda)) then

            ier=9999
            call xplasma_errmsg_append(s,'  cannot construct f(R,Z): '//trim(fname))
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_rzprof: cannot use scrapeoff length for function')
            call xplasma_errmsg_append(s, &
                 '  extrapolation: no base function specified.')
            return
         endif
      endif

      zlamda=zlarge
      isource=0
      if(present(lamda)) then
         if(zlamda.le.ZERO) then
            zlamda=zlarge
         else
            zlamda=max(lamda_min,lamda)
         endif
         isource=isource+1
      endif

      if(present(data)) then
         need_dmap=.FALSE.
         isource=isource+1
         if((size(data,1).ne.nR).or.(size(data,2).ne.nZ)) then
            ier=510
            call xplasma_errmsg_append(s,'  cannot construct f(R,Z): '//trim(fname))
            call xplasma_errmsg_append(s, &
                 '  ?xplasma_rzprof: input data array dimensions do not match')
            call xplasma_errmsg_append(s, &
                 '   grid sizes of implicitly or explicitly chosen R&Z grids.')
         endif
      else
         need_dmap=.TRUE.
      endif

      if(isource.gt.1) then
         ier=9999
         call xplasma_errmsg_append(s,'  cannot construct f(R,Z): '//trim(fname))
         call xplasma_errmsg_append(s, &
              '  ?xplasma_rzprof: input data array and scrape-off extrapolation (lamda)')
         call xplasma_errmsg_append(s, &
              '   cannot both be present.')
         call xplasma_errmsg_append(s, &
              '   choose one or the other to define region beyond plasma boundary')
         return

      endif

      !----------------------------------------------
      outside_only = need_dmap.and.standard_RZ.and.(zsm_edge.eq.ZERO)

      allocate(zdata(nR,nZ),zR(nR),zZ(nZ),zZtmp(nR),zth(nR*nZ))
      allocate(zwk1(nR*nZ),zwk2(nR*nZ),zwk3(nR*nZ),zwk4(nR*nZ),zrho(nR*nZ))
      allocate(zdist(nR*nZ),inside(nR*nZ))

      call xplasma_grid(s,id_Rg,zR,ier); if (ier.ne.0) return
      call xplasma_grid(s,id_Zg,zZ,ier); if (ier.ne.0) return

      if(idf.eq.0) then

         ! just use the data provided to create a profile...

         zdata = data

      else

         if(standard_RZ) then
            maptype=3

            !  make sure rhomap is available

            call xplasma_ctrans(sp,iertmp, &
                 R_in=zR(nR/2), Z_in=zZ(nZ/2), &
                 rho_out=zrho(1), theta_out=zth(1), maptype=3)

            call xplasma_find_item(s,'__FASTMAP',id_map,iertmp)
            if(iertmp.ne.0) then
               call xplasma_errmsg_append(s, &
                ' %xplasma_rzprof: __FASTMAP not found: no scrapeoff region?')
               maptype=2
            else
               call xplasma_blackbox_retrieve(s, id_map, iertmp, &
                    ia_ptr=idata, r8a_ptr=eqbuf)
               lrhomap=idata(1)
               lchimap=idata(2)
               jj=lrhomap-nR-1
               kk=lchimap-nR-1
            endif
         else
            maptype=2
         endif
      
         jvec=0
         ii=-nR
         do iZ=1,nZ
            ii=ii+nR

            if(maptype.eq.2) then
               zZtmp=zZ(iZ)
               call xplasma_ctrans(s,.TRUE.,ier, R_in=zR,Z_in=zZtmp, &
                    rho_out=zrho(ii+1:ii+nR),theta_out=zth(ii+1:ii+nR), &
                    maptype=maptype)
               if(ier.ne.0) exit
            else
               jj=jj+nR
               kk=kk+nR
               zrho(ii+1:ii+nR)=eqbuf(jj+1:jj+nR)
               zth(ii+1:ii+nR)=eqbuf(kk+1:kk+nR)
            endif

            if(need_dmap) then
               ! using lamda extrapolation for pts beyond plasma bdy

               do i=1,nR
                  if(zrho(ii+i).gt.ONE) then
                     jvec=jvec+1
                     zwk1(jvec)=zZ(iZ)
                     zwk2(jvec)=zR(i)
                     inside(ii+i)=.FALSE.
                     zrho(ii+i)=ONE
                  else
                     inside(ii+i)=.TRUE.
                     zdist(ii+i)=ZERO
                  endif
               enddo
            else
               do i=1,nR
                  if(zrho(ii+i).le.ONE) then
                     jvec=jvec+1
                     zwk1(jvec)=zR(i)
                     zwk2(jvec)=zrho(ii+i)
                     zwk3(jvec)=zth(ii+i)
                     inside(ii+i)=.TRUE.
                  else
                     inside(ii+i)=.FALSE.
                  endif
               enddo

            endif
         enddo

         if(need_dmap) then
            call xplasma_bdfind(s,jvec,zwk2(1:jvec),zwk1(1:jvec),ier, &
                 maptype=maptype,outside_only=outside_only, &
                 theta_out=zwk3(1:jvec), dist=zwk4(1:jvec))
         else

            ! use idf & iRflag

            call xplasma_eval_prof(s,idf, &
                 xplasma_rho_coord,zwk2(1:jvec), &
                 xplasma_theta_coord,zwk3(1:jvec), &
                 zwk4(1:jvec),ier)

            if(iRflag) then
               zwk4(1:jvec)=bphi_ccw*zwk4(1:jvec)/zwk1(1:jvec)
            endif

         endif

         if(ier.eq.0) then
            if(need_dmap) then
               jvec=0
               ii=-nR
               do iZ=1,nZ
                  ii=ii+nR
                  do i=1,nR
                     if(.not.inside(ii+i)) then
                        jvec=jvec+1
                        zth(ii+i)=zwk3(jvec)
                        zdist(ii+i)=zwk4(jvec)
                     endif
                  enddo
               enddo

               call xplasma_eval_prof(s,idf, &
                    xplasma_rho_coord,zrho, xplasma_theta_coord,zth, &
                    zwk2,ier)

               if(ier.eq.0) then
                  ii=-nR
                  do iZ=1,nZ
                     ii=ii+nR
                     do i=1,nR
                        zdata(i,iZ)=zwk2(ii+i)
                        if((zlamda.lt.zlarge).and.(zdist(ii+i).gt.ZERO)) then
                           zdata(i,iZ)=zdata(i,iZ)*exp(-zdist(ii+i)/zlamda)
                        endif
                     enddo
                  enddo
               endif

            else
               ! using data for pts beyond plasma bdy

               jvec=0
               ii=-nR
               do iZ=1,nZ
                  ii=ii+nR
                  do i=1,nR
                     if(inside(ii+i)) then
                        jvec=jvec+1
                        zdata(i,iZ)=zwk4(jvec)  ! interior-- use function eval
                        data(i,iZ)=zwk4(jvec)   ! pass back also...
                     else
                        zdata(i,iZ)=data(i,iZ)  ! exterior-- use passed data
                     endif
                  enddo
               enddo

            endif

         endif
      endif

      if(ier.eq.0) then

         if(zsm_edge.gt.ZERO) then
            if(idf.eq.0) then
               sp => s
               call eqi_rzsmedg(zR,nR,zZ,nZ,zdata,zsm_edge,iertmp)
            else
               ! use module-internal blending routine; mesh with core data
               call xpblend(s,zR,nR,zZ,nZ,zdata,zsm_edge,idf,iRflag,bphi_ccw, &
                    maptype, iertmp)
            endif
         endif

         call xplasma_create_2dprof(s,fname,id_Rg,id_Zg,zdata,id_out,ier, &
              ispline=jspline,assoc_id=idf,label=label,units=units)

      endif

      deallocate(inside,zdata,zR,zZ,zZtmp,zwk1,zwk2,zwk3,zwk4,zrho,zth,zdist)

      nullify(eqbuf,idata)

    end subroutine xplasma_rzprof

    !-------------------------------------------------------
    subroutine xplasma_rzprof_fun(s,fname,user_fcn,user_iarg,id_out,ier, &
         id_Rgrid,id_Zgrid,ispline,sm_edge,id_fun_in,no_extrap, &
         label,units)

      ! create a profile f(R,Z) from a user provided callable function
      ! possibly combined with existing data.

      ! the callable function covers the entire region inside the limiters
      ! if id_fun_in is omitted; it covers only the region outside the plasma
      ! and inside the limiters if id_fun_in is included.

      ! modes of use:
      !   (id_fun_in omitted) -- just use the given callable function.
      !   (id_fun_in provided)-- use existing profile
      !       at id_fun_in to give variation inside plasma; use the given
      !       function to define the variation beyond the plasma boundary.

      ! the region beyond the limiters is defined by a simple extrapolation.
      ! unless no_extrap is set, in which case the function is used even
      ! beyond the limiters.

      !   sm_edge -- Meters -- gives smoothing in vicinity of boundary
      !       a hat function convolution of half width sm_edge is applied
      !       in the vicinity of the boundary

      type(xplasma), pointer :: s

      character*(*), intent(in) :: fname  ! name of profile to create...

      real*8, external :: user_fcn        ! user provided function...
      integer, intent(in) :: user_iarg    ! argument for user_fcn:

      !  real*8 function user_fcn(user_iarg,R,Z,Phi,ier)

      integer, intent(out) :: id_out   ! ID of function just created
      integer, intent(out) :: ier      ! completion code 0=OK

      !------

      integer, intent(in), optional :: id_Rgrid,id_Zgrid  ! R & Z grids to use
      !  (default: __RGRID & __ZGRID)

      integer, intent(in), optional :: ispline  ! interpolation order
      ! 0 (default): piecewise linear; 1: C1 Hermite; 2: C2 Spline

      real*8, intent(in), optional :: sm_edge   ! edge smoothing control
      ! default: no smoothing

      integer,intent(in), optional :: id_fun_in ! f(rho) from which to 
      ! form f(R,Z) (default: 0, none)
      ! note: if id_fun_in = -1, use the formula bphi_ccw*g(rho)/R(rho,theta)

      logical, intent(in), optional :: no_extrap ! T to suppress extrapolation
      ! beyond limiter.  In this case user_fcn calls are used for (R,Z) even
      ! for points outside the limiter.

      character*(*), intent(in), optional :: label,units  ! labeling info

      !--------------------------------
      integer :: id_Rg,id_Zg,nR,nZ,i,iZ,jvec
      integer :: jspline,idf,maptype,iertmp,bphi_ccw
      real*8 :: zsm_edge,zphi0
      real*8, dimension(:,:), allocatable :: zdata
      real*8, dimension(:), allocatable :: zR,zZ,zwkZ,zrho,zth
      real*8, dimension(:), allocatable :: zdist_lim,zfeval,zwkR
      logical, dimension(:), allocatable :: inside
      logical, dimension(:,:), allocatable :: inside_lim
      logical :: standard_RZ,extrap,axisymm,scrapeoff,iRflag

      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter :: ONE = 1.0d0
      !--------------------------------
      !  dmc: this routine might be sped up by making calls covering
      !   all RxZ instead of just one row at a time.  Other routines
      !   have been optimized in this way... but not yet this one.

      id_out=0

      id_Rg=0; id_Zg=0
      if(present(id_Rgrid)) id_Rg=id_Rgrid
      if(present(id_Zgrid)) id_Zg=id_Zgrid

      call xplasma_ck_rzgrid(s,id_Rg,nR,id_Zg,nZ,standard_RZ,ier)
      if(ier.ne.0) return
               
      call xplasma_global_info(s,ier, axisymm=axisymm,scrapeoff=scrapeoff, &
           bphi_ccw=bphi_ccw)
      if(ier.ne.0) return

      jspline=0
      if(present(ispline)) jspline=ispline

      zsm_edge=0
      if(present(sm_edge)) zsm_edge=sm_edge

      idf=0
      iRflag=.FALSE.
      if(present(id_fun_in)) then
         if(id_fun_in.gt.0) then
            call xplasma_ck_fun(s,id_fun_in,ier)
            if(ier.ne.0) return
            idf=id_fun_in
         else if(id_fun_in.eq.-1) then
            iRflag=.TRUE.  ! sign*g/R
            call xplasma_common_ids(s,ier, id_g=idf)
         endif
      endif

      extrap=.TRUE.
      if(present(no_extrap)) extrap = .not.no_extrap

      !------------------------------------------
      allocate(zdata(nR,nZ),zR(nR),zZ(nZ),zwkZ(nR))
      allocate(zrho(nR),zth(nR),zdist_lim(nR),zfeval(nR),zwkR(nR))
      allocate(inside(nR),inside_lim(nR,nZ))

      call xplasma_grid(s,id_Rg,zR,ier); if (ier.ne.0) return
      call xplasma_grid(s,id_Zg,zZ,ier); if (ier.ne.0) return

      if(standard_RZ) then
         maptype=3
      else
         maptype=2
      endif

      zphi0 = ZERO

      do iZ=1,nZ
         zwkZ=zZ(iZ)

         if(idf.eq.0) then
            if(extrap) then
               call xplasma_lim_distance(s,zR,zwkZ,zdist_lim,ier, &
                    maptype=maptype)
               if(ier.ne.0) exit
            else
               zdist_lim=ZERO
            endif

            do i=1,nR
               if(zdist_lim(i).gt.ZERO) then
                  zdata(i,iZ)=ZERO
                  inside_lim(i,iZ)=.FALSE.
               else
                  zdata(i,iZ) = user_fcn(user_iarg,zR(i),zwkZ(1),ZERO,ier)
                  if(ier.ne.0) exit
                  inside_lim(i,iZ)=.TRUE.
               endif
            enddo
            if(ier.ne.0) exit

         else

            call xplasma_ctrans(s,.TRUE.,ier, R_in=zR,Z_in=zwkZ, &
                 rho_out=zrho,theta_out=zth, maptype=maptype)
            if(ier.ne.0) exit

            if(extrap) then
               call xplasma_lim_distance(s,zR,zwkZ,zdist_lim,ier, &
                    maptype=maptype)
               if(ier.ne.0) exit
            else
               zdist_lim=ZERO
            endif

            jvec=0
            do i=1,nR
               inside_lim(i,iZ)=.TRUE.
               if(zdist_lim(i).gt.ZERO) then
                  zdata(i,iZ)=ZERO
                  inside(i)=.FALSE.
                  inside_lim(i,iZ)=.FALSE.
               else if(zrho(i).gt.ONE) then
                  zdata(i,iZ) = user_fcn(user_iarg,zR(i),zwkZ(1),ZERO,ier)
                  if(ier.ne.0) exit
                  inside(i)=.FALSE.
               else
                  jvec=jvec+1
                  inside(i)=.TRUE.
                  zwkR(jvec)=zR(i)
                  zrho(jvec)=zrho(i)
                  zth(jvec)=zth(i)
               endif
            enddo
            if(ier.ne.0) exit

            if(jvec.gt.0) then
               call xplasma_eval_prof(s,idf, &
                    xplasma_rho_coord,zrho(1:jvec), &
                    xplasma_theta_coord,zth(1:jvec), &
                    zfeval(1:jvec),ier)
               if(ier.ne.0) exit

               if(iRflag) then
                  zfeval(1:jvec) = bphi_ccw*zfeval(1:jvec)/zwkR(1:jvec)
               endif

               jvec=0
               do i=1,nR
                  if(inside(i)) then
                     jvec=jvec+1
                     zdata(i,iZ)=zfeval(jvec)
                  endif
               enddo
            endif
         endif

      enddo

      if(ier.eq.0) then

         sp => s

         if(extrap) then
            ! compute simple extrapolation for beyond-limiter data points
            call eqi_rzfx(zR,nR,zZ,nZ,inside_lim,zdata)
         endif

         if(zsm_edge.gt.ZERO) then
            ! optional smoothing...
            if(idf.eq.0) then
               call eqi_rzsmedg(zR,nR,zZ,nZ,zdata,zsm_edge,iertmp)
            else
               ! use module-internal blending routine; mesh with core data
               call xpblend(s,zR,nR,zZ,nZ,zdata,zsm_edge,idf,iRflag,bphi_ccw, &
                    maptype, iertmp)
            endif
         endif

         call xplasma_create_2dprof(s,fname,id_Rg,id_Zg,zdata,id_out,ier, &
              ispline=jspline,assoc_id=idf,label=label,units=units)

      endif

      deallocate(zdata,zR,zZ,zwkZ,zrho,zth,zdist_lim,zfeval,zwkR)
      deallocate(inside,inside_lim)

    end subroutine xplasma_rzprof_fun

    !-------------------------------------------------------
    subroutine xplasma_ck_fun(s,idf,ier)

      !  ** private **
      !  verify function id -- must be fcn of (rho,theta)

      type(xplasma), pointer :: s

      integer, intent(in) :: idf
      integer, intent(out) :: ier

      !------------------
      integer :: id_grid1,id_grid2,icoord1,icoord2
      character*32 zname1,zname2
      !------------------

      do
         call xplasma_prof_info(s,idf,ier, gridId1=id_grid1, gridId2=id_grid2)
         if(ier.ne.0) exit

         call xplasma_grid_info(s,id_grid1,ier, coord=icoord1, name=zname1)
         if(ier.ne.0) exit

         if(id_grid2.gt.0) then
            call xplasma_grid_info(s,id_grid2,ier, coord=icoord2, name=zname2)
            if(ier.ne.0) exit

         else
            icoord2=0
         endif

         if((icoord1.eq.xplasma_rho_coord).or. &
              (icoord1.eq.xplasma_theta_coord)) then

            continue

         else

            ier=9999
            call xplasma_errmsg_append(s, &
                 '   passed profile needs to be defined vs. flux coordinates.')
            call xplasma_errmsg_append(s, &
                 '   instead it is a function of: '//trim(zname1))
            exit
         endif

         if((icoord2.eq.0).or.(icoord2.eq.xplasma_rho_coord).or. &
              (icoord2.eq.xplasma_theta_coord)) then

            continue

         else

            ier=9999
            call xplasma_errmsg_append(s, &
                 '   passed profile needs to be defined vs. flux coordinates.')
            call xplasma_errmsg_append(s, &
                 '   instead it is a function of: '//trim(zname2))
            exit
         endif

         exit
      enddo

      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ...error in xplasma_profs module.')
      endif

    end subroutine xplasma_ck_fun

    !-------------------------------------------------------
    subroutine xplasma_ck_rzgrid(s,id_Rgrid,nR,id_Zgrid,nZ,standard_RZ,ier)

      !  ** private ** 
      !  summary info on R & Z grids

      type(xplasma), pointer :: s

      integer, intent(inout) :: id_Rgrid,id_Zgrid ! grid IDs: use default if 0
      integer, intent(out) :: nR,nZ               ! grid sizes
      logical, intent(out) :: standard_RZ         ! T if standard grid is used
      integer, intent(out) :: ier                 ! error code, 0=OK

      !-------------------------
      integer :: icoord,isize,id_Rgrid0,id_Zgrid0
      character*32 :: zname
      !-------------------------

      standard_RZ=.FALSE.
      nR=0
      nZ=0

      call xplasma_find_item(s,'__Rgrid',id_Rgrid0,ier)
      if(ier.ne.0) return

      call xplasma_find_item(s,'__Zgrid',id_Zgrid0,ier)
      if(ier.ne.0) return

      if(id_Rgrid.eq.0) then
         id_Rgrid=id_Rgrid0
      endif

      if(id_Zgrid.eq.0) then
         id_Zgrid=id_Zgrid0
      endif

      if((id_Rgrid.eq.id_Rgrid0).and.(id_Zgrid.eq.id_Zgrid0)) then
         standard_RZ=.TRUE.
      endif

      call xplasma_grid_info(s,id_Rgrid,ier, coord=icoord, size=nR, name=zname)
      if(ier.eq.0) then
         if(icoord.ne.xplasma_R_coord) then
            ier=701
            call xplasma_errmsg_append(s,' passed id grid name: '//trim(zname))
         endif
      endif
      if(ier.ne.0) return

      call xplasma_grid_info(s,id_Zgrid,ier, coord=icoord, size=nZ, name=zname)
      if(ier.eq.0) then
         if(icoord.ne.xplasma_Z_coord) then
            ier=701
            call xplasma_errmsg_append(s,' passed id grid name: '//trim(zname))
         endif
      endif
      if(ier.ne.0) return

    end subroutine xplasma_ck_rzgrid

    !=====================================

    subroutine xpblend(s,zR,nR,zZ,nZ,zdata,zsm,idf,iRflag,ibccw,maptype,ier)

      !  ** PRIVATE **

      !  smooth zdata(R,Z) beyond edge region by matching point and slope
      !  of interior function (idf), in extrapolated (rho,theta) space and
      !  then mapping back to (R,Z).

      !  if iRflag is set, expect idf=id_g and use bphi_ccw*f(rho)/R(rho,theta)
      !      (i.e. the toroidal field) for the interior function.

      type (xplasma), pointer :: s

      integer :: nR,nZ        ! grid sizes
      real*8 :: zR(nR),zZ(nZ) ! grids
      real*8 :: zdata(nR,nZ)  ! data

      real*8 :: zsm           ! smoothing distance

      integer :: idf          ! interior matching function ID

      logical :: iRflag       ! T to use ibccw*f(rho)/R; idf points to f(rho)
      integer :: ibccw        ! +/- 1 for use if (iRflag) true

      integer :: maptype      ! mapping selector for (R,Z) -> (rho,theta)

      integer, intent(out) :: ier    ! status code; 0=OK

      !-----------------------------------------------------------
      integer :: ii,jvec,jvec2,iR,iZ
      real*8 :: del,del1,del2,dR1,dZ1
      real*8 :: dat0,dat1,dat2,datx,b,denom,dxtrap,xx,ff

      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter :: HALF = 0.5d0
      real*8, parameter :: ONE  = 1.0d0
      real*8, parameter :: TWO  = 2.0d0
      real*8, parameter :: THREE  = 3.0d0

      ! data for 1 row @ fixed Z:
      real*8, dimension(:), allocatable :: zzwk,zrho,zth,zRb,zZb,z1vec

      ! collected data for edge blending:
      integer :: isizx,isizf
      real*8, dimension(:), allocatable :: delx,rhox,thx,zRx,datbx, &
           Rinx,Zinx,rho_inx,th_inx,datinx
      integer, dimension(:), allocatable :: iRsave,iZsave
      !-----------------------------------------------------------

      del2 = min( (zR(nR)-zR(1))/max(40,nR), (zZ(nZ)-zZ(1))/max(40,nZ) )
      del1 = HALF*del2

      ier=0

      allocate(zzwk(nR),zrho(nR),zth(nR),zRb(nR),zZb(nR),z1vec(nR))
      z1vec = ONE      ! vector of length nR of 1s...

      isizf=20
10    continue
      if(isizf.eq.20) then
         isizf=10
      else if(isizf.eq.10) then
         isizf=5
      else 
         isizf=isizf-1
         if(isizf.eq.0) then
            call xplasma_errmsg_append(s,' ?xpblend: too many points in range for smoothed blending.')
            ier=9999
            return
         endif
      endif

      isizx = (nR*nZ)/isizf
      if(allocated(delx)) then
         deallocate(delx,rhox,thx,datbx,iRsave,iZsave)
         deallocate(Rinx,Zinx,rho_inx,th_inx,datinx)
      endif
      allocate(delx(isizx),rhox(isizx),thx(isizx),zRx(isizx))
      allocate(datbx(isizx),iRsave(isizx),iZsave(isizx))
      allocate(Rinx(2*isizx),Zinx(2*isizx),rho_inx(2*isizx),th_inx(2*isizx))
      allocate(datinx(2*isizx))

      jvec = 0
      jvec2 = 0

      do iZ=1,nZ
         zzwk = zZ(iZ) ! expand this Z to vector of length nR...

         call xplasma_ctrans(s,.TRUE.,ier, R_in=zR, Z_in=zzwk, &
              rho_out=zrho, theta_out=zth, maptype=maptype)
         if(ier.ne.0) exit

         call xplasma_ctrans(s,.TRUE.,ier, rho_in=z1vec, theta_in=zth, &
              R_out=zRb, Z_out=zZb)
         if(ier.ne.0) exit

         do iR=1,nR
            if(zrho(iR).gt.ONE) then
               del = sqrt((zR(iR)-zRb(iR))**2 + (zZ(iZ)-zZb(iR))**2)
               dR1 = (zRb(iR)-zR(iR))/del
               dZ1 = (zZb(iR)-zZ(iZ))/del
               if((ZERO.lt.del).AND.(del.lt.zsm)) then
                  ! point in range for smoothing

                  jvec = jvec + 1
                  if(jvec.gt.isizx) go to 10  ! check bound

                  delx(jvec)=del
                  iRsave(jvec)=iR
                  iZsave(jvec)=iZ
                  rhox(jvec)=ONE
                  thx(jvec)=zth(iR)
                  zRx(jvec)=zRb(iR)

                  Rinx(jvec2+1)=zRb(iR) + del1*dR1
                  Rinx(jvec2+2)=zRb(iR) + del2*dR1

                  Zinx(jvec2+1)=zZb(iR) + del1*dZ1
                  Zinx(jvec2+2)=zZb(iR) + del2*dZ1

                  jvec2 = jvec2 + 2
               endif
            endif
         enddo

      enddo
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ?error 1 in xpblend internal routine')
         return
      endif

      ! evaluate data at boundary

      call xplasma_eval_prof(s,idf, &
                 xplasma_rho_coord,rhox(1:jvec), &
                 xplasma_theta_coord,thx(1:jvec), &
                 datbx(1:jvec),ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ?error 2 in xpblend internal routine')
         return
      endif

      if(iRflag) then
         datbx(1:jvec) = ibccw*datbx(1:jvec)/zRx(1:jvec)
      endif

      ! evaluate interior points to allow matching assessment
      ! 1st convert (R,Z) to (rho,theta); use maptype=2 for reasonable accuracy

      call xplasma_ctrans(s,.TRUE.,ier, &
           R_in=Rinx(1:jvec2),Z_in=Zinx(1:jvec2), &
           rho_out=rho_inx(1:jvec2),theta_out=th_inx(1:jvec2), &
           maptype=2)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ?error 3 in xpblend internal routine')
         return
      endif

      ! evaluate interior data 

      call xplasma_eval_prof(s,idf, &
                 xplasma_rho_coord,rho_inx(1:jvec2), &
                 xplasma_theta_coord,th_inx(1:jvec2), &
                 datinx(1:jvec2),ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ?error 2 in xpblend internal routine')
         return
      endif

      if(iRflag) then
         datinx(1:jvec2) = ibccw*datinx(1:jvec2)/Rinx(1:jvec2)
      endif

      ! find linear extrapolation along fixed exterior theta lines
      ! i.e. a profile with continuous value and 1st derivative; 2nd
      ! derivative definitely not continuous.

      denom = del1*del2*(del2-del1)

      jvec2=0
      do ii=1,jvec
         dat0=datbx(ii)
         dat1=datinx(jvec2+1)
         dat2=datinx(jvec2+2)
         jvec2 = jvec2 + 2

         !  3 pts: (0,dat0), (del1,dat1), (del2,dat2) -- parabolic fit
         !          bdy       1st int.     2nd int.
         !  d(x) = A*x**2 + B*x + C    d(0)=C=dat0; d'(0)=B
         !    denom = del1*del2*(del2-del1)
         !        A = -[(dat1-dat0)*del2 + (dat2-dat0)*del1]/denom
         !        B =  [(dat1-dat0)*del2**2 + (dat2-dat0)*del1**2]/denom

         B = ((dat1-dat0)*del2**2 - (dat2-dat0)*del1**2)/denom

         del=delx(ii)
         xx=del/zsm             ! normalized extrapolation distance

         dxtrap = dat0 - B*del  ! linear extrapolation of interior data

         iR = iRsave(ii)
         iZ = iZsave(ii)

         ! blend using Hermite cubic f*data + (1-f)*dxtrap
         ! f(0)=0, f'(0)=0, f(1)=1, f'(1)=0
         !     => f(x) = x**2*(-2*x + 3)

         ff = xx*xx*(THREE-TWO*xx)

         zdata(iR,iZ) = ff*zdata(iR,iZ) + (ONE-ff)*dxtrap
      enddo

    end subroutine xpblend

    !------------------------------------------
    subroutine xplasma_irhofun(s,id_axis,zlbl,inprof,zprof,iflag,id,ierr, &
         smdelx)

      !  make XPLASMA profile -- integrated quantity; smooth by 1/2 zone width
      !  to insure smooth derivative from spline

      type (xplasma), pointer :: s

      integer, intent(in) :: id_axis    ! axis id -- must be rho or akin to rho

      character*(*), intent(in) :: zlbl ! name for profile function to create

      integer, intent(in) :: inprof     ! size of the integrated data profile

      real*8, intent(in) :: zprof(inprof) ! the integrated data provided
              ! if inprof = size(x_axis) zprof(1)=0 is expected
              ! if inprof = size(x_axis)-1 the axial data point
              !    is presumed to be omitted.

      integer, intent(in) :: iflag      ! =1 -- volume normalization;
              ! derivative evaluations -> W/m^3, #/sec/m^3, etc.
              ! =2 -- area normalization;
              ! derivative evaluations -> A/m^2 (current density).

              ! >2 -- ID of normalization profile, should be akin to vol(rho)
              !       or area(rho)

      integer, intent(out) :: id        ! id of stored profile (if successful)
      integer, intent(out) :: ierr      ! completion code, 0=OK.

      real*8, intent(in), optional :: smdelx   ! smoothing width (option)

      !  a single zone width hat function smooth (effectively, half a zone
      !  width) is the minimum used by eqi_irhofun.  For more smoothing
      !  than this, give smdelx = [a number greater than x(2)-x(1)] where
      !  x(1:2) are the 1st two grid points identified via id_axis.

      !  Note, the actual smoothing width is proportional to the grid
      !  resolution and so could vary across the grid; a fixed multiplier
      !  of smdelx/(x(2)-x(1)) would be applied in this case.  In the
      !  evaluation routine (eqi_irhofun) there is an upper limit imposed
      !  which corresponds a multiplier that yields (1/4) the profile width
      !  at the grid spacing at rho=0.

      !  if an error occurs and ierr is set, id=0 will be returned.
      !-----------------------------------------------------------------
      integer :: ierr_save
      real*8 :: zsm
      !-----------------------------------------------------------------

      sp => s

      zsm = 0.0d0
      if(present(smdelx)) zsm = smdelx

      call eqi_irhofun(id_axis,zlbl,inprof,zprof,iflag,zsm,id,ierr)

    end subroutine xplasma_irhofun

    !-------------------------------------------------------
    subroutine xplasma_geqdsk_rewrite(filename,ier)

      !  rewrite g-eqdsk file e.g. under new filename, after reading
      !  (e.g. inside an eqi_fromgeqdsk call).

      !  WARNING: the information written is based on the most recent
      !  read of g-eqdsk information.  Information to reconstruct the
      !  g-eqdsk from the data as originally read is NOT saved with each
      !  xplasma object.  Hence: no s pointer argument.

      character*(*), intent(in) :: filename
      integer, intent(out) :: ier

      call eqm_wr_geqdsk(filename,ier)
      if(ier.ne.0) ier=9999

    end subroutine xplasma_geqdsk_rewrite

    !-------------------------------------------------------
    subroutine xplasma_wr_geqdsk(s,ierr, &
         lun_geqdsk, filename, label, Rmin, Rmax, Zmin, Zmax, &
         cur, id_qprof, id_pprof, id_psi_in, &
         nh, nv, nbdy)

      !  Build and write a G-eqdsk file (disk file) from current xplasma
      !  contents.  This differs from "xplasma_geqdsk_rewrite" as the latter
      !  just echoes data read in from another G-eqdsk file or MDSplus data
      !  structure.  This actual constructs the information from current
      !  xplasma contents -- interpolation is involved.

      !  The xplasma object must contain a complete free boundary or 
      !  extrapolated equilibrium, so that Psi(R,Z) covering a rectangle
      !  enclosing the plasma is defined.

      !-----------------
      !  required:

      type (xplasma), pointer :: s  ! XPLASMA object containing equilibrium

      integer, intent(out) :: ierr  ! status code on exit: 0=OK

      !-----------------
      !  optional:

      integer, intent(in), optional :: lun_geqdsk  
      ! FORTRAN LUN on which to write file.  If omitted, the LUN used for
      ! reading G-eqdsk files, in the geqdsk_mds library, is used.
      !   (call geq_getlun(ilun)) (geqdsk_mod default value: 77 as of 7/2006).

      character*(*), intent(in), optional :: filename
      ! default " "; if non-blank, it is the
      ! name of the file to write.  If blank or omitted, the code simply 
      ! writes to lun_geqdsk; it would be up to the caller to OPEN a file.

      character*(*), intent(in), optional :: label
      ! default " "; default means: use xplasma global label; 
      ! if non-blank, the 1st 48 characters are used as a label string in the
      ! G-eqdsk file being written.

      real*8, intent(in), optional :: Rmin,Rmax, Zmin,Zmax
      ! [R,Z] domain over which Psi(R,Z) is written in the G-eqdsk file.
      ! default: use the [R,Z] grid limits already stored in xplasma.  If
      ! explicit limits are provided, overriding the defaults, these must
      ! not extend beyond the grid limits.

      real*8, intent(in), optional :: cur
      ! total plasma current.  Default: use value implied by equilibrium data.

      integer, intent(in), optional :: id_qprof
      ! ID of profile f(rho) defining q(rho).  Default: use value derived from
      ! equilibrium-- calculated here, it will be named Q_EQDSK

      integer, intent(in), optional :: id_pprof
      ! ID of profile f(rho) defining equilibrium pressure.  Default: use value
      ! previously tagged as equilibrium pressure.  If this is not available,
      ! the code will attempt to compute a pressure using the surface averaged
      ! GS equation JxB=grad(P).  This could lead to a non-physical result if
      ! there are errors in the equilibrium data already provided, or if the
      ! assumption of a scalar P is inappropriate.  If a pressure profile is
      ! calculated, it will be named P_EQDSK-- done using an xplasma_gen_p
      ! call.

      integer, intent(in), optional :: id_psi_in
      ! ID of associated pair of Psi(rho) or Psi(R,Z) profiles
      ! If defaulted, the code looks for standard names w/in xplasma object.

      integer, intent(in), optional :: nh,nv
      ! number of horizontal and vertical grid points, respectively.  If
      ! defaulted, the sizes of the [R,Z] grids are used.  NOTE that nh also
      ! controls the size of the 1d profiles f,ff',P,P',and q written in the
      ! G-eqdsk profiles.  These 1d profiles are written over an implied 
      ! evenly spaced Psi grid going from Psi(min) at the axis to Psi(max)
      ! at the boundary.
      
      integer, intent(in), optional :: nbdy
      ! number of points to use to described plasma boundary and limiter.
      ! default: 200.

      !-------------------------------------
      integer :: ilun,istat,ifile,ilen,id_Rc,id_Zc,idp,idq,ii,idpsi
      integer :: inh,inv,inbdy,inum
      integer :: idwk,iertmp

      character*48 glabel
      real*8 :: zRmin,zRmax,zZmin,zZmax,zcur

      !  for plasma current estimate:
      integer, parameter :: inxi=21
      real*8 :: zxi(inxi),zii(inxi)

      integer :: inbdy_dflt=201
      !-------------------------------------

      if(present(lun_geqdsk)) then
         ilun = lun_geqdsk
      else
         call eqi_geq_getlun(ilun)
      endif

      ierr = 0
      ifile= 0
      if(present(filename)) then
         open(unit=ilun,file=trim(filename),status='unknown',iostat=istat)
         if(istat.ne.0) then
            ierr = 9999
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_profs(xplasma_wr_geqdsk): file open failure: '// &
                 trim(filename))
         else
            ifile=1
         endif
      endif
      if(ierr.ne.0) return

      do

         ! label...

         if(present(label)) then
            glabel=label
         else
            glabel='xplasma_wr_geqdsk:'
            ilen=len(trim(glabel))
            call xplasma_global_info(s,ierr,initLabel=glabel(ilen+1:))
         endif
         if(ierr.ne.0) exit
            
         ! Rmin&max, Zmin&max...

         call xplasma_find_item(s,'__R_coord',id_Rc,ierr); if(ierr.ne.0) exit
         call xplasma_find_item(s,'__Z_coord',id_Zc,ierr); if(ierr.ne.0) exit

         call xplasma_coord_info(s,id_Rc,ierr, ngrids=inum); if(ierr.ne.0) exit
         if(inum.eq.0) then
            ierr=615  ! need R & Z rectangle; Psi(R,Z)
            exit
         endif

         call xplasma_coord_info(s,id_Zc,ierr, ngrids=inum); if(ierr.ne.0) exit
         if(inum.eq.0) then
            ierr=615  ! need R & Z rectangle; Psi(R,Z)
            exit
         endif

         call xplasma_coord_info(s,id_Rc,ierr, xmin=zRmin, xmax=zRmax)
         if(ierr.ne.0) exit

         call xplasma_coord_gridId(s,id_Rc,1,idwk,ierr)
         if(ierr.ne.0) exit

         call xplasma_grid_size(s,idwk,inh,ierr)
         if(ierr.ne.0) exit

         call xplasma_coord_info(s,id_Zc,ierr, xmin=zZmin, xmax=zZmax)
         if(ierr.ne.0) exit

         call xplasma_coord_gridId(s,id_Zc,1,idwk,ierr)
         if(ierr.ne.0) exit

         call xplasma_grid_size(s,idwk,inv,ierr)
         if(ierr.ne.0) exit

         if(present(Rmin)) then
            zRmin=max(zRmin,Rmin)
         endif

         if(present(Rmax)) then
            zRmax=min(zRmax,Rmax)
         endif

         if(present(Zmin)) then
            zZmin=max(zZmin,Zmin)
         endif

         if(present(Zmax)) then
            zZmax=min(zZmax,Zmax)
         endif

         !  total current 

         if(present(cur)) then
            zcur=cur
         else
            !  estimate the current

            call xplasma_author_set(s,'xplasma_profs',iertmp)
            do ii=1,inxi
               zxi(ii)=(ii-1)*1.0d0/(inxi-1)
            enddo

            call xplasma_create_integ(s,'__xpasma_profs_tmp_integ', &
                 zxi,idwk,iertmp, cache_enable=.FALSE.)
                 
            call xplasma_rho_zonint(s,idwk,'Itor',zii,ierr)
            if(ierr.ne.0) exit

            zcur=zii(inxi)

            call xplasma_remove_item(s,idwk,iertmp)

            call xplasma_author_clear(s,'xplasma_profs',iertmp)
         endif

         !  P & q

         if(present(id_qprof)) then
            idq=id_qprof
         else
            call xplasma_gen_q(s,'Q_EQDSK',2,idq,ierr)
            if(ierr.ne.0) exit
         endif

         if(present(id_pprof)) then
            idp=id_pprof
         else
            call xplasma_common_ids(s,ierr,id_p=idp)
            if(ierr.ne.0) exit
            if(idp.eq.0) then
               call xplasma_gen_p(s,'P_EQDSK',2,idp,ierr)
               if(ierr.ne.0) exit
            endif
         endif

         if((idp.eq.0).or.(idq.eq.0)) then
            call xplasma_errmsg_append(s, &
                 ' xplasma_profs: zero ids for Pressure or q profiles')
            ierr=9999
            exit
         endif

         !  Psi (id can be zero)

         idpsi=0
         if(present(id_psi_in)) then
            idpsi = id_psi_in
         endif

         !  output grid sizes (default is from R & Z coordinates' 1st grids).

         if(present(nh)) then
            inh=nh
         endif

         if(present(nv)) then
            inv=nv
         endif

         if(present(nbdy)) then
            inbdy=max(61,nbdy)
         else
            inbdy=inbdy_dflt
         endif

         !===============================
         sp => s
         call eqi_geqdsk(ilun,glabel,zRmin,zRmax,zZmin,zZmax,zcur, &
              idp,idq,idpsi,inh,inv,inbdy,ierr)
         !===============================
         exit
      enddo

      if(ierr.ne.0) call xplasma_errmsg_append(s,' error in xplasma_wr_geqdsk')
      if(ifile.eq.1) then
         if(ierr.ne.0) then
            close(unit=ilun,status='delete')
         else
            close(unit=ilun)
         endif
      endif

    end subroutine xplasma_wr_geqdsk

    !-------------------------------------------------------
    subroutine xplasma_rhopsi_gen(s,npsi,ierr, tol, psivals, rhovals)

      !  generate an evenly spaced Psi vector (poloidal flux, Wb/rad) that
      !  spans the plasma from axis to edge.  Psi=0 corresponds to the
      !  magnetic axis; Psibdy is taken to be a positive number; the sign
      !  is omitted.  (The sign is available separately-- it corresponds 
      !  to the direction of the toroidal plasma current, jphi_ccw,
      !  available as an optional argument in xplasma_global_info).

      !  also, find the corresponding rho = sqrt(Phi_tor/Phi_tor_bdy)
      !  values that map to the Psi values w/in tolerance (tol).

      !---------------
      !  required arguments:

      type (xplasma), pointer :: s
      integer, intent(in) :: npsi   ! number of Psi values wanted (min. 2)
      integer, intent(out) :: ierr  ! completion status, 0=OK

      !---------------
      !  optional arguments:
      real*8, intent(in), optional :: tol  ! mapping tolerance
      !  default: something close to real*8 machine precision

      real*8, dimension(:), intent(out), optional :: Psivals
      !  The evenly spaced Psi vector generated, Wb/rad, (1:npsi)

      real*8, dimension(:), intent(out), optional :: Rhovals
      !  Rho values satisfying 
      !     abs(Psi(Rhovals(i))-Psivals(i)) <= tol*[Psi at bdy]
      !  (1:npsi) -- sqrt(toroidal_flux/toroidal_flux_at_bdy)
      !  0 on axis, 1 at the edge.

      !-------------------------------------
      real*8 :: psimin,psimax
      real*8, dimension(:), allocatable :: zpsis,zrhos
      integer :: ipsi
      !-------------------------------------

      ierr=0
      if(npsi.lt.2) then
         ierr=9999
         call xplasma_errmsg_append(s, &
              ' ?xplasma_rhopsi_gen: argument error: npsi >= 2 required.')
         return
      endif

      if(present(psivals)) then
         psivals=0
         if(size(psivals).lt.npsi) then
            ierr=612
            call xplasma_errmsg_append(s, &
                ' ?xplasma_rhopsi_gen: "psivals" array size >= npsi required.')
         endif
      endif

      if(present(rhovals)) then
         rhovals=0
         if(size(rhovals).lt.npsi) then
            ierr=612
            call xplasma_errmsg_append(s, &
                ' ?xplasma_rhopsi_gen: "rhovals" array size >= npsi required.')
         endif
      endif

      if(ierr.ne.0) return

      !-----------------------
      !  error checks passed

      call xplasma_psi_range(s,psimin,psimax)

      allocate(zpsis(npsi),zrhos(npsi))

      zpsis(1)=psimin
      zpsis(npsi)=psimax
      do ipsi=2,npsi-1
         zpsis(ipsi)=((npsi-ipsi)*psimin + (ipsi-1)*psimax)/(npsi-1)
      enddo

      if(present(rhovals)) then
         call xplasma_rhopsi_find(s,zpsis,zrhos,ierr, tol=tol)
      endif

      if(ierr.eq.0) then
         if(present(psivals)) psivals(1:npsi)=zpsis(1:npsi)
         if(present(rhovals)) rhovals(1:npsi)=zrhos(1:npsi)
      endif

      deallocate(zpsis,zrhos)

    end subroutine xplasma_rhopsi_gen

    !-------------------------------------------------------
    subroutine xplasma_rhopsi_find(s,psivals,rhovals,ierr, tol, iwarn)

      !  find rho values corresponding to a specified set of Psi values
      !    Psi -- Poloidal flux, Wb/rad
      !    rho -- sqrt(Tor_flux/Tor_flux_at_bdy)

      !  all input Psi values should be in the range [Psimin,Psimax] where
      !  Psimin corresponds to the magnetic axis and Psimax-Psimin
      !  corresponds to the (unsigned) poloidal flux, Wb/rad, enclosed 
      !  within the core plasma.  Usually Psimin=0 is set.

      type (xplasma), pointer :: s
      real*8, dimension(:), intent(in) :: Psivals  ! Psi values in any order
      real*8, dimension(:), intent(out) :: rhovals ! rho values output
      !  sizes of Psivals and rhovals must match

      integer, intent(out) :: ierr   ! status code, =0 on normal exit
      !  error occurs if xplasma is unitialized or contains no MHD equilibrium;
      !  Psi-out-of-range is handled (see iwarn, below).

      real*8, intent(in), optional :: tol     ! accuracy tolerance
      !  on output, rho values satisfy
      !     abs(Psi(rhovals(i))-Psivals(i)) <= tol*[Psi at bdy]
      !  (1:npsi) -- sqrt(toroidal_flux/toroidal_flux_at_bdy)
      !  0 on axis, 1 at the edge.

      integer, intent(out), optional :: iwarn  ! #of Psi values out of range

      !  if Psi <= Psimin, rho=0 is returned; if Psi >= Psimax rho=1 is
      !  returned.

      !------------------------------------
      real*8 :: psimin,psimax,rhomin,rhomax,ztol,psitol
      real*8, dimension(:), allocatable :: zpsi,zrho,zdum,zmina,zmaxa
      integer :: ipsi,inpsi,jwarn,insize
      logical, dimension(:), allocatable :: iok

      external eqi_xpsi_fun
      !------------------------------------

      ierr=0
      if(size(psivals).ne.size(rhovals)) then
         ierr=612
         call xplasma_errmsg_append(s, &
              ' ?xplasma_rhopsi_find: "psivals" and "rhovals" array argument sizes inconsistent.')
         return
      else
         insize=size(psivals)
      endif

      rhovals = 0
      if(present(iwarn)) iwarn=insize

      rhomin = 0
      rhomax = 1

      !-----------------
      !  get bdy & search error tolerance

      if(present(tol)) then
         ztol=tol
      else
         call xplasma_global_info(s,ierr, bdytol=ztol)
         if(ierr.ne.0) return
      endif

      !-----------------
      !  get available range of Psi values
      
      call xplasma_psi_range(s,psimin,psimax)

      psitol = ztol*max(abs(psimin),abs(psimax))

      !-----------------
      !  get the list of interior values for which search is required

      inpsi=0
      allocate(zpsi(insize),zrho(insize),zdum(insize),iok(insize))
      iok = .FALSE.

      jwarn=0

      do ipsi=1,insize
         if(psivals(ipsi).lt.(psimin+psitol)) then
            rhovals(ipsi)=rhomin
            if(psivals(ipsi).lt.(psimin-psitol)) jwarn = jwarn + 1

         else if(psivals(ipsi).gt.(psimax-psitol)) then
            rhovals(ipsi)=rhomax
            if(psivals(ipsi).gt.(psimax+psitol)) jwarn = jwarn + 1

         else
            inpsi=inpsi+1
            zpsi(inpsi)=psivals(ipsi)
            zrho(inpsi)=(psivals(ipsi)-psimin)/(psimax-psimin) ! crude guess
         endif
      enddo

      !  root finder...

      if(inpsi.gt.0) then

         sp => s

         allocate(zmina(inpsi),zmaxa(inpsi)); zmina=rhomin; zmaxa=rhomax

         call zridderx(inpsi, iok, zmina, zmaxa, ztol, psitol, &
              eqi_xpsi_fun, zrho, ierr, inpsi, zpsi, 1, zdum, 1)

         deallocate(zmina,zmaxa)

         if(ierr.ne.0) then

            ierr=9999
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_rhopsi_find: root finder failure.')
            jwarn=insize

         else

            inpsi = 0
            do ipsi=1,insize
               if(psivals(ipsi).lt.(psimin+psitol)) then
                  continue

               else if(psivals(ipsi).gt.(psimax-psitol)) then
                  continue

               else
                  inpsi = inpsi + 1
                  rhovals(ipsi) = zrho(inpsi)
               endif
            enddo

         endif

      endif
      if(present(iwarn)) iwarn = jwarn

      deallocate(zpsi,zrho,zdum,iok)

    end subroutine xplasma_rhopsi_find

end module xplasma_profs
