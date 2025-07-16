module xplasma_flxint

  !  methods for numerical integration
  !
  !  the following produce 1d results:
  !
  !  xplasma_rho_zonint -- f(rho) as a result of integration or average
  !     given a known integration id name, e.g. "Vol" or "dVol" or "<R^2>S".
  !
  !  xplasma_2d_zonint -- f(theta,rho) as a result of integration or average
  !     given a known integration id name, e.g. "dVol" or "<R^2>S".

  !  Added dmc 13 Dec 2006 -- fetch dVol or dArea elements for a rho grid
  !    or... for a 2d theta x rho grid.
  
  use xplasma_obj
  use eqi_intg_module
  implicit NONE

  private

  public :: xplasma_rho_zonint,xplasma_2d_zonint,xplasma_integ_eligible

  public :: xplasma_rhogrid_dva
  public :: xplasma_rhoth_grid_dva

  real*8, parameter :: czero = 0.0d0

  !---------------------------
  contains

    subroutine xplasma_rho_zonint(s,idi,int_name,result,ier, &
         dvol_out,iwarn_dvol,dvdrho_out,iwarn_dvdrho)

      ! integrator -- over rho zones result of form f(rho) indexed by zone
      !                      1 result per zone (#zones)
      !               or rho surfaces, also result of form f(rho)
      !                      1 result per surface (#zones + 1)

      type (xplasma), pointer :: s
      integer, intent(in) :: idi ! integrator data structure ("black box") id
      character*(*),intent(in) :: int_name ! name of desired integration
      real*8, dimension(:), intent(out) :: result
      integer, intent(out) :: ier  ! status code 0=OK

      !  ** optional outputs **

      real*8, dimension(:), optional :: dvol_out   ! zone volumes
      integer, intent(out), optional :: iwarn_dvol ! =0 if zone volumes OK
      ! iwarn_dvol=1: zone volumes not available; =2: array dimension error.

      real*8, dimension(:), optional :: dvdrho_out ! dV/drho @ surfaces
      integer, intent(out), optional :: iwarn_dvdrho ! =0 if dV/drho OK
      ! iwarn_dvdrho=1: dV/drho not available; =2: array dimension error.

      !-----------------------------------
      real*8, dimension(:), allocatable :: rho_int,rho_brk
      real*8, dimension(:), allocatable :: dvol,dvdrho,bmaxa
      integer, dimension(:), allocatable :: indexb

      real*8 :: rhomin
      integer :: inzons,inzonb,insurfs,isizer,isizei,ilen,iertmp,i,j,ii,inum
      integer :: inrho_brk,infine,inumu
      logical :: dvol_avail,dvdrho_avail,dvol_need,dvdrho_need,ilfine
      logical :: bmaxa_avail,bmaxa_need

      character*32 zname,iiname
      character*128 zmsg
      real*8, parameter :: c2pi = 6.2831853071795862D+00
      real*8, dimension(:), allocatable :: zwk
      real*8, dimension(:), allocatable :: rhofine,wrhof
      real*8, dimension(:), pointer :: thvec,wth,cache
      real*8 :: thzon1(2)
      integer :: cache_size
      logical :: cache_init(eqi_intg_nfi),ilsurf,ilocal
      logical :: iret_vol,iret_dvdrho
      !-----------------------------------
      !  info on the integrator dataset

      result=0

      thzon1(1)=0
      thzon1(2)=c2pi

      call xplasma_integ_info(s,idi,ier, &
           dvol_avail=dvol_avail, dvdrho_avail=dvdrho_avail, &
           bmaxa_avail=bmaxa_avail, &
           n_rho_int = insurfs, &
           n_rho_brk = inrho_brk)
      
      if(ier.ne.0) return

      inzons=insurfs-1
      inzonb=inrho_brk-1

      iret_vol=.FALSE.
      if(present(iwarn_dvol)) iwarn_dvol=1  ! for now...
      if(present(dvol_out)) then
         dvol_out=0
         if(size(dvol_out).eq.inzons) then
            iret_vol=.TRUE.
         else
            if(present(iwarn_dvol)) iwarn_dvol=2  ! array size error
         endif
      endif

      iret_dvdrho=.FALSE.
      if(present(iwarn_dvdrho)) iwarn_dvdrho=1  ! for now...
      if(present(dvdrho_out)) then
         dvdrho_out=0
         if(size(dvdrho_out).eq.insurfs) then
            iret_dvdrho=.TRUE.
         else
            if(present(iwarn_dvdrho)) iwarn_dvdrho=2  ! array size error
         endif
      endif

      allocate(rho_int(insurfs),rho_brk(inrho_brk), &
           dvol(inzons),dvdrho(insurfs),bmaxa(insurfs))

      do 

         ilocal = .FALSE.
         ilfine = .FALSE.  ! rho eval at surfaces only if(ilsurf)...

         call xplasma_integ_info(s,idi,ier, &
              rho_int=rho_int,rho_brk=rho_brk, &
              dvol=dvol, dvdrho=dvdrho, bmaxa=bmaxa, &
              thvec_ptr=thvec, thwts_ptr=wth, cache_ptr=cache, &
              rhomin=rhomin,cache_init=cache_init)
      
         if(ier.ne.0) exit

         if(dvol_avail.and.iret_vol) then
            dvol_out = dvol
            if(present(iwarn_dvol)) iwarn_dvol=0
         endif

         if(dvdrho_avail.and.iret_dvdrho) then
            dvdrho_out = dvdrho
            if(present(iwarn_dvol)) iwarn_dvol=0
         endif

         zname = int_name
         call uupper(zname)
         ilen = len(trim(zname))

         !-----------------------
         sp => s   ! *** so integrand routines can see s *** cf eqi_intg_module
         !-----------------------
         call xplasma_integ_num(zname,inum,ilsurf, &
              dvol_need,dvdrho_need,bmaxa_need)

         if(inum.eq.0) then
            call xplasma_errmsg_append(s, &
                 '  unrecognized integration name:  "'//trim(int_name)//'".')
            ier=611
            exit
         endif

         if(.not.ilsurf) then
            nullify(cache)   ! retained cache for surface integrals only...
         endif

         if(.not.associated(cache)) then
            if(ilsurf) then
               cache_size=size(thvec)*insurfs*eqi_intg_nfi
            else
               ilfine=.TRUE.  ! rho eval at many points inside each rho zone
               infine=10*inzonb
               allocate(rhofine(infine),wrhof(infine))
               call xplasma_integ_mkvec(inrho_brk,rho_brk,rhofine,wrhof)
               cache_size=size(thvec)*infine*eqi_intg_nfi
               allocate(indexb(infine))
               call xplasma_integ_mkindex(insurfs,rho_int, &
                    infine,rhofine,indexb)
            endif
            allocate(cache(cache_size))
            ilocal = .TRUE.   ! flags locally allocated CACHE
            cache_init = .FALSE.
         endif

         isizer=size(result)
         if(ilsurf) then
            isizei=insurfs
         else
            isizei=inzons
         endif

         !  check expected size of result against size of array provided...

         if((inum.eq.eqi_intg_dvol).or.(inum.eq.eqi_intg_darea)) then
            if(isizer.ne.inzons) then
               ier=610
               isizei=inzons
            endif
         else if(isizer.ne.isizei) then
            ier=610
         endif

         if(ier.ne.0) then
            zmsg=' '
            call xplasma_get_item_info(s,idi,iertmp, name=iiname)
            write(zmsg,*) '  integrator dataset name: ',trim(iiname)
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) '  expected size of result of "',trim(zname), &
                 '" integral/avg is: ',isizei
            call xplasma_errmsg_append(s,zmsg)
            zmsg=' '
            write(zmsg,*) '  size of result array passed is: ',isizer
            exit
         endif

         if(ilfine) then
            allocate(zwk(infine))
            deallocate(bmaxa,dvdrho)
            allocate(bmaxa(infine)); bmaxa=czero
            allocate(dvdrho(infine))
            !  these to be recalculated on the finer scale:
            dvdrho_avail=.FALSE.
         else
            allocate(zwk(insurfs))
         endif

         !  pre-fetch volume if needed...

         if(ilfine.and.dvol_need) then
            call xplasma_integ_work(eqi_intg_dvdrho,size(thvec),thvec,wth, &
                 2,thzon1,infine,rhofine,rhomin, &
                 cache_init,cache,zwk,bmaxa,dvdrho,ier)  ! zwk not used
            if(ier.ne.0) exit
         endif

         if(dvol_need.and.(.not.dvol_avail)) then

            if(ilfine) then
               dvol=czero
               do i=1,infine
                  ii = indexb(i)
                  dvol(ii)=dvol(ii)+dvdrho(i)*wrhof(i)
               enddo
            else
               !
               call xplasma_integ_work(eqi_intg_vol,size(thvec),thvec,wth, &
                    2,thzon1,insurfs,rho_int,rhomin, &
                    cache_init,cache,dvdrho,bmaxa,zwk,ier) ! dvdrho not used
               dvol(1:inzons)=zwk(2:insurfs)-zwk(1:inzons)
            endif

            if(ier.eq.0) then
               call xplasma_integ_dva(s,idi,iertmp, dvol=dvol)
               if(iret_vol) then
                  dvol_out = dvol
                  if(present(iwarn_dvol)) iwarn_dvol=0
               endif
            endif

         endif
         if(ier.ne.0) exit

         !  pre-fetch dV/drho if needed...
         !  note-- if no cache with integrator, then the dVdrho integral
         !  must be evaluated in order to fill the local cache with local
         !  dV/drho variation information-- hence ".or.ilocal".

         if(dvdrho_need.and.((.not.dvdrho_avail).or.ilocal)) then

            if(ilfine) then
               ier=9999
               call xplasma_errmsg_append(s, &
                    ' ? xplasma_flxint internal error: ilfine.and.dvdrho_need')
               exit
            endif

            call xplasma_integ_work(eqi_intg_dvdrho,size(thvec),thvec,wth, &
                 2,thzon1,insurfs,rho_int,rhomin, &
                 cache_init,cache,zwk,bmaxa,dvdrho,ier)  ! zwk not used..

            if((ier.eq.0).and.(.not.dvdrho_avail)) then
               call xplasma_integ_dva(s,idi,iertmp, dvdrho=dvdrho)
               if(iret_dvdrho) then
                  dvdrho_out = dvdrho
                  if(present(iwarn_dvdrho)) iwarn_dvdrho=0
               endif
            endif

         endif
         if(ier.ne.0) exit

         !  now evaluate the requested integral or average...

         if(inum.eq.eqi_intg_vol) then

            result(1)=0
            do i=2,insurfs
               result(i)=result(i-1)+dvol(i-1)
            enddo

         else if(inum.eq.eqi_intg_dvol) then

            result = dvol

         else if(inum.eq.eqi_intg_dvdrho) then

            result = dvdrho

         else if(ilsurf)then

            if(ilfine) then
               ier=9999
               call xplasma_errmsg_append(s, &
                    ' ? xplasma_flxint internal error: ilfine.and.ilsurf')
               exit
            endif

            inumu=inum
            if(inumu.eq.eqi_intg_darea) inumu=eqi_intg_area ! use Gauss form.

            call xplasma_integ_work(inumu,size(thvec),thvec,wth, &
                 2,thzon1,insurfs,rho_int,rhomin, &
                 cache_init,cache,dvdrho,bmaxa,zwk,ier)

            if(inum.eq.eqi_intg_darea) then
               result(1:inzons)=zwk(2:insurfs)-zwk(1:inzons)
            else
               result(1:insurfs) = zwk(1:insurfs)
            endif

         else

            if(.not.ilfine) then
               ier=9999
               call xplasma_errmsg_append(s, &
                    ' ? xplasma_flxint internal error: .not.(ilfine.or.ilsurf)')
               exit
            endif

            call xplasma_integ_work(inum,size(thvec),thvec,wth, &
                 2,thzon1,infine,rhofine,rhomin, &
                 cache_init,cache,dvdrho,bmaxa,zwk,ier)

            result=czero
            do i=1,infine
               ii=indexb(i)
               if(dvol_need) then
                  result(ii)=result(ii)+dvdrho(i)*zwk(i)*wrhof(i)
               else
                  result(ii)=result(ii)+zwk(i)*wrhof(i)
               endif
            enddo
            if(dvol_need) then
               do i=1,inzons
                  result(i)=result(i)/dvol(i)
               enddo
            endif

         endif
         exit
      enddo

      if(.not.ilfine) then
         if((ier.eq.0).and.bmaxa_need.and.(.not.bmaxa_avail)) then
            ! Bmax(rho) was computed
            call xplasma_integ_dva(s,idi,iertmp, bmaxa=bmaxa)
         endif
      endif

      deallocate(rho_brk,rho_int,dvol,dvdrho,bmaxa)
      if(allocated(zwk)) deallocate(zwk)
      if(ilocal) then
         deallocate(cache)
         if(ilfine) deallocate(rhofine,wrhof,indexb)
      else
         !  record cache markers (saved with integrator dataset).
         if(ier.eq.0) call xplasma_integ_scache(s,idi,cache_init,ier)
      endif

      nullify(thvec,wth,cache)

    end subroutine xplasma_rho_zonint

    subroutine xplasma_2d_zonint(s,idi,int_name,result,ier, &
         dvol_out,iwarn_dvol,dvdrho_out,iwarn_dvdrho)

      ! integrator -- 2d domain; return integrations or averages
      !   over theta^rho zones or across theta zones at rho surfaces.

      type (xplasma), pointer :: s
      integer, intent(in) :: idi ! integrator data structure ("black box") id
      character*(*),intent(in) :: int_name ! name of desired integration
      real*8, dimension(:,:), intent(out) :: result
      integer, intent(out) :: ier  ! status code 0=OK

      ! result is dimension [#theta-zones,#rho-zones] or
      !        [#theta-zones,#rho-surfaces] according as the integrand
      !        specified is zone oriented or surface oriented.
      !          The #zones are defined when the integrator dataset (idi)
      !          is defined.

      !  ** optional outputs **

      real*8, dimension(:,:), optional :: dvol_out   ! zone volumes
      integer, intent(out), optional :: iwarn_dvol   ! =0 if zone volumes OK
      ! iwarn_dvol=1: zone volumes not available; =2: array dimension error.

      real*8, dimension(:,:), optional :: dvdrho_out ! dV/drho @ surfaces
      integer, intent(out), optional :: iwarn_dvdrho ! =0 if dV/drho OK
      ! iwarn_dvdrho=1: dV/drho not available; =2: array dimension error.

      ! Options for integrands is the same as for xplasma_rho_zonint;
      ! Surface oriented results e.g. LPOL are on rho surfaces at each theta 
      ! zone;
      ! Zone orented results e.g. dVol are on rho^theta zones.

      !-----------------------------------
      real*8, dimension(:), allocatable :: rho_int,rho_brk,th_int
      real*8, dimension(:,:), allocatable :: dvol_2d,dvdrho_2d
      real*8, dimension(:), allocatable :: bmaxa
      real*8 :: rhomin

      integer :: inzons,insurfs,isizer,isizei,ilen,iertmp,i,j,ii,inum
      integer :: inzonb,inthzons,inthbdys,isizeth,ith
      integer :: inrho_brk,infine
      logical :: dvol_need,dvdrho_need,ilfine
      logical :: bmaxa_need

      character*32 zname,iiname
      character*128 zmsg
      real*8, parameter :: c2pi = 6.2831853071795862D+00
      real*8, dimension(:,:), allocatable :: zwk
      real*8, dimension(:), allocatable :: rhofine,wrhof
      real*8, dimension(:), pointer :: thvec,wth,cache
      integer, dimension(:), allocatable :: indexb
      integer :: cache_size
      logical :: cache_init(eqi_intg_nfi),ilsurf,ilocal
      logical :: iret_vol,iret_dvdrho
      !-----------------------------------
      !  info on the integrator dataset

      result=0

      call xplasma_integ_info(s,idi,ier, &
           n_rho_int = insurfs, &
           n_rho_brk = inrho_brk, &
           n_th_int = inthbdys)
      
      if(ier.ne.0) return

      inzons=insurfs-1
      inzonb=inrho_brk-1
      inthzons=inthbdys-1

      iret_vol=.FALSE.
      if(present(iwarn_dvol)) iwarn_dvol=1  ! for now...
      if(present(dvol_out)) then
         dvol_out=0
         if((size(dvol_out,1).eq.inthzons).and. &
              (size(dvol_out,2).eq.inzons)) then
            iret_vol=.TRUE.
         else
            if(present(iwarn_dvol)) iwarn_dvol=2  ! array size error
         endif
      endif

      iret_dvdrho=.FALSE.
      if(present(iwarn_dvdrho)) iwarn_dvdrho=1  ! for now...
      if(present(dvdrho_out)) then
         dvdrho_out=0
         if((size(dvdrho_out,1).eq.inthzons).and. &
              (size(dvdrho_out,2).eq.insurfs)) then
            iret_dvdrho=.TRUE.
         else
            if(present(iwarn_dvdrho)) iwarn_dvdrho=2  ! array size error
         endif
      endif

      allocate(rho_int(insurfs),rho_brk(inrho_brk),th_int(inthbdys), &
           dvol_2d(inthzons,inzons),dvdrho_2d(inthzons,insurfs),bmaxa(insurfs))
      dvol_2d=0; dvdrho_2d=0; bmaxa=0

      do 

         ! note local cache only for 2d integration...
         ilocal = .FALSE.   ! set .TRUE. after allocation

         call xplasma_integ_info(s,idi,ier, &
              rho_int=rho_int,rho_brk=rho_brk,th_int=th_int, &
              thvec_ptr=thvec, thwts_ptr=wth, rhomin=rhomin)
      
         if(ier.ne.0) exit

         zname = int_name
         call uupper(zname)
         ilen = len(trim(zname))

         !-----------------------
         sp => s   ! *** so integrand routines can see s *** cf eqi_intg_module
         !-----------------------

         call xplasma_integ_num(zname,inum,ilsurf, &
              dvol_need,dvdrho_need,bmaxa_need)

         if(inum.eq.0) then
            call xplasma_errmsg_append(s, &
                 '  unrecognized integration name:  "'//trim(int_name)//'".')
            ier=611
            exit
         endif

         if((inum.eq.eqi_intg_vol).or.(inum.eq.eqi_intg_area)) then
            call xplasma_errmsg_append(s,' error in xplasma_2d_zonint:')
            call xplasma_errmsg_append(s, &
                 '  integration name invalid for 2d integration:  "'//trim(int_name)//'".')
            ier=611
            exit
         endif

         nullify(cache)

         if(ilsurf) then
            ilfine=.FALSE.
            infine=0
            cache_size=size(thvec)*insurfs*eqi_intg_nfi
         else
            ilfine=.TRUE.
            infine=10*inzonb
            cache_size=size(thvec)*infine*eqi_intg_nfi
         endif

         allocate(cache(cache_size))
         ilocal = .TRUE.   ! flags locally allocated CACHE
         cache_init = .FALSE.

         if((.not.ilsurf).or.(inum.eq.eqi_intg_darea).or.(dvol_need)) then
            infine=10*inzonb
            allocate(rhofine(infine),wrhof(infine))
            call xplasma_integ_mkvec(inrho_brk,rho_brk,rhofine,wrhof)
            allocate(indexb(infine))
            call xplasma_integ_mkindex(insurfs,rho_int, &
                 infine,rhofine,indexb)
         endif

         isizeth=size(result,1)
         isizer=size(result,2)
         if(ilsurf) then
            isizei=insurfs
         else
            isizei=inzons
         endif

         !  check expected size of result against size of array provided...

         if((inum.eq.eqi_intg_dvol).or.(inum.eq.eqi_intg_darea)) then
            if(isizer.ne.inzons) then
               ier=610
               isizei=inzons
            endif
         else if(isizer.ne.isizei) then
            ier=610
         endif

         if(isizeth.ne.inthzons) ier=610

         if(ier.ne.0) then

            zmsg=' '
            call xplasma_get_item_info(s,idi,iertmp, name=iiname)
            write(zmsg,*) '  integrator dataset name: ',trim(iiname)
            call xplasma_errmsg_append(s,zmsg)

            call xplasma_errmsg_append(s, &
                 '  integrand name: "'//trim(int_name)//'".')

            zmsg=' '
            write(zmsg,*) '  expected size of 2d result (dim1): ',inthzons
            call xplasma_errmsg_append(s,zmsg)

            zmsg=' '
            write(zmsg,*) '  result(dim1) size is: ',isizeth
            call xplasma_errmsg_append(s,zmsg)
            
            zmsg=' '
            write(zmsg,*) '  expected size of 2d result (dim2): ',isizei
            call xplasma_errmsg_append(s,zmsg)

            zmsg=' '
            write(zmsg,*) '  size of result(dim2) is: ',isizer

            exit
         endif

         if(ilfine) then
            allocate(zwk(inthzons,infine))
            deallocate(bmaxa,dvdrho_2d)
            allocate(bmaxa(infine)); bmaxa=czero
            allocate(dvdrho_2d(inthzons,infine))
         else
            allocate(zwk(inthzons,insurfs))
         endif

         if(ilfine.and.dvol_need) then
            call xplasma_integ_work(eqi_intg_dvdrho,size(thvec),thvec,wth, &
                 inthbdys,th_int,infine,rhofine,rhomin, &
                 cache_init,cache,zwk,bmaxa,dvdrho_2d,ier)  ! zwk not used
            if(ier.ne.0) exit
         endif

         if(dvol_need) then

            if(ilfine) then
               dvol_2d=czero
               do i=1,infine
                  ii = indexb(i)
                  do ith=1,inthzons
                     dvol_2d(ith,ii)=dvol_2d(ith,ii)+dvdrho_2d(ith,i)*wrhof(i)
                  enddo
               enddo
            else
               ! dvdrho_2d not used...:
               call xplasma_integ_work(eqi_intg_vol,size(thvec),thvec,wth, &
                    inthbdys,th_int,insurfs,rho_int,rhomin, &
                    cache_init,cache,dvdrho_2d,bmaxa,zwk,ier)
               if(ier.eq.0) then
                  call xplasma_integ_gwork(eqi_intg_dvol,inthbdys,th_int, &
                       insurfs,rho_int,infine,rhofine,wrhof,indexb, &
                       zwk,dvol_2d,ier)
               endif
            endif

            if(iret_vol) then
               dvol_out=dvol_2d
               if(present(iwarn_dvol)) iwarn_dvol=0
            endif

         endif
         if(ier.ne.0) exit

         !  pre-fetch dV/drho if needed...
         !  note-- if no cache with integrator, then the dVdrho integral
         !  must be evaluated in order to fill the local cache with local

         if(dvdrho_need) then

            if(ilfine) then
               ier=9999
               call xplasma_errmsg_append(s, &
                    ' ? xplasma_flxint internal error: ilfine.and.dvdrho_need')
               exit
            endif

            call xplasma_integ_work(eqi_intg_dvdrho,size(thvec),thvec,wth, &
                 inthbdys,th_int,insurfs,rho_int,rhomin, &
                 cache_init,cache,zwk,bmaxa,dvdrho_2d,ier)  ! zwk not used..

            if(iret_dvdrho) then
               dvdrho_out=dvdrho_2d
               if(present(iwarn_dvdrho)) iwarn_dvdrho=0
            endif

         endif
         if(ier.ne.0) exit

         !  now evaluate the requested integral or average...

         if(inum.eq.eqi_intg_dvol) then

            result = dvol_2d

         else if(inum.eq.eqi_intg_dvdrho) then

            result = dvdrho_2d

         else if(ilsurf) then

            if(ilfine) then
               ier=9999
               call xplasma_errmsg_append(s, &
                    ' ? xplasma_flxint internal error: ilfine.and.ilsurf')
               exit
            endif

            if(inum.eq.eqi_intg_darea) then
               !  compute rho=const. Gaussian segments...
               call xplasma_integ_work(eqi_intg_area,size(thvec),thvec,wth, &
                    inthbdys,th_int,insurfs,rho_int,rhomin, &
                    cache_init,cache,dvdrho_2d,bmaxa,zwk,ier)
               if(ier.ne.0) exit
               !  compute theta=const. Gaussian segments and combine...
               call xplasma_integ_gwork(eqi_intg_darea,inthbdys,th_int, &
                    insurfs,rho_int,infine,rhofine,wrhof,indexb, &
                    zwk,result,ier)
               if(ier.ne.0) exit

            else
               call xplasma_integ_work(inum,size(thvec),thvec,wth, &
                    inthbdys,th_int,insurfs,rho_int,rhomin, &
                    cache_init,cache,dvdrho_2d,bmaxa,zwk,ier)
               if(ier.ne.0) exit
               result(1:inthzons,1:insurfs) = zwk(1:inthzons,1:insurfs)

            endif

         else

            if(.not.ilfine) then
               ier=9999
               call xplasma_errmsg_append(s, &
                    ' ? xplasma_flxint internal error: .not.(ilfine.or.ilsurf)')
               exit
            endif

            call xplasma_integ_work(inum,size(thvec),thvec,wth, &
                 inthbdys,th_int,infine,rhofine,rhomin, &
                 cache_init,cache,dvdrho_2d,bmaxa,zwk,ier)

            result=czero
            do i=1,infine
               ii=indexb(i)
               if(dvol_need) then
                  do ith=1,inthzons
                     result(ith,ii)=result(ith,ii)+ &
                          dvdrho_2d(ith,i)*zwk(ith,i)*wrhof(i)
                  enddo
               else
                  do ith=1,inthzons
                     result(ith,ii)=result(ith,ii)+zwk(ith,i)*wrhof(i)
                  enddo
               endif
            enddo
            if(dvol_need) then
               do i=1,inzons
                  do ith=1,inthzons
                     result(ith,i)=result(ith,i)/dvol_2d(ith,i)
                  enddo
               enddo
            endif

         endif
         exit
      enddo

      deallocate(rho_brk,rho_int,th_int,dvol_2d,dvdrho_2d,bmaxa)
      if(allocated(zwk)) deallocate(zwk)

      if(ilocal) then
         deallocate(cache)
      else
         !  record cache markers (saved with integrator dataset).
         if(ier.eq.0) call xplasma_integ_scache(s,idi,cache_init,ier)
      endif

      if(infine.gt.0) then
         deallocate(rhofine,wrhof,indexb)
      endif

      nullify(thvec,wth,cache)

    end subroutine xplasma_2d_zonint

    subroutine xplasma_integ_eligible(sp,integ_name,id,ier)

      ! determine if a named function (integ_name) is eligible for 
      ! integration.  At present, this means a function f(rho,theta)
      ! or f(theta,rho).  Eligibility may be extended beyond this at
      ! some point in the future, if the need arises.

      type (xplasma), pointer :: sp
      character*(*), intent(in) :: integ_name  ! name of profile
      integer, intent(out) :: id               ! id of profile if eligible
      integer, intent(out) :: ier              ! status code

      ! if the function is not eligible, id=0 and ier=611 are returned.
      !------------------------------

      character*32 iname

      integer :: iertmp,itype,idx1,idx2,idx3,idcf,icoord,id_R,idcR
      logical :: iok
      !------------------------------

      iname = integ_name
      call uupper(iname)

      call xplasma_find_item(sp,iname,id,iertmp)

      if(id.gt.0) then
         call xplasma_get_item_info(sp,id,iertmp, itype=itype)

         if(itype.eq.xplasma_profType) then
            ! is is a profile-- see if dependents coordinates are (rho,th)

            call xplasma_prof_info(sp,id,iertmp, & 
                 gridId1=idx1, gridId2=idx2, gridId3=idx3, counter=idcf)

            iok=.TRUE.
            call xplasma_grid_info(sp,idx1,iertmp, coord=icoord)
            if(.not.((icoord.eq.xplasma_rho_coord).or. &
                 (icoord.eq.xplasma_theta_coord))) iok=.FALSE.

            call xplasma_grid_info(sp,idx2,iertmp, coord=icoord)
            if(.not.((icoord.eq.xplasma_rho_coord).or. &
                 (icoord.eq.xplasma_theta_coord))) iok=.FALSE.

            if(idx3.gt.0) iok=.FALSE.

            if(.not.iok) then
               call xplasma_errmsg_append(sp, &
                    ' cannot flux surface average profile <'// &
                    iname//'> -- not defined vs. (rho,theta)')
               id = 0
            else
               call xplasma_common_ids(sp,iertmp, id_R=id_R)
               call xplasma_prof_info(sp,id_R,iertmp, counter=idcR)

               if(idcR.gt.idcf) then
                  call xplasma_errmsg_append(sp, &
                       ' cannot flux surface average profile <'// &
                       iname// &
                       '> -- stale, not updated since equilibrium change.')

                  id = 0

               else

                  continue

               endif
            endif

         else
            call xplasma_errmsg_append(sp,'  '//trim(iname)//' not a profile.')
            id=0
         endif
      endif

      if(id.eq.0) ier=611

    end subroutine xplasma_integ_eligible

    subroutine xplasma_rhogrid_dva(s,idrho,ierr, dvol, darea)

      !  retrieve volume elements associated with a grid

      type(xplasma), pointer :: s
      integer, intent(in) :: idrho                ! grid ID
      integer, intent(out) :: ierr                ! completion code, 0=OK

      real*8, dimension(:), intent(out), optional :: dvol   ! volume elements
      real*8, dimension(:), intent(out), optional :: darea  ! area elements
      !---------------------------------------
      character*32 gname,bbname,iname
      character*128 zmsg
      integer :: idbb,idi,inx,iz,inz,icoord,iertmp,idum(4),itype
      real*8, dimension(:), pointer :: r8buf
      real*8, dimension(:), allocatable :: rhogrid
      !---------------------------------------

      if(present(dvol)) dvol=0
      if(present(darea)) darea=0

      call xplasma_grid_info(s,idrho,ierr, size=inx, coord=icoord, name=gname)

      !  the idrho argument must point to a rho grid...
      if((ierr.ne.0).or.(icoord.ne.xplasma_rho_coord)) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_rhogrid_dva: subroutine argument idrho: not a rho grid.')
         ierr=9999
         return
      endif

      !  the number of zones in the rhogrid (size of grid - 1) must match
      !  the passed array size

      inz=inx-1
      if(present(dvol)) then
         if(size(dvol).ne.inz) then
            zmsg=' '
            write(zmsg,'(" ?xplasma_rhogrid_dva: mismatch: #zones in grid = ",i4,"; dVol array size = ",i4)') &
                 inz,size(dvol)
            call xplasma_errmsg_append(s,zmsg)
            ierr=9999
            return
         endif
      endif

      if(present(darea)) then
         if(size(darea).ne.inz) then
            zmsg=' '
            write(zmsg,'(" ?xplasma_rhogrid_dva: mismatch: #zones in grid = ",i4,"; darea array size = ",i4)') &
                 inz,size(darea)
            call xplasma_errmsg_append(s,zmsg)
            ierr=9999
            return
         endif
      endif


      !  OK:  look up blackbox.  If not found, compute it.

      bbname = '__DVA_'//trim(gname)
      call xplasma_blkbxId(s,bbname,idbb)

      if(idbb.eq.0) then

         !  compute the DVA structure

         iname = '__TMP_DVA'

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         allocate(rhogrid(inx))
         call xplasma_grid(s,idrho,rhogrid,iertmp)
         call xplasma_create_integ(s,iname,rhogrid,idi,ierr, &
              cache_enable=.TRUE.)
         deallocate(rhogrid)

         if(ierr.ne.0) then
            call xplasma_errmsg_append(s, &
                 '?integrator creation failed in xplasma_rhogrid_dva')
            call xplasma_author_clear(s,xplasma_xmhd,iertmp)
            return
         endif

         idum=0
         itype=15  ! arb. type code

         allocate(r8buf(2*inz))

         call xplasma_rho_zonint(s,idi,'dVol',r8buf(1:inz),iertmp)
         call xplasma_rho_zonint(s,idi,'dArea',r8buf(inz+1:2*inz),ierr)

         ierr=max(ierr,iertmp)
         if(ierr.ne.0) then
            call xplasma_errmsg_append(s, &
                 '?integrator execution failed in xplasma_rhogrid_dva')
         endif

         if(ierr.eq.0) then
            call xplasma_create_blackBox(s,bbname,itype,idbb,ierr, &
                 iarray=idum,r8array=r8buf, &
                 label='dVol and dArea for '//trim(gname),units='mixed MKS')
         endif

         ! cleanup

         deallocate(r8buf)
         call xplasma_remove_item(s,idi,iertmp)
         call xplasma_author_clear(s,xplasma_xmhd,iertmp)

         if(ierr.ne.0) then
            call xplasma_errmsg_append(s,'error in xplasma_rhogrid_dva')
            return
         endif

      endif

      call xplasma_blackBox_retrieve(s,idbb,iertmp, r8a_ptr = r8buf)

      do iz=1,inz
         if(present(dvol)) dvol(iz)=r8buf(iz)
         if(present(darea)) darea(iz)=r8buf(inz+iz)
      enddo

    end subroutine xplasma_rhogrid_dva

    subroutine xplasma_rhoth_grid_dva(s,idg1,idg2,ierr, ccwflag, dvol, darea)

      !  retrieve volume elements associated with combined (rho,th) grids
      !    idg1 can be rho or th; idg2 can be th or rho
      !    to reverse th ordering in output, set ccwflag=.FALSE.

      type(xplasma), pointer :: s
      integer, intent(in) :: idg1,idg2            ! grid IDs
      integer, intent(out) :: ierr                ! completion code, 0=OK

      logical, intent(in), optional :: ccwflag    ! .TRUE. means 
      ! theta is oriented counter-clockwise, as in xplasma2 internal variables;
      ! .FALSE. means clockwise orientation; th index reversal required.

      real*8, dimension(:,:), intent(out), optional :: dvol   ! volume elements
      real*8, dimension(:,:), intent(out), optional :: darea  ! area elements

      !------------------------------------------------------------
      character*32 gname1,gname2,rhoname,thname,bbname,iname,iname2
      character*128 zmsg
      integer :: idbb,idi,idi2,inx1,inx2,icoord1,icoord2,idum(4),itype
      integer :: idrho,idth,inrho,inth,ith,irho,ibuf,ithx
      integer :: ier1,ier2,iertmp
      real*8, dimension(:), pointer :: r8buf
      real*8, dimension(:), allocatable :: rhogrid,thgrid
      real*8, dimension(:,:), allocatable :: zdva
      logical :: iccwflag,rhofirst
      !------------------------------------------------------------

      if(present(dVol)) dVol=0
      if(present(dArea)) dArea=0

      iccwflag = .TRUE.
      if(present(ccwflag)) iccwflag = ccwflag

      call xplasma_grid_info(s,idg1,ier1, size=inx1, coord=icoord1, &
           name=gname1)

      call xplasma_grid_info(s,idg2,ier2, size=inx2, coord=icoord2, &
           name=gname2)

      ierr=max(ier1,ier2)
      if(ierr.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_rhoth_grid_dva: invalid grid IDs received.')
         return
      endif

      if((icoord1.eq.icoord2).or. &
           ((icoord1.ne.xplasma_rho_coord).and.(icoord1.ne.xplasma_theta_coord)).or. &
           ((icoord2.ne.xplasma_rho_coord).and.(icoord2.ne.xplasma_theta_coord))) then
         zmsg = ' ?xplasma_rhoth_grid_dva: "rho" and/or "theta" grid missing: '//trim(gname1)//' '//trim(gname2)
         call xplasma_errmsg_append(s,zmsg)
         ierr=9999
         return
      endif

      if(present(dVol)) call checksize(dVol,'dVol')
      if(present(dArea)) call checksize(dArea,'dArea')
      if(ierr.ne.0) return

      !  OK:  look up blackbox.  If not found, compute it.

      if(icoord1.eq.xplasma_rho_coord) then
         rhoname=gname1
         inrho=inx1
         idrho=idg1
         thname=gname2
         inth=inx2
         idth=idg2
         rhofirst = .TRUE.
      else
         rhoname=gname2
         inrho=inx2
         idrho=idg2
         thname=gname1
         inth=inx1
         idth=idg1
         rhofirst = .FALSE.
      endif

      bbname = '__DVA_'//trim(thname)//'_'//trim(rhoname)
      call xplasma_blkbxId(s,bbname,idbb)

      if(idbb.eq.0) then

         !  compute the DVA structure

         iname = '__TMP_DVA'
         iname2= '__TMP_DVA2'

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         allocate(rhogrid(inrho))
         call xplasma_grid(s,idrho,rhogrid,iertmp)

         allocate(thgrid(inth))
         call xplasma_grid(s,idth,thgrid,iertmp)

         call xplasma_create_integ(s,iname,rhogrid,idi,ierr, &
              cache_enable=.TRUE.)

         if(ierr.ne.0) then
            call xplasma_errmsg_append(s, &
                 '?integrator creation failed in xplasma_rhoth_grid_dva')
            call xplasma_author_clear(s,xplasma_xmhd,iertmp)
            deallocate(rhogrid,thgrid)
            return
         endif

         call xplasma_augment_integ(s,iname2,idi,thgrid,idi2,ierr)

         if(ierr.ne.0) then
            call xplasma_errmsg_append(s, &
                 '?integrator augmentation failed in xplasma_rhoth_grid_dva')
            call xplasma_author_clear(s,xplasma_xmhd,iertmp)
            deallocate(rhogrid,thgrid)
            return
         endif

         idum=0
         itype=14  ! arb. type code

         allocate(r8buf(2*(inrho-1)*(inth-1)),zdva(inth-1,inrho-1))

         ibuf=0
         call xplasma_2d_zonint(s,idi2,'dVol',zdva,iertmp)
         do irho=1,inrho-1
            do ith=1,inth-1
               ibuf=ibuf+1
               r8buf(ibuf)=zdva(ith,irho)
            enddo
         enddo

         call xplasma_2d_zonint(s,idi2,'dArea',zdva,ierr)
         do irho=1,inrho-1
            do ith=1,inth-1
               ibuf=ibuf+1
               r8buf(ibuf)=zdva(ith,irho)
            enddo
         enddo

         ierr=max(ierr,iertmp)
         if(ierr.ne.0) then
            call xplasma_errmsg_append(s, &
                 '?integrator execution failed in xplasma_rhoth_grid_dva')
         endif

         if(ierr.eq.0) then
            call xplasma_create_blackBox(s,bbname,itype,idbb,ierr, &
                 iarray=idum,r8array=r8buf, &
                 label='dVol and dArea for '//trim(thname)//'x'//trim(rhoname), &
                 units='mixed MKS')
         endif

         ! cleanup

         deallocate(r8buf,zdva,rhogrid,thgrid)
         call xplasma_remove_item(s,idi2,iertmp)
         call xplasma_remove_item(s,idi,iertmp)
         call xplasma_author_clear(s,xplasma_xmhd,iertmp)

         if(ierr.ne.0) then
            call xplasma_errmsg_append(s,'error in xplasma_rhogrid_dva')
            return
         endif

      endif

      call xplasma_blackBox_retrieve(s,idbb,iertmp, r8a_ptr = r8buf)

      if(present(dVol)) then
         ibuf=0
         do irho=1,inrho-1
            do ith=1,inth-1
               ithx=ith
               if(.not.ccwflag) ithx=inth-ith
               ibuf=ibuf+1
               if(rhofirst) then
                  dvol(irho,ithx)=r8buf(ibuf)
               else
                  dvol(ithx,irho)=r8buf(ibuf)
               endif
            enddo
         enddo
      endif

      if(present(dArea)) then
         ibuf=(inth-1)*(inrho-1)
         do irho=1,inrho-1
            do ith=1,inth-1
               ithx=ith
               if(.not.ccwflag) ithx=inth-ith
               ibuf=ibuf+1
               if(rhofirst) then
                  darea(irho,ithx)=r8buf(ibuf)
               else
                  darea(ithx,irho)=r8buf(ibuf)
               endif
            enddo
         enddo
      endif

      contains
        subroutine checksize(arr2d,aname)

          ! verify that arr2d dimensions are correct

          real*8, dimension(:,:), intent(in) :: arr2d
          character*(*), intent(in) :: aname

          if((size(arr2d,1).ne.(inx1-1)).or.(size(arr2d,2).ne.(inx2-1))) then

             call xplasma_errmsg_append(s, &
                  '?xplasma_rhoth_grid_dva: '//aname//' array size error:')
             zmsg=' '
             write(zmsg,1001) '1st',size(arr2d,1),(inx1-1)
             call xplasma_errmsg_append(s,zmsg)
             zmsg=' '
             write(zmsg,1001) '2nd',size(arr2d,2),(inx2-1)
             call xplasma_errmsg_append(s,zmsg)
             ierr=9999
          endif
1001      format(3x,'2d array ',a,' dimension size: ',i4,'; expected: ',i4)

        end subroutine checksize

    end subroutine xplasma_rhoth_grid_dva

end module xplasma_flxint
