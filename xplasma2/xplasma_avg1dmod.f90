module xplasma_avg1dmod

  !  Added dmc 16 Oct 2006 -- 1d averager "xplasma_avg1d"
  !   compute int[f*dV]/int[dV] or int[f*w*dV]/int[w*dV]

  use xplasma_obj
  use xplasma_profs

  implicit NONE

  private

  public :: xplasma_avg1d,xplasma_avg1d_single,xplasma_avg1d_vec

  interface xplasma_avg1d
     module procedure xplasma_avg1d_single,xplasma_avg1d_vec
  end interface

  contains

    subroutine xplasma_avg1d_single(s,xtarg,idd,idf,answr,ierr,idw,ismoo)

      !  1d average
      !   unweighted:
      !               int[xtarg(j):xtarg(j+1)]{dx*(d(f[idd](x))/dx)*f[idf](x)}
      !    answr(j) = ------------------------------------------------------
      !                   int[xtarg(j):xtarg(j+1)]{dx*(d(f[idd](x))/dx)}
      !
      !   weighted -- multiply numerator and denominator integrands by
      !               w[idw](x).

      type(xplasma), pointer :: s

      real*8, dimension(:), intent(in) :: xtarg  ! integration boundaries
      !  xtarg(1:2) -- 1st zone; xtarg(2:3) -- 2nd zone, etc.

      integer, intent(in) :: idd    ! ID of differential function
      ! typically Volume(rho) or Area(rho); integral will then be 
      ! weighted by dVolume/drho or dArea/drho.

      integer, intent(in) :: idf    ! ID of function to be averaged

      real*8, dimension(:), intent(out) :: answr ! result of averaging

      integer, intent(out) :: ierr  ! completion code (0=OK)

      integer, intent(in), optional :: idw   ! weighting function ID (optional)

      logical, intent(in), optional :: ismoo ! optional smoothing (default: F)

      ! if SMOOTHING is selected, all functions and weighting functions must
      ! be defined over the same GRID; if not, it is sufficient if they are
      ! defined over the same COORDINATE.

      !---------------------------------------------------
      !  local:
      integer :: idfa(1),idwa(1)
      logical :: jsmoo
      real*8 :: ansa(size(answr),1)
      real*8, parameter :: ZERO = 0.0d0
      !---------------------------------------------------

      idfa = idf
      if(present(idw)) then
         idwa = idw
      else
         idwa = 0
      endif

      jsmoo=.FALSE.
      if(present(ismoo)) jsmoo=ismoo

      call xplasma_avg1d_vec(s,xtarg,idd,idfa,ansa,ierr,idw=idwa,ismoo=jsmoo)

      if(ierr.ne.0) then
         answr = ZERO
      else
         answr = ansa(:,1)
      endif

    end subroutine xplasma_avg1d_single

    subroutine xplasma_avg1d_vec(s,xtarg,idd,idf,answr,ierr,idw,ismoo)

      !  1d averager -- same as xplasma_avg1d_single but apply averaging
      !        operator to a list of profile functions idf(...).

      type(xplasma), pointer :: s

      real*8, dimension(:), intent(in) :: xtarg  ! integration boundaries
      !  xtarg(1:2) -- 1st zone; xtarg(2:3) -- 2nd zone, etc.

      integer, intent(in) :: idd    ! ID of differential function
      ! typically Volume(rho) or Area(rho); integral will then be 
      ! weighted by dVolume/drho or dArea/drho.

      integer, intent(in), dimension(:) :: idf ! ID of functions to be averaged

      real*8, dimension(:,:), intent(out) :: answr ! results of averaging

      integer, intent(out) :: ierr  ! completion code (0=OK)

      integer, intent(in), dimension(:), optional :: idw  ! weighting functions

      logical, intent(in), optional :: ismoo ! optional smoothing (default: F)

      ! if SMOOTHING is selected, all functions and weighting functions must
      ! be defined over the same GRID of "rho"; if smoothing is not selected,
      ! it is sufficient if all functions are defined over the same COORDINATE.

      !---------------------------------------------------
      ! local:
      integer :: ixid(100),isizes(100),inxid,icoord,jcoord,idx,jdx,irank,isize
      integer :: ii,ix,inx,inf,id,isizex,inbrk,inbrk_new,inint,iertmp,invec
      integer :: jj,isize_total,iaoff,ixtarg,inmap,imap,iadr,inxu,idarea,iflag
      integer :: idtmp,idvol
      logical :: perio
      logical :: jsmoo
      real*8 :: xmin,xmax,bdytol,ztol,zwgt
      real*8, dimension(:), allocatable :: x,xu,xgrid,xbrk,xbrk_new
      real*8, dimension(:), allocatable :: xvec,wvec,wkf
      real*8, dimension(:,:), allocatable :: eval_arr,fu,denoma
      integer, dimension(:), allocatable :: idint,ideriv,idwa,idwmap

      real*8 :: wkx(size(xtarg)),wkd(size(xtarg))

      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter :: c2pi = 6.2831853071795862D+00
      character*32 profname
      !---------------------------------------------------

      ierr=0
      answr=0

      jsmoo=.FALSE.
      if(present(ismoo)) jsmoo=ismoo

      !----------------------------
      ! check rank of differential element

      call xplasma_prof_info(s,idd,ierr, gridId1=idx, rank=irank, &
           name=profname)

      if(irank.gt.1) then
         ierr=512

         call xplasma_errmsg_append(s, &
              ' ?xplasma_avg1d: rank of profile >1: '//trim(profname))
      endif
      if(ierr.ne.0) return

      !-----------------------------
      !  check argument array sizes

      inx=size(xtarg)
      if(inx.ne.(size(answr,1)+1)) then
         ierr=510
         call xplasma_errmsg_append(s, &
              ' ?xplasma_avg1d: size of result vector must be 1 less than size of the "xtarg".')
      endif

      inf=size(idf)
      if(inf.ne.size(answr,2)) then
         ierr=510
         call xplasma_errmsg_append(s, &
              ' ?xplasma_avg1d: number of result vectors must match number of profile input IDs.')
      endif
      if(ierr.ne.0) return
         
      if(present(idw)) then
         if(inf.ne.size(idw)) then
            ierr=510
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_avg1d: the number of weight profiles must match')
            call xplasma_errmsg_append(s, &
                 '  the number of integrand profiles.')
         endif
      endif
      if(ierr.ne.0) return

      call xplasma_global_info(s,ierr, bdytol=bdytol)
      if(ierr.ne.0) return

      !------------------------------
      !  check grid zone bdy inputs
      !    make reference to grid of "idd" profile

      call xplasma_grid_info(s,idx,ierr, &
           size=isizex, xmin=xmin, xmax=xmax, perio=perio, &
           coord=icoord)

      ztol = bdytol*max(abs(xmin),abs(xmax))

      do ix=1,inx-1
         if(xtarg(ix+1).le.xtarg(ix)+2*ztol) then
            ierr=603
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_avg1d: "xtarg" values must be strict ascending.')
            exit
         endif
      enddo
      if(ierr.ne.0) return

      if(.not.perio) then
         if((xtarg(1).lt.xmin-ztol).or.(xtarg(inx).gt.xmax+ztol)) then
            ierr=609
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_avg1d: out-of-range "xtarg" values detected.')
         endif

      else
         if((xtarg(inx)-xtarg(1)).gt.c2pi+ztol) then
            ierr=609
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_avg1d: periodic "xtarg" coordinate range exceeds 2pi.')
         endif
      endif

      if(ierr.ne.0) return

      allocate(x(inx))
      x = xtarg

      if(.not.perio) then
         x(1)=max(xmin,x(1))
         x(inx)=min(xmax,x(inx))
      else
         x(inx)=min(x(inx),x(1)+c2pi)
      endif

      !------------------------------
      !  OK, form list of grids; all must be over the same coordinate
      !  Check that all profiles are rank 1 over the same coordinate

      inxid=1
      ixid(1)=idx
      isizes(1)=isizex
      jdx=0

      if(present(idw)) then
         do ii=1,inf
            id=idw(ii)
            if(id.eq.0) cycle
            call xplasma_prof_info(s,id,ierr, rank=irank, gridId1=idx)
            if(ierr.ne.0) then
               call xplasma_errmsg_append(s, &
                    ' ?xplasma_avg1d: error looking up weight profile info.')
            else
               if(irank.gt.1) then
                  ierr=512
                  call xplasma_errmsg_append(s, &
                       ' ?xplasma_avg1d: weight profile rank exceeds 1')
               endif
            endif
            if(ierr.eq.0) then
               call xplasma_grid_info(s,idx,ierr, coord=jcoord, size=isizex)
               if(icoord.ne.jcoord) then
                  ierr=511
                  call xplasma_errmsg_append(s, &
                       ' ?xplasma_avg1d: weight and differential profiles (idd) and (idw)')
                  call xplasma_errmsg_append(s, &
                       '  must be defined over a common coordinate.')
               endif
               if(jsmoo) then
                  if(jdx.eq.0) jdx=idx
                  if(jcoord.ne.xplasma_rho_coord) then
                     ierr=511
                     call xplasma_errmsg_append(s, &
                          ' ?xplasma_avg1d: only functions of form f(rho)')
                     call xplasma_errmsg_append(s, &
                          '  can be smoothed with this software.')
                  endif
                  if(idx.ne.jdx) then
                     ierr=511
                     call xplasma_errmsg_append(s, &
                          ' ?xplasma_avg1d: when smoothing is activated,')
                     call xplasma_errmsg_append(s, &
                          '  all profiles must be defined over the same grid.')
                  endif
               endif
            endif
            if(ierr.eq.0) call addx(ixid,isizes,inxid,idx,isizex)
            if(ierr.ne.0) exit
         enddo

         if(ierr.ne.0) then
            deallocate(x)
            return
         endif
      endif

      do ii=1,inf
         id=idf(ii)
         call xplasma_prof_info(s,id,ierr, rank=irank, gridId1=idx)
         if(ierr.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_avg1d: error looking up profile information.')
            exit
         else if(irank.gt.1) then
            ierr=512
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_avg1d: an integrand profile rank exceeds 1')
            exit
         endif
         call xplasma_grid_info(s,idx,ierr, coord=jcoord, size=isizex)
         if(icoord.ne.jcoord) then
            ierr=511
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_avg1d: all integrand profiles must be defined over a common coordinate.')
            exit
         endif
         if(jsmoo) then
            if(jdx.eq.0) jdx=idx
            if(jcoord.ne.xplasma_rho_coord) then
               ierr=511
               call xplasma_errmsg_append(s, &
                    ' ?xplasma_avg1d: only functions of form f(rho)')
               call xplasma_errmsg_append(s, &
                    '  can be smoothed with this software.')
            endif
            if(idx.ne.jdx) then
               ierr=511
               call xplasma_errmsg_append(s, &
                    ' ?xplasma_avg1d: when smoothing is activated,')
               call xplasma_errmsg_append(s, &
                    '  all profiles must be defined over the same grid.')
            endif
            if(ierr.ne.0) exit
         endif

         call addx(ixid,isizes,inxid,idx,isizex)
      enddo

      if(ierr.ne.0) then
         deallocate(x)
         return
      endif

      !------------------------------------
      ! having gotten this far, all grids are over the same coordinate, and
      ! a list of distinct grids has been formed.
      !
      ! make a compressed list of weighting functions (if any)

      inmap=0
      allocate(idwa(inf),idwmap(inf))

      if(present(idw)) then
         do ii=1,inf
            if(idw(ii).eq.0) then
               idwmap(ii)=0
            else
               imap=0
               do jj=1,inmap
                  if(idwmap(jj).eq.idw(ii)) then
                     imap=jj
                     exit
                  endif
               enddo
               if(imap.gt.0) then
                  idwmap(ii)=imap
               else
                  inmap=inmap+1
                  idwa(inmap)=idw(ii)
                  idwmap(ii)=inmap
               endif
            endif
         enddo
      else
         idwa=0
         idwmap=0
      endif

      ! now "inmap" is the number of unique weighting profiles; could be zero

      ! now form a single grid which spans the user requested integration
      ! range but also has all interior grid boundaries added to it...

      !--------------
      ! choose a safe upper limit on the number of x break points in the
      ! integration

      if(jsmoo) then
         inxu=isizex
         allocate(xu(inxu))
         call xplasma_grid(s,idx,xu,ierr)
         call xplasma_profId(s,'XPLASMA_AREA',idarea)
         call xplasma_profId(s,'XPLASMA_VOLUME',idvol)
         if(idd.eq.idarea) then
            iflag=2
         else if(idd.eq.idvol) then
            iflag=1
         else
            iflag=idd  ! assumed >2 and will be...
         endif
      else
         inxu=inx
         allocate(xu(inxu))
         xu=x
      endif

      isize_total=inxu
      do ii=1,inxid
         isize_total=isize_total+isizes(ii)
      enddo

      allocate(xbrk(isize_total),xbrk_new(isize_total))

      inbrk=inxu
      xbrk(1:inxu)=xu(1:inxu)

      !--------------
      ! use library routine to insert break points without duplicates

      do ii=1,inxid
         allocate(xgrid(isizes(ii)))
         call xplasma_grid(s,ixid(ii),xgrid,iertmp)
         if(perio) then
            call fluxav_brk_addperio(xbrk,inbrk,xgrid,isizes(ii), &
                 xbrk_new,inbrk_new)
         else
            call fluxav_brk_addsub(xbrk,inbrk,xgrid,isizes(ii), &
                 xbrk_new,inbrk_new)
         endif
         call fluxav_brktol(xbrk_new,inbrk_new,ztol)
         if(inbrk_new.gt.inbrk) then
            inbrk = inbrk_new
            xbrk(1:inbrk)=xbrk_new(1:inbrk)
         endif
         deallocate(xgrid)
      enddo

      !  gather information needed to interpolate profile to necessary
      !  locations

      inint = inf + 1 + inmap

      allocate(idint(inint),ideriv(inint))
      ideriv(1)=1
      ideriv(2:inint)=0
      idint(1)=idd
      idint(2:inf+1) = idf
      if(inmap.gt.0) then
         idint(inf+2:inf+2+inmap-1) = idwa(1:inmap)
      endif

      !  x vector -- 10 gauss points per break interval

      invec = 10*(inbrk-1)
      allocate(xvec(invec),wvec(invec),eval_arr(invec,inint))

      !  get Gauss integration vector and weights

      call xplasma_integ_mkvec(inbrk,xbrk,xvec,wvec)

      !  interpolate the data needed for the integration
      !  this includes the differential profile, integrand profiles, and
      !  whatever weighting profiles (if any) that were specified...

      do
         call xplasma_eval_prof(s,idint,xvec,eval_arr,ierr, ideriv1s=ideriv)
         if(ierr.ne.0) exit

         if(.not.jsmoo) then
            !  compute the integration
            !   (answr array already was cleared to zero

            allocate(denoma(inx-1,inf)); denoma=ZERO

            ixtarg = 1

            do ix=1,invec
               if(xvec(ix).gt.xu(ixtarg+1)) then
                  ! start next zone
                  ixtarg = ixtarg + 1
               endif

               do ii=1,inf
                  if(idwmap(ii).gt.0) then
                     iadr=inf+1+idwmap(ii)
                     zwgt = eval_arr(ix,1)*eval_arr(ix,iadr)*wvec(ix)
                  else
                     zwgt = eval_arr(ix,1)*wvec(ix)
                  endif
                  zwgt=abs(zwgt)
                  denoma(ixtarg,ii)=denoma(ixtarg,ii)+zwgt
                  answr(ixtarg,ii)=answr(ixtarg,ii)+zwgt*eval_arr(ix,ii+1)
               enddo
            enddo

            do ii=1,inf
               call zclean(denoma(:,ii),1)
               call zclean(answr(:,ii),2)
               do ixtarg=1,inx-1
                  if(denoma(ixtarg,ii) == ZERO) then
                     answr(ixtarg,ii) = ZERO
                  else
                   !  write(*,'(I2,1x,I2,4(1x,E13.6))') ixtarg,ii,answr(ixtarg,ii),denoma(ixtarg,ii) &
                   !       ,answr(ixtarg,ii)/denoma(ixtarg,ii),eval_arr((ixtarg-1)*10+5,ii+1)

                     answr(ixtarg,ii)=answr(ixtarg,ii)/denoma(ixtarg,ii)
                  endif
               enddo
            enddo

         else
            ! calculate on original grid...

            allocate(fu(inxu-1,inf)); fu=ZERO
            allocate(denoma(inxu-1,inf)); denoma=ZERO

            ixtarg = 1

            do ix=1,invec
               if(xvec(ix).gt.xu(ixtarg+1)) then
                  ! start next zone; inherit sum from prior zone
                  ixtarg = ixtarg + 1
               endif

               do ii=1,inf
                  if(idwmap(ii).gt.0) then
                     iadr=inf+1+idwmap(ii)
                     zwgt = eval_arr(ix,1)*eval_arr(ix,iadr)*wvec(ix)
                  else
                     zwgt = eval_arr(ix,1)*wvec(ix)
                  endif
                  zwgt=abs(zwgt)
                  fu(ixtarg,ii)=fu(ixtarg,ii)+eval_arr(ix,ii+1)*zwgt
                  denoma(ixtarg,ii)=denoma(ixtarg,ii)+zwgt
               enddo
            enddo

            do ii=1,inf
               call zclean(denoma(:,ii),1)
               call zclean(fu(:,ii),2)
               do ixtarg=1,inxu-1
                  if(denoma(ixtarg,ii).eq.ZERO) then
                     fu(ixtarg,ii)=ZERO
                  else
                     fu(ixtarg,ii)=fu(ixtarg,ii)/denoma(ixtarg,ii)
                  endif
               enddo
            enddo

            ! compute smoothed integrals of functions
            call xplasma_author_set(s,'xplasma_profs',iertmp)
            do ii=1,inf
               if(idwmap(ii).eq.0) then
                  allocate(wkf(inxu))
                  call xplasma_eval_prof(s,idd,xu,wkf,ierr)
                  if(ierr.ne.0) exit
                  fu(1,ii)=fu(1,ii)*(wkf(2)-wkf(1))  ! multiply in *dV
                  do ix=2,inxu-1
                     fu(ix,ii)=fu(ix-1,ii)+fu(ix,ii)*(wkf(ix+1)-wkf(ix))
                  enddo
                  deallocate(wkf)
                  call xplasma_irhofun(s,idx,'___TMP',inxu-1,fu(1:inxu-1,ii), &
                       iflag,idtmp,ierr)
                  if(ierr.ne.0) exit
                  call xplasma_eval_prof(s,idtmp,x,wkx,ierr)
                  if(ierr.ne.0) exit
                  call xplasma_eval_prof(s,idd,x,wkd,ierr)
                  if(ierr.ne.0) exit
                  answr(1:inx-1,ii)=(wkx(2:inx)-wkx(1:inx-1))/ &
                       (wkd(2:inx)-wkd(1:inx-1))
               else
                  allocate(wkf(inxu))
                  wkf(2:inxu-1)=(fu(1:inxu-2,ii)+fu(2:inxu-1,ii))/2
                  call fextrap(xu,fu(:,ii),wkf(1),wkf(inxu))
                  call xplasma_create_1dprof(s,'___TMP',idx,wkf,idtmp,ierr, &
                       ispline=1,ibca=1,zbca=ZERO)   ! Hermite interpolation
                  deallocate(wkf)
                  if(ierr.ne.0) exit
                  call xplasma_eval_prof(s,idtmp, &
                       (x(1:inx-1)+x(2:inx))/2, answr(1:inx-1,ii), ierr)
                  if(ierr.ne.0) exit
               endif
            enddo
            call xplasma_remove_item(s,idtmp,iertmp)
            call xplasma_author_clear(s,'xplasma_profs',iertmp)

            deallocate(fu)
         endif
         exit
      enddo

      if(allocated(xu)) deallocate(xu)
      deallocate(xvec,wvec,eval_arr,xbrk,xbrk_new,x)
      deallocate(idint,ideriv,idwa,idwmap)
      if(allocated(denoma)) deallocate(denoma)

    contains
      subroutine fextrap(xbdy,fzon,f0,f1)
        real*8, dimension(:), intent(in) :: xbdy
        real*8, dimension(:), intent(in) :: fzon
        real*8, intent(out) :: f0,f1

        ! quadratic extrap to center; linear extrap to edge; constrain
        ! against sign change

        !---------------------------------
        real*8 :: dx1,dx2
        integer :: inz,inb
        real*8, parameter :: tenth = 0.1d0
        !---------------------------------

        inz=size(fzon)
        inb=size(xbdy)

        dx1=xbdy(2)/2
        dx2=(xbdy(2)+xbdy(3))/2
        f0 = (dx2*dx2*fzon(1)-dx1*dx1*fzon(2))/(dx2*dx2-dx1*dx1)
        if((f0-tenth*fzon(1))*fzon(1).le.0.0d0) then
           f0=tenth*fzon(1)
        endif

        dx1=(xbdy(inb)-xbdy(inb-1))/2
        dx2=xbdy(inb) - (xbdy(inb-1)+xbdy(inb-2))/2
        f1 = (dx2*fzon(inz)-dx1*fzon(inz-1))/(dx2-dx1)
        if((f1-tenth*fzon(inz))*fzon(inz).le.0.0d0) then
           f1=tenth*fzon(inz)
        endif
      end subroutine fextrap

      subroutine zclean(arr,itest)
        ! any value less than 10**-10 * max(abs(arr)) set to zero

        real*8, dimension(:) :: arr
        integer :: itest

        !-------------------------
        integer :: ixi,isizi
        real*8 :: amaxa,ztest
        real*8, parameter :: small_rat = 1.0d-10
        !-------------------------

        isizi=size(arr)
        amaxa=ZERO

        do ixi=1,isizi
           amaxa=max(amaxa,abs(arr(ixi)))
        enddo

        do ixi=1,isizi
           if(itest.eq.1) then
              ztest=arr(ixi)
           else if(itest.eq.2) then
              ztest=abs(arr(ixi))/itest
           endif
           if(ztest.lt.small_rat*amaxa) then
              arr(ixi)=ZERO
           endif
        enddo

      end subroutine zclean

    end subroutine xplasma_avg1d_vec
         
    subroutine addx(ixid,isizes,inxid,idx,isizex)
      ! **private**
      ! add to list of x axis IDs

      integer, dimension(:) :: ixid    ! the list of grid IDs
      integer, dimension(:) :: isizes  ! the list of grid sizes
      integer :: inxid  ! list size (possibly incremented on output
      integer :: idx    ! grid to add IF it is not there already.
      integer :: isizex ! grid size

      integer :: ii,imatch

      imatch=0
      do ii=1,inxid
         if(idx.eq.ixid(ii)) then
            imatch=ii
            exit
         endif
      enddo

      if(imatch.eq.0) then
         inxid=inxid+1
         ixid(inxid)=idx
         isizes(inxid)=isizex
      endif

    end subroutine addx

end module xplasma_avg1dmod
