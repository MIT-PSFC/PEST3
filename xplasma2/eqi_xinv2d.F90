subroutine eqi_xinv2d(ivec,zrho,zth,zphi,zRtarg,zZtarg,nregion, &
     tol,ztol,zrho_bdy,ierr)
  
  !  Newton inverse map routine (internal) -- axisymmetric case
  !  magnetic coordinates

  use xplasma_definitions
  use eqi_rzbox_module
  IMPLICIT NONE

  integer ivec                      ! vector dimension
  REAL*8 zrho(ivec),zth(ivec)       ! in:  init guess, out:  answer
  REAL*8 zphi(ivec)                 ! in:  phi plane (unmodified)
  REAL*8 zRtarg(ivec),zZtarg(ivec)  ! in:  target point
  integer nregion(ivec)             ! out:  region code
  REAL*8 tol                        ! in:  relative error tolerancey
  REAL*8 ztol                       ! in:  error tolerance in (R,Z) units
  real*8 zrho_bdy                   ! last allowed rho value
  
  integer ierr                      ! out:  completion code, 0=OK

  !-------------------------------------
  !  working vectors...

  !  book-keeping, main loop
  integer, dimension(:), allocatable :: istat,indxo,ireset,ibcros,ibchkd,indxb
  real*8, dimension(:), allocatable :: zrho_use,zth_use,zphi_use
  real*8, dimension(:), allocatable :: zrtarg_use,zztarg_use,zdum,zdist
  real*8, dimension(:), allocatable :: zthb,zrtargb,zztargb,zphib
  real*8, dimension(:), allocatable :: zRdiff,zZdiff,zdeti,zdrho,zdth
  real*8, dimension(:,:), allocatable :: zRZ

  logical :: iskipi,iskipo,iexact
  !------------------

  integer itcount,itck,ibck   ! max iterations; iteration control
  integer ibck2,iertmp
  integer i,j,jvec,nc,iaxis

  REAL*8 zdthmax,zdrhomax

  REAL*8 zrhop,zrhomin

  !------------------
  !  more working vectors...

  integer iflist(6),idths(6),idrhos(6),lRnum,lZnum
  real*8 zstep

  REAL*8 zRaxis,zZaxis              ! 2d axisymm:  axis location

  real*8, parameter :: ceps15 = 1.0d-15
  real*8, parameter :: ZERO =0.0d0
  real*8, parameter :: ONE = 1.0d0
  real*8, parameter :: HALF = 0.5d0
  real*8, parameter :: cpi = 3.1415926535897931D+00

  integer, parameter :: iR=1
  integer, parameter :: iRth=2
  integer, parameter :: iRrho=3

  integer, parameter :: iZ=4
  integer, parameter :: iZth=5
  integer, parameter :: iZrho=6

  integer, parameter :: itmax=40

  integer, parameter :: ilun_dbg=6
  !-------------------------------------

  zdrhomax=0.25d0
  zdthmax=0.60d0

  ierr=0

  call xplasma_common_ids(sp,ierr, id_R=lRnum, id_Z=lZnum)
  if(ierr.ne.0) return

  allocate(zRZ(ivec,6),zRdiff(ivec),zZdiff(ivec),zdeti(ivec))
  allocate(zdrho(ivec),zdth(ivec),zthb(ivec),zrtargb(ivec),zztargb(ivec))
  allocate(zphib(ivec),zdist(ivec),zdum(ivec))
  allocate(zrtarg_use(ivec),zztarg_use(ivec),zphi_use(ivec))
  allocate(zrho_use(ivec),zth_use(ivec))
  allocate(istat(ivec),indxo(ivec),ireset(ivec),ibcros(ivec))
  allocate(ibchkd(ivec),indxb(ivec))

  itcount=0
  iaxis=0

  zrhomin=max(ceps15,tol/10)

  jvec=ivec
  do i=1,ivec
     istat(i)=0
     indxo(i)=i
     ireset(i)=0
     ibcros(i)=0
     ibchkd(i)=0
     nregion(i)=0
     zrho_use(i)=max(zrhomin,zrho(i))
     zth_use(i)=zth(i)
     zphi_use(i)=zphi(i)
     zrtarg_use(i)=zrtarg(i)
     zztarg_use(i)=zztarg(i)
  enddo

  iflist(1:3)=lRnum
  iflist(4:6)=lZnum
  idths=0; idths(2)=1; idths(5)=1
  idrhos=0; idrhos(3)=1; idrhos(6)=1

  zstep = ONE

  !  iterative search loop
  do

     if(itcount.eq.(itmax/2)) zstep = HALF

     !  get value & rho,th derivatives

     call xplasma_RZeval_2d(sp,iflist, &
          xplasma_theta_coord,zth_use, xplasma_rho_coord,zrho_use, &
          zRZ, ierr, ideriv1s=idths,ideriv2s=idrhos)

     if(ierr.ne.0) exit

     !  convergence test

     nc=0
     itck=0
     ibck=0
     ibck2=0
     do i=1,jvec
        zRdiff(i)=zRtarg_use(i)-zRZ(i,iR)
        zZdiff(i)=zZtarg_use(i)-zRZ(i,iZ)

#ifdef __DEBUG
        if(itcount.ge.itmax-10) then
           write(ilun_dbg,1001) &
                indxo(i),zrho_use(i),zth_use(i),zRdiff(i),zZdiff(i)
1001       format(1x,i4,' (rho,th) = (',2(1x,1pd11.4),') (dR,dZ) = (', &
                2(1x,1pd11.4))
        endif
#endif

        if(max(abs(zRdiff(i)),abs(zZdiff(i))).lt.ztol) then
           !  OK
           itck=itck+1
           istat(i)=1
           j=indxo(i)
           zrho(j)=zrho_use(i)
           zth(j)=zth_use(i)
           if(zrho_use(i).gt.ONE) then
              nregion(j)=2
           else
              nregion(j)=1
           endif
        else
           nc=nc+1
           !  another iteration

           zdeti(i)=ONE/(zRZ(i,iRrho)*zRZ(i,iZth)-zRZ(i,iZrho)*zRZ(i,iRth))

           zdrho(i)=zdeti(i)*(zRdiff(i)*zRZ(i,iZth)-zZdiff(i)*zRZ(i,iRth))
           if(zdrho(i).gt. zdrhomax) then
              zdrho(i)=zdrhomax
           else if(zdrho(i).lt.-zdrhomax) then
              zdrho(i)=-zdrhomax
           endif

           zdth(i)=zdeti(i)*(-zRdiff(i)*zRZ(i,iZrho)+zZdiff(i)*zRZ(i,iRrho))
           if(zdth(i).gt.(zdthmax)) then
              zdth(i)=zdthmax
           else if(zdth(i).lt.-zdthmax) then
              zdth(i)=-zdthmax
           endif

           zrhop=zrho_use(i)
           zrho_use(i)=zrho_use(i)+zstep*zdrho(i)
           if((zrho_use(i).eq.ONE).and.(zrhop.gt.ONE).and. &
                (zdrho(i).eq.-zdrhomax)) then
              zrho_use(i)=HALF*(zrho_use(i)+zrhop)
           endif

           !  if interior bdy check was done: use this information.

           if(ibchkd(i).eq.1) then
              zrho_use(i)=min(zrho_use(i),(ONE-ceps15))
           endif
           if(ibchkd(i).eq.2) then
              zrho_use(i)=max(zrho_use(i),(ONE+ceps15))
           endif

           !  check for interior bdy crossing...

           if((zrhop-ONE)*(zrho_use(i)-ONE).le.ZERO) then
              ibcros(i)=ibcros(i)+1
           else
              ibcros(i)=0
           endif

           if(zrho_use(i).gt.zrho_bdy) then
              itck=itck+1
              ibck=ibck+1
              istat(i)=-1
              if(iaxis.eq.0) iaxis=-1

           else if(ibcros(i).eq.4) then

              !  solution is oscillating across the (ONE) boundary-- use bdfind

              itck=itck+1
              ibck2=ibck2+1
              istat(i)=-2
              if(iaxis.eq.0) iaxis=-1

           else
              if(zrho_use(i).lt.ZERO) then

                 !  if cross axis:  for next step, go to other side;
                 !  delta(angle)=0.9*pi
                 !  make it hunt around the axis...

                 zrho_use(i)=-zrho_use(i)
                 if(zth_use(i).gt.ZERO) then
                    zth_use(i)=zth_use(i)-0.9d0*cpi
                 else
  !  .8 instead of .9 here, prevents possible oscillation
                    zth_use(i)=zth_use(i)+0.8d0*cpi
                 endif
              else

  !  normal next step

                 zth_use(i)=zth_use(i)+zstep*zdth(i)
              endif
              if(zrho_use(i).lt.zrhomin) then
                 zrho_use(i)=zrhomin
              endif
           endif

#ifdef __DEBUG
           if(itcount.ge.itmax-10) then
              write(ilun_dbg,*) indxo(i),' status = ',istat(i), &
                   ' new rho,th = ',zrho_use(i),zth_use(i)
              write(ilun_dbg,*) '  Rtarg,Ztarg = ',zRtarg_use(i),zZtarg_use(i)
           endif
#endif
        endif
     enddo

     !  get axis location if needed

     if(iaxis.eq.-1) then
        iaxis=1
        call xplasma_mag_axis(sp,iertmp, zRaxis,zZaxis)

        !  check if any points are close to axis -- if so, declare them
        !  converged...

        do i=1,jvec
           if(istat(i).ne.1) then
              if(max(abs(zRtarg_use(i)-zRaxis), &
                   abs(zZtarg_use(i)-zZaxis)).lt.ztol) then
                 if(istat(i).eq.-1) ibck=ibck-1
                 istat(i)=1
                 nc=nc-1               ! decrement non-converged ctr
                 j=indxo(i)
                 zrho(j)=ZERO
                 zth(j)=ZERO
                 nregion(j)=1
              endif
           endif
        enddo

     endif

     if(ibck.gt.0) then

        !  pts strayed out of bounds... use eqi_bdfind
        !  process resets...

        j=0
        do i=1,jvec
           if(istat(i).eq.-1) then
              if(ireset(i).eq.1) then
                 nregion(i)=4
                 ierr=2511
                 call xplasma_errmsg_append(sp, &
                      ' search crossed boundary more than once.')
                 exit
              endif
              ireset(i)=1
              j=j+1
              zRtargb(j)=zRtarg_use(i)
              zZtargb(j)=zZtarg_use(i)
              zphib(j)=zphi_use(i)
              indxb(j)=i
           endif
        enddo
        if(j.ne.ibck) call errmsg_exit(' ?eqi_xinv2d book-keeping bug!')
        if(ierr.gt.0) exit

        iskipi=.FALSE.
        iskipo=.FALSE.
        iexact=.TRUE.
        call eqi_bdfind(ibck,zRtargb,zZtargb,zphib,zrho_bdy,iskipi,iskipo, &
             iexact,zthb,zdum,zdist,ierr)
        if(ierr.ne.0) exit

        do j=1,ibck
           i=indxb(j)
           zrho_use(i)=zrho_bdy
           zth_use(i)=zthb(j)
           if(zdist(j).gt.ZERO) then
              !  this code left over from original xplasma, but, since
              !  eqi_bdfind screening was already done, it should not be
              !  reached...
              istat(i)=1
              nregion(indxo(i))=3
              zrho(indxo(i))=zrho_use(i)
              zth(indxo(i))=zth_use(i)
              nc=nc-1
           else
              istat(i)=0
              itck=itck-1
           endif
        enddo
     endif

     if(ibck2.gt.0) then

        !  pts oscillating across (ONE) interior boundary

        j=0
        do i=1,jvec
           if(istat(i).eq.-2) then
              if(ibchkd(i).ne.0) then
                 ierr=2511
                 call xplasma_errmsg_append(sp, &
                      '?eqi_xinv2d:  convergence failure (bdy oscillation)')
                 nregion(i)=4
                 exit
              endif
              j=j+1
              zRtargb(j)=zRtarg_use(i)
              zZtargb(j)=zZtarg_use(i)
              zphib(j)=zphi_use(i)
              indxb(j)=i
           endif
        enddo
        if(ierr.ne.0) exit

        if(j.ne.ibck2) call errmsg_exit(' ?eqi_xinv2d book-keeping bug!')

        iskipi=.FALSE.
        iskipo=.FALSE.
        iexact=.TRUE.
        call eqi_bdfind(ibck2,zRtargb,zZtargb,zphib,ONE,iskipi,iskipo,iexact, &
             zthb,zdum,zdist,ierr)
        if(ierr.ne.0) exit

        do j=1,ibck2
           i=indxb(j)
           if(zdist(j).gt.ZERO) then
              zrho_use(i)=ONE+ceps15
              ibchkd(i)=2
           else
              zrho_use(i)=ONE-ceps15
              ibchkd(i)=1
           endif
           zth_use(i)=zthb(j)
           if(abs(zdist(j)).lt.ztol) then
              istat(i)=1                    ! close enough: treat as solution
              if(zdist(j).gt.ZERO) then
                 nregion(indxo(i))=2
              else
                 nregion(indxo(i))=1
              endif
              zrho(indxo(i))=zrho_use(i)
              zth(indxo(i))=zth_use(i)
              nc=nc-1
           else
              istat(i)=0                    ! keep looking...
              itck=itck-1
           endif
        enddo
     endif

     !--------------------------
      if(nc.gt.0) then
         !  another iteration is in store...
         itcount=itcount+1
#ifdef __DEBUG
         if(itcount.eq.itmax-10) then
            write(ilun_dbg,*) '%eqi_xinv2d:  convergence warning!'
         endif
#endif
         if(itcount.gt.itmax) then
            call xplasma_errmsg_append(sp, &
                 '?eqi_xinv2d:  convergence failure (itmax)')
            ierr=2511
            exit
         endif

         if(itck.eq.0) cycle         ! next iteration now

         !  status change -- one of the points converged, or, one of the
         !  points went out of bounds...

         !  first reduce vector to unconverged points oinly

         itck=0
         j=0
         do i=1,jvec
            if(istat(i).ne.1) then
               j=j+1
               if(j.lt.i) then
                  istat(j)=istat(i)
                  indxo(j)=indxo(i)
                  ireset(j)=ireset(i)
                  ibcros(j)=ibcros(i)
                  ibchkd(j)=ibchkd(i)
                  zrho_use(j)=zrho_use(i)
                  zth_use(j)=zth_use(i)
                  zphi_use(j)=zphi_use(i)
                  zrtarg_use(j)=zrtarg_use(i)
                  zztarg_use(j)=zztarg_use(i)
               endif
            endif
         enddo
         jvec=j
         if(jvec.eq.0) call errmsg_exit('eqi_xinv2d: algorithm error!')

  !  ok now do the next iteration

         cycle
      endif

      exit
   enddo

   deallocate(zRZ,zRdiff,zZdiff,zdeti,zdrho,zdth,zthb,zrtargb,zztargb)
   deallocate(zphib,zdist,zdum,zrtarg_use,zztarg_use,zphi_use)
   deallocate(zrho_use,zth_use)
   deallocate(istat,indxo,ireset,ibcros,ibchkd,indxb)

 end subroutine eqi_xinv2d
