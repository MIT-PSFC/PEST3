subroutine eqi_bprobe(inum,itype,imaptype,ra,za,da,rb,zb,db,ierr)

  !  generate a set of probe segments (Ra,Za)--(Rb,Zb) such that
  !  the (Ra,Za) are all inside, and the (Rb,Zb) all outside, the
  !  given surface:

  !    itype=1:  plasma boundary
  !    itype=2:  limiter

  !  also make sure to get (inum) segments with start and stop point
  !  monotonically increasing in poloidal angle, and, a sequence whose
  !  bdy intersections will procede in an orderly manner around the
  !  boundary contour.

  use xplasma_definitions
  use eqi_rzbox_module
  implicit NONE

  !  input:
  integer inum                      ! number of segments wanted
  integer itype                     ! bdy to straddle
  integer imaptype                  ! mapping accuracy control
  ! for computing distance: 1=slow&exact, 2=medium&accurate, 3=fast&approximate

  !  output:
  real*8 ra(inum),za(inum)          ! segment points just inside bdy
  real*8 da(inum)                   ! distance inside bdy (.lt.0)
  real*8 rb(inum),zb(inum)          ! segment points just outside bdy
  real*8 db(inum)                   ! distance outside bdy (.gt.0)

  integer ierr                      ! completion code, 0=OK

  !---------------------------------

  integer :: nrhox,nchix
  real*8, dimension(:), allocatable :: zrho,zchi,zphi,zdist,zchix,zrhox
  real*8, dimension(:), allocatable :: zchiv
  real*8, dimension(:,:), allocatable :: zbuf,zbuf2

  integer i,j,id_rhox,id_chix,ilist(4),idchis(4)
  integer ifnd,imiss,indx

  real*8 zdx,znorm,Rfastmin,Rfastmax,Zfastmin,Zfastmax,zindx,zxchi

  logical :: sol,axisymm

  real*8, parameter :: ZERO=0.0d0
  real*8, parameter :: ONE=1.0d0

  !---------------------------------

  ierr=0

  call xplasma_global_info(sp,ierr, scrapeoff=sol, axisymm=axisymm)
  if(ierr.ne.0) return

  if(.not.axisymm) then
     call xplasma_errmsg_append(sp,' ?eqi_bprobe: axisymmetry required.')
     ierr=107
  endif

  if(.not.sol) then
     call xplasma_errmsg_append(sp, &
          ' ?eqi_bprobe: scrapeoff region must be defined before this routine can beused.')
     ierr=109
  endif
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'__RHOX',id_rhox,ierr)
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'__THETAX',id_chix,ierr)
  if(ierr.ne.0) return

  call xplasma_grid_size(sp,id_rhox,nrhox,ierr)
  if(ierr.ne.0) return

  call xplasma_grid_size(sp,id_chix,nchix,ierr)
  if(ierr.ne.0) return

  allocate(zchix(nchix))
  call xplasma_grid(sp,id_chix,zchix,ierr)

  allocate(zrhox(nrhox))
  call xplasma_grid(sp,id_rhox,zrhox,ierr)

  call xplasma_RZminmax_extended(sp, &
       Rfastmin,Rfastmax,Zfastmin,Zfastmax, ierr)
  if(ierr.ne.0) return

  allocate(zrho(inum),zchi(inum))
  allocate(zchiv(nrhox))
  allocate(zbuf(inum,4),zbuf2(0:nrhox,2),zdist(0:nrhox))

  !  (R,Z) function nos.

  call xplasma_common_ids(sp,ierr, id_R=ilist(1), id_Z=ilist(2))
  ilist(3)=ilist(1)
  ilist(4)=ilist(2)
  idchis(1)=0; idchis(2)=0
  idchis(3)=1; idchis(4)=1

  !  straddle plasma bdy
  !  this is also done initially for the limiter calculation (itype=2)
  !  in case the limiter is on the plasma bdy

  zdx=min((Rfastmax-Rfastmin),(Zfastmax-Zfastmin))/(10*nchix) ! small step

  do i=1,inum
     zindx=ONE+(i-1)*(nchix-ONE)/(inum)   ! last pt short of 1st pt
     indx=max(1,min((nchix-1),int(zindx)))
     zxchi=max(ZERO,min(ONE,(zindx-indx)))

     zrho(i)=ONE
     zchi(i)=zchix(indx)+zxchi*(zchix(indx+1)-zchix(indx))

  enddo

  call xplasma_RZeval_2d(sp,ilist, &
       xplasma_theta_coord,zchi(1:inum),xplasma_rho_coord,zrho(1:inum), &
       zbuf,ierr, ideriv1s=idchis)

  if(ierr.ne.0) return

  !  zbuf(1:ivec,1):  Rbdy
  !  zbuf(1:ivec,2):  Zbdy
  !  zbuf(1:ivec,3):  dR/dchi
  !  zbuf(1:ivec,4):  dZ/dchi

  do i=1,inum
     znorm=ONE/sqrt(zbuf(i,2)**2+zbuf(i,4)**2)
     ra(i)=zbuf(i,1)-zdx*zbuf(i,4)*znorm
     za(i)=zbuf(i,2)+zdx*zbuf(i,3)*znorm
     da(i)=-zdx
     rb(i)=zbuf(i,1)+zdx*zbuf(i,4)*znorm
     zb(i)=zbuf(i,2)-zdx*zbuf(i,3)*znorm
     db(i)=zdx
  enddo

  if(itype.eq.2) then

     !  limiter spec is known to exist, as it is a prerequisite for 
     !  setting up a scrapeoff region

     imiss=0
     do i=1,inum

        !  find point, along each theta line, to straddle limiter.

        zchiv=zchi(i)
        call xplasma_RZeval_2d(sp,ilist(1:2), &
             xplasma_theta_coord,zchiv,xplasma_rho_coord,zrhox, &
             zbuf2(1:nrhox,1:2),ierr)

        if(ierr.ne.0) return

        !  zbuf2(1:nrhox,1):  R along external theta line
        !  zbuf2(1:nrhox,2):  Z along external theta line

        !  provide an extra point on the inside
        zbuf2(0,1) = 2*zbuf2(1,1) - zbuf2(2,1)
        zbuf2(0,2) = 2*zbuf2(1,2) - zbuf2(2,2)

        call xplasma_lim_distance(sp, zbuf2(0:nrhox,1),zbuf2(0:nrhox,2), &
             zdist(0:nrhox), ierr, maptype=imaptype)

        if(ierr.ne.0) return

        !  find straddling segment-- error if not found.
        !  take first one if more than one.

        ifnd=-1
        do j=0,nrhox-1
           if((zdist(j).le.ZERO).and.(zdist(j+1).gt.ZERO)) then
              if(ifnd.eq.-1) then
                 ifnd=j
                 ra(i)=zbuf2(j,1)
                 za(i)=zbuf2(j,2)
                 da(i)=zdist(j)
                 rb(i)=zbuf2(j+1,1)
                 zb(i)=zbuf2(j+1,2)
                 db(i)=zdist(j+1)
              endif
           endif
        enddo
        if(ifnd.eq.-1) then
           !  use precomputed bdy straddle -- but correct the distances
           call xplasma_lim_distance(sp,ra(i),za(i),da(i),ierr, &
                maptype=imaptype)
           if(ierr.ne.0) return

           call xplasma_lim_distance(sp,rb(i),zb(i),db(i),ierr, &
                maptype=imaptype)
           if(ierr.ne.0) return

           !  verify the straddle

           if(da(i)*db(i).gt.ZERO) then
              imiss=imiss+1         ! record miss
           endif
        endif

     enddo                          ! chi loop

     if(imiss.gt.0) then
        call xplasma_errmsg_append(sp,' ?eqi_bprobe: algorithm failure.')
        call xplasma_errmsg_append(sp,'  failed to straddle limiter everywhere')
        call xplasma_errmsg_append(sp,'  plasma boundary might be crossing the limiter')
        ierr=9999
     endif

  endif                             ! itype

end subroutine eqi_bprobe
