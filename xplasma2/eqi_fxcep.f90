subroutine eqi_fxcep(inum,zbdy_info,imaptype,idinit, &
     ra,za,phia,da,rb,zb,phib,db,tol, &
     rx,zx,phix,istat,ierr)

  !  given a set of line segments which straddle a boundary, find the points
  !  lying on the segments which intersect the boundary.

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE

  !  input:

  integer inum                      ! number of segments

  real*8 :: zbdy_info               ! if >= 0: rho of surface to seek
                                    ! if  < 0: look for limiter

  integer :: imaptype               ! map accuracy control
  ! for computing distance: 1=slow&exact, 2=medium&accurate, 3=fast&approximate
  integer idinit                    ! =1:  da & db already initialized

  real*8 ra(inum),za(inum),phia(inum) ! end pts A of segments
  real*8 rb(inum),zb(inum),phib(inum) ! end pts B of segments

  !  input if idinit=1, output otherwise:

  real*8 da(inum)                   ! distance of A pts from bdy
  real*8 db(inum)                   ! distance of B pts from bdy

  !  input:

  real*8 tol                        ! relative accuracy tolerance

  !  approx. accuracy of intercept position:  tol*[system size]
  !  system size = 0.5*(Rmax-Rmin)

  !  output:

  real*8 rx(inum),zx(inum),phix(inum) ! the intercepts returned

  !  if (from da,db) a segment does not straddle the boundary,
  !                  for that segment rx,zx contains the nearer endpoint

  integer istat(inum)               ! segment status

  !  istat(j)=2 -- both endpoints outside bdy
  !  istat(j)=1 -- 1st endpoint outside, 2nd endpoint inside
  !  istat(j)=-1 -- 1st endpoint inside, 2nd endpoint outside
  !  istat(j)=-2 -- both endpoints inside bdy
  !    (rx(j),zx(j)) contains bdy intercept if abs(istat(j))=1.

  integer ierr                      ! completion code, 0=OK

  !  ierr only set in case of serious error
  !-----------------------------------------------------------------

  integer i,isave,iersum,itype
  logical imask(inum)
  real*8 mindab,maxdab

  real*8 zinput(inum,9),zdist(inum)
  real*8 zoutput(inum,3)
  real*8 z0(inum),z1(inum),zxx(inum)

  real*8 ztol,Rfastmax,Rfastmin,rhosurf

  real*8, parameter :: ZERO=0.0d0
  real*8, parameter :: ONE= 1.0d0

  !-----------------------------------------------------------------
  !  external subroutines -- for zriddery, root finder

  external eqi_fxcep_bdy
  external eqi_fxcep_lim

  !-----------------------------------------------------------------
  !  surface or limiter search?

  if(zbdy_info.lt.0.0d0) then
     itype=2
     rhosurf=1.0d0
  else
     itype=1
     rhosurf=zbdy_info
  endif
  
  !-----------------------------------------------------------------
  !  get end pt distances if necessary

  iersum = 0
  if(idinit.ne.1) then

     if(itype.eq.1) then
        call xplasma_bdfind(sp,inum,ra,za,ierr, &
             maptype=imaptype, rho_in=rhosurf, phi_in=phia, dist=da)
        iersum=max(ierr,iersum)
        call xplasma_bdfind(sp,inum,rb,zb,ierr, &
             maptype=imaptype, rho_in=rhosurf, phi_in=phia, dist=db)
        iersum=max(ierr,iersum)

     else
        call xplasma_lim_distance(sp,ra,za,da,ierr, maptype=imaptype)
        iersum=max(ierr,iersum)
        call xplasma_lim_distance(sp,rb,zb,db,ierr, maptype=imaptype)
        iersum=max(ierr,iersum)

     endif

  endif
  ierr=iersum
  if(ierr.gt.0) return

  !  determine status of each segment; endpt=0 is considered "outside".

  do i=1,inum
     mindab=min(da(i),db(i))
     maxdab=max(da(i),db(i))
     if(mindab.ge.ZERO) then
        istat(i)=2                  ! all outside
        imask(i)=.TRUE.
     else if(maxdab.lt.ZERO) then
        istat(i)=-2                 ! all inside
        imask(i)=.TRUE.

        !  if we get past the above two clauses, then, either da(i).lt.0 or
        !  db(i).lt.0, and the other end point is .ge.0

     else if(db(i).lt.ZERO) then
        istat(i)=1                  ! A outside, B inside
        imask(i)=.FALSE.
     else
        istat(i)=-1                 ! A outside, B inside
        imask(i)=.FALSE.
     endif

     z0(i)=ZERO
     z1(i)=ONE

     zinput(i,1)=ra(i)
     zinput(i,2)=za(i)
     zinput(i,3)=phia(i)
     zinput(i,4)=rb(i)
     zinput(i,5)=zb(i)
     zinput(i,6)=phib(i)
     zinput(i,7)=tol
     zinput(i,8)=imaptype
     zinput(i,9)=rhosurf

     if(imask(i)) then
        if(abs(da(i)).lt.abs(db(i))) then
           rx(i)=ra(i)
           zx(i)=za(i)
           phix(i)=phia(i)
        else
           rx(i)=rb(i)
           zx(i)=zb(i)
           phix(i)=phib(i)
        endif
     endif

  enddo

  !  each segment is parametrized with x=0 ~ pt A, x=1 ~ pt B
  !  call the vectorized root finder

  call xplasma_RZminmax(sp,ONE,ierr, rmin=Rfastmin, rmax=Rfastmax)
  if(ierr.ne.0) return
  ztol=tol*(Rfastmax-Rfastmin)/2

  !  note using zriddery instead of zridderx because...
  !  passed subroutines can themselves make zridderx calls!
  !  (bugfix dmc July 2003)

  if(itype.eq.1) then
     call zriddery(inum,imask,z0,z1,ZERO,ztol, &
          eqi_fxcep_bdy,zxx,ierr,inum,zinput,9,zoutput,3)
     if(ierr.ne.0) then
        call xplasma_errmsg_append(sp,' ?eqi_fxcep(1) zriddery failure.')
     endif

  else if(itype.eq.2) then
     call zriddery(inum,imask,z0,z1,ZERO,ztol, &
          eqi_fxcep_lim,zxx,ierr,inum,zinput,9,zoutput,3)
     if(ierr.ne.0) then
        call xplasma_errmsg_append(sp,' ?eqi_fxcep(2) zriddery failure.')
     endif

  else
     call xplasma_errmsg_append(sp,' ?eqi_fxcep:  invalid itype (not 1 or 2)')
     ierr=9999
  endif

  if(ierr.ne.0) return

  !  copy results

  do i=1,inum
     if(.not.imask(i)) then
        rx(i)=zoutput(i,1)
        zx(i)=zoutput(i,2)
        phix(i)=zoutput(i,3)
     endif
  enddo

end subroutine eqi_fxcep

!-----------------------------------------------------
subroutine eqi_fxcep_bdy(ivec,imask,zxx,zdist, &
     ivecd,zinput,ninput,zoutput,noutput)

  !  distance lookup routine (for root finder)
  !  search is for plasma bdy

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE

  integer ivec                      ! vector dimension
  logical imask(ivec)               ! .TRUE. to skip pts
  real*8 zxx(ivec)                   ! position parameter (in)
  real*8 zdist(ivec)                ! distance to bdy (out)

  integer ivecd
  integer ninput,noutput
  real*8 zinput(ivecd,ninput)       ! R,Z,phi data (in)
  real*8 zoutput(ivecd,noutput)     ! R,Z,phi data (out)

  !------------------------

  call eqi_fxcep_b0(1,ivec,imask,zxx,zdist, &
       ivecd,zinput,ninput,zoutput,noutput)

end subroutine eqi_fxcep_bdy
 
!-----------------------------------------------------
subroutine eqi_fxcep_lim(ivec,imask,zxx,zdist, &
     ivecd,zinput,ninput,zoutput,noutput)

  !  distance lookup routine (for root finder)
  !  search is for plasma limiter

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE

  integer ivec                      ! vector dimension
  logical imask(ivec)               ! .TRUE. to skip pts
  real*8 zxx(ivec)                   ! position parameter (in)
  real*8 zdist(ivec)                ! distance to bdy (out)

  integer ivecd
  integer ninput,noutput
  real*8 zinput(ivecd,ninput)       ! R,Z,phi data (in)
  real*8 zoutput(ivecd,noutput)     ! R,Z,phi data (out)

  !------------------------

  call eqi_fxcep_b0(2,ivec,imask,zxx,zdist, &
       ivecd,zinput,ninput,zoutput,noutput)

end subroutine eqi_fxcep_lim
 
!-----------------------------------------------------
subroutine eqi_fxcep_b0(itype,ivec,imask,zxx,zdist, &
     ivecd,zinput,ninput,zoutput,noutput)

  !  distance lookup routine (for root finder)
  !  search is for plasma bdy

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE

  integer itype                     ! =1:  bdy, =2:  limiter
  integer ivec                      ! vector dimension
  logical imask(ivec)               ! .TRUE. to skip pts
  real*8 zxx(ivec)                  ! position parameter (in)
  real*8 zdist(ivec)                ! distance to bdy (out)
  
  integer ivecd
  integer ninput,noutput
  real*8 zinput(ivecd,ninput)       ! R,Z,phi data (in)
  real*8 zoutput(ivecd,noutput)     ! R,Z,phi data (out)

  !------------------------

  integer i,jvec,ierr
  integer :: imaptype
  real*8 r(ivec),z(ivec),phi(ivec),zdistwk(ivec),zdumv(ivec),ztol
  real*8 rhotarg(ivec)
  integer idumv(ivec)

  integer, parameter :: ilun_dbg=6
  real*8, parameter :: ONE= 1.0d0
  !------------------------

  jvec=0
  do i=1,ivec
     if(.not.imask(i)) then
        jvec=jvec+1
        r(jvec)=(ONE-zxx(i))*zinput(i,1)+zxx(i)*zinput(i,4)
        z(jvec)=(ONE-zxx(i))*zinput(i,2)+zxx(i)*zinput(i,5)
        phi(jvec)=(ONE-zxx(i))*zinput(i,3)+zxx(i)*zinput(i,6)
        rhotarg(jvec)=zinput(i,9)
     endif
  enddo

  ztol=min(1.0d-5,zinput(jvec,7)) ! all zinput(:,7) values equal
  imaptype=zinput(jvec,8)+0.5d0   ! all zinput(:,8) values equal

  if(itype.eq.1) then
     call xplasma_ctrans(sp,.TRUE.,ierr, &
          r_in=r(1:jvec),z_in=z(1:jvec),tol=ztol,maptype=imaptype, &
          rho_out=zdistwk(1:jvec))
     zdistwk(1:jvec)=zdistwk(1:jvec)-rhotarg(1:jvec)
  else
     call xplasma_lim_distance(sp,r(1:jvec),z(1:jvec),zdistwk(1:jvec),ierr, &
          maptype=imaptype)
     if(ierr.ne.0) then
        write(ilun_dbg,*) &
             ' ?eqi_fxcep: unexpected error inside root finder call:'
        call xplasma_error(sp,ierr,6)
     endif
  endif

  jvec=0
  do i=1,ivec
     if(.not.imask(i)) then
        jvec=jvec+1
        zdist(i)=zdistwk(jvec)
        zoutput(i,1)=r(jvec)
        zoutput(i,2)=z(jvec)
        zoutput(i,3)=phi(jvec)
     endif
  enddo

end subroutine eqi_fxcep_b0
