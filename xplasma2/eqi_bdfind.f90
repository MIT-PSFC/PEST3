subroutine eqi_bdfind(ivec,zR,zZ,zphi,zrho_in,iskipi,iskipo,iexact, &
     zchi_out,zphi_out,zdist,ierr)

  !  **vectorized** dmc Apr 2000

  !  accurately find the point on a boundary surface which is
  !  closest to a target point (zR,zZ,zphi)

  !  to begin the search, evaluate a coarse grid of points and look
  !  for the distance minima; use a root finder to zero in exactly
  !  (to machine precision).

  !  input:

  use xplasma_definitions
  use eqi_rzbox_module
  IMPLICIT NONE

  integer ivec                      ! vector length
  REAL*8 zR(ivec),zZ(ivec),zphi(ivec) ! target point
  REAL*8 zrho_in                    ! the (1) surface to be searched

  logical :: iskipi                 ! .TRUE. to skip accurate solution for
  !  points clearly *inside* the boundary; an approximate solution is
  !  returned.

  logical :: iskipo                 ! .TRUE. to skip accurate solution for
  !  points clearly *outside* the boundary; an approximate solution is
  !  returned.
  
  logical :: iexact                 ! .TRUE. for exact solution using
  !  root finder (for points not affected by iskip* switches); .FALSE.
  !  for distance-circles approximation method (quite accurate, much
  !  faster).

  !  points that are more than about (1/2) a segment length away are
  !  deemed clearly inside or outside.

  !  output:

  REAL*8 zchi_out(ivec)             ! chi of nearest approach
  REAL*8 zphi_out(ivec)             ! phi of nearest approach

  REAL*8 zdist(ivec)                ! distance outside surface;

  integer ierr                      ! completion code, 0=OK

  !  if zdist.lt.0.0 then point is *inside* the surface;
  !  if zdist.eq.0.0 then the point is *on* the surface.

  !------------------------------------------------------------

  real*8, dimension(:), allocatable :: zchi,zrho,zchimin,zdistmin,zdistcut, &
       zsignmin,zchix,zchi1,zchi2,zRsurf,zZsurf,zchi_aug,zdwk
  real*8, dimension(:,:), allocatable :: zdista,zderv,zsign, &
       zrzbuf,zconst1,zout1
  integer, dimension(:), allocatable :: ichimin,iseg

  real*8, dimension(:), pointer :: rzrad
  integer, dimension(:), pointer :: krzrad

  logical imask(nsrch1),jmask(ivec)

  REAL*8 zchi_inc,zdr,zdz,z2seg2,zrho_extrap,zd2,zbogus,zrho_ii

  external feq_bdy2

  REAL*8 zeps,zeta

  integer i,j,nc,ifnd(ivec),iflag,ichi,iertmp,iexec,idum
  integer :: iflist(4),idchis(4)

  real*8 :: Rfastmin,Rfastmax,Zfastmin,Zfastmax
  real*8 :: ztesta,zdervp,zchip

  integer,parameter :: iR=1  ! R index
  integer,parameter :: iZ=2  ! Z index
  integer,parameter :: iRth=3  ! dR/dtheta index
  integer,parameter :: iZth=4  ! dZ/dtheta index

  real*8, parameter :: czero= 0.0d0
  real*8, parameter :: cone = 1.0d0
  real*8, parameter :: c2pi = 6.2831853071795862D+00
  real*8, parameter :: ceps10 = 1.0d-10

  integer :: ibb,id_bb
  integer, parameter :: bbtype = 97
  character*32 bbnames(3)
  data bbnames/'__BDY_CIRC_DATA','__EXTRAP_BDY_CIRC_DATA', &
       '__GEN_CIRC_DATA'/

  !------------------------------------------------------------
  !  only axisymmetric magnetic coordinates based equilibria -- for now

  allocate(zchi(nsrch1),zchimin(ivec),zdistmin(ivec),zsignmin(ivec), &
       zchix(ivec),zchi1(ivec),zchi2(ivec),ichimin(ivec),zrho(nsrch1))
  allocate(zRsurf(0:nsrch1+1),zZsurf(0:nsrch1+1),zchi_aug(0:nsrch1+1))
  allocate(zdista(nsrch1,ivec),zderv(nsrch1,ivec),zsign(nsrch1,ivec), &
       zconst1(ivec,3),zout1(ivec,2),zrzbuf(nsrch1,4))
  allocate(iseg(ivec),zdwk(ivec),zdistcut(ivec))

  ierr=0

  !  ...assuming axisymmetry here; negative zR allowed because of
  !  2d extrapolated surfaces pushed through the axis...

  !  initial search -- nsrch0 pts

  imask=.FALSE.
  jmask=.FALSE.

  zchi(1)=czero
  zchi_aug(1)=czero

  zchi_inc=c2pi/nsrch0
  zrho_ii=max(zrho_in,ceps10)
  do i=1,nsrch0
     zchi(i+1)=zchi(i)+zchi_inc
     zchi_aug(i+1)=zchi(i)+zchi_inc
  enddo

  zchi(nsrch1)=c2pi
  zchi_aug(nsrch1)=c2pi

  zchi_aug(nsrch1+1)=zchi_aug(2)+c2pi
  zchi_aug(0)=zchi(nsrch0)-c2pi

  ibb=3
  iflag=1
  if(abs(zrho_ii-cone).lt.2.0d-15) then
     iflag=-1  ! plasma bdy
     ibb=1
  endif

  call eqi_extrap_rho_bdy(zrho_extrap)
  if(zrho_extrap.gt.cone) then
     if(abs(zrho_ii-zrho_extrap).lt.2.0d-15) then
        iflag=-1  ! extrap. space bdy.
        ibb=2
     endif
  endif

  zbogus = -1.0d10  ! place holder

  call xplasma_common_ids(sp,ierr, id_R=iflist(1), id_Z=iflist(2))
  iflist(3:4)=iflist(1:2)
  idchis(1:2)=0
  idchis(3:4)=1

  call eqi_rzbox_binfo_get(zrho_ii,zrzbuf,ierr)
  if(ierr.ne.0) then
     zrho=zrho_ii   ! vector
     call xplasma_RZeval_2d(sp,iflist, &
             xplasma_rho_coord,zrho, &
             xplasma_theta_coord,zchi, &
             zrzbuf,ierr, ideriv2s=idchis)
     if(ierr.eq.0) call eqi_rzbox_binfo_set(zrho_ii,zrzbuf,idum,ierr)
  endif

  do i=1,nsrch1
     zRsurf(i)=zrzbuf(i,iR)
     zZsurf(i)=zrzbuf(i,iZ)
  enddo
  zRsurf(0)=zRsurf(nsrch1-1)
  zRsurf(nsrch1+1)=zRsurf(2)
  zZsurf(0)=zZsurf(nsrch1-1)
  zZsurf(nsrch1+1)=zZsurf(2)

  do j=1,ivec

     zconst1(j,1)=zR(j)
     zconst1(j,2)=zZ(j)
     zconst1(j,3)=zrho_ii

     ichimin(j)=0
     zdistmin(j)=1.0d30
     do i=1,nsrch1
        zdista(i,j)=(zR(j)-zRZbuf(i,iR))**2+(zZ(j)-zRZbuf(i,iZ))**2
        zderv(i,j)=(zR(j)-zRZbuf(i,iR))*zRZbuf(i,iRth)+ &
             (zZ(j)-zRZbuf(i,iZ))*zRZbuf(i,iZth)
        zsign(i,j)=sign(cone, &
             ((zR(j)-zRZbuf(i,iR))*zRZbuf(i,iZth)- &
             (zZ(j)-zRZbuf(i,iZ))*zRZbuf(i,iRth)))
        if(zdista(i,j).le.zdistmin(j)) then
           ichimin(j)=i
           zchimin(j)=zchi(i)
           zdistmin(j)=zdista(i,j)
           zsignmin(j)=zsign(i,j)
        endif
     enddo
     zdist(j)=zbogus
     zdistcut(j)=2*zdistmin(j)

     ifnd(j)=1
     if(iskipi.or.iskipo) then
        ichi=ichimin(j)
        zdr=zRsurf(ichi+1)-zRsurf(ichi-1)
        zdz=zZsurf(ichi+1)-zZsurf(ichi-1)
        z2seg2=zdr*zdr + zdz*zdz  ! square of distance subtended by 2 segments
     endif

     !  1/16 = 0.0625d0 = ((1/2)/2)**2
     !  if distance**2 is greater than (1/16)*distance**2 subtended by 2
     !  segments, it is deemed not close... this corresponds to (1/2) the
     !  length of a single segment.

     if(iskipi) then
        if(zsignmin(j).lt.0.0d0) then
           if((zdistmin(j)/z2seg2).gt.0.0625d0) then
              jmask(j)=.TRUE. ! definitely inside
              ifnd(j)=nsrch1
           endif
        endif
     endif
     
     if(iskipo) then
        if(zsignmin(j).gt.0.0d0) then
           if((zdistmin(j)/z2seg2).gt.0.0625d0) then
              jmask(j)=.TRUE. ! definitely outside
              ifnd(j)=nsrch1
           endif
        endif
     endif

  enddo

  ! set error tolerances for root finder...

  if(iexact) then
     call xplasma_RZminmax(sp,cone,ierr, &
          rmin=Rfastmin, rmax=Rfastmax, zmin=Zfastmin, zmax=Zfastmax)

     zeps=c2pi*2.5d-14
     zeta=2.5d-15*(Rfastmax-Rfastmin)*(Zfastmax-Zfastmin)

  else

     ! retrieve pretabulated circle fit data if it exists...
     ! otherwise just allocate a workspace for the data which will be computed.

     call xplasma_find_item(sp,bbnames(ibb),id_bb,iertmp,nf_noerr=.TRUE.)

     if(id_bb.eq.0) then

        call xplasma_author_set(sp,xplasma_xmhd,iertmp)

        call xplasma_create_blackBox(sp,bbnames(ibb),bbtype, &
             id_bb,iertmp, &
             iasize=nsrch1, r8asize=3*nsrch1+1, label='circle fit data')

        call xplasma_author_clear(sp,xplasma_xmhd,iertmp)

     endif

     !  retrieve pointers; data may be updated...

     call xplasma_blackBox_retrieve(sp,id_bb,iertmp, &
          ia_ptr=krzrad, r8a_ptr=rzrad)

     !  if just created-- data values are preset to zero...

     if(rzrad(nsrch1*3+1).eq.czero) rzrad(nsrch1*3+1)=zrho_ii

     if(iflag.gt.0) then
        !  not a bdy surface... 1 struct shared by all other surfaces...
        if(zrho_ii.ne.rzrad(nsrch1*3+1)) then
           krzrad=0  ! recompute: new rho value
        endif
     endif

  endif

  ! the search incorporates an overlapping end zone to accomodate
  ! imperfections in periodicity of derivatives...

10 continue
  nc=0
  do j=1,ivec
     do i=ifnd(j),nsrch0
        ztesta=abs(zderv(i,j))
        zdervp=zderv(i+1,j)
        zchip=zchi(i+1)
        iexec=0
        if(iexact) then
           if(ztesta.le.zeta) then
              if(zdista(i,j).lt.zdistmin(j)) then
                 zchimin(j)=zchi(i)
                 zdistmin(j)=zdista(i,j)
                 zsignmin(j)=zsign(i,j)
              endif
           else if((zderv(i,j).gt.czero).and.(zdervp.lt.czero)) then
              iexec=1
           endif
        else
           if((zderv(i,j).ge.czero).and.(zdervp.le.czero)) then
              iexec=1
           endif
        endif
        if(min(zdista(i,j),zdista(i+1,j)).gt.zdistcut(j)) then
           iexec=0
        endif
        if(iexec.eq.1) then
           nc=nc+1
           ifnd(j)=i+1
           iseg(j)=i
           zchi1(j)=zchi(i)
           zchi2(j)=zchip
           go to 20
        endif
     enddo
     ifnd(j)=nsrch1+1
     jmask(j)=.true.
20   continue
  enddo

  if(nc.eq.0) go to 100

  if(iexact) then

     !  root finder

     call zridderx(ivec,jmask,zchi1,zchi2, &
          zeps,zeta,feq_bdy2,zchix,ierr, &
          ivec,zconst1,3,zout1,2)

  else

     !  circles approximation

     call eqi_bdfind_circ(ivec,jmask,zR,zZ,iseg,zchix,zdwk, &
          nsrch1,zRsurf,zZsurf,zchi_aug,krzrad,rzrad,ierr)

  endif

  if(ierr.ne.0) go to 990

  do j=1,ivec
     if(.not.jmask(j)) then

        if(iexact) then
           if(zout1(j,1).lt.zdistmin(j)) then
              zchimin(j)=zchix(j)
              zdistmin(j)=zout1(j,1)
              zsignmin(j)=zout1(j,2)
           endif

        else
           zd2=zdwk(j)*zdwk(j)
           if(zd2.lt.zdistmin(j)) then
              zchimin(j)=zchix(j)
              zdistmin(j)=zd2
              zdist(j)=zdwk(j)
              if(zdist(j).lt.0.0) then
                 zsignmin(j)=-1.0d0
              else
                 zsignmin(j)=1.0d0
              endif
           endif
        endif

     endif
  enddo

  go to 10

  !  ok go with this

100 continue
  do j=1,ivec
     zphi_out(j)=zphi(j)
     zchi_out(j)=zchimin(j)
     if(iexact.or.(zdist(j).eq.zbogus)) then
        zdist(j)=zsignmin(j)*sqrt(zdistmin(j))
     endif
  enddo

990 continue
  deallocate(zchi,zrho,zchimin,ichimin,zdistmin,zsignmin,zchix,zchi1,zchi2)
  deallocate(zdista,zderv,zsign,zrzbuf,zconst1,zout1,zRsurf,zZsurf)
  deallocate(iseg,zdwk,zdistcut)

  nullify(rzrad,krzrad)

end subroutine eqi_bdfind

!-----------------------------------------------------------
subroutine feq_bdy2(iveci,imask,zchi,zanswer,ivecd,zinput,idimi,zoutput,idimo)

  use xplasma_definitions
  use eqi_rzbox_module

  !  given (rho,chi) return d dot d/dchi[(R,Z) along zrho bdy]
  !    where d = displacement vector, [(R,Z)target - (R,Z)(zrho,zchi)]

  !  the quantity is positive if points along the zrho bdy, closer to the
  !  target, are in the +chi direction.

  IMPLICIT NONE
  !  input:
  integer iveci                     ! input vector dimension

  logical imask(abs(iveci))         ! =.FALSE.:  evaluate...
  REAL*8 zchi(abs(iveci))           ! zchi=angle param. on this surface

  !  output:
  real*8 zanswer(abs(iveci))
  !  zanswer:  d dot d/dchi(R,Z)(zchi)|rho=zrho -- function value

  !  more input:
  integer ivecd                     ! auxilliary vector dimension
  integer idimi                     ! aux. input dimension
  real*8 zinput(ivecd,idimi)        ! (R(...),Z(...),rho(...))

  integer idimo                     ! aux. output dimension

  !  output:
  real*8 zoutput(ivecd,idimo)       ! (d**2(...), sign(...))

  !-------------------------------------

  real*8, dimension(:,:), allocatable :: zRZbuf
  integer ierr

  integer i,ivec,jvec,iflist(4),idchis(4),lRnum,lZnum,idtmp
  real*8 zR,zZ
  real*8, dimension(:), allocatable :: zchi_use,zrho_use

  integer,parameter :: iR=1  ! R index
  integer,parameter :: iZ=2  ! Z index
  integer,parameter :: iRth=3  ! dR/dtheta index
  integer,parameter :: iZth=4  ! dZ/dtheta index

  real*8, parameter :: cone = 1.0d0

  !-------------------------------------

  ivec=abs(iveci)

  allocate(zchi_use(ivec),zrho_use(ivec),zRZbuf(ivec,4))

  jvec=0
  do i=1,ivec
     if(.not.imask(i)) then
        jvec=jvec+1
        zrho_use(jvec)=zinput(i,3)
        zchi_use(jvec)=zchi(i)
     endif
  enddo

  call xplasma_common_ids(sp,ierr, id_R=lRnum,id_Z=lZnum)
  if((lRnum.eq.0).or.(lZnum.eq.0)) then
     zanswer=0
     zoutput=0
     go to 900   ! this shouldn't happen!
  endif

  iflist(1)=lRnum
  iflist(2)=lZnum
  iflist(3)=lRnum
  iflist(4)=lZnum
  idchis(1:2)=0
  idchis(3:4)=1

  call xplasma_RZeval_2d(sp,iflist, &
       xplasma_rho_coord,zrho_use(1:jvec), &
       xplasma_theta_coord,zchi_use(1:jvec), &
       zrzbuf(1:jvec,1:4),ierr, ideriv2s=idchis)
  if(ierr.ne.0) go to 900

  jvec=0
  do i=1,ivec
     if(.not.imask(i)) then
        jvec=jvec+1
        zR=zinput(i,1)
        zZ=zinput(i,2)
        zoutput(i,1)= (zR-zRZbuf(jvec,iR))**2 + (zZ-zRZbuf(jvec,iZ))**2
        zanswer(i)=(zR-zRZbuf(jvec,iR))*zRZbuf(jvec,iRth)+ &
             (zZ-zRZbuf(jvec,iZ))*zRZbuf(jvec,iZth)

        !  when near closest approach:  +1 if outside, -1 if inside...

        zoutput(i,2)=sign(cone, &
             ((zR-zRZbuf(jvec,iR))*zRZbuf(jvec,iZth)- &
             (zZ-zRZbuf(jvec,iZ))*zRZbuf(jvec,iRth)))

     endif
  enddo

900 continue

  deallocate(zchi_use,zrho_use,zRZbuf)

end subroutine feq_bdy2

subroutine eqi_extrap_rho_bdy(zrho_out)

  !  return max rho of extrapolated space, or -1 if none is defined

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE

  real*8, intent(out) :: zrho_out

  !------------
  integer :: id,ier
  !------------

  zrho_out = -1

  call xplasma_find_item(sp,"__RHOX",id,ier,nf_noerr=.TRUE.)

  if((ier.ne.0).or.(id.eq.0)) return

  call xplasma_grid_info(sp,id,ier, xmax=zrho_out)

end subroutine eqi_extrap_rho_bdy

