subroutine equigs_write(idi_xpsi,idi_q,ilun,filepath,ierr)

  ! write the ascii equilibrium data file used as input by TORIC4
  ! mod dmc Sept 2005 -- require at least 80 radial grid points
  ! MOD DMC March 2006 -- write on normalized Psi grid not normalized Phi grid.

  implicit NONE

  integer, intent(in) :: idi_xpsi        ! id of x_psi(rho) profile
  !  if zero, construct this locally

  integer, intent(in) :: idi_q           ! id of q(rho) profile
  !  if zero, construct this locally

  integer, intent(in) :: ilun            ! i/o unit number
  character*(*), intent(in) :: filepath  ! path to file to be written
  integer, intent(out) :: ierr           ! exit status, 0=OK

  integer :: inx0,ingrid0,iper,iwarn
  integer, parameter :: inx_min = 120
  integer :: id_rho=0,id_g=0,id_q=0,ingrid,inx,kmom,iopen,idum,i,j,m
  integer :: id_xpsi=0

  real*8 :: zdum,zraxis,zr0,zrmin,zrmax,zb0

  real*8, dimension(:), allocatable :: zxgrid,zxgrid0,zxpsi,zicur,zg,zq,mombuf
  real*8, dimension(:,:), allocatable :: zrx,zzx          ! R & Z expanded
  real*8, dimension(:,:), allocatable :: zrc,zrs,zzc,zzs  ! R & Z moments

  real*8, dimension(:), allocatable :: zrho,ztheta,zphi
  integer, dimension(:), allocatable :: iregion

  integer :: inthx
  integer :: ith,ix,ibrk,im
  real*8, parameter :: cpi = 3.1415926535897931D+00

  logical, parameter :: ipatch = .TRUE.
  real*8, parameter :: zxpatch = 0.15D0
  real*8, parameter :: zxrtol = 1.0d-6

  integer :: id_locrho
  character*32 :: loc_rho_name,loc_xpsi_name
  !---------------------------

  ierr=1
  if(filepath.eq.' ') return

  ! fetch from xplasma: grid for data
  call eq_ganum('__RHO',id_rho)
  if(id_rho.eq.0) return
  call eq_ngrid(id_rho,ingrid0)
  inx0=ingrid0-1  ! xplasma native gridsize

  inx=0          ! compute grid size to use...
  iper=0
  do
     inx=inx+inx0
     iper=iper+1
     if(inx.ge.inx_min) exit
  enddo
  ingrid=inx+1

  ! also get "G" fcn id (for R*Bphi)
  call eq_gfnum('G',id_g)
  if(id_g.eq.0) return

  ! also get "Q" fcn id (for q profile)
  id_q = idi_q
  if(id_q.eq.0) call eqm_gen_q('q_equigs_write',id_q,ierr)
  if(id_q.eq.0) return

  ! and get xpsi id
  id_xpsi = idi_xpsi  ! could be zero; deal with this below...

  allocate(zxgrid0(0:inx0))
  call eq_grid(id_rho,zxgrid0,ingrid0,idum,ierr)

  allocate(zxgrid(0:inx),zicur(0:inx),zg(0:inx),zq(0:inx),zxpsi(0:inx))
  zxgrid(0)=zxgrid0(0)
  m=0
  do i=1,inx0
     do j=1,iper
        m=m+1
        zxgrid(m)=(zxgrid0(i-1)*(iper-j)+zxgrid0(i)*j)/iper
     enddo
  enddo

  if(id_xpsi.eq.0) then
     loc_xpsi_name = 'xpsi_equigs_write'
     loc_rho_name = 'rho_equigs_write'
     call eq_gfnum(loc_xpsi_name,id_xpsi)

     if(id_xpsi.eq.0) then
        !  build xpsi = sqrt(Psipol/Psipol(a)) vs. rho

        call eq_ganum(loc_rho_name,id_locrho)
        if(id_locrho.eq.0) then
           call eqm_uaxis(loc_rho_name,id_rho,0,zxgrid,ingrid,zxrtol, &
                id_locrho,ierr)
           if(ierr.ne.0) return
        endif

        call equigs_xpsi(id_locrho,id_q,loc_xpsi_name,id_xpsi,ierr)
        if(ierr.ne.0) return
     endif
  endif

  call eq_rgetf(ingrid,zxgrid,id_xpsi,0,zxpsi,ierr)
  call eq_rgetf(ingrid,zxgrid,id_g,0,zg,ierr)
  call eq_rgetf(ingrid,zxgrid,id_q,0,zq,ierr)

  zdum=1.0d-7*zxgrid(inx)
  call eq_flxint_init(0,ingrid,zxgrid,zdum,ierr)
  call eq_flxint('ITOR',1,zicur,1,ingrid,ierr)
  zicur = zicur*(4*cpi)*1.0d-7   ! "I2MEX" units, (T/m), from (A)

  call xmoments_kmom_get(kmom)
  
  !  note constraints on moment count:  5.le.kmom.le.16

  kmom=min(16,max(5,kmom))

  inthx=max(80,10*kmom)+1

  allocate(zrho(inthx),ztheta(inthx),zphi(inthx),iregion(inthx))
  allocate(zrx(inthx,0:inx),zzx(inthx,0:inx))
  allocate(zrc(0:kmom,0:inx),zrs(0:kmom,0:inx),mombuf(0:inx))
  allocate(zzc(0:kmom,0:inx),zzs(0:kmom,0:inx))

  do
     zphi = 0
     do ith=1,inthx
        ztheta(ith)=((ith-1)*2*cpi)/(inthx-1)
     enddo

     ix=inx
     zrho = zxgrid(ix)
     call eq_rzget(inthx,zrho,ztheta,zphi,zrx(1:inthx,ix),zzx(1:inthx,ix), &
          iregion,ierr)
     if(ierr.ne.0) exit

     zrmin=minval(zrx(1:inthx,ix))
     zrmax=maxval(zrx(1:inthx,ix))
     zr0=(zrmin+zrmax)/2
     zb0=zg(inx)/zr0

     zrs(0,0:inx)=0
     zzs(0,0:inx)=0
     call eqmom_rcos(zxgrid,ingrid,0,0,mombuf,iwarn); zrc(0,0:inx)=mombuf
     call eqmom_zcos(zxgrid,ingrid,0,0,mombuf,iwarn); zzc(0,0:inx)=mombuf
     do im=1,kmom
        call eqmom_rcos(zxgrid,ingrid,0,im,mombuf,iwarn); zrc(im,0:inx)=mombuf
        call eqmom_rsin(zxgrid,ingrid,0,im,mombuf,iwarn); zrs(im,0:inx)=mombuf
        call eqmom_zcos(zxgrid,ingrid,0,im,mombuf,iwarn); zzc(im,0:inx)=mombuf
        call eqmom_zsin(zxgrid,ingrid,0,im,mombuf,iwarn); zzs(im,0:inx)=mombuf
     enddo

     if(ipatch) then
        !
        !  patch the 0'th moment behavior towards the mag. axis:
        !  assure dR0/dx -->0, dZ0/dx -->0 as x -->0
        !
        !  also make moments smaller zxrtol*R0(0) go away...

        call scrunch_cleanup(zxgrid,ingrid,zrc,zrs,zzc,zzs, kmom, &
             .TRUE., zxpatch,zxrtol, ierr)

        if(ierr.ne.0) exit

     endif
     zRaxis=zrc(0,0)

     iopen = 0

     open(unit=ilun, file=trim(filepath), status='unknown', iostat=ierr)
     if(ierr.ne.0) exit

     iopen = 1

     write(ilun,'(A)')  'Major radius (central)'
     write(ilun,'(E18.9)') zr0

     write(ilun,'(A)')  'Major radius (axis)'
     write(ilun,'(E18.9)') zraxis

     write(ilun,'(A)')  'Magnetic field at major radius'
     write(ilun,'(E18.9)') zb0

     write(ilun,'(A)')  'Total toroidal current'
     write(ilun,'(E18.9)') zicur(inx)

     write(ilun,'(A)')  'Number of poloidal modes'
     write(ilun,'(I5)') kmom

     write(ilun,'(A)')  'Number of radial mesh points'
     write(ilun,'(I5)') ingrid

     write(ilun,'(A)')  'Radial mesh: sqrt(Psi/Psilim)'
     write(ilun,'(4E18.9)') (zxpsi(i), i=0,inx)

     write(ilun,'(A)')  'Fourier equilibrium coefficients'
     write(ilun,'(4E18.9)') (zrc(0,i),i=0,inx)
     write(ilun,'(4E18.9)') (zzc(0,i),i=0,inx)
     do  m=1,kmom
        write(ilun,'(4E18.9)') (zrc(m,i),i=0,inx)
        write(ilun,'(4E18.9)') (zzs(m,i),i=0,inx)
        write(ilun,'(4E18.9)') (zrs(m,i),i=0,inx)
        write(ilun,'(4E18.9)') (zzc(m,i),i=0,inx)
     enddo

     write(ilun,'(A)')  'Safety factor'
     write(ilun,'(4E18.9)') (zq(i), i=0,inx)

     write(ilun,'(A)')  'Current profile'
     write(ilun,'(4E18.9)') (zicur(i), i=0,inx)

     write(ilun,'(A)')  'R*B_phi (m*T)'
     write(ilun,'(4E18.9)') (zg(i), i=0,inx)

     write(ilun,'(A)')  'Radial mesh: sqrt(Phi/Philim)'
     write(ilun,'(4E18.9)') (zxgrid(i), i=0,inx)

     write(ilun,'(A)')  'Last surface (set to 1: TRANSP)'
     write(ilun,'(E18.9)')  zxpsi(inx)

     exit
  enddo

  if(iopen.eq.1) close(unit=ilun)
  deallocate(zxgrid,zxgrid0,zicur,zg,zq,zxpsi,mombuf)
  deallocate(zrho,ztheta,zphi,iregion)
  deallocate(zrx,zzx)
  deallocate(zrc,zrs,zzc,zzs)

end subroutine equigs_write
