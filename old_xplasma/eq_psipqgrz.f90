subroutine eq_psipqgrz(fpath,ixtra, &
     id_pmhd,id_q, &
     ierr)
 
  !  write J. Menard - style - psipqgRZ file
  !  the radial grid used is the xplasma grid, unless ixtra specifies
  !  extra radial points on a per-zone basis
 
  !  dmc Jan 2002 -- xplasma module not used; using data access routines.
 
  use ezcdf
  implicit NONE
 
  character*(*), intent(in) :: fpath    ! ezcdf filename
 
  integer, intent(in) :: ixtra    ! extra points per zone (can be zero)
 
  integer, intent(in) :: id_pmhd  ! id for pressure profile
  integer, intent(in) :: id_q     ! id for q profile
 
  integer, intent(out) :: ierr    ! completion code, 0=OK
 
  !-------------------
 
  real*8, dimension(:), allocatable :: rhotmp,phitmp
  real*8, dimension(:), allocatable :: rho,chi,psi,pmhd,g,q
  real*8, dimension(:,:), allocatable :: r,z
  integer, dimension(:), allocatable :: nregion
 
  integer :: id_psi=0,id_g=0,id_rho=0,id_chi=0
  integer isnccwb,isnccwi
 
  integer lunerr
 
  integer nrho,nchi,igot,i,j,k,istat
  integer dimlens(3),ncid
 
  real*8 zdchi,ztest
 
  real*8, parameter :: z2pi = 6.2831853071795862D+00
 
  character*20 errstr
 
  !-------------------
 
  ierr=0
 
  call eq_get_lunerr(lunerr)
 
  errstr = ' '
 
  !-------------------
  !  check arguments...
  !-------------------
 
  if((ixtra.lt.0).or.(ixtra.gt.100)) then
     write(lunerr,*) ' ?eq_psipqgrz: "ixtra" argument not in range 0 - 100:', &
          ixtra
     ierr=1
  endif
 
  if(id_pmhd.eq.0) then
     write(lunerr,*) ' ?eq_psipqgrz: pressure profile id (id_pmhd) is zero.'
     ierr=1
  endif
 
  if(id_q.eq.0) then
     write(lunerr,*) ' ?eq_psipqgrz: q profile id (id_q) is zero.'
     ierr=1
  endif
 
  if(ierr.ne.0) return
 
  !-------------------
  !  check xplasma...
  !-------------------
 
  call eq_gfnum('g',id_g)
  call eq_gfnum('psi',id_psi)
 
  if(id_g.eq.0) then
     write(lunerr,*) ' ??eq_psipqgrz:  "g" not in xplasma data.'
     ierr=1
  endif
 
  if(id_psi.eq.0) then
     write(lunerr,*) ' ??eq_psipqgrz:  "psi" not in xplasma data.'
     ierr=1
  endif
 
  call eq_bchk_sign(isnccwb,isnccwi)
 
  if(isnccwb.eq.0) then
     write(lunerr,*) ' ??eq_psipqgrz:  incomplete B field setup in xplasma.'
     ierr=1
  endif
 
  if(ierr.ne.0) return
 
  !--------------------------------------------------------
  ! get grids
  !--------------------------------------------------------
 
  call eq_ganum('__RHO',id_rho)
  call eq_ganum('__CHI',id_chi)
 
  call eq_ngrid(id_rho,nrho)
  call eq_ngrid(id_chi,nchi)
 
  ! allocate extra space for chi; if it is unevenly spaced, a finer grid
  ! will be used.
  ! ***after allocation: go to 999 for cleanup if there is an error.
 
  allocate(rhotmp(nrho),chi(2*nchi+1))
 
  call eq_grid(id_rho,rhotmp,nrho,igot,ierr)
  if(ierr.ne.0) go to 999
 
  call eq_grid(id_chi,chi,nchi,igot,ierr)
  if(ierr.ne.0) go to 999
 
  ! check chi spacing
 
  zdchi = z2pi/(nchi-1)
  igot = nchi
  do i=1,nchi-1
     ztest = ((chi(i+1)-chi(i))-zdchi)/zdchi
     if(abs(ztest).gt.1.0d-3) igot=2*nchi-1
  enddo
 
  ! if chi not evenly spaced, form an evenly spaced grid with more points.
 
  if(igot.gt.nchi) then
     nchi=igot
     do i=1,nchi
        chi(i)=(i-1)*z2pi/(nchi-1)
     enddo
  endif
 
  ! implement ixtra feature-- add extra points to rho grid if requested.
 
  if(ixtra.eq.0) then
     allocate(rho(nrho))
     rho=rhotmp
  else
     igot=nrho
     nrho=1+(ixtra+1)*(igot-1)
     allocate(rho(nrho))
     rho(nrho)=rhotmp(igot)
     k=0
     do i=1,igot-1
        do j=1,ixtra+1
           k=k+1
           rho(k)=rhotmp(i)+(j-1)*(rhotmp(i+1)-rhotmp(i))/(ixtra+1)
        enddo
     enddo
  endif
 
  deallocate(rhotmp)
  allocate(rhotmp(nchi),phitmp(nchi))
 
  !-------------------------------------------------------
  ! get data
  !-------------------------------------------------------
 
  allocate(psi(nrho),pmhd(nrho),g(nrho),q(nrho))
  allocate(r(nchi,nrho),z(nchi,nrho))
 
  errstr = 'eq_rgetf(psi)'
  call eq_rgetf(nrho,rho,id_psi,0,psi,ierr)
  if(ierr.ne.0) go to 999
 
  errstr = 'eq_rgetf(pmhd)'
  call eq_rgetf(nrho,rho,id_pmhd,0,pmhd,ierr)
  if(ierr.ne.0) go to 999
 
  errstr = 'eq_rgetf(g)'
  call eq_rgetf(nrho,rho,id_g,0,g,ierr)
  if(ierr.ne.0) go to 999
 
  errstr = 'eq_rgetf(q)'
  call eq_rgetf(nrho,rho,id_q,0,q,ierr)
  if(ierr.ne.0) go to 999
 
  errstr = 'eq_rzget'
 
  allocate(nregion(nchi))
  phitmp=0
  do i=1,nrho
     rhotmp=rho(i)
     call eq_rzget(nchi,rhotmp,chi,phitmp,r(1:nchi,i),z(1:nchi,i),nregion, &
          ierr)
     if(ierr.ne.0) go to 999
  enddo
 
  !-------------------
  !  try to open file...
  !-------------------
 
  call ezcdf_open(ncid,fpath,'w',ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_psipqgrz: file creation failure: ', fpath
     go to 999
  endif
 
  !-------------------
  !  write file...
  !-------------------
 
  ! define variables
 
  dimlens=0
  call cdfDefVar(ncid,'isnccwb',dimlens,'INT')
  call cdfDefVar(ncid,'isnccwi',dimlens,'INT')
  call cdfDefVar(ncid,'ns',dimlens,'INT')
  call cdfDefVar(ncid,'nt1',dimlens,'INT')
 
  dimlens(1)=nrho
  call cdfDefVar(ncid,'rho',dimlens,'R8')
  call cdfDefVar(ncid,'psi',dimlens,'R8')
  call cdfDefVar(ncid,'p',dimlens,'R8')
  call cdfDefVar(ncid,'q',dimlens,'R8')
  call cdfDefVar(ncid,'g',dimlens,'R8')
 
  dimlens(1)=nchi
  dimlens(2)=nrho
  call cdfDefVar(ncid,'R',dimlens,'R8')
  call cdfDefVar(ncid,'Z',dimlens,'R8')
 
  ! put variables
 
  call cdfPutVar(ncid,'isnccwb',isnccwb)
  call cdfPutVar(ncid,'isnccwi',isnccwi)
  call cdfPutVar(ncid,'ns',nrho)
  call cdfPutVar(ncid,'nt1',nchi)
 
  call cdfPutVar(ncid,'rho',rho)
  psi = -isnccwi * psi
  call cdfPutVar(ncid,'psi',psi)
  call cdfPutVar(ncid,'p',pmhd)
  call cdfPutVar(ncid,'q',q)
  g = isnccwb * g
  call cdfPutVar(ncid,'g',g)
 
  call cdfPutVar(ncid,'R',r)
  call cdfPutVar(ncid,'Z',z)
 
  call ezcdf_close(ncid)
 
  !------------------------------------
  ! cleanup and exit
 
999 continue
 
  deallocate(rhotmp,phitmp,stat=istat)
  deallocate(rho,chi,stat=istat)
  deallocate(psi,pmhd,q,g,stat=istat)
  deallocate(nregion,stat=istat)
  deallocate(r,z,stat=istat)
 
end subroutine eq_psipqgrz
