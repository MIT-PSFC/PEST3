subroutine eqm_bbin(zbb,iauto,inum,ztol,id,ierr)

  use xplasma_obj_instance
  use eq_module

  !  set up mod(B) grid -- must be strict ascending

  IMPLICIT NONE

  !  input:

  integer inum                      ! #of grid bdy pts
  integer iauto                     ! .ge.1 for automatic grid generation
  REAL*8 zbb(inum)                  ! the B grid (zone bdys, inum-1 zones).
  REAL*8 ztol                       ! even spacing tolerance

  !  zbb is output if iauto is set .gt.0.
  !    if iauto=1, every surface is searched for Bmin,Bmax.
  !    if iauto.gt.1, every iauto'th surface is searched, starting with
  !    the outermost surface.

  !  output:

  integer id                        ! axis id code (returned)
  integer ierr                      ! =0: OK

  !-----------------

  integer i,iwarn,inrho,id_rho,iertmp
  real*8 zbmin,zbmax,zbmin1,zbmax1
  real*8, dimension(:), allocatable :: zrho
  real*8, parameter :: ZERO = 0.0d0
  !-------------------------------------------------------------
  !  error checks:

  ierr=0
  iwarn=0

  call xplasma_find_item(s,'__RHO',id_rho,ierr)
  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     write(lunerr,*) ' ?eqm_bbin: radial flux coordinate not found.'
     return
  endif

  call xplasma_grid_size(s,id_rho,inrho,ierr)
  if(ierr.ne.0) return

  allocate(zrho(inrho))
  call xplasma_grid(s,id_rho,zrho,iertmp)

  if(iauto.gt.0) then

     !  automatic generation of B grid -- search every iauto'th surface
     !  for Bmin, Bmax; create grid to span the range.

     do i=inrho,1,-iauto
        call eq_glimb(zrho(i),zbmin1,zbmax1,ierr)
        if(ierr.ne.0) then
           write(lunerr,*) ' ?eqm_bbin -- Bmin/Bmax search failed.'
           return
        endif
        if(i.eq.inrho) then
           zbmin=zbmin1
           zbmax=zbmax1
        else
           zbmin=min(zbmin,zbmin1)
           zbmax=max(zbmax,zbmax1)
        endif
     enddo
     
     if((zbmin.eq.ZERO).and.(zbmax.eq.ZERO)) then
        zbmin=0.1d0
        zbmax=0.1d0
     endif

     !  spread the range a bit...

     zbmin=0.95d0*zbmin
     zbmax=1.05d0*zbmax

     !  generate the grid

     do i=1,inum
        zbb(i)=zbmin+(i-1)*(zbmax-zbmin)/(inum-1)
     enddo

  endif

  deallocate(zrho)

  !  monotonicity check

  do i=1,inum-1
     if(zbb(i+1).le.zbb(i)) then
        iwarn=iwarn+1
        ierr=ierr+1
        if(iwarn.le.3) then
           call eq_errmsg(' ?eqm_bbin:  mod(B) grid not strict ascending.')
           write(lunerr,*) ' zbb(',i,')=',zbb(i),' zbb(',i,'+1)=',zbb(i+1)
        endif
     endif
  enddo

  if(inum.le.4) then
     ierr=ierr+1
     call eq_errmsg(' ?eqm_bbin:  less than 4 points in mod(B) grid.')
     write(lunerr,*) ' inum=',inum
  endif

  if(ierr.gt.0) then
     ierr=9999
     return
  endif

  !  OK, create grid

  call xoi_author_set(iertmp)
  call xplasma_create_grid(s,'Bgrid',xplasma_B_coord,zbb,id,ierr, &
       label='Grid of B values based on plasma Bmin & Bmax')
  call xoi_author_clear(iertmp)

  if(ierr.gt.0) return

end subroutine eqm_bbin
