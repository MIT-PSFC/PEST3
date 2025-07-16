subroutine eqm_rzfunda(zname,ifun,zdata,inr,inz,iorder,zsm,ierr)

  !  set up function "zname" on (R,Z) grid, Akima-Hermite interpolation
  !  data values from user supplied array zdata(1:nR,1:nZ)

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  input:
  character*(*) zname               ! user supplied name of function
  integer inR,inZ                   ! array dimensions
  real*8 zdata(inR,inZ)             ! the function data

  integer iorder                    ! order of interpolating function

  !  iorder=0  -- piecewise bilinear function
  !  iorder=1  -- piecewise bicubic hermite function
  !  iorder=2  -- piecewise bicubic spline function

  !  iorder=99 -- redefinition of existing f(R,Z) function

  !  iorder=100 -- bilinear data (may be replacement of same sized object)
  !  iorder=101 -- Akima Hermite (may be replacement of same sized object)
  !  iorder=102 or 103 -- spline (may be replacement of same sized object)

  real*8 zsm                        ! smoothing parameter for edge (m)

  !  output:
  integer ifun                      ! xplasma function id (returned)
  integer ierr                      ! completion code, 0=OK

  !---------------------------------------------
  integer :: iorderi,id_Rg,id_Zg,iertmp,nR,nZ
  !---------------------------------------------
  ierr=0

  iorderi=iorder
  if(iorderi.lt.99) iorderi=max(0,min(2,iorderi))
  if(iorderi.ge.100) iorderi=max(0,min(2,(iorderi-100)))

  !  iorderi=99 gets special treatment below...

  call xplasma_find_item(s,'__Rgrid',id_Rg,iertmp)
  ierr=max(ierr,iertmp)
  call xplasma_find_item(s,'__Zgrid',id_Zg,iertmp)
  ierr=max(ierr,iertmp)

  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  call xplasma_grid_info(s,id_Rg,iertmp, size=nR)
  call xplasma_grid_info(s,id_Zg,iertmp, size=nZ)

  if(inR.lt.nR) then
     ierr=1
     write(lunerr,*) ' ?eqm_rzfunda: R dimension inR=',inR
     write(lunerr,*) '  not compatible with "__Rgrid" size: ',nR
  endif

  if(inZ.lt.nZ) then
     ierr=1
     write(lunerr,*) ' ?eqm_rzfunda: Z dimension inZ=',inZ
     write(lunerr,*) '  not compatible with "__Zgrid" size: ',nZ
  endif

  if(ierr.ne.0) return

  !-----------------
  !  iorderi=99 means request to use prior splining method for this call;
  !  this implies profile must already exist:  check...

  if(iorderi.eq.99) then
     call xplasma_find_item(s,zname,ifun,iertmp)
     if(iertmp.eq.0) then
        call xplasma_prof_info(s,ifun,ierr, splineType=iorderi)
     else
        ierr=iertmp
     endif
     if(ierr.ne.0) then
        call xplasma_error(s,ierr,lunerr)
        call eq_errmsg(' ?eqm_rhofun: iorder=99 call failed: '//trim(zname)// &
             ' does not exist.')
        ierr=1
     endif
  endif
  if(ierr.ne.0) return

  !-----------------
  !  OK...

  call xoi_author_set(iertmp)

  call xplasma_rzprof(s,zname,ifun,ierr, &
       id_Rgrid=id_Rg,id_Zgrid=id_Zg, &
       ispline=iorderi,sm_edge=zsm, &
       data=zdata)

  call xoi_author_clear(iertmp)

  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     ierr=1
  endif

end subroutine eqm_rzfunda
