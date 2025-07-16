subroutine eqm_rzfun(ifun,lamda,iorder,zsm,ierr)

  !  set up function "ifun" on (R,Z) grid
  !         "ifun" is already in form f(rho) or f(rho,theta).

  !  recommend:  Akima-Hermite interpolation (iorder=1)

  !  lamda -- scrapeoff width parameter, exp(-d/lamda) decay beyond the
  !  last flux surface; if lamda.le.0.0, use flat extrapolation instead

  !  pass the function number (ifun):  e.g.
  !     call eq_gfnum('TE',ifun)       ! get TE function number
  !     call eqm_rzfun(ifun,czero,iorder,zsm,ierr) ! flat extrapolation for TE

  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  input:
  integer ifun                      ! function number (existing f(rho))
  real*8 lamda                      ! scrape-off parameter (m)

  integer iorder                    ! order of interpolating function

  !  iorder=0  -- piecewise bilinear function
  !  iorder=1  -- piecewise bicubic hermite function
  !  iorder=2  -- piecewise bicubic spline function

  !  iorder=99 -- redefine existing f(R,Z) function -- fit order as before.

  !  iorder=100 -- bilinear data (may be replacement of same sized object)
  !  iorder=101 -- Akima Hermite (may be replacement of same sized object)
  !  iorder=102 or 103 -- spline (may be replacement of same sized object)

  real*8 zsm                        ! smoothing parameter for edge (m)

  !  output:
  integer ierr                      ! completion code, 0=OK

  !---------------------------------------------
  character*32 pname,units
  character*80 label
  integer :: itype,iorderi,iertmp,idnew
  !---------------------------------------------

  call xplasma_prof_info(s,ifun,ierr, name=pname, gridType=itype, &
       label=label, units=units)

  if(label.eq.' ') label=pname

  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     ierr=1
     return
  endif

  if(itype.ne.xplasma_flux_coords) then
     write(lunerr,*) ' ?eqm_rzfun: cannot do scrapeoff region "lamda" extrapolation on: ',trim(pname)
     write(lunerr,*) '  Profile is not defined vs. flux coordinates.'
     ierr=1
     return
  endif

  pname = trim(pname)//'_FRZ'   ! name for f(R,Z) derived from f(rho)...

  iorderi=iorder
  if(iorderi.lt.99) iorderi=max(0,min(2,iorderi))
  if(iorderi.ge.100) iorderi=max(0,min(2,(iorderi-100)))

  !-----------------
  !  iorderi=99 means request to use prior splining method for this call;
  !  this implies profile must already exist:  check...

  if(iorderi.eq.99) then
     call xplasma_find_item(s,pname,ifun,iertmp)
     if(iertmp.eq.0) then
        call xplasma_prof_info(s,ifun,ierr, splineType=iorderi)
     else
        ierr=iertmp
     endif
     if(ierr.ne.0) then
        call xplasma_error(s,ierr,lunerr)
        call eq_errmsg(' ?eqm_rhofun: iorder=99 call failed: '//trim(pname)// &
             ' does not exist.')
        ierr=1
     endif
  endif
  if(ierr.ne.0) return

  !-----------------
  !  OK...

  call xoi_author_set(iertmp)

  call xplasma_rzprof(s,pname,idnew,ierr, &
       ispline=iorderi, sm_edge=zsm, id_fun_in=ifun, lamda=lamda, &
       label=trim(label)//' mapped to (R,Z).', units=units)

  call xoi_author_clear(iertmp)

  if(ierr.ne.0) then
     call xplasma_error(s,ierr,lunerr)
     ierr=1
  endif

end subroutine eqm_rzfun
