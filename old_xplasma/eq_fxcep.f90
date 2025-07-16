subroutine eq_fxcep(inum,itype,idinit, &
     ra,za,phia,da,rb,zb,phib,db,tol, &
     rx,zx,phix,istat,ierr)

  !  given a set of line segments which straddle a boundary, find the points
  !  lying on the segments which intersect the boundary.

  use xplasma_obj_instance
  use eq_module

  implicit NONE

  !  input:

  integer inum                      ! number of segments
  integer itype                     ! =1:  plasma bdy; =2:  limiter
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
  integer :: maptype
  !-----------------------------------------------------------------
  
  call xplasma_max_maptype(s,maptype,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' xplasma_max_maptype call error in eq_fxcep:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  if(idinit.eq.1) then

      call xplasma_fxcep(s,itype, &
           ra,za,phia,rb,zb,phib, &
           rx,zx,phix,istat,ierr, &
           maptype=maptype, tol=tol, da_input=da, db_input=db)

  else

      call xplasma_fxcep(s,itype, &
           ra,za,phia,rb,zb,phib, &
           rx,zx,phix,istat,ierr, &
           maptype=maptype, tol=tol, da_output=da, db_output=db)

  endif

end subroutine eq_fxcep
