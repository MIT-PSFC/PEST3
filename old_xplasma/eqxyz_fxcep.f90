subroutine eqxyz_fxcep(inum,itype,idinit, &
     xa,ya,za,da,xb,yb,zb,db,tol, &
     xx,yx,zx,istat,ierr)

  !  Cartesian interface to eq_fxcep

  implicit NONE

  !  input:

  integer inum                      ! number of segments
  integer itype                     ! =1:  plasma bdy; =2:  limiter
  integer idinit                    ! =1:  da & db already initialized

  real*8 xa(inum),ya(inum),za(inum) ! end pts A of segments
  real*8 xb(inum),yb(inum),zb(inum) ! end pts B of segments

  !  input if idinit=1, output otherwise:

  real*8 da(inum)                   ! distance of A pts from bdy
  real*8 db(inum)                   ! distance of B pts from bdy

  !  input:

  real*8 tol                        ! relative accuracy tolerance

  !  approx. accuracy of intercept position:  tol*[system size]
  !  system size = 0.5*(Rmax-Rmin)

  !  output:

  real*8 xx(inum),yx(inum),zx(inum) ! the intercepts returned

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

  real*8 ra(inum),phia(inum),rb(inum),phib(inum)
  real*8 rx(inum),phix(inum)

  !-----------------------------------------

  call eq_rcyl(inum,xa,ya,ra,phia)
  call eq_rcyl(inum,xb,yb,rb,phib)

  call eq_fxcep(inum,itype,idinit, &
       ra,za,phia,da,rb,zb,phib,db,tol, &
       rx,zx,phix,istat,ierr)

  call eq_xcyl(inum,rx,phix,xx,yx)

end subroutine eqxyz_fxcep
