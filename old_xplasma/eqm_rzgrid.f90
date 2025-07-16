subroutine eqm_rzgrid(zR,zZ,iautoR,iautoZ,inumR,inumZ,ztol,id_R,id_Z,ierr)

  !  set up (R,Z) grid to cover plasma and scrape-off layer
  !  *** mechanical boundary must be specified first (eqm_cbdy or eqm_tbdy) ***

  use xplasma_obj_instance
  use eq_module

  !  for R & Z... (iautoR & iautoZ referred to as "iauto" here...)

  !  if iauto=0, use the user specified grid, after order & coverage tests
  !    two rows of "guard points" are added at either end of the R and Z grids
  !  if iauto=-1 the user supplied grids are used verbatim

  !  if iauto=1, use an automatically generated, evenly spaced grid

  !  note:  on each end of each grid, two extra points are tacked on; this
  !  is for the benefit of Akima/Hermite interpolation:  avoidance of grid
  !  boundary effects.

  IMPLICIT NONE

  !  input:

  integer inumR,inumZ               ! #of grid bdy pts for R, Z
  integer iautoR,iautoZ             ! =1 to generate R, Z grids
  REAL*8 zR(inumR),zZ(inumZ)        ! the grids; output if iauto=1.
  REAL*8 ztol                       ! tolerance **defunct -- now ignored**

  integer id_R,id_Z                 ! axis id codes (returned)
  integer ierr                      ! =0: OK

  !-----------------

  integer i,iwarn,itype,iertmp
  integer ista,iinc,istaz,iincz

  real*8 zRuse(-1:inumR+2),zZuse(-1:inumZ+2)
  real*8 :: bdy_Rmin,bdy_Rmax,bdy_Zmin,bdy_Zmax

  integer :: id_lim
  !-------------------------------------------------------------

  !  error checks:

  if(min(inumR,inumZ).lt.10) then
     write(lunerr,*) ' eqm_rzgrid:  inumR=',inumR,' inumZ=',inumZ
     write(lunerr,*) '   a minimum value of  10  is expected for each grid.'
     ierr=1
     return
  endif

  call xplasma_find_item(s,'__LIMITER',id_lim,ierr, nf_noerr=.TRUE.)

  if((id_lim.le.0).and.((iautoR.gt.0).or.(iautoZ.gt.0))) then
     call eq_errmsg(' ?eqm_rzgrid:  specify limiter/wall boundary first:')
     write(lunerr,*) &
          '  ...automatic (R,Z) grid requires prior bdy specification.'
     write(lunerr,*) '  ...see e.g.:  eqm_cbdy or eqm_tbdy subroutines.'
     ierr=1
     return
  endif

  if(max(iautoR,iautoZ).gt.0) then
     call eq_bdlims(itype,bdy_Rmin,bdy_Rmax,bdy_Zmin,bdy_Zmax,iwarn)
  endif

  if(iautoR.gt.0) then
     do i=1,inumR
        zR(i)=bdy_Rmin+(i-1)*(bdy_Rmax-bdy_Rmin)/(inumR-1)
     enddo
  endif

  if(iautoZ.gt.0) then
     do i=1,inumZ
        zZ(i)=bdy_Zmin+(i-1)*(bdy_Zmax-bdy_Zmin)/(inumZ-1)
     enddo
  endif

  ierr=0
  iwarn=0

  if(ierr.gt.0) return

  zRuse(1:inumR)=zR
  zZuse(1:inumZ)=zZ

  zRuse(0)=zRuse(1)-(zRuse(2)-zRuse(1))
  zRuse(-1)=zRuse(1)-ctwo*(zRuse(2)-zRuse(1))
  zZuse(0)=zZuse(1)-(zZuse(2)-zZuse(1))
  zZuse(-1)=zZuse(1)-ctwo*(zZuse(2)-zZuse(1))

  zRuse(inumR+1)=zRuse(inumR)+(zRuse(inumR)-zRuse(inumR-1))
  zRuse(inumR+2)=zRuse(inumR)+ctwo*(zRuse(inumR)-zRuse(inumR-1))
  zZuse(inumZ+1)=zZuse(inumZ)+(zZuse(inumZ)-zZuse(inumZ-1))
  zZuse(inumZ+2)=zZuse(inumZ)+ctwo*(zZuse(inumZ)-zZuse(inumZ-1))

  if(iautoR.ge.0) then

     !  don't use small R guard points if the R value is too low...

     ista=-1
     iinc=2
     if(zRuse(0).lt.0.5*zRuse(1)) then
        ista=1
     else if(zRuse(-1).lt.0.5*zRuse(1)) then
        ista=0
     endif
  else
     !  use the user R vector verbatim

     ista=1
     iinc=0

  endif

  if(iautoZ.ge.0) then

     !  extended Z vector

     istaz=-1
     iincz=2

  else

     !  user Z vector, verbatim

     istaz=1
     iincz=0

  endif

  !  OK, create grids...

  call xoi_author_set(iertmp)

  call xplasma_create_RZgrid(s, &
       zRuse(ista:inumR+iinc),zZuse(istaz:inumZ+iincz),id_R,id_Z,ierr)

  call xoi_author_clear(iertmp)

  if(ierr.ne.0) then

     id_R=0
     id_Z=0

     write(lunerr,*) ' ?eqm_rzgrid: error creating R and/or Z grids: '
     call xplasma_error(s,ierr,lunerr)

     ierr=1
  endif

end subroutine eqm_rzgrid
