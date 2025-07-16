subroutine mcgrid_find(id)

  !  look for a unique mcgrid ID in current "s"
  !  return the ID which will be greater than zero.

  !  if none is found, return 0.

  !  if there are multiple mcgrids defined return -N, where N are defined.

  use xplasma_obj_instance

  !-----------------------
  integer :: ierr
  !-----------------------

  id=0
  call xoi_init(ierr)
  if(ierr.ne.0) return

  call xplasma_mcgrid_find(s,ierr, id_mcgrid1=id)
  if(ierr.ne.0) id=0

end subroutine mcgrid_find

subroutine mcgrid_define(inphi,iudsym,inth0,inznbmri,inznbmr,id)
 
  !  define a 2d irregular "MC" polar grid
  !  fewer zones near the axis, more zones towards the edge, each
  !  zone has roughly the same volume... see mcgrid_mod.f90

  implicit NONE
 
  integer, intent(in) :: inphi      ! no. of phi zones
  integer, intent(in) :: iudsym     ! =1: updown symmetry; =2: no symmetry
  integer, intent(in) :: inth0      ! no. of zones spanning theta=[0,pi] @axis
  integer, intent(in) :: inznbmri   ! no. of radial zones spanning [axis,bdy]
  integer, intent(out) :: inznbmr   ! total no. of radial zones
  !  inznbmr = inznbmri + inznbmri/3 -- extension of grid beyond plasma
  !    boundary; usage is optional.
 
  integer, intent(out) :: id        ! grid id, positive integer; 0 means error

  !----------------------

  CALL mcgrid_named_define(' ',inphi,iudsym,inth0,inznbmri,inznbmr,id)

end subroutine mcgrid_define

subroutine mcgrid_named_define(zname,inphi,iudsym,inth0,inznbmri,inznbmr,id)
 
  !  define a 2d irregular "MC" polar grid
  !  fewer zones near the axis, more zones towards the edge, each
  !  zone has roughly the same volume... see mcgrid_mod.f90

  use xplasma_obj_instance
  use eq_module
 
  implicit NONE
 
  character*(*), intent(in) :: zname  ! author name -- or blank for default

  integer, intent(in) :: inphi      ! no. of phi zones
  integer, intent(in) :: iudsym     ! =1: updown symmetry; =2: no symmetry
  integer, intent(in) :: inth0      ! no. of zones spanning theta=[0,pi] @axis
  integer, intent(in) :: inznbmri   ! no. of radial zones spanning [axis,bdy]
  integer, intent(out) :: inznbmr   ! total no. of radial zones
  !  inznbmr = inznbmri + inznbmri/3 -- extension of grid beyond plasma
  !    boundary; usage is optional.
 
  integer, intent(out) :: id        ! grid id, positive integer; 0 means error
 
  !-----------------------------------------
  integer :: ierr
  !-----------------------------------------
 
  inznbmr = inznbmri+max(3,inznbmri/3)

  if(zname.eq.' ') then
     call xoi_author_set(ierr)
  else
     call xplasma_author_set(s,zname,ierr)
  endif

  call xplasma_mcgrid_define(s,'F77_MCGRID', &
       inth0,inznbmri,id,ierr, &
       nphi = inphi, &
       udsym = (iudsym.eq.1), &
       nrow_ext_in = inznbmr-inznbmri)
 
  if(zname.eq.' ') then
     call xoi_author_clear(ierr)
  else
     call xplasma_author_clear(s,zname,ierr)
  endif

end subroutine mcgrid_named_define
 
subroutine mcgrid_getsym(id,iudsym)
 
  !  for specified mcgrid, return
  !        iudsym=2  if the grid is updown Asymmetric
  !        iudsym=1  if the grid is updown symmetric
 
  use xplasma_obj_instance
  use eq_module
 
  implicit NONE
 
  integer, intent(in) :: id      ! mcgrid id
  integer, intent(out) :: iudsym
 
  !-----------------------------------------------------
  logical :: klsym
  integer :: ierr
  !-----------------------------------------------------

  call xplasma_mcgrid_info(s,id,ierr, udsym = klsym)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getsym error:'
     call xplasma_error(s,ierr,lunerr)
     iudsym=0
  else
     if(klsym) then
        iudsym=1
     else
        iudsym=2
     endif
  endif
 
end subroutine mcgrid_getsym
 
subroutine mcgrid_getnumr(id,iznbmr,iznbmri)
 
  !  for specified mcgrid, return
  !        iznbmr= # of zone rows
  !        iznbmri= # of zone rows *inside plasma*
 
  use xplasma_obj_instance
  use eq_module
 
  implicit NONE
 
  integer, intent(in) :: id      ! mcgrid id
  integer, intent(out) :: iznbmr,iznbmri
 
  !-----------------------------------------------------
  integer :: ierr
  !-----------------------------------------------------

  call xplasma_mcgrid_info(s,id,ierr, nzrow=iznbmri, nzrow_ext=iznbmr)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getnumr error:'
     call xplasma_error(s,ierr,lunerr)
     iznbmr=0
     iznbmri=0
  endif
 
end subroutine mcgrid_getnumr
 
subroutine mcgrid_getnumz(id,ifbzns,ifbznsi)
 
  !  for specified mcgrid, return
  !        ifbzns= total # of zones
  !        ifbznsi= total # of zones *inside plasma*
 
  use xplasma_obj_instance
  use eq_module
 
  implicit NONE
 
  integer, intent(in) :: id      ! mcgrid id
  integer, intent(out) :: ifbzns,ifbznsi
 
  !-----------------------------------------------------
  integer :: ierr
  !-----------------------------------------------------

  call xplasma_mcgrid_info(s,id,ierr, nzons=ifbznsi, nzons_ext=ifbzns)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getnumz error:'
     call xplasma_error(s,ierr,lunerr)
     ifbzns=0
     ifbznsi=0
  endif
  
end subroutine mcgrid_getnumz
 
subroutine mcgrid_getnuma(id,inth0,inphi)
 
  !  for specified mcgrid, return
  !        inth0 = # of theta zones covering [0,pi] in first row at axis
  !        inphi = # of phi zones (usually: 1)
 
  use xplasma_obj_instance
  use eq_module
 
  implicit NONE
 
  integer, intent(in) :: id      ! mcgrid id
  integer, intent(out) :: inth0,inphi
 
  !-----------------------------------------------------
  integer :: ierr
  !-----------------------------------------------------

  call xplasma_mcgrid_info(s,id,ierr, nth0=inth0, nphi=inphi)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getnuma error:'
     call xplasma_error(s,ierr,lunerr)
     inth0=0
     inphi=0
  endif
  
end subroutine mcgrid_getnuma
 
subroutine mcgrid_getnumzns(id,inthzns)
 
  !  id<0:  return information on zone rows inside plasma only
  !  id>0:  return information on all zone rows
  ! 
  !  for specified mcgrid, return
  !        inthzns(1) = # of theta zones in first row at axis
  !        inthzns(2) = # of theta zones in second row from axis
  !            ...
  !        inthzns(N) = # of theta zones in last row at edge of region

  !  see mcgrid_getnumr to find number of zone rows
 
  use xplasma_obj_instance
  use eq_module
 
  implicit NONE
 
  integer, intent(in) :: id      ! mcgrid id
  integer, intent(out) :: inthzns(*)
 
  ! -----------------------------------------------------
  integer inum,inumx,ida,ierr
  ! -----------------------------------------------------

  ida=abs(id)

  call mcgrid_getnumr(ida,inumx,inum)
  if(inum.eq.0) then
     inthzns(1)=0
     return
  endif

  if(id.eq.ida) then
     ! all zone rows
     call xplasma_mcgrid_info(s,ida,ierr, nths=inthzns(1:inumx))
  else
     ! interior zone rows only
     call xplasma_mcgrid_info(s,ida,ierr, nths=inthzns(1:inum))
  endif

  if(ierr.ne.0) then
     write(lunerr,*) ' ?mcgrid_getnumzns error:'
     call xplasma_error(s,ierr,lunerr)
     if(id.eq.ida) then
        inthzns(1:inumx)=0
     else
        inthzns(1:inum)=0
     endif
  endif
 
end subroutine mcgrid_getnumzns
