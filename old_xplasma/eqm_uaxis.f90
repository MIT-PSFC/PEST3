subroutine eqm_uaxis(aname,ikin,iper,zaxis,inum,ztol,id,ierr)

  use xplasma_obj_instance
  use eq_module

!  set up user defined grid -- must be strict ascending

  IMPLICIT NONE

!  input:

  character*(*), intent(in):: aname ! name of axis (max 20 chars)
  integer, intent(in) :: ikin       ! id of "related" axis (0 if none)
  integer, intent(in) :: iper       ! =1 if a periodic dimension, 0 o.w.
                                    ! =-1 periodic, CW instead of CCW
                                    !  -1 ignored unless ikin specifies
                                    !     theta or phi coordinate

  integer, intent(in) :: inum       ! #of grid bdy pts
  REAL*8, intent(in) :: zaxis(inum) ! the grid (zone bdys, inum-1 zones).
  REAL*8, intent(in) :: ztol        ! even spacing tolerance

!  the user defined grid can be "related" to another grid, whose id is
!  specified by the argument "ikin".  When two "related" grids exist,
!  a data "rebinning" operation can be enabled.

!  output:

  integer id                        ! axis id code (returned)
  integer ierr                      ! =0: OK

!-----------------

  real*8 zaxis_use(inum)
  real*8 :: ztolx
  real*8, dimension(:), allocatable :: zaxis_prev

  integer icoord,itype,ilen,jper,iertmp,idgrid,ingrid,ix,ixdiff
  character*32 zcoord
  logical :: lccw,lper

!-------------------------------------------------------------
!  error checks:

  ierr=0
  id=0

  call xoi_author_set(ierr)
  if(ierr.ne.0) return

  !  see if coordinate is already defined...
  !  if so, do not try to redefine it (avoid "grid in use" error).

  call xplasma_gridId(s,aname,idgrid)
  if(idgrid.gt.0) then

     call xplasma_grid_size(s,idgrid,ingrid,iertmp)
     if(ingrid.eq.inum) then
        call xplasma_global_info(s,iertmp,bdytol=ztolx)
        allocate(zaxis_prev(ingrid))
        call xplasma_grid(s,idgrid,zaxis_prev,iertmp)
        if(iertmp.eq.0) then
           ixdiff=0
           do ix=1,ingrid
              if(abs(zaxis(ix)-zaxis_prev(ix)).gt.ztolx) then
                 ixdiff=ix
                 exit
              endif
           enddo

           deallocate(zaxis_prev)
           if(ixdiff.eq.0) then
              ! infer that axis is already defined.
              id = idgrid
              return
           endif
        endif
     endif
  endif

!  get coordinate corresponding to ikin argument (if non-zero)
!  or create coordinate to contain this grid (if ikin is zero)

  jper=min(1,abs(iper))

  if(ikin.eq.0) then
     icoord=0
     ilen=len(trim(aname))
     if(ilen.le.len(zcoord)-6) then
        zcoord = trim(aname)//'_COORD'
     else
        zcoord = aname
        ilen=len(zcoord)
        zcoord(ilen-5:ilen)='_COORD'
     endif

     call xplasma_create_coord(s,zcoord,(jper.ne.0),icoord,ierr)
     if(ierr.ne.0) call xplasma_error(s,ierr,lunerr)

  else

     !  find coordinate implied by "ikin" grid id

     call xplasma_get_item_info(s,ikin,ierr, itype=itype)
     if(ierr.ne.0) then
        call xplasma_error(s,ierr,lunerr)
        
     else

        if(itype.eq.xplasma_gridType) then
           call xplasma_grid_info(s,ikin,ierr, coord=icoord)
           if(ierr.ne.0) call xplasma_error(s,ierr,lunerr) ! no error expected.
        else if(itype.eq.xplasma_coordType) then
           icoord = ikin
        else
           write(lunerr,*) ' ?eqm_uaxis: argument "ikin" invalid: ',ikin
           ierr=1
        endif

        if(ierr.eq.0) then

           !  verify consistency of "iper" argument with associated coordinate

           call xplasma_coord_isPeriodic(s,icoord,lper,ierr)
           if(lper) then
              if(jper.eq.0) then
                 ierr=1
                 write(lunerr,*) &
                      ' ?eqm_uaxis: iper.eq.0 but ikin indicates periodic coordinate.'
              endif
           else
              if(jper.ne.0) then
                 ierr=1
                 write(lunerr,*) &
                      ' ?eqm_uaxis: iper.ne.0 but ikin indicates non-periodic coordinate.'
              endif
           endif

        endif
     endif
  endif

  if(ierr.eq.0) then

!  OK, see if ccw reversal is needed
!  this could be, if the axis is kin to iaxis_chi or iaxis_phi
!  first-- find ancestral kin (the one with no further kin)

     if((iper.eq.-1).and. &
          ((icoord.eq.xplasma_theta_coord).or.(icoord.eq.xplasma_phi_coord))) then

        lccw = .FALSE.

     else

        lccw = .TRUE.
        
     endif

!  OK, copy to module; call r8genxpkg

     zaxis_use = zaxis
     call eqi_evenx(zaxis_use,inum,ztol)

     call xplasma_create_grid(s,aname,icoord,zaxis_use, id,ierr, lccw)

  endif

  call xoi_author_clear(iertmp)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_uaxis: error detected.'
     call xplasma_error(s,ierr,lunerr)
  endif

  return

end subroutine eqm_uaxis
