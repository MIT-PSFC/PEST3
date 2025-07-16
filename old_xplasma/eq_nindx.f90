subroutine eq_nindx(iaxis,ivec,zval,kindx,iwarn)

  use xplasma_obj_instance
  use eq_module

  !  **vectorized**
  !  return indices kindx corresponding to passed values zval, on axis
  !  identified by iaxis

  IMPLICIT NONE

  integer iaxis                     ! xplasma axis id code
  integer ivec                      ! vector size
  REAL*8 zval(abs(ivec))                 ! (iaxis) coordinate value (input)
  integer kindx(abs(ivec))               ! zone index (output)
  integer iwarn                     ! set .gt.0 if any pts out of range

  !  the subroutine returns
  !     kindx = 0  ... if zval.lt. (1st axis value)
  !     kindx = size(axis) if rho_bdy .gt. (last axis value)
  !        but only if the axis is not periodic.  periodic values get
  !        shifted into range automatically

  !     if 1.le.kindx.lt.size(iaxis) then
  !        kindx satisfies eqbuf(lrho+kindx-1).le.zval.le.eqbuf(lrho+kindx)

  !----------------------------------------

  INTEGER i,id_axis,inum_outside,inumx,ierr,iertmp,inum
  integer, dimension(:), pointer :: ix
  real*8 :: xmin,xmax

  type(xpeval) :: xinfo

  logical :: ccwflag

  !----------------------------------------

  id_axis = iaxis
  inum = abs(ivec)
  ccwflag = ivec.gt.0

  call xplasma_x_lookup(s,id_axis,zval(1:inum),xinfo,ierr, &
       ccwflag=ccwflag, force_bounds=.TRUE., n_out_of_bounds=inum_outside)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_nindx: x lookup error:'
     call xplasma_error(s,ierr,lunerr)
     iwarn=inum+1
     kindx=0
     call xpeval_free(xinfo)
     return
  endif

  iwarn = inum_outside

  call xpeval_info(xinfo,ierr, ix=ix)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_nindx: unexpected xpeval_info error!'
     call xplasma_error(s,ierr,lunerr)
     iwarn=inum+1
     kindx=0
     nullify(ix)
     call xpeval_free(xinfo)
     return
  endif

  if(iwarn.gt.0) then
     call xplasma_grid_info(s,id_axis,iertmp, &
          size=inumx, xmin=xmin, xmax=xmax)
  endif

  do i=1,inum
     kindx(i)=ix(i)
     if(iwarn.gt.0) then
        if(zval(i).gt.xmax) then
           kindx(i)=inumx
        else if(zval(i).lt.xmin) then
           kindx(i)=0
        endif
     endif
  enddo

  nullify(ix)
  call xpeval_free(xinfo)

end subroutine eq_nindx
