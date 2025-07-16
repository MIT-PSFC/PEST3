subroutine eq_hrhochi(ivec,zrho,zchi,nlist,iflist,ictwant,ivecd,zans,ierr)

  !  evaluate 1 or more functions f(rho,chi) derivatives df/drho, df/dchi
  !    f(rho) also allowed, in which case df/dchi=0 is returned.

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  input:

  integer ivec                      ! input vector dimension
  !  input negative number for clockwise chi specification

  REAL*8 zrho(abs(ivec))                 ! argument -- where to evaluate
  REAL*8 zchi(abs(ivec))                 ! argument -- where to evaluate
  integer nlist                     ! number of functions
  integer iflist(nlist)             ! functions to evaluate

  integer ictwant(6)                ! output selector (pspline style)

  !  each element of ictwant(1:6) should be 0 or 1
  !  use 1 to indicate a desired item:

  !   ictwant(1) -- function values            ...1 if wanted, 0 otherwise
  !   ictwant(2) -- functions' df/drho values
  !   ictwant(3) -- functions' df/dchi values
  !   ictwant(4) -- functions' d2f/drho2 values
  !   ictwant(5) -- functions' d2f/dchi2 values
  !   ictwant(6) -- functions' d2f/drhodchi values

  integer ivecd                     ! output vector dimension

  !  output:

  REAL*8 zans(ivecd,*)              ! evaluation results
  integer ierr                      ! completion code (0=OK)

  !  ** let ictnum = # of non-zero elements of ictwant.
  !         (ictnum=sum(ictwant)).  Then...
  !  zans(1:ivec,1)   = value or derivative for function id iflist(1)
  !                     indicated by first non-zero value of ictwant(...)
  !  zans(1:ivec,2)   = value or derivative for function id iflist(1)
  !                     indicated by second non-zero value of ictwant(...)
  !    ...
  !  zans(1:ivec,ictnum)   = value or derivative for function id iflist(1)
  !                     indicated by last non-zero value of ictwant(...)

  !  or generally, for 1 .le. k .le. ictnum,
  !    ...
  !  zans(1:ivec,(j-1)*ictnum+k) = value or derivative for function
  !                                id iflist(j), corresponding to k'th
  !                                non-zero element of ictwant(...)

  !--------------------------------------------
  logical :: ccw_chi,rzeval
  integer :: iveca,i,id_R,id_Z,ii,igrp,inumd,itot,ideriv1f(6),ideriv2f(6)
  integer, dimension(:), allocatable :: ids,idrhos,idchis
  !--------------------------------------------

  iveca=abs(ivec)
  ccw_chi = (ivec.gt.0)

  if(ivecd.lt.iveca) then
     ierr=1
     write(lunerr,*) ' ?eq_frhochi: answer array dimension ivecd too small:'
     write(lunerr,*) '  ivec = ',abs(ivec),'; ivecd = ',ivecd
     return
  endif

  inumd=0
  if(ictwant(1).eq.1) then
     inumd=inumd+1
     ideriv1f(inumd)=0
     ideriv2f(inumd)=0
  endif
  if(ictwant(2).eq.1) then
     inumd=inumd+1
     ideriv1f(inumd)=1
     ideriv2f(inumd)=0
  endif
  if(ictwant(3).eq.1) then
     inumd=inumd+1
     ideriv1f(inumd)=0
     ideriv2f(inumd)=1
  endif
  if(ictwant(4).eq.1) then
     inumd=inumd+1
     ideriv1f(inumd)=2
     ideriv2f(inumd)=0
  endif
  if(ictwant(5).eq.1) then
     inumd=inumd+1
     ideriv1f(inumd)=0
     ideriv2f(inumd)=2
  endif
  if(ictwant(6).eq.1) then
     inumd=inumd+1
     ideriv1f(inumd)=1
     ideriv2f(inumd)=1
  endif

  itot=inumd*nlist
  allocate(ids(itot),idrhos(itot),idchis(itot))

  igrp=-inumd
  do ii=1,nlist
     igrp=igrp+inumd
     ids(igrp+1:igrp+inumd)=iflist(ii)
     idrhos(igrp+1:igrp+inumd)=ideriv1f(1:inumd)
     idchis(igrp+1:igrp+inumd)=ideriv2f(1:inumd)
  enddo

  call xplasma_common_ids(s,ierr, id_R=id_R,id_Z=id_Z)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_hrhochi: interpolation call unsuccessful, ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
     return
  endif
  rzeval=.TRUE.
  do i=1,nlist
     if((iflist(i).ne.id_R).and.(iflist(i).ne.id_Z)) then
        rzeval=.FALSE.
        exit
     endif
  enddo

  if(rzeval) then
     call xplasma_RZeval_2d(s,ids, &
          xplasma_rho_coord,zrho, xplasma_theta_coord,zchi, &
          zans(1:iveca,1:itot),ierr, &
          ccwflag2=ccw_chi, ideriv1s=idrhos, ideriv2s=idchis)
  else
     call xplasma_eval_prof(s,ids, &
          xplasma_rho_coord,zrho, xplasma_theta_coord,zchi, &
          zans(1:iveca,1:itot),ierr, &
          ccwflag2=ccw_chi, ideriv1s=idrhos, ideriv2s=idchis)
  endif

  deallocate(ids,idrhos,idchis)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_frhochi: interpolation call unsuccessful, ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_hrhochi
