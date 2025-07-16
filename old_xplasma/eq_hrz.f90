subroutine eq_hRZ(ivec,zR,zZ,nlist,iflist,ictwant,ivecd,zans,ierr)

  !  evaluate 1 or more functions f(R,Z) and derivatives

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  input:

  integer ivec                      ! input vector dimension

  REAL*8 zR(ivec)                   ! argument -- where to evaluate
  REAL*8 zZ(ivec)                   ! argument -- where to evaluate
  integer nlist                     ! number of functions
  integer iflist(nlist)             ! functions to evaluate

  integer ictwant(6)                ! output selector (pspline style)

  !  each element of ictwant(1:6) should be 0 or 1
  !  use 1 to indicate a desired item:

  !   ictwant(1) -- function values            ...1 if wanted, 0 otherwise
  !   ictwant(2) -- functions' df/dR values
  !   ictwant(3) -- functions' df/dZ values
  !   ictwant(4) -- functions' d2f/dR2 values
  !   ictwant(5) -- functions' d2f/dZ2 values
  !   ictwant(6) -- functions' d2f/dRdZ values

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
  integer :: iveca,i,id_R,id_Z,ii,igrp,inumd,itot,ideriv1f(6),ideriv2f(6)
  integer, dimension(:), allocatable :: ids,idRs,idZs
  !--------------------------------------------

  iveca=ivec

  if(ivecd.lt.iveca) then
     ierr=1
     write(lunerr,*) ' ?eq_fRZ: answer array dimension ivecd too small:'
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
  allocate(ids(itot),idRs(itot),idZs(itot))

  igrp=-inumd
  do ii=1,nlist
     igrp=igrp+inumd
     ids(igrp+1:igrp+inumd)=iflist(ii)
     idRs(igrp+1:igrp+inumd)=ideriv1f(1:inumd)
     idZs(igrp+1:igrp+inumd)=ideriv2f(1:inumd)
  enddo

  call xplasma_eval_prof(s,ids, &
       xplasma_R_coord,zR, xplasma_Z_coord,zZ, &
       zans(1:iveca,1:itot),ierr, &
       ideriv1s=idRs, ideriv2s=idZs)

  deallocate(ids,idRs,idZs)

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_fRZ: interpolation call unsuccessful, ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_hRZ
