subroutine eq_frhochi(ivec,zrho,zchi,nlist,iflist,ivecd,zans,ierr)

  !  evaluate 1 or more functions f(rho,chi) (f(rho) also allowed).
  !  see also eq_grhochi (below) to evaluate 1st derivatives

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  input:

  integer ivec                      ! input vector dimension
  !  input negative number for clockwise chi specification

  REAL*8 zrho(abs(ivec))            ! argument -- where to evaluate
  REAL*8 zchi(abs(ivec))            ! argument -- where to evaluate
  integer nlist                     ! number of functions
  integer iflist(nlist)             ! functions to evaluate

  integer ivecd                     ! output vector dimension

  !  output:

  REAL*8 zans(ivecd,nlist)          ! evaluation results
  integer ierr                      ! completion code (0=OK)

  !--------------------------------------------
  logical :: ccw_chi,rzeval
  integer :: iveca,id_R,id_Z,i
  !--------------------------------------------

  iveca=abs(ivec)
  ccw_chi = (ivec.gt.0)

  if(ivecd.lt.iveca) then
     ierr=1
     write(lunerr,*) ' ?eq_frhochi: answer array dimension ivecd too small:'
     write(lunerr,*) '  ivec = ',abs(ivec),'; ivecd = ',ivecd
     return
  endif

  call xplasma_common_ids(s,ierr, id_R=id_R,id_Z=id_Z)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_frhochi: interpolation call unsuccessful, ierr=',ierr
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
     call xplasma_RZeval_2d(s,iflist(1:nlist), &
          xplasma_rho_coord,zrho, xplasma_theta_coord,zchi, &
          zans(1:iveca,1:nlist),ierr, ccwflag2=ccw_chi)
  else
     call xplasma_eval_prof(s,iflist(1:nlist), &
          xplasma_rho_coord,zrho, xplasma_theta_coord,zchi, &
          zans(1:iveca,1:nlist),ierr, ccwflag2=ccw_chi)
  endif


  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_frhochi: interpolation call unsuccessful, ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  endif

end subroutine eq_frhochi

!===================================================================
subroutine eq_grhochi(ivec,zrho,zchi,nlist,iflist,ivecd,zans,ierr)

  !  evaluate 1 or more functions f(rho,chi) derivatives df/drho, df/dchi
  !    f(rho) also allowed, in which case df/dchi=0 is returned.

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  input:

  integer ivec                      ! input vector dimension
  REAL*8 zrho(abs(ivec))            ! argument -- where to evaluate
  REAL*8 zchi(abs(ivec))            ! argument -- where to evaluate
  integer nlist                     ! number of functions
  integer iflist(nlist)             ! functions to evaluate

  integer ivecd                     ! output vector dimension

  !  output:

  REAL*8 zans(ivecd,2,nlist)        ! evaluation results
  integer ierr                      ! completion code (0=OK)

  !  zans(1:ivec,1,j) = df/drho for function id iflist(j)
  !  zans(1:ivec,2,j) = df/dchi for function id iflist(j)

  !--------------------------------------------
  logical :: ccw_chi,rzeval
  integer :: iveca,i,id_R,id_Z
  integer, dimension(:), allocatable :: iflist2,iderivs1,iderivs2
  real*8,dimension(:,:), allocatable :: zabuf
  !--------------------------------------------

  iveca=abs(ivec)
  ccw_chi = (ivec.gt.0)

  if(ivecd.lt.iveca) then
     ierr=1
     write(lunerr,*) ' ?eq_frhochi: answer array dimension ivecd too small:'
     write(lunerr,*) '  ivec = ',abs(ivec),'; ivecd = ',ivecd
     return
  endif

  allocate(iflist2(2*nlist),zabuf(1:iveca,2*nlist))
  allocate(iderivs1(2*nlist),iderivs2(2*nlist))

  iflist2(1:nlist)=iflist(1:nlist)
  iflist2(nlist+1:2*nlist)=iflist(1:nlist)

  iderivs1(1:nlist)=1
  iderivs1(nlist+1:2*nlist)=0

  iderivs2(1:nlist)=0
  iderivs2(nlist+1:2*nlist)=1

  call xplasma_common_ids(s,ierr, id_R=id_R,id_Z=id_Z)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_grhochi: interpolation call unsuccessful, ierr=',ierr
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
     call xplasma_RZeval_2d(s,iflist2, &
       xplasma_rho_coord,zrho, xplasma_theta_coord,zchi, &
       zabuf,ierr, ccwflag2=ccw_chi, &
       ideriv1s=iderivs1, ideriv2s=iderivs2)

  else
     call xplasma_eval_prof(s,iflist2, &
          xplasma_rho_coord,zrho, xplasma_theta_coord,zchi, &
          zabuf,ierr, ccwflag2=ccw_chi, &
          ideriv1s=iderivs1, ideriv2s=iderivs2)
  endif

  if(ierr.ne.0) then
     write(lunerr,*) ' ?eq_frhochi: interpolation call unsuccessful, ierr=',ierr
     call xplasma_error(s,ierr,lunerr)
  else
     do i=1,nlist
        zans(1:iveca,1,i)=zabuf(1:iveca,i)
        zans(1:iveca,2,i)=zabuf(1:iveca,i+nlist)
     enddo
  endif

end subroutine eq_grhochi
