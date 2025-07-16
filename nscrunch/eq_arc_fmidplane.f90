subroutine eq_arc_fmidplane(ivec,iok,th,ansr,ivecd, &
     zinput,ninput,zoutput,noutput)

  ! subroutine for comput/zridderx root finder

  ! find midplane intercept; zinput contains midplane Z value

  use eq_arc_scrunch_mod
  implicit NONE

  integer :: ivec
  logical :: iok(ivec)  ! mask
  real*8 :: th(ivec)    ! theta
  real*8 :: ansr(ivec)  ! Z(th)-Zmidplane

  integer :: ivecd
  integer :: ninput,noutput
  real*8 :: zinput(ivecd,ninput)
  real*8 :: zoutput(ivecd,noutput)  ! not used here

  !-----------------
  real*8 :: Zmidplane
  integer, dimension(3) :: ict = (/ 1, 0, 0 /)
  integer :: idum,ier,jvec,ii

  real*8 :: thbuf(ivec),abuf(ivec)
  !-----------------

  Zmidplane = zinput(1,1)

  jvec=0
  do ii=1,ivec
     if(.not.iok(ii)) then
        jvec=jvec+1
        thbuf(jvec)=th(ii)
     endif
  enddo

  call r8spvec(ict,jvec,thbuf,jvec,abuf,ntha,tha_pkg,Zspl,idum,ier)

  jvec=0
  do ii=1,ivec
     if(.not.iok(ii)) then
        jvec=jvec+1
        ansr(ii)=abuf(jvec)-Zmidplane
     endif
  enddo

end subroutine eq_arc_fmidplane
