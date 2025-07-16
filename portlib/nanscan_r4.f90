subroutine nanscan_r4(ilun,zexpr,inum,zdata_r4,icount)

  ! look for NaNs or other floating invalids, by scanning ascii encoding

  implicit NONE

  integer, intent(in) :: ilun        ! I/O unit number
  character*(*), intent(in) :: zexpr ! expression or label identifying the data
  integer, intent(in) :: inum        ! number of data points
  real, intent(inout) :: zdata_r4(inum) ! data to be scanned
  integer, intent(inout) :: icount   ! NaN count -- incremented

  !--------------------------
  character*30 :: nanscan
  integer :: iwarn,ii
  !--------------------------

  iwarn=0
  do ii=1,inum
     nanscan=' '
     write(nanscan,*) zdata_r4(ii)
     if(index(nanscan,'.').eq.0) then

        icount = icount + 1
        zdata_r4(ii) = 0.0

        iwarn=iwarn+1
        if(iwarn.le.5) then
           write(ilun,*) ' ***NaNscan: at index ',ii,': ',trim(nanscan), &
                ' detected in ',trim(zexpr)
        else if(iwarn.eq.6) then
           write(ilun,*) ' [more NaNs were detected...]'
        endif

     endif
  enddo

end subroutine nanscan_r4
