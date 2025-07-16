logical function lmds_errstat(lunmsg,istat)
 
  ! return .FALSE. if status is OK
  ! return .TRUE. & write a message if there is an error.
 
  implicit NONE
 
  integer, intent(in) :: lunmsg  ! fortran LUN for message output
  integer, intent(in) :: istat   ! status code from preceding MDSplus operation
 
  character*100 zmsg
  integer isto2,idum,ilz,str_length,mds_errstr
 
  !------------------------------------
 
  isto2=istat/2
  lmds_errstat = ((isto2*2).eq.istat).and.(istat.ne.0)
  if(lmds_errstat) then
 
     zmsg='(no message)'
     idum = mds_errstr(istat,zmsg)
     ilz = str_length(zmsg)
 
     write(lunmsg,*) ' '
     write(lunmsg,*) ' MDSplus error status:  ',istat,', message:'
     write(lunmsg,*) zmsg(1:ilz)
     write(lunmsg,*) ' '
 
  else if(istat.eq.0) then
 
     write(lunmsg,*) ' (lmds_errstat warning: status code = ZERO)'
 
  endif
 
end function lmds_errstat
