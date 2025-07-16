subroutine trx_tlims(stime,ftime,iok)
  print *,'-------------------------------------------------'
  print *,'trx_tlims is not installed on your system'
  print *,'-------------------------------------------------'
  call exit(1)
end subroutine trx_tlims
subroutine trx_time(time,dtime,iwarn,iok)
  print *,'-------------------------------------------------'
  print *,'trx_time is not installed on your system'
  print *,'-------------------------------------------------'
  call exit(1)
end subroutine trx_time
subroutine trx_connect(path, iok)
  character*(*) path
  print *,'-------------------------------------------------'
  print *,'trx_connect is not installed on your system'
  print *,'-------------------------------------------------'
  call exit(1)
end subroutine trx_connect
subroutine trx_ready(label,iok)
  character*(*) label
  print *,'-------------------------------------------------'
  print *,'trx_ready is not installed on your system'
  print *,'-------------------------------------------------'
  call exit(1)
end subroutine trx_ready
 
subroutine trx_msgs(iunit, mes)
  character*(*) mes
   print *,'-------------------------------------------------'
  print *,'trx_msgs is not installed on your system'
  print *,'-------------------------------------------------'
  call exit(1)
end subroutine trx_msgs
 
