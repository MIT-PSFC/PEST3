subroutine rpdims(istype,irank,idims,zxabb,iok)
  print *,'-------------------------------------------------'
  print *,' rpdims is not installed on your system'
  print *,'-------------------------------------------------'
  call exit(1)
end subroutine rpdims
subroutine rplabel(name,label,units,imulti,istype)
  character*(*) name,label,units
  print *,'-------------------------------------------------'
  print *,' rplabel is not installed on your system'
  print *,'-------------------------------------------------'
  call exit(1)
  end subroutine rplabel
subroutine t1mhdeq(time, dtime, ns_in, nt1, igot, units,the, &
       & xin, zin, rho, psi, p, q, g, toroidalFlux, plasmaCurrent, iok)
  character*(*) units
  print *,'-------------------------------------------------'
  print *,'t1mhdeq is not installed on your system'
  print *,'-------------------------------------------------'
  call exit(1)
end subroutine t1mhdeq
