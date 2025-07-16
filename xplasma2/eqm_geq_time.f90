subroutine eqm_geq_numt(nt)

  ! return number of time points available, if geqdsk module contains time
  ! dependent MDS+ data.  If not, nt=0 is returned.

  ! reference is to most recently read geqdsk data!

  use eqi_geq_mod
  implicit NONE

  integer, intent(out) :: nt
 
  call geq_gtime_num(geq,nt)

end subroutine eqm_geq_numt

subroutine eqm_geq_time(ztime,idim,nt)

  ! return time points, if geqdsk module contains time dependent MDS+ data.
  ! If not, ztime(:)=0.0 and nt=0 is returned.

  ! Also, if ztime(:) as specified by dimension idim is too small, then
  ! the routine returns ztime(:)=0.0 and nt=-<the number of time pts needed>.

  ! reference is to most recently read geqdsk data!

  use eqi_geq_mod
  implicit NONE

  integer, intent(in) :: idim
  real*8 :: ztime(idim)
  integer, intent(out) :: nt

  call geq_gtime(geq,ztime,idim,nt)

  if(idim.lt.nt) nt=-nt

end subroutine eqm_geq_time

subroutine eqm_geq_curtime(ztime1)

  ! return the currently selected G-eqdsk time from MDS+ data
  ! if a file was read (i.e. no MDS+, no time dependence), 
  !   ztime1=0.00 will be returned.

  use eqi_geq_mod
  implicit NONE

  real*8, intent(out) :: ztime1

  call geq_curtime(geq,ztime1)

end subroutine eqm_geq_curtime
