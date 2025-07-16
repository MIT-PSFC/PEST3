subroutine xplasma_checkName_sub(zname,ier)
  !
  !  check for legal name; no module dependence...
  !
  !  ier=1 means there is an illegal character--
  !   non-trailing blank, non-alphanumeric, not "$" or "_", or,
  !   the first character is a numeric digit.
  !   or, the name is all blank
  !
  !  ier=2 means the non-blank length exceeds the limit (32 characters)
  !
  implicit NONE
  character*(*), intent(in) :: zname
  integer, intent(out) :: ier

  !----------
  ! local:

  character*32 znamu  ! uppercase version

  character*10 :: zdigits = "0123456789"
  character*28 :: zalphpp = "ABCDEFGHIJKLMNOPQRSTUVWXYZ_$"

  integer :: i,idig,ialph,ilen
  
  !----------

  ier=0

  ilen=len(trim(zname))
  if(ilen.gt.len(znamu)) then
     ier=2
     return
  endif

  if(ilen.eq.0) then
     ier=1
     return
  endif

  znamu=zname
  call uupper(znamu)

  ialph=index(zalphpp,znamu(1:1))
  if(ialph.eq.0) then
     ier=1
     return
  endif

  do i=2,ilen
     ialph=index(zalphpp,znamu(i:i))
     idig =index(zdigits,znamu(i:i))
     if(max(ialph,idig).eq.0) then
        ier=1
        exit
     endif
  enddo

end subroutine xplasma_checkName_sub
