subroutine scrsetlun(ilun)
!
! set lun for messages
!
  use scrunch_inc1
 
  lunmsg = ilun
 
  return
 
end subroutine scrsetlun
 
subroutine scrgetlun(ilun)
!
! set lun for messages
!
  use scrunch_inc1
 
  ilun = lunmsg
 
  return
 
end subroutine scrgetlun
