!define _TEST
  subroutine cdump( cvar, cval, istat )
  !.. cdump                            Dump Char string in Octal ( desperate.. )
  !                                    codesys/source/portlib/cdump.
  !
  !
  !  Author                            Jim.Conboy@ccfe.ac.uk
  !  Version                           1.00,  06Jan2012
  !  Modifications
  !
  !  1.00  06Jan2012                   Jim.Conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  character(len=*), intent(in)       :: cvar         & ! variable to dump
                                       ,cval           !  & its value
  integer, intent(out)               :: istat          ! return status, 0=OK
  !
  !..local
  integer, parameter                 :: ilun=6       & ! io unit
                                       ,ncln=30        ! nchar / line
  integer                            :: i,j            ! loop index 
  integer                            :: l,k            ! length of cval
  integer,dimension(ncln)            :: ival           ! integer value
  !
  character(len=*),parameter         :: cSr='portlib/cdump'
  !----------------------------------::================!.........................|
  !
  istat = 0
  l = min(len(cval),len_trim(cval)+1)
  write(ilun,'(4a,i9)')  cSr,' -I- of ',cvar,' 1:',l
  l_ln: do i=0,l-1,ncln
     k = min(ncln,l-i)
     l_char: do j=1,k
        ival(j) = ichar( cval(i+j:i+j) )
     enddo l_char
     write(ilun,100) i+1,i+k,(cval(i+j:i+j),j=1,k)
     write(ilun,101)         (ival(j)  ,j=1,k)
  enddo l_ln
  !
100 format(2i6,t15,30(a1,2x))
101 format(    t15,30(Z2,1x))
  !
  end subroutine cdump
#ifdef _TEST
  program test_cdump
  !.. test_cdump                       as it says on the can
  !                                    codesys/source/portlib/test_cdump.
  !  Use
  !            >cpp --traditional $CS/source/portlib/cdump.f90 ./cdump.F
  !            >lf95 --nfix cdump.F -o cdump
  !            >./cdump
  !
  !  Author                            <author>@<place>
  !  Version                           1.00,  06Jan2012
  !  Modifications
  !
  !  1.00  06Jan2012                   <author>@<place>
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  !use..
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..local
  integer                            :: i              ! loop index 
  character(len=10)                  :: cval
  !
  character(len=*),parameter         :: cSr='<dir>/test_cdump'
  !----------------------------------::================!.........................|
  !
  call cdump( 'Test1', 'abcde/fghij/klmno/pqrst/uvwxy/z    ', i )
  !
  call cdump( 'Test2', '0123456789 !"£$%^&*()_+', i )
  !
  cval = 'ABCDE'
  cval(3:3) = char(0)
  call cdump( 'Test3', cval, i )
  !
  end program test_cdump

#endif
