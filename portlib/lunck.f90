  subroutine lunck( lunmsg, lun, lun2 )
  !.. lunck                            print lun / file association 
  !                                    codesys/source/portlib/lunck.
  !
  !
  !  Author                            Jim.Conboy@ccfe.ac.uk
  !  Version                           1.00,  11Nov2011
  !  Modifications
  !
  !  1.00  11Nov2011                   Jim.Conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  use logmod
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  integer, intent(in)                :: lunmsg      &  ! write to , -1 => infolog
                                       ,lun            ! unit to examine
  integer, intent(in), optional      :: lun2           ! .. lun to lun2, if present
  !
  !..local
  integer                            :: ilun        &  ! loop index 
                                       ,ilun1,ilun2    !      limits
  logical                            :: connected      ! obvious
  character(len=512)                 :: cbuf           ! output
  !
  character(len=*),parameter         :: cSr='portlib/lunck'
  !----------------------------------::================!.........................|
  !
  ilun1 = max(0, lun )
  if( present(lun2) )                  then
     ilun2 = max( lun2, ilun1 )
                                       else
     ilun2 = ilun1
                                       endif
!
  write(cbuf,'(i4,a,i4)') ilun1,' to ',ilun2
  if( lunmsg > -1 )                    then
     write(lunmsg,*) cSr,trim(cbuf)
                                       else
     call infolog( cSr, trim(cbuf) )
                                       endif
!
  do ilun=ilun1,ilun2
     inquire( ilun, opened=connected )
     if( connected )                   then
       write( cbuf,'(1x,a,i4,a)') 'Lun ',ilun,':'
       inquire(ilun, name=cbuf(11:) )
       if( lunmsg > -1 )            then
         write(lunmsg, * )  trim(cbuf)
                                    else
         call infolog( trim(cbuf) )
                                    endif
                                        endif
  enddo
!
  end subroutine lunck
