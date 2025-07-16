      SUBROUTINE scraborter(iu,message)
      implicit none
      integer iu
      character*(*) message
 
      write(iu,*) '  nscrunch(scraborter) abnormal exit message:'
      write(iu,*) trim(message)
 
      call bad_exit
 
      END SUBROUTINE scraborter
 
      SUBROUTINE scrabwarn(iu,message)
      implicit none
      integer iu
      character*(*) message
 
      write(iu,*) '  nscrunch(scrabwarn) warning message:'
      write(iu,*) trim(message)
  
    END SUBROUTINE scrabwarn
 
