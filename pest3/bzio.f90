MODULE bzio
  INTEGER ::  NREAL
END MODULE bzio
!...........................................................................
!     source code for a zio emulator. uses fortran input/output extensions 
!     available under unicos. see unicos i/o technical note sec.  4.4_r8  pp. 55-58
!     j. manickam  9.16_r8 .91
!     includes entries for 
!     zop
!     zcl
!     zwr
!     zrd
!.........................
!     usage notes:
!     a) the word addressable (wa) i/o routines are in:
!     libf     on the cray-2
!     libio    on cx/cea systems.
!     b) since the wa  routines of unicos identify
!     disk file by only one parameter, either a unit number, n, which 
!     implies a file called fort.n, or a filename, eg. "test.out", we
!     have to choose one of these alternatives. we choose fort.n. this 
!     means that the variable 'name' in the zop CALL pstis irrelevant.
!
!...........................................................................
!     
      SUBROUTINE pstzop(ioc,name,nsize,idisk,icode,ilab)
!     
      USE bzio
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER :: IOC, NSIZE, IDISK, ICODE, ILAB, IERR
      character(6) :: name
      LOGICAL	 :: unitok
      LOGICAL	 :: unitop
      REAL*8    :: b
!
      inquire(unit=ioc,exist=unitok,opened=unitop)
!
! 						record length determined by 
!						real test variable b to 
!						ensure portability
!
      inquire( iolength=NREAL ) b
      if(unitok  .AND.  unitop)then
!     unit number already assigned .
         write(0,'(''error , file already opened to ioc '',i3)')ioc
         stop 'ERROR: in bzio'
      endif
      open(unit=ioc,file=name,access='direct',recl=NREAL,err=99)
      return
 99   write(*,*)'error opening file', ioc
      END SUBROUTINE pstzop
!...................................................................
         SUBROUTINE pstzcl(ioc,ierr)
         IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
         INTEGER :: ioc, ierr
!...................................................................
!
	close(ioc,err=99)
	return
 99	write(*,*)'error closing file', ioc
	stop 'ERROR: in bzio'
        END SUBROUTINE pstzcl
!...................................................................
        SUBROUTINE pstzwr(ioc,a,length,nadres,lgivup,irr)
!...................................................................
!     
        USE bzio
        IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
        INTEGER IOC, LENGTH, NADRES, LGIVUP, IRR, j
!
!
        REAL*8, DIMENSION(*) 	 :: a
!
        do j=1,length
	   write(ioc,rec=(nadres-1)*NREAL+j,err=99) a(j)
        end do
!
	return
 99	write(*,*)'error writing onto unit = ',ioc, &
 	' at address = ',nadres
	stop 'ERROR: in bzio'
        END SUBROUTINE pstzwr
!...................................................................
         SUBROUTINE pstzrd(ioc,a,length,nadres,lgivup,irr)
!...................................................................
!     
 USE bzio
         IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
         INTEGER :: ioc, length,nadres,lgivup, irr,j
!
 REAL*8, DIMENSION(*) 	 :: a
!
 do j=1,length
	read(ioc,rec=(nadres-1)*NREAL+j,err=99) a(j)
 end do
!
	return
 99	write(*,*)'error reading unit = ',ioc, &
 	'at address = ',nadres
	stop 'ERROR: in bzio'
        END SUBROUTINE pstzrd
