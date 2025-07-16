!
!
   subroutine pstsecond(pt,jyear,jmonth,jday,jhour,jmin)
   IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 PT
      INTEGER JYEAR
      INTEGER JMONTH
      INTEGER JDAY
      INTEGER JHOUR
      INTEGER JMIN
      INTEGER JSEC
      INTEGER JMS

!
!  give time in seconds
!
 INTEGER, dimension(8) :: it
     character*10 :: date, time, zone


 print*,'CALL date_and_time'
 CALL date_and_time(date,time,zone,it)
     jyear = it(1)
     jmonth= it(2)
     jday  = it(3)    
     jhour = it(5)
     jmin  = it(6)
     jsec  = it(7)
     jms   = it(8)


     pt = jhour* 3600._r8  + jmin* 60._r8  + jsec + jms/ 1000._r8 

    return
   end subroutine pstsecond
      subroutine pstsetlnk ( i )
!----------------------------
 USE pstcom
 USE comggf

!
      if ( i  >  1 )  go to 10
      return
   10 continue
      open( unit=itty, file='pest3.log', form='formatted')
      open( unit=outmod, file='pest3_l.log', form='formatted')
      open( unit=outilc, file='p34ilc', form='formatted')
      open( unit=outdel, file='p34del', form='formatted')
      open( unit=outequ, file='p34equ', form='formatted')
!
      return
      end


