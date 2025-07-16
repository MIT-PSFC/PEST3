         SUBROUTINE pstalva(sl,v,u,mj)
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

      INTEGER MJ
      INTEGER I

!
! up-down asym. check ap  17.10_r8 .95
!
!     lcm (alva)
!***         dimension v(1),u(1)
!
!     multiply sl=v*conjg(u)
!
 COMPLEX*16, DIMENSION(*) 	 :: v
 COMPLEX*16, DIMENSION(*) 	 :: u
 COMPLEX*16	 	 	 :: sl
         sl= 0._r8 
         do  10  i=1,mj
!***         sl=sl+v(i)*u(i)
         sl=sl+v(i)*conjg( u(i) )
  10     continue
!
         return
         end




