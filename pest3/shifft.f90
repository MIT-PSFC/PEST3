!----------------------------------------------------------------------
!      4.1_r8 .2.1   subroutine shifft
!.......................................................................
!     shifft moves t"s up in storage.
!
      subroutine pstshifft
!.......................................................................
 USE pstcom

 USE l22com

 USE r33com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER MM1
      INTEGER MTOT
      INTEGER IS


!
      if ( m1  <=  2 )   stop ' ERROR: m1<=2 in pstshifft'
!                       ......
      mm1 = m1 - 2
!
!      number of surfaces to be moved..
!
      mtot = mdiv + 1
      do 20 surf = 1, mtot
      is = surf + mdiv
      tw(surf,:,:) = tw(is,:,:)
      ty(surf,:,:) = ty(is,:,:)
      tz(surf,:,:) = tz(is,:,:)
      td(surf,:,:) = td(is,:,:)
   20 continue
      return
      end


