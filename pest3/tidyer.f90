!
      subroutine psttidyer
!
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IS


!!$ LOGICAL	 	 	 :: lquit
!!$ COMMON /c2f90lquit/ lquit 
      !!CALL pstzcl(outmp1,is, 7000)
      CALL pstzcl(inpot,is, 7000)
      if(.not.lmarg)CALL pstzcl(inkin,is, 7000)
      CALL pstempty ( outmod )
  200 continue
      CALL pstempty(outpst)
 7000 return
      end


