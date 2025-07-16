SUBROUTINE scrScrunchMe (rmb,zmb,nmb,mpolvar,lasym,lzsm1only, &
     &                    rmbout,zmbout,nmbout,fsq)
!
!  dmc Sep 2004 --- believed to be dead code, this routine has been spiked.
!   however, improvement with scrdesreg may help in case this is ever revived.
!
!***  Driver program for DESCUR
!***  Dick Wieland / PPPL
!***
!***  SCRUNCHes a closed boundary, returning Fourier moments
!***  Cf. file DESCUR.F for comments
!***
!dob  HACKER: DOB to run descur scruncher from EFIT for each of the NW efit
!dob          surfaces
!***  Edited - PMS & JPJ
 
  Implicit none
 
!***  Usage:
!      Input:
!	rmb(1:nmb),zmb(1:nmb) : {R,Z} parameterization of contour
!       mpolvar               : 0:mpolvar-1 moments desired
!       lasym                 : TRUE if asymmetric analysis desired
!       lzsm1only             : TRUE if Zsin(m=1) term only desired for Z
!      Output:
!	rmbout(0:nmbout,2),zmbout(0:nmbout,2)
!                             : Moments returned
!       fsq                   : chi-squared fit
!***
 
  Integer nmb,nmbout
  REAL*8 rmb(nmb),zmb(nmb)
  REAL*8 rmbout(0:nmbout,2),zmbout(0:nmbout,2)
  Integer mpolvar
  REAL*8 fsq
  Logical lasym,lzsm1only
  
  Integer icon,ier
  Logical ioflag
 
  call errmsg_exit(' ??scrScrunchme: defunct routine (dmc Sep 2004)!')
  icon = 0
  ioflag = .true.
  ier = 1

!!!  CALL scrdesreg (rmb,zmb,nmb,mpolvar,icon,ioflag, &
!!!       &             fsq,lasym,lzsm1only, &
!!!       &             nmbout,rmbout,zmbout,ier)
 
  return
 
End SUBROUTINE scrScrunchMe
 
