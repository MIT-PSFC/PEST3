      REAL*8 function r8_shafnl(Delprime)
! Returns the nonlinear Shafranov shift derivative "ShafNL",
! given the linear Shafranov shift derivative "Delprime".
! ShafNL is related to Delprime by the equation:
!	ShafNL = Delprime * (1-ShafNL**2)**(3/2)
! This insures that -1 <= ShafNL <= 1.  Delprime may lie outside this
! interval, in which case the flux surfaces would nonsensically overlap
! if the nonlinear correction were not used.  SMP was used to solve the
! above cubic equation.
 
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      REAL*8 delprime,b,xprs
!============
      real*8 :: third = (1.D0/3.D0)
 
      logical shafnl_flag
      data shafnl_flag /.true./ ! .false. = always use linear shift!
 
! For extremely small values of Delprime, the linear formula is perfectly
! O.K., while the nonlinear formula suffers from round-off errors.
 
      if(abs(Delprime) .le. 0.001D0 .or. .not. shafnl_flag) then
      r8_shafnl=Delprime
      return
      endif
 
      b=abs(Delprime)
      xprs=(-27.0D0/b**2+(729.0D0/b**4+108.0D0/b**6)**0.5D0)**third
      r8_shafnl=1.0D0-2.0D0**(1.D0/3.D0)/b**2/xprs + xprs/3/2**third
 
      r8_shafnl=Sqrt(r8_shafnl)*SIGN (1.0D0, Delprime)
      return
      end
 
 
