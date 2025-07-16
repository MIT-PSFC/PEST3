 REAL*8 function pythag(a,b)
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

!
!     finds sqrt(a**2+b**2) without overflow or destructive underflow
!
 REAL*8	 	 	 :: a
 REAL*8	 	 	 :: b
 REAL*8	 	 	 :: p
 REAL*8	 	 	 :: r
 REAL*8	 	 	 :: s
 REAL*8	 	 	 :: t
 REAL*8	 	 	 :: u
      p = max(abs(a),abs(b))
      if (p  ==   0.0E0_r8  ) go to 20
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t =  4.0E0_r8   + r
         if (t  ==   4.0E0_r8  ) go to 20
         s = r/t
         u =  1.0E0_r8   +  2.0E0_r8  *s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end


