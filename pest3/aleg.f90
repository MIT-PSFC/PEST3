      SUBROUTINE pstaleg(x,nloc,pm,pn,pp, aleg0,aleg1 )
!
!      subroutine to calculate half integral legendre functions.
!      uses upwards recurrence relations starting from elliptic
!      integrals,evaluated from hasting's formulae.
!      these expressions are very bad for large values of nloc.
!
      IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

      real*8 X
      INTEGER NLOC
      real*8 PM
      real*8 PN
      real*8 PP
      real*8 ALEG0
      real*8 ALEG1
      real*8 PYE
      INTEGER INIT
      real*8 AK0
      real*8 AK1
      real*8 AK2
      real*8 AK3
      real*8 AK4
      real*8 BK0
      real*8 BK1
      real*8 BK2
      real*8 BK3
      real*8 BK4
      real*8 AE1
      real*8 AE2
      real*8 AE3
      real*8 AE4
      real*8 BE1
      real*8 BE2
      real*8 BE3
      real*8 BE4
      real*8 PII
      real*8 SQPI
      real*8 GAMG
      real*8 XXQ
      real*8 S
      real*8 YSQ
      real*8 Y
      real*8 W
      real*8 V
      real*8 X1
      real*8 X2
      real*8 X3
      real*8 X4
      real*8 ELIPK
      real*8 ELIPE
      INTEGER KLOC
      real*8 AK
      real*8 AK02

      data pye/ 3.141592653589793_r8 / , init/0/, &
     ak0 / 1.38629436112_r8 / , &
     ak1 / 0.09666344259_r8 / , &
     ak2 / 0.03590092383_r8 / , &
     ak3 / 0.03742563713_r8 / , &
     ak4 / 0.01451196212_r8 / , &
     bk0 / 0.5_r8 / , &
     bk1 / 0.12498593597_r8 / , &
     bk2 / 0.06880248576_r8 / , &
     bk3 / 0.03328355346_r8 / , &
     bk4 / 0.00441787012_r8 / , &
     ae1 / 0.44325141463_r8 / , &
     ae2 / 0.0626060122_r8 / , &
     ae3 / 0.04757383546_r8 / , &
     ae4 / 0.01736506451_r8 / , &
     be1 / 0.2499836831_r8 / , &
     be2 / 0.09200180037_r8 / , &
     be3 / 0.04069697526_r8 / , &
     be4 / 0.00526449639_r8 / 
!
!!$      if ( init  >   0 )  go to 110
      pii =  2.0_r8  / pye
      sqpi = sqrt ( pye )
      init = init + 1
!!$  110 continue
!
      gamg = sqpi
      xxq = x*x
      s = (x- 1.0_r8 )/(x+ 1.0_r8 )
      ysq = xxq- 1.0_r8 
      y = sqrt(ysq)
      w = x+y
      v =  2.0_r8 *y/w
      x1 =  2.0_r8  / (x+ 1.0_r8 )
      x2 = x1*x1
      x3 = x2*x1
      x4 = x3*x1
!
      elipk = ak0+ak1*x1+ak2*x2+ak3*x3+ak4*x4 &
              - (bk0+bk1*x1+bk2*x2+bk3*x3+bk4*x4)*log(x1)
!
      pn = pii*sqrt( 2.0_r8 /(x+ 1.0_r8 ))*elipk
!
      aleg0 = pn
!
      x1 =  1.0_r8  / w**2
      x2=x1*x1
      x3 = x2*x1
      x4 = x3*x1
!
      elipe= 1.0_r8 
      if(abs(x1)  >    1.0E-6_r8  ) &
  elipe= 1.0_r8 +ae1*x1+ae2*x2+ae3*x3+ae4*x4 &
          - (be1*x1+be2*x2+be3*x3+be4*x4)*log(x1)
!
      pp = (pii*sqrt(w)*elipe-x*pn)/( 2.0_r8 *y)
!
      aleg1 = pp
!
      kloc=0
      ak =  0.0_r8 
      if ( nloc  ==  0 )  go to 10
    5 kloc=kloc+1
      ak = float(kloc)
      ak02 =  0.5_r8  - ak
      pm = pn
      pn = pp
      pp = - 2.0_r8 *ak*x*pn/y - ak02*ak02*pm
      gamg = gamg / ak02
      if ( kloc  /=  nloc )  go to 5
   10 continue
!
      return
      end


