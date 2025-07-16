!------------------------------------------------------
      SUBROUTINE pstdelwbb
!------------------------------------------------------
! computes the big-big contribution to the energy...
!
! a. pletzer october  91
! modified for up-down asymmetric plasas on 19/1/94
!
 USE pstcom

 USE r33com

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER I
      INTEGER ISLEFT
      INTEGER ISRIGH
      real*8 ZLEFT
      real*8 ZRIGH
      INTEGER J1
      INTEGER JBEG1
      INTEGER JEND1
      real*8 ZDPSL
      real*8 ZDPSR
      real*8 ZDPST
      INTEGER JMS


!
!
 COMPLEX*16	 	 	 :: zf0oo
 COMPLEX*16	 	 	 :: zf0oe
 COMPLEX*16	 	 	 :: zf0eo
 COMPLEX*16	 	 	 :: zf0ee
 COMPLEX*16	 	 	 :: zf1oo
 COMPLEX*16	 	 	 :: zf1oe
 COMPLEX*16	 	 	 :: zf1eo
 COMPLEX*16	 	 	 :: zf1ee
 COMPLEX*16	 	 	 :: zf2oo
 COMPLEX*16	 	 	 :: zf2oe
 COMPLEX*16	 	 	 :: zf2eo
 COMPLEX*16	 	 	 :: zf2ee
 COMPLEX*16	 	 	 :: zaoo
 COMPLEX*16	 	 	 :: zaoe
 COMPLEX*16	 	 	 :: zaeo
 COMPLEX*16	 	 	 :: zaee
 COMPLEX*16	 	 	 :: zboo
 COMPLEX*16	 	 	 :: zboe
 COMPLEX*16	 	 	 :: zbeo
 COMPLEX*16	 	 	 :: zbee
 COMPLEX*16	 	 	 :: zcoo
 COMPLEX*16	 	 	 :: zcoe
 COMPLEX*16	 	 	 :: zceo
 COMPLEX*16	 	 	 :: zcee
 INTEGER	 	 	 :: surfl
 INTEGER	 	 	 :: surfc
 INTEGER	 	 	 :: surfr
      do 20 i=1,ms
      w1l1oo(i,ms) =  0.0_r8 
      w1l1oe(i,ms) =  0.0_r8 
      w1l1eo(i,ms) =  0.0_r8 
      w1l1ee(i,ms) =  0.0_r8 
 20   continue
!
! test whether big prescribed solution vanishes at adjacent
!  rational surface.
!
      isleft = 1
      do i=1,msing(ms) + 1
         isleft = isleft + msub(i)
      end do
!
      isrigh = 1
      do i=1,msing(ms+2) - 1
         isrigh = isrigh + msub(i)
      end do
!
      zleft =  0._r8 
      zrigh =  0._r8 
      jmax1 = lmax(1) - lmin(1) + 1
      do 13 j1 = 1,jmax1
      zleft = zleft + x1frbo(isleft,j1,ms) &
               + x1frbe(isleft,j1,ms) &
               + c1frbo(isleft,j1,ms) &
               + c1frbe(isleft,j1,ms)
      zrigh = zrigh + x1frbo(isrigh,j1,ms) &
               + x1frbe(isrigh,j1,ms) &
               + c1frbo(isrigh,j1,ms) &
               + c1frbe(isrigh,j1,ms)
 13   continue
!
      if(abs(zleft) >   1.E-10_r8 ) &
  write(itty,*)'delwbb: big solution overlaps prev. rat. surf.'
      if(abs(zrigh) >   1.E-10_r8 ) &
  write(itty,*)'delwbb: big solution overlaps next rat. surf.'
!
      jbeg1 = mbeg + 1
      jend1 = msing(ms+2) - 1
!
! run from rational surface ms-1 to rational surface
! ms+ 1...
!
      do m1=jbeg1,jend1
         jmax1 = lmax(1) - lmin(1) + 1
         surfl = 1
         !
         do i=1,m1-2
            surfl = surfl + msub(i)
         end do
         !
         surfc = surfl + 1
         surfr = surfl + 2
         if(m1 == mp) surfc = surfl
         zdpsl = psinew(surfc) - psinew(surfl)
         zdpsr = psinew(surfr) - psinew(surfc)
         zdpst = zdpsl + zdpsr
         !
         do jms = 1,ms
            zf0oo =  0.0_r8 
            zf0oe =  0.0_r8 
            zf0eo =  0.0_r8 
            zf0ee =  0.0_r8 
            zf1oo =  0.0_r8 
            zf1oe =  0.0_r8 
            zf1eo =  0.0_r8 
            zf1ee =  0.0_r8 
            zf2oo =  0.0_r8 
            zf2oe =  0.0_r8 
            zf2eo =  0.0_r8 
            zf2ee =  0.0_r8 
            !
            do 31 j1 = 1,jmax1
               zf0oo = zf0oo &
                    + conjg( x1frbo(surfl,j1,jms) )*alxi1o(surfl,j1) &
                    - conjg( x1frbo(surfl,j1,jms) )*qdchio(surfl,j1) &
                    - conjg( x1fpro(surfl,j1,jms) )*ddchio(surfl,j1)
               zf0oe = zf0oe &
                    + conjg( x1frbo(surfl,j1,jms) )*alxi1e(surfl,j1) &
                    - conjg( x1frbo(surfl,j1,jms) )*qdchie(surfl,j1) &
                    - conjg( x1fpro(surfl,j1,jms) )*ddchie(surfl,j1)
               zf0eo = zf0eo &
                    + conjg( x1frbe(surfl,j1,jms) )*alxi1o(surfl,j1) &
                    - conjg( x1frbe(surfl,j1,jms) )*qdchio(surfl,j1) &
                    - conjg( x1fpre(surfl,j1,jms) )*ddchio(surfl,j1)
               zf0ee = zf0ee &
                    + conjg( x1frbe(surfl,j1,jms) )*alxi1e(surfl,j1) &
                    - conjg( x1frbe(surfl,j1,jms) )*qdchie(surfl,j1) &
                    - conjg( x1fpre(surfl,j1,jms) )*ddchie(surfl,j1)
               zf1oo = zf1oo &
                    + conjg( x1frbo(surfc,j1,jms) )*alxi1o(surfc,j1) &
                    - conjg( x1frbo(surfc,j1,jms) )*qdchio(surfc,j1) &
                    - conjg( x1fpro(surfc,j1,jms) )*ddchio(surfc,j1)
               zf1oe = zf1oe &
                    + conjg( x1frbo(surfc,j1,jms) )*alxi1e(surfc,j1) &
                    - conjg( x1frbo(surfc,j1,jms) )*qdchie(surfc,j1) &
                    - conjg( x1fpro(surfc,j1,jms) )*ddchie(surfc,j1)
               zf1eo = zf1eo &
                    + conjg( x1frbe(surfc,j1,jms) )*alxi1o(surfc,j1) &
                    - conjg( x1frbe(surfc,j1,jms) )*qdchio(surfc,j1) &
                    - conjg( x1fpre(surfc,j1,jms) )*ddchio(surfc,j1)
               zf1ee = zf1ee &
                    + conjg( x1frbe(surfc,j1,jms) )*alxi1e(surfc,j1) &
                    - conjg( x1frbe(surfc,j1,jms) )*qdchie(surfc,j1) &
                    - conjg( x1fpre(surfc,j1,jms) )*ddchie(surfc,j1)
               zf2oo = zf2oo &
                    + conjg( x1frbo(surfr,j1,jms) )*alxi1o(surfr,j1) &
                    - conjg( x1frbo(surfr,j1,jms) )*qdchio(surfr,j1) &
                    - conjg( x1fpro(surfr,j1,jms) )*ddchio(surfr,j1)
               zf2oe = zf2oe &
                    + conjg( x1frbo(surfr,j1,jms) )*alxi1e(surfr,j1) &
                    - conjg( x1frbo(surfr,j1,jms) )*qdchie(surfr,j1) &
                    - conjg( x1fpro(surfr,j1,jms) )*ddchie(surfr,j1)
               zf2eo = zf2eo &
                    + conjg( x1frbe(surfr,j1,jms) )*alxi1o(surfr,j1) &
                    - conjg( x1frbe(surfr,j1,jms) )*qdchio(surfr,j1) &
                    - conjg( x1fpre(surfr,j1,jms) )*ddchio(surfr,j1)
               zf2ee = zf2ee &
                    + conjg( x1frbe(surfr,j1,jms) )*alxi1e(surfr,j1) &
                    - conjg( x1frbe(surfr,j1,jms) )*qdchie(surfr,j1) &
                    - conjg( x1fpre(surfr,j1,jms) )*ddchie(surfr,j1)
31             continue
               !
               zaoo = (zdpsl*zf2oo - zdpst*zf1oo + zdpsr*zf0oo) &
                    /(zdpst*zdpsl*zdpsr)
               zaoe = (zdpsl*zf2oe - zdpst*zf1oe + zdpsr*zf0oe) &
                    /(zdpst*zdpsl*zdpsr)
               zaeo = (zdpsl*zf2eo - zdpst*zf1eo + zdpsr*zf0eo) &
                    /(zdpst*zdpsl*zdpsr)
               zaee = (zdpsl*zf2ee - zdpst*zf1ee + zdpsr*zf0ee) &
                    /(zdpst*zdpsl*zdpsr)
               zboo = ( - zdpsl**2 *zf2oo + zdpst**2 *zf1oo &
                    - (zdpst**2 - zdpsl**2) *zf0oo ) &
                    /(zdpst*zdpsl*zdpsr)
               zboe = ( - zdpsl**2 *zf2oe + zdpst**2 *zf1oe &
                    - (zdpst**2 - zdpsl**2) *zf0oe ) &
                    /(zdpst*zdpsl*zdpsr)
               zbeo = ( - zdpsl**2 *zf2eo + zdpst**2 *zf1eo &
                    - (zdpst**2 - zdpsl**2) *zf0eo ) &
                    /(zdpst*zdpsl*zdpsr)
               zbee = ( - zdpsl**2 *zf2ee + zdpst**2 *zf1ee &
                    - (zdpst**2 - zdpsl**2) *zf0ee ) &
                    /(zdpst*zdpsl*zdpsr)
               zcoo = zf0oo
               zcoe = zf0oe
               zceo = zf0eo
               zcee = zf0ee
               !
               w1l1oo(jms,ms) =  w1l1oo(jms,ms) + &
                    twopi2*zdpst*( zaoo*zdpst**2/ 3.0_r8  &
                    + zboo*zdpst/ 2.0_r8  &
                    + zcoo )
               w1l1oe(jms,ms) =  w1l1oe(jms,ms) + &
                    twopi2*zdpst*( zaoe*zdpst**2/ 3.0_r8  &
                    + zboe*zdpst/ 2.0_r8  &
                    + zcoe )
               w1l1eo(jms,ms) =  w1l1eo(jms,ms) + &
                    twopi2*zdpst*( zaeo*zdpst**2/ 3.0_r8  &
                    + zbeo*zdpst/ 2.0_r8  &
                    + zceo )
               w1l1ee(jms,ms) =  w1l1ee(jms,ms) + &
                    twopi2*zdpst*( zaee*zdpst**2/ 3.0_r8  &
                    + zbee*zdpst/ 2.0_r8  &
                    + zcee )
               !
            end do
         end do
!
      w1l1oo(ms,1:ms) = conjg( w1l1oo(1:ms,ms) )
      w1l1oe(ms,1:ms) = conjg( w1l1eo(1:ms,ms) )
      w1l1eo(ms,1:ms) = conjg( w1l1oe(1:ms,ms) )
      w1l1ee(ms,1:ms) = conjg( w1l1ee(1:ms,ms) )
!
      return
      end


