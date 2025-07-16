!
       SUBROUTINE pstconcom
!
! concomitant and outer matching matrix calculation...
!
 USE pstcom

 USE combla

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER JNPAR
      INTEGER JMSS
      INTEGER JPAR
      INTEGER JMS
      INTEGER JMAX3
      INTEGER J1
      INTEGER J2



!
 COMPLEX*16	 	 	 :: zwoo
 COMPLEX*16	 	 	 :: zweo
 COMPLEX*16	 	 	 :: zwoe
 COMPLEX*16	 	 	 :: zwee
      ibas = 1
!
      jnpar = 2
      if ( lcirc ) jnpar = 1
!
        apr =  0.0_r8 
        bpr =  0.0_r8 
      gampr =  0.0_r8 
      delpr =  0.0_r8 
!
      do 200 jmss = 1,nosing
      do 200 jpar = 1,jnpar
!
       if( jpar == 1 ) then
!
! odd responses...
!
! delpr and gampr...
!
          do 10 jms = 1,nosing
             zwoo =  0.0_r8 
             zwoe =  0.0_r8 
             do m1 = 1,mp
                jmax1 = lmax(1) - lmin(1) + 1
                jmax3 = ibas*jmax1
                do j1 = 1,jmax3
                   zwoo = zwoo + &
                        conjg( xisolo(m1,j1,jmss) )*w0l1ou(m1,j1,jms)
                   zwoe = zwoe + &
                        conjg( xisolo(m1,j1,jmss) )*w0l1eu(m1,j1,jms)
                end do
             end do
             !
             !      write(itty  ,104) jms,jmss,
             !    .    real(zwoe),aimag(zwoe), real(zwoo),aimag(zwoo)
             !    .  , real(w1l1eo(jms,jmss)),aimag(w1l1eo(jms,jmss))
             !    .  , real(w1l1oo(jms,jmss)),aimag(w1l1oo(jms,jmss))
             write(outmod,104) jms,jmss, &
                  real(zwoe),aimag(zwoe), real(zwoo),aimag(zwoo) &
                  , real(w1l1eo(jms,jmss)),aimag(w1l1eo(jms,jmss)) &
                  , real(w1l1oo(jms,jmss)),aimag(w1l1oo(jms,jmss))
             !
             cjump(jms,nosing+jmss) =  zwoe -  w1l1eo(jms,jmss)
             gampr(jmss,jms) = cjump(jms,nosing+jmss)
             cjump(nosing+jms,nosing+jmss) =  zwoo - w1l1oo(jms,jmss)
             delpr(jmss,jms) = cjump(nosing+jms,nosing+jmss)
             !
10       continue
!
      else
!
! even responses...
!
!
! apr and bpr...
!
         do 20 jms = 1,nosing
            zweo =  0.0_r8 
            zwee =  0.0_r8 
            do m1 = 1,mp
               jmax1 = lmax(1) - lmin(1) + 1
               jmax3 = ibas*jmax1
               do j1 = 1,jmax3
                  zweo = zweo + &
                       conjg( xisole(m1,j1,jmss) )*w0l1ou(m1,j1,jms)
                  zwee = zwee + &
                       conjg( xisole(m1,j1,jmss) )*w0l1eu(m1,j1,jms)
               end do
            end do
            !
            !      write(itty  ,103) jms,jmss
            !    .  , real(zwee),aimag(zwee), real(zweo),aimag(zweo)
            !    .  , real(w1l1ee(jms,jmss)),aimag(w1l1ee(jms,jmss))
            !    .  , real(w1l1oe(jms,jmss)),aimag(w1l1oe(jms,jmss))
            write(outmod,103) jms,jmss &
                 , real(zwee),aimag(zwee), real(zweo),aimag(zweo) &
                 , real(w1l1ee(jms,jmss)),aimag(w1l1ee(jms,jmss)) &
                 , real(w1l1oe(jms,jmss)),aimag(w1l1oe(jms,jmss))
            !
            cjump(nosing+jms,jmss) = zweo - w1l1oe(jms,jmss)
            bpr(jmss,jms) = cjump(nosing+jms,jmss)
            cjump(jms, jmss) = zwee - w1l1ee(jms,jmss)
            apr(jmss,jms) = cjump(jms, jmss)
            !
20       continue
!
       end if
!
 200   continue
!
       write(itty,101)
       write(outmod,101)
       do 31 j1 = 1,nosing
       write(itty  ,106)  &
 ( real(  apr(j1,j2)),j2=1,nosing), &
 ( real(  bpr(j1,j2)),j2=1,nosing)
       write(outmod,105)  &
 ( real(  apr(j1,j2)),j2=1,nosing), &
 ( real(  bpr(j1,j2)),j2=1,nosing)
 31   continue
!
      do 32 j1 = 1,nosing
       write(itty  ,106) &
 ( real(gampr(j1,j2)),j2=1,nosing), &
 ( real(delpr(j1,j2)),j2=1,nosing)
       write(outmod,105) &
 ( real(gampr(j1,j2)),j2=1,nosing), &
 ( real(delpr(j1,j2)),j2=1,nosing)
 32    continue
!
      write(itty  ,107)
      write(outmod,107)
      do j1 = 1,nosing
      write(itty  ,106) &
 (aimag(  apr(j1,j2)),j2=1,nosing),  &
 (aimag(  bpr(j1,j2)),j2=1,nosing)
      write(outmod,105) &
 (aimag(  apr(j1,j2)),j2=1,nosing), &
 (aimag(  bpr(j1,j2)),j2=1,nosing)
      end do
!
      do j1 = 1, nosing
      write(itty  ,106) &
 (aimag(gampr(j1,j2)),j2=1,nosing), &
 (aimag(delpr(j1,j2)),j2=1,nosing)
      write(outmod,105) &
 (aimag(gampr(j1,j2)),j2=1,nosing), &
 (aimag(delpr(j1,j2)),j2=1,nosing)
      end do
!
       return
 101  format(1x/'2 times matching matrix (re part)'/)
 107  format(1x/'2 times matching matrix (im part)'/)
 102  format(10e15.7)
 105  format(6(1x,f10.4,' ',f10.4))
 106  format(6(1x,e10.2,' ',e10.2))
!
 103  format(/1x,'energies:'/1x, &
 'jms = ',i3,' jmss = ',i3/  &
 'zwee = ',e15.7,' + i*',e15.7/  &
 ' zweo = ',e15.7,' + i*',e15.7/ &
 'w1l1ee(jms,jmss) = ',e15.7,' + i*',e15.7/ &
 ' w1l1oe(jms,jmss) = ',e15.7,' + i*',e15.7/)
!
!
 104  format(/1x,'energies:'/ &
  1x,       'jms = ',i3,' jmss = ',i3/  &
 'zwoe = ',e15.7,' + i*',e15.7/      ' zwoo = ',e15.7,' + i*',e15.7/ &
 'w1l1eo(jms,jmss) = ',e15.7,' + i*',e15.7/ &
 ' w1l1oo(jms,jmss) = ',e15.7,' + i*',e15.7/)
       end



