!-----------------------------------------------------------
      SUBROUTINE pstdelwrb
!
! this SUBROUTINE pstcomputes the right-hand side contribution to the energy
! \int d\psi d\theta d\zeta \xi_0 l \xi_1.
!
! a.pletzer aug.90
!
!.............................................................
!
 USE pstcom

 USE r33com

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER JBEG1
      INTEGER JEND1
      INTEGER M11
      INTEGER J1
      INTEGER LM1
      INTEGER J



!
      if ( ibas /=  1 ) then
      write(itty,*) 'ibas /=  not implemented in delwrb'
      stop 'ERROR: Only ibas=1 basis functions are supported in delwrb'
      end if
!
      jbeg1 = mbeg
      jend1 = msing(ms+2)
!
!      run from rational surface ms-1 to rational surface ms+ 1._r8 ..
!
      do 90 m1 = mbeg, jend1
!
!        limit of fourier modes.
!
      jmax1  = lmax(1) - lmin(1) + 1
         m11 = m1
         if ( m1  <  mp  )   m11 = m1 + 1
!
!      do integrals around node m1...
!
         call pststorrb
!
      do  13   j1 = 1, jmax1
!
      w0l1ou(m1,j1,ms) = twopi2*w0l1ou(m1,j1,ms)
      w0l1eu(m1,j1,ms) = twopi2*w0l1eu(m1,j1,ms)
   13 continue
!
   90 continue
!
      if(lchkap) then
!
      lm1 = lmin(1) - 1
!
      write(outmod,1004)
      write(outmod,1002) (j, j=1+lm1, jmax1+lm1)
      do m1=mbeg,jend1
      write(outmod,1003)m1,psinod(m1),( real(w0l1ou(m1,j1,ms)) &
                                 ,j1=1,jmax1)
      end do
!
      write(outmod,1005)
      write(outmod,1002) (j, j=1+lm1, jmax1+lm1)
      do m1=mbeg,jend1
      write(outmod,1003)m1,psinod(m1),(aimag(w0l1ou(m1,j1,ms)) &
                                 ,j1=1,jmax1)
      end do
!
      end if
!
      return
 1002 format(1x,'m1','     psinod  ',9(10x,i2,1x))
 1003 format(i3,1x,f10.6,9(1x,e12.4))
 1004 format(1x,'         re w0l1ou')
 1005 format(1x,'         im w0l1ou')
!
      end



