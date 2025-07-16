!-----------------------------------------------------------------------
      SUBROUTINE pstreset
!
! set quantities to zero
!
! a. pletzer 21 april  95._r8 
!------------------------------------------------------------------------
!
 USE pstcom

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER JS
      INTEGER JSS
      INTEGER JM
      INTEGER J

      ensq = n*n
      no2pi = n / twopi
      no2pi2 = no2pi * no2pi
!
      write(outmod,8020) imode
      write(itty,8020)   imode
 8020 format(" run imode",20x,"=  ",i5/)

!
      do js  = 1, nsg
      do jss = 1, nsg
      w1l1oo(js ,jss) =  0.0_r8 
      w1l1oe(js ,jss) =  0.0_r8 
      w1l1eo(js ,jss) =  0.0_r8 
      w1l1ee(js ,jss) =  0.0_r8 
      end do
      do jm = 1, nfm
      do j  = 1, nfe
      xisolo(j ,jm,js) =  0.0_r8 
      xisole(j ,jm,js) =  0.0_r8 
      w0l1ou(j ,jm,js) =  0.0_r8 
      w0l1eu(j ,jm,js) =  0.0_r8 
      end do
      end do
      end do
!
      return
      end


