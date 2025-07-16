!-----------------------------------------------------------------------
SUBROUTINE pstrepet2
  !
  ! set up new mesh node number and extrapolate matching data from previous
  ! runs.
  !
  ! a. pletzer nov 2000
  !------------------------------------------------------------------------
  !
  USE pstcom
  USE comggf
  IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

  COMPLEX*16, DIMENSION(2*nsg,2*nsg), save 	 :: sy
  COMPLEX*16, DIMENSION(2*nsg,2*nsg), save	 :: sxy
  REAL*8, save :: s, sx, sxx
  complex*16, dimension(2*nosing, 2*nosing) :: a_lsq, b_lsq
  real*8    , dimension(2*nosing, 2*nosing) :: sig_lsq
  complex*16, dimension(2*nsg, 2*nsg, nrun), save :: old_data

  integer i, j, jmm, jr, jrate
  real*8 zmm2, zd

  !
  if(imode==1) then
     ! initialize
     sx =  0.0_r8 
     sxx =  0.0_r8 
     sy =  (0.0_r8, 0.0_r8)
     sxy =  (0.0_r8, 0.0_r8)
  endif

  if( .not. lsing ) return

  lcont = .true.
  if(imode+1 == nrun .or. mm(imode+1) == 0) lcont = .false.

  ! store old data

  old_data(       1:  nosing,        1:  nosing, imode) =   apr(1:nosing, 1:nosing)
  old_data(       1:  nosing, nosing+1:2*nosing, imode) =   bpr(1:nosing, 1:nosing)
  old_data(nosing+1:2*nosing,        1:  nosing, imode) = gampr(1:nosing, 1:nosing)
  old_data(nosing+1:2*nosing, nosing+1:2*nosing, imode) = delpr(1:nosing, 1:nosing)

  jrate = 2
  if(lcub) jrate = 6
  jr = -jrate

  jmm = mm(imode)
  zmm2 = (real(jmm, r8))**jr
  s = s +  1._r8 
  sx = sx + zmm2
  sxx = sxx + zmm2**2
  sy(1:2*nosing, 1:2*nosing)   =  sy(1:2*nosing, 1:2*nosing)  &
       & + old_data(1:2*nosing, 1:2*nosing,imode)
  sxy(1:2*nosing, 1:2*nosing)   = sxy(1:2*nosing, 1:2*nosing)  &
       & + old_data(1:2*nosing, 1:2*nosing,imode)*zmm2

  m = mm(imode + 1)


  write(itty,9005)imode
  write(outmod,9005)imode
9005 format(//,20x,' end of case no',i3,/)


  if(m > 0) then
     mp = m + 1
     w1l1oo =  0.0_r8 
     w1l1oe =  0.0_r8 
     w1l1eo =  0.0_r8 
     w1l1ee =  0.0_r8 
     xisolo =  (0.0_r8, 0.0_r8)
     xisole =  (0.0_r8, 0.0_r8) 
     w0l1ou =  0.0_r8 
     w0l1eu =  0.0_r8 
!!$         nosurf = ( mp - 1 ) * mdiv + 1
!!$         nosurf = nusurf
!!$         nsf0 = nusurf
!!$         nsf   = nsf0
!!$         nsf2 = 2*nsf
!!$         nfe = 1 + nsf/2
!!$         nfe1=nfe+1
     return
  endif
  !
  ! last iteration
  !
  jmm = mm(imode)
  !
  ! convergence rate
  !
  jrate = 2
  if(lcub) jrate = 6
  jr = -jrate
  zmm2 = (real(jmm))**jr
  !
  !
  ! use least square regression a_lsq + b_lsq*zmm2
  ! formula to extrapolate to infinite mm... 
  !
  zd = s * sxx - sx**2
  if(imode==1) zd = 0.0_r8
  if(zd/=0.0_r8) then
     a_lsq = (sxx* sy(1:2*nosing,1:2*nosing) - &
          & sx*sxy(1:2*nosing,1:2*nosing))/zd
     b_lsq = (  s*sxy(1:2*nosing,1:2*nosing) - &
          & sx* sy(1:2*nosing,1:2*nosing))/zd
  else
     a_lsq = sy(1:2*nosing,1:2*nosing)
     b_lsq = (0.0_r8, 0.0_r8)
  endif

  sig_lsq = 0.0_r8
  
  do i = 1, imode
     sig_lsq = sig_lsq + ( &
          & abs( a_lsq + b_lsq/real(mm(i)**2, r8) - &
          & old_data(1:2*nosing, 1:2*nosing,i) ) &
          &              )**2
  enddo
  ! most mean square variation use 1/(N-2)
  ! for simplicity use 1/N
  sig_lsq = sqrt( sig_lsq/real(imode, r8) )

  ! values extrapolated to infinite no of finite elements

  apr(1:nosing, 1:nosing) = a_lsq(       1:  nosing,        1:  nosing)
  bpr(1:nosing, 1:nosing) = a_lsq(       1:  nosing, nosing+1:2*nosing)
  gampr(1:nosing, 1:nosing) = a_lsq(nosing+1:2*nosing,        1:  nosing)
  delpr(1:nosing, 1:nosing) = a_lsq(nosing+1:2*nosing, nosing+1:2*nosing)

  ! error estimates

  error_apr(1:nosing,1:nosing)   = sig_lsq(       1:  nosing,        1:  nosing)
  error_bpr(1:nosing,1:nosing)   = sig_lsq(       1:  nosing, nosing+1:2*nosing)
  error_gampr(1:nosing,1:nosing) = sig_lsq(nosing+1:2*nosing,        1:  nosing)
  error_delpr(1:nosing,1:nosing) = sig_lsq(nosing+1:2*nosing, nosing+1:2*nosing)

  write(itty  ,108)
  write(outmod,108)
  !
108 format(//1x,'2 times matching matrix extrapolated to infinite m' , &
       ' (re part)'/)
  do i = 1, 2*nosing
     write(itty  ,111)( real(a_lsq(i,j)), j = 1,2*nosing)
     write(outmod,111)( real(a_lsq(i,j)), j = 1,2*nosing)
  end do
111 format(12(1x,f10.5))
  !
  write(itty  ,109)
  write(outmod,109)
  !
109 format(//1x,'2 times matching matrix extrapolated to infinite m', &
       ' (im part)'/)
  do i = 1, 2*nosing
     write(itty  ,111)(aimag(a_lsq(i,j)), j = 1,2*nosing)
     write(outmod,111)(aimag(a_lsq(i,j)), j = 1,2*nosing)
  end do

  call pstWresult

end SUBROUTINE pstrepet2



