!---------------------------------------------------------------------------
!
      SUBROUTINE pstdeltap
!
!
 USE pstcom

 USE combla

 USE comggf

 USE comivi

 use l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 EPSMAC
      real*8 EPSCON
      real*8 DTRY
      INTEGER NTYPE
      INTEGER NSTEPS
      INTEGER I
      INTEGER IOE
      INTEGER MSS
      INTEGER IP


 LOGICAL	 	 	 :: lquit
 INTEGER, DIMENSION(100) 	 :: ineg
 COMMON /c2f90ineg/ ineg 
!!$      namelist /eigdat/ alam,epsmac,epscon,dtry,nda,ndb,ntype,nitmax, &
!!$  nsteps, lquit,liter
!
real*8 lambda
!!$ 
      ineg = 0
!
!      set constants...
!
      epsmac =  1.0E-12_r8  
      nda = inpot
      ng = mp - 1
      mf = jsub(1) ! number of fourier modes
      ml = mf
      ndlt = outpst
      m12 = mf + ml ! number of upper diagonals
      nlong = m12 * ( m12+1 ) / 2
!     print *,'nlong d',mp,mf,ml,m12,nlong

 if (isolver > 0) then
    call pstlaeigen(wlambda)
 endif
!
!!$      if (imode == 1) read(inmode,eigdat)
!
!      here we read directly from
!      potdsk
!
      nds = nda
!
!      set c(x) for use in decomposition..
!
      CALL pstfxfu
!
!     cholesky decomposition of a, result in a...
!
!...aplet
      neg = 0
!...aplet
!
      CALL pstconalr(epsmac)
!
!      check for singular a, ideal instabilities...
!
      if ( nsing  ==  -1 )  go to 120
!aplet  if ( neg  >  0 )     go to 130
!  13.8_r8 .91
!
!      write(itty,*)' neg = ',neg
!      write(outmod,*)' neg = ',neg
!
      nIdealInstabilities=0
      if (neg > 0) nIdealInstabilities = neg
      if (neg > 0) then
         idealInstabilityNode(1:neg) = ineg(1:neg)
         write(itty,201)neg
         write(outmod,201)neg
!
!
 201  format(//1x,'plasma is ideally unstable: ',i3,' unstable' ,' mode(s)'/ &
        1x,'---------------------------')
      write(itty,203) (ineg(i), i=1,neg)
      write(outmod,203) (ineg(i), i=1,neg)
!
 203  format(/1x,'nodes where instability was detected:',    &
    /20i4)
      else
      write(itty,202)
      write(outmod,202)
!
 202  format(//1x,'plasma is ideally stable: '&
    /1x,'-------------------------')
      end if
!  13.8_r8 .91
!
! aplet  26.11_r8 .92
      if (.not. lsing) return
! aplet  26.11_r8 .92
!
!     solve for delta-primes at each rational surface...
!
      ioe = 2
      if ( lcirc )  ioe = 1
!
      do 100 mss = 1, nosing
!
!       do for both odd/even solutions...
!
      do 90 ip = 1, ioe
!
!      solve equations ; put rhs into xt, result in ut...
!
      call pststrhs2(mss,ip)


!
      CALL pstcbxlu(1)


      CALL pstcdlhxv

!
!      find delta-prime...
!
      CALL pstexsolx(mss,ip)
!
!      check the inversion by iterating
!!$      if (.not. lquit) go to 90
!!$      if(ip  /=  1) go to 90
!!$      CALL pstchkinv(mss)
!!$!
   90 continue
!
  100 continue
!
! compute concomitant matrix and matching matrix...
!
      CALL pstconcom
!
 101  continue
!
      return
!
  120 continue
!
      write(outmod,9000)
      write(itty,9000)
 9000 format(1x,"  *****  delta-w singular  *****" )
      CALL pstempty ( outmod )
      stop 'ERROR: in deltap'
!
  130 continue
!
      write(outmod,9001) neg
      write(itty,9001)   neg
 9001 format(1x,"  ideal instability, ",i4," unstable modes" )
      CALL pstempty ( outmod )
      return
!
      end


