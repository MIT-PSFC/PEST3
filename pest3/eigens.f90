!--------------------------------------------------------------------------
      SUBROUTINE psteigens
!----------------------
! up-down asymmetry ap  17.10_r8 .95 (no modifs required)
!
 USE pstcom
 use eigens

!     dimension  al(3), npin(14), npout(4)
 USE combla
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      integer :: nstep = 1
      real*8 FACT
      INTEGER NOVLAP
      real*8 ALAML
      real*8 ALAMU
      INTEGER ICNT
      INTEGER I
      real*8 ALAMT
      INTEGER NTIPE
      INTEGER NEG
      INTEGER NIT
      INTEGER NEGL
      INTEGER NEGU
      real*8 EIGVL
      INTEGER NEGS


 LOGICAL	 	 	 :: lquit
 LOGICAL	 	 	 :: ltry
 REAL*8, DIMENSION(10) 	 :: eigvag
 REAL*8, DIMENSION(3) 	 :: al
 REAL*8	 	 	 :: epscon
 REAL*8	 	 	 :: epsmac
 INTEGER, DIMENSION(14) 	 :: npin
 INTEGER, DIMENSION(4) 	 :: npout
 COMMON /c2f90eigvag/ eigvag 
 COMMON /c2f90al/ al 
 COMMON /c2f90epscon/ epscon 
 COMMON /c2f90epsmac/ epsmac 
 COMMON /c2f90npin/ npin 
 COMMON /c2f90npout/ npout 
!!$      namelist /eigdat/ alam,epsmac,epscon,dtry,nda,ndb,ntype,nitmax, &
!!$  nsteps, leigens_quit,liter
!
!        read constants for eigenvalue problem
!
      eigens_lquit = .false.
      fact =  1.0_r8 
      write(outpst,920) fact
      if ( liter )  go to 20
      alam=zero
      epsmac= 1.0E-12_r8  
      epscon= 1.0E-3_r8  
      dtry= 0.1_r8 
      nda=inpot
      ndb=inkin
      ntype=0
      nitmax=10
      nsteps = 1
      novlap = ipol * jsub(1)
      lquit = eigens_lquit
      epscon= eigens_epscon
      epsmac= eigens_epsmac
      go to 30
   20 continue
      alam =  0.0_r8 
   30 continue
!
!      quit, and continue with goeigen ??
!
      if(lquit) stop 'ERROR in eigens'
!
      alam = alam * fact
      dtry = dtry * fact
      alaml = alam
      alamu = alaml - nsteps * dtry
      nsteps = nsteps + 1
      ltry = .false.
      icnt = 0
      nz = 0
!
!          1.0_r8   empty buffers and reset for read.
!
!        solve eigenvalue problem
!       --------------------------
!
!       set constants
!
      al(1)=alam
      npin(1)=ntype
      npin(2)=mp-1
      npin(3)=novlap
      npin(4)=novlap
      npin(5)=nitmax
      npin(7)=outmod
      npin(8)=nda
      npin(9)=ndb
      npin(10)=outpst
print *,'eigenvalue pstsearch launched'
print *,'nstep = ',nstep
      do 400 i = 1, nsteps
print *,'i = ',i
print *,'alam = ', alam

!
!      define alam
!
  100 alam = alaml - (i-1) * dtry
      if(ltry) alam = alamt
      ltry = .false.
      ntipe = 0
!      special for alam = 0, do shift only do not iterate
      if(alam ==  0.0_r8 ) ntipe = 1
!
!      shift and iterate
!
      CALL pstivi(ntipe,alam,neg,nit)
print *,'after ivi'
print *,'ntipe = ',ntipe
print *,'alam = ',alam
print *,'neg = ',neg
print *,'nit = ',nit
!
!      first time in count the number of eigenvalues less than alaml
      if ( i ==  1  .AND.  icnt == 0) negl = neg
      if(i == nsteps) negu = neg
!
!      no eigenvalues less than zero , then quit
      if(negl  ==  0  .AND.  alam  ==   0.0_r8 ) go to 200
!      possible overwriting error
      if(neg  >  50) go to 900
!      no iterations done then redefine alam and proceed
      if(ntipe  ==  1) go to 400
!      does the eigenvalue lie in the range of interest ?
      eigvl = al(2)
print *,'al(2) = ',al(2)
!      eigenvalue lies above the window
      if(eigvl  >  alaml)go to 700
!      eigenvalue lies below the window
      if(eigvl  <  alamu) go to 750
!      within the window limits
!       is it a converged value ?
      if(nit  <=  nitmax) go to 500
      if( i  ==  nsteps  .OR.  neg  ==  0 ) go to 160
      go to 400
  160 continue
      icnt = icnt + 1
      ltry = .true.
      if( icnt  ==  3 ) go to 800
      alamt = al(2)
      go to 100
!      yes. then output value etc.
  500 nz = nz + 1
      eigvag(nz) = eigvl /fact
      call pstsaver(nz)
!      continue or quit ?
      if (neg  ==  0) go to 150
      if( neg  ==  1  .AND.  alam  >  eigvl) go to 150
      if( negu  ==  negl) go to 200
      go to 400
  700 continue
!      eigenvalue lies above window
!      if eigenvalue is positive reset alam and try again
      if(eigvl  <   0.0_r8 )go to 750
      if(alam  >  alamu) go to 400
      icnt = icnt + 1
      ltry = .true.
      if(icnt  ==  3) go to 800
      alamt = alam /  10.0_r8 
      go to 100
  750 continue
!      converged ? yes, then output. no, ignore.
      if(nit  <=  nitmax) go to 500
  300 continue
!
!       has not converged use rayleigh quotient
!
      alam = al(3)
      go to 160
!
  400 continue
!
!      output eigenvalues etc.
!
  150 continue
      CALL pstoutput(nz)
      negs = negl - neg
      write(itty,930) negs,nz
      write(itty,940) ( eigvag(i),i=1,nz )
      wsq = eigvag(1)
      write(outmod, 930) negs,nz
      write(outmod, 940)(eigvag(i), i = 1, nz)
      write(outmod,950)
      CALL pstempty(outmod)
      return
!
!       no negative eigenvalues
!
  200 continue
!
      if ( negl == 0  .AND.  alam ==  0._r8  )  wsq =  0.0_r8 
!
      alamu = alamu / fact
      alaml = alaml / fact
      write( outmod, 970) alaml, alamu
      write(outmod,950)
      go to 150
!
  800 continue
      write(itty,990)
      write(outpst,990)
      write(outmod,990)
      go to 150
!
  900 continue
      write(itty,980)neg
      write(outpst,980)neg
      write(outmod,980)neg
!!$      stop
!
  920 format(/,' fact= ',e14.7/)
!
  930 format(1x,1x," there are ", i3, &
 " eigenvalues within the window,", i3," were found",&
 /," the eigenvalues are : ",/)
  940 format(15x,e14.7,/)
  950 format( 1x," ********************************************** ",/)
!
  960 format(" there are ",i3," eigenvalues in the specified range", &
 /,  2x," only two will be found",/)
!
  970 format(1x," there are no eigenvalues in the range from",1x,e14.7,&
 " to", 1x,e14.7/)
  980 format(" neg is very large =",i5," probable error***",/)
  990 format(" pstsearch for eigenvalues seems to be in a loop***",/)
      return
!
      end subroutine


