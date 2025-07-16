      SUBROUTINE pstdefglo
!......................................................................
!      this SUBROUTINE pstdefaults logical unit numbers
!      and logical control switches for package.
!......................................................................
!
 USE pstcom
 USE l22com
 USE comggf
 use eigens
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 RSAVE
      real*8 BETAPS
      INTEGER J
 LOGICAL	 	 	 :: ldelta
 LOGICAL	 	 	 :: nodchk
 LOGICAL	 	 	 :: lremap
 LOGICAL	 	 	 :: lmapck
 LOGICAL	 	 	 :: lvacdsk
 CHARACTER*8, DIMENSION(2) 	 :: under
 REAL*8, DIMENSION(nsg1,nsg1) 	 :: zonew
 REAL*8	 	 	 :: fracem
 REAL*8, DIMENSION(nsg1) 	 :: faclay
 INTEGER, DIMENSION(3,nsg1) 	 :: mzoned
 INTEGER	 	 	 :: maxis


!     include 'comgg2.inc'
!
!      1.1_r8 .1   constants
!
      zero   =  0.0E0_r8  
      pt1    =  1.E-1_r8 
      half   =  0.5E0_r8  
      one    =  1.0E0_r8  
      two    =  2.0E0_r8  
      three  =  3.0E0_r8  
      four   =  4.0E0_r8  
      five   =  5.0E0_r8  
      seven  =  7.0E0_r8  
!
      pye    =  3.141592653589E0_r8  
      twopi  = two * pye
      twopi2 = twopi * twopi

      wlambda = 0.0_r8

      alx     = 1
      alz     = 1
      n       = 1.0_r8
!!$      lmax1   = 4
!!$      jmax2   = 2*lmax1 + 1
      ensq = n * n
      m       = 2
      mp     = 3
      minc   = 10
      mth    = 10
      mth1   = mth + 1
      mth2   = mth1 + 1
      mthh   = mth2/2
      mdiv   = 2
      nosurf = ( mp - 1 ) * mdiv + 1
      n1surf = 4
      npsurf = 6
!aplet  30.1_r8 .97      dpsi(1)= one / mdiv / m
      bit    =  1.0E-12_r8  
      amu0   = four * pye *  1.E-7_r8 
      dth    = twopi / mth
      dthinc = dth  / minc
      gamma  = five / three
      p0     = zero
      upsiln = one
      r      =  20.0E0_r8  
      r2     = r * r
      r4     = r2 * r2
      r6     = r4 *r2
      isymz  = 2
      wsqold =  0.0_r8 
      wsq =  0.0_r8 
!
!      default to force i/o into separate files for equilibrium,
!      mapping and normal modes.
!
      itty = 12
      inpest = 20
      outpst = 21
      inpot = 7
      inkin = 8

      iomode = 30
      outmp1 = 32

!!$      outmap = 25
      outd2 = 41
      inmode = 26
      outmod = 27
     outdel = 29
     outilc = 24
!!$      outgrd = 31
!!$      outnod = 33
!!$      outlay = 35
!
!!$      outso = 45
!!$      outse = 46
!!$      outbo = 47
!!$      outbe = 48
!!$      outnd = 49
!!$      outla = 51
!!$      outgd = 52
!      outev = 53
!
!!$      outneq = 37
      outvec = 44
      over = 50
!
      neqdsk = 10000
      nmod2 = 10000
!
!      logicals
!
      lcont = .false.
      lsymz  = .false.
      check1 = .false.
      fast   = .false.
      ke     = .false.
      laddmo = .false.
      wall   = .false.
      liter = .false.
!
!      the following are pointers to tell how many times the code has
!      gone through equilibrium, mapping and normal modes calculations.
!
      iequil = 0
      imap = 0
      imode = 0
! default input 
      infwal = .true.
      ldelta = .true.
      lvacdsk = .false.
!
! (zero pressure switch)
      lzerop = .false.
!
! defaults

      isolver = 0 ! solver (0= built-in, 1=LAPACK..)
!
      betmsh =  1.0_r8 
      lsymhi = .true.
!
      alfmsh =  1.0_r8 
      widmsh =  0.01_r8 
      zonew(1,:) =  1.0_r8 
      mzoned(1,:) = 1

      rsave  = r
      betaps = betap


      msin(1) = 1
      msin(2:nsg) = 0
!
      scale =  1.0_r8 

! 
! for simpson's quadrature rule pquad=1/3
!
      pquad =  1._r8 / 3._r8 
      lsymhi = .false.
      lsmoth = .false.
      lchkap = .false.
      lplot = .false.
      lplot2 = .false.
      dlayb =  0.9_r8 
      betmsh =  1.0_r8 
      lsmnus = .false.
      lpack  = .true. !.false.
      ltransf = .false.

 m=100 
 mdiv=2 
 mth = 32
 lsymz=.false.  
 lfunin=.false.
 lcirc=.false.   

 lpmax=1   
 phi= 0.0_r8   

! cplots
 lpmax=1  
 phi= 0.0_r8   

! debugs
 checkv=.true. 
 checkd=.false.  
 checki=.false.
  check1 = .false. 
 infwal=.TRUE. 
 wall=.FALSE.  
 ke=.false.
  fast=.true.   
 symvac=.true.  
 lebc=.false.   
  lmapck = .true.

! vacdat
 aw= 100._r8   
 bw= 0._r8  
 nsing=500 
 epsq= 1.E-5_r8  
 nout=50  

! shape
 a = - 10.0_r8 
 b =  0.30_r8 
 gext= 1.0_r8   
 f0= 0.02_r8  
 r= 20.0_r8     

! cprofl
 alphar= 0.0_r8  
 betar= 0.0_r8   
 delr= 1.0_r8   
 psi0r= 1.0_r8 
 alphap= 2.5_r8  
 betap= 4.0_r8   
 dlp= 1.0_r8    
 psi0p= 0.0_r8 
 betap =  4.0_r8   
 p0= 0.0_r8   
 gamma= 1.66666666667_r8   
 il=1 

! vardat
 scale=+ 1.00_r8  
 varmin= 1.0_r8   
 varmax= 4.5_r8   
 nvar=1 
 
! resis
 lmsin=.true.  
 msin = 0
 msin(1) = 1
 mm = 0
 mm(1) = 20
 mm(2) = 102
 dlays = 0.0_r8
 dlays(1) = 0.005_r8
 nodchk=.false.
 maxis=10 
 fracem= 0.25_r8  
 faclay= .1_r8 
 lextra=.true.
 lsing=.true.
 lcub=.false. 
 lmarg = .true.  
 lpstps=.true.   
 lsub=.true.  

! aplet
 lchkap = .false. 
 lsmoth = .false.
 lplot = .false.
 lplot2=.false. 
 dlayb =  0.9_r8 
 lsymhi = .false.
 alfmsh= 1.0_r8
 alfmsh(1) = 2.0_r8 

 widmsh = 0.01_r8
 widmsh(1)= 0.05_r8
 widmsh(2)= 0.02_r8
 widmsh(3)= 0.01_r8
 widmsh(4)= 0.005_r8
 widmsh(5)= 0.0010_r8

 betmsh= 1.0_r8 
 xsplus =  0._r8
 xsmnus =  0._r8 
 lsmnus = .true.
 lpack = .true.

! growth
      do j = 1,nsg
      etalay(j) =  0._r8 
      rholay(j) =  1._r8 
      end do
 slundq =  1.E+6_r8

! eigdat
 eigens_liter=.false.
 alam= 0.1
 dtry=0.01 
 nsteps=3 
 eigens_epscon=1.e-5 
 eigens_epsmac=1.0e-12  

      dw =  0.0_r8               ! typically 0 <= triangularity < 1
      sw =  0.0_r8               ! typically 0 <= squareness <=1

      msin = 0
      msin(1:2) = (/1, 0/)

!
 100  format(/8x,50('%')/8x,'% outer region resistive stability', ' code pest-3.4  %'/8x,50('%')/)
!
!
      dlay = dlays(1)
      if ( dlay  <=   0._r8  ) dlay =  0.005_r8 
      m = mm(1)
!
!!$      lmax(1:nfe) = lmax1
!!$      lmin(1:nfe) = -lmax1
!
      mp     = m + 1
!
!     only good for up-down symmetric cases?...
!
      isymz = 1
      mth1   = mth + 1
      mth2   = mth1 + 1
      mthh   = mth2/2
      no2pi  = n / twopi
      no2pi2 = no2pi * no2pi
      ensq = n * n
!
      dth    = twopi / mth
      r2     = r * r
      r4     = r2 * r2
      r6     = r4 * r2
!
      if(abs(scale-one)  >    1.E-6_r8 )lscale= .true.
!
      r      = rsave
      betap  = betaps
!
!    set  dlay , the layer width for singular cases, if not set by input.
!
       if ((dlay  <   1.E-3_r8 )   .OR.   (dlay >   0.1_r8 ) ) dlay= .005_r8 
       if (dlay  >    0.1_r8 ) dlay= .005_r8 

      end subroutine 


