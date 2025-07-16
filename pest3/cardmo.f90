      SUBROUTINE pstcardmo
!.......................................................................
 USE pstcom

 USE l22com

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

      real*8 RSAVE
      real*8 BETAPS
      INTEGER I
      INTEGER J
      INTEGER JM
      INTEGER JM2
!
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

 COMMON /c2f90ldelta/ ldelta 
 COMMON /c2f90lmapck/ lmapck 
 COMMON /c2f90lvacdsk/ lvacdsk 
 COMMON /c2f90nodchk/ nodchk 
 COMMON /c2f90zonew/ zonew 
 COMMON /c2f90mzoned/ mzoned 
 COMMON /c2f90maxis/ maxis 
 COMMON /c2f90fracem/ fracem 
 COMMON /c2f90faclay/ faclay 
!!$      namelist /modes/ lmax1,lmin1,m,mth,n,mdiv,lsymz,lfunin,lcirc,jmax2
!!$      namelist / cplots / lpmax, phi
!!$      namelist / debugs / checkd, checki, checkv,check1, ke, fast, &
!!$     lebc,symvac, wall, infwal, ldelta, lmapck, lso, lvacdsk
!!$      namelist / vacdat / aw, bw, dw, sw, nsing, epsq, nout, delg
!!$      namelist / shape  / a, b, r,  gext, f0
!!$      namelist / cprofl / alphar, betar, delr, psi0r, &
!!$                     alphap, betap, dlp, psi0p, p0, &
!!$                     gamma, il
!!$      namelist / vardat / varmin, varmax, nvar, scale
!!$      namelist / resis/ mm,lcub,lmarg,lsing,lsub,lpstps &
!!$      ,lremap,dlays,zonew,mzoned,nodchk,lmsin,msin,lrego,lextra &
!!$      , maxis, fracem, faclay
!!$      namelist/aplet/ laplet,lchkap,lsmoth,lplot,lplot2 &
!!$     ,dlayb, lsymhi, alfmsh, widmsh, betmsh &
!!$     ,xsplus,xsmnus, pquad, lsmnus, lpack, ltransf
!!$      namelist/growth/ slundq, etalay, rholay
!
      data under / "-----","-----" /
!
!      there are m zones, m+1 elements..
!
!          checkd.... dumps vacuum matrices
!          checki....writes out indices in indint.
!          checkv....sets checks on vacuum,finds eigenvalues and
!                    aborts.
!          check1....writes out pot,kin in delpla and minim.
!          lmapck....will plot the remapped profiles for comparison
!          ldelta....assumes that the pot and kin disks are avaibale
!          symvac....forces vacuum matrix to be symmetric
!          fast......in this version fast=.true. eliminates the element
!                    at the origin.
!          lsub......if true  msub(i) from input data will determine
!                    the subdivision of zones between rational surfaces
!                    otherwise,  msub=mm.
!          lsing.....if true  includes singular expansion functions at
!                    rational surfaces.
!
!aplet: 14/jun/ 91._r8  eleminate sing. expansion but keep lsing as switch.
!
!          lmarg.....if true  solves for marginal stability...
!                    otherwise, finds growth rate for model ke
!                    normalization.
!          lpstps....if true,reads old mapping files...
!          lsymz.....if true, will not adjust equally spaced mesh..
!          ldelta....if false, will skip matelt (calc of pot and kin)
!
!        scale.....allows rescaling of equilibrium if  /=   1._r8  0._r8  if
!                  scale <  0._r8  then -scale represents the new q value 
!                  on axis. otherwise, if scale >  1.0_r8  then scale
!                  represents new q value at the edge.
!
!        dlay......sets the layer width for the element adjacent to 
!                  the rational surfaces. allowed range
!                  is  0.001_r8  to  0.1_r8  if dlay is less than or equal to zero
!                  it is set by the dpsi of the first zone.
!        mzoned.... .3_r8  element array giving the number of finite elements
!                   in each zone.
!
!        msin.......binary array which selects (if =1) or discards 
!                   (if = 0) rational surfaces of index i,  where i runs 
!                   from the magnetic axis to the edge. alternatively,
!                   if msin > 1 then the msin represent the min and max
!                   values of the resonant poloidal modes to be treated.
!
!        lextra.....run with various number of mesh nodes
!
!        mm.........array giving the number of mesh nodes fort each run
!                   (last value must be zero).
!
!        alfmsh.....array containing the exponent for maesh packing near
!                   rational surfaces.
!
!        widmsh.....corresponding mesh width in psi
!
!        betmsh.....global mesh packinkg exponent of distortion function
!
!        dlayb......width of the prescribed solution support normalized
!                   to the distance between two rational surfaces
!                   ( 0 < dlayb < 1 ).
!
!        lsymhi.....if true then symmetric prescribed shape function.
!
!        xsplus, xsmnus..coeff. to insert small solution in
!                        the prescribed solution.
!
!        lsmnus.....sets appropriate xsmnus coeff. for  .5_r8  < mu <  1._r8 
!
!        lpack......automatic mesh packing (finite beta) if true. if zero beta
!                   plasma then set alfmsh= 1._r8 
!
!        ltransf....use transformed operator version if TRUE.
!
!        slundq.....lundquist number.
!
!        nqscan.....number of growth rates to determine roots of
!                   dispersion relation.
!
!        etalay......resistivity in the various layers (in kev)
!                   check if slundq /= 0, if not then compute lundquist
!                   number using etalay.
!
!        rholay.....density in the various layers (in kg/m3/ mu0 )
!
!       1.2_r8 .1   job title
!             ... .....
!
      isolver = 0

      infwal = .true.
      ldelta = .true.
      lvacdsk = .false.
!
! (zero pressure switch)
      lzerop = .false.
!
! defaults
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
      lpack  = .true. ! .false.
      ltransf = .false.

!
!!$      rewind inmode
!
!       1.2_r8 .2   read data from card.
!     ............................
!
!
! modes

!!$ lmax1= +4
!!$ lmin1= -4
!!$ n =  1.0_r8 
!!$ jmax2 = 9
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
 mm(1) = 50
 mm(2) = 103
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
!!$ liter=.false.
!!$ alam=  0.1_r8 
!!$ dtry= 0.01_r8  
!!$ nsteps=3 
!!$ epscon= 1.E-5_r8  
!!$ epsmac= 1.0E-12_r8    


!!$ 1221 continue
!!$      print *,'read modes'
!!$      read(inmode,modes)
!
!!$      if(jmax2 /= 1 ) jmax2 = 2*max(lmax1, -lmin1) + 1


      nfn = jmax2
      if( nfn > nths/2 ) then
         write(itty, '(a,i3,a,i3)') &
              & ' Mapper resolution appears insufficient NFN=',nfn, &
              & ' NTHS0 = ',nths0
         write(itty, *) 'Rerun mapper with higher NTHS0 (MTH) or decrease LMAX1-LMIN1'
         stop 'ERROR in cardmo'
      endif

!
!!$      print *,'read cplots'
!!$      read(inmode,cplots)
!!$      print *,'read debugs'
!!$      read(inmode,debugs)
      dw =  0.0_r8               ! typically 0 <= triangularity < 1
      sw =  0.0_r8               ! typically 0 <= squareness <=1
!!$      print *,'read vacdat'
!!$      read(inmode,vacdat)
!!$      print *,'read shape'
!!$      read(inmode,shape)
!!$      print *,'read cprofl'
!!$      read(inmode,cprofl)
!!$      print *,'read vardat'
!!$      read(inmode,vardat)
!!$      print *,'read resis'
!!$      read(inmode,resis)
!
      do 10 i = 2, nsg
      if( msin(i)  ==  0 ) go to 10
      if( msin(i)  >   lmax1   .OR.   &
     &          msin(i)  <  -lmax1 ) msin(i) = 0
 10   continue


!
!!$      print *,'read aplet'
!!$      read(inmode,aplet)
!
      if(betmsh <=  0._r8 ) betmsh =  1._r8 
!
!!$      slundq = 1.0e6_r8
!
!!$      print *,'read growth'
!!$      read(inmode,growth)
!
!!$! set dims here
!!$
!!$      nfn = jmax2
!!$      nfm = nfn
!!$      nsg3=ipolp*2*(nfm+2)
!!$      nthsfm=4 * nthss * nfm
!!$      nkd=4*nfm
!!$      nke=2*nfm*nfe
!!$      nfv0=nfm * nfm
!!$      nmt=ipolp* nfm * ( nfe-2 ) + nfm
!!$      nmt2=nmt*nmt
!!$      nmtg = nfm*nfe
!!$      nfmgbx = 2*(nfm+1)
!!$      nfm21=2*nfm + 1
!!$      nfm2=nfm21*nfm21
!!$      nfm41=2*nfm21
!!$      nfm8=nfm41*nfm41
!!$      nkb=ipolp*nfm + 1
!!$      nkc = ipolp * 2 * nfm + 2
!!$      nad = (nkc)*(nkc+1)/2
!!$      nkc2 = nkc*nkc
!!$
!!$      CALL pstmemory(3)

!
!     output data just read in
!     ...... .... .... .... ...
!
!
!!$      write ( outmod, modes )
!!$      write(outmod,cplots)
!
!!$      j = lmax1 - lmin1 + 1
!!$      jm = max( abs(lmax1), abs(lmin1) )
!!$      jm2 = (jmax2-1)/2
!!$      write (itty   , 9004 ) n,lmax1,lmin1,jmax2
!!$      write (outmod , 9004 ) n,lmax1,lmin1,jmax2
!!$!
!!$      if( jm >  jmax2   .AND.   j /= 1 ) then
!!$      write (itty   , 9005 )
!!$      stop
!!$      end if
!!$      if( (lmax1 /= jm2   .OR.   abs(lmin1) /= jm2)  &
!!$       .and. (j /= 1)   ) then
!!$       write (itty   , 9006 )
!!$      end if 
!!$!
!!$      if ( jmax2 > nfn ) then
!!$       write(itty,*)' INPUT ERROR: jmax2 should be <= nfn '
!!$       stop
!!$      end if      
!!$!
!!$      write ( outpst, 9003 )lmax(1),lmin(1),m,mdiv,n
!!$      write ( outmod, shape )
!!$      write ( outmod, cprofl )
!!$      write ( outmod, debugs )
!!$      write(outmod,vacdat)
!!$      write(outmod,vardat)
!!$      write(outmod,resis)
!!$      write(outmod,aplet)
!!$      write(outmod,growth)
!
!      1.2_r8  .3_r8    set any additional constants.
!     .....................................
!
      dlay = dlays(1)
      if ( dlay  <=   0._r8  ) dlay =  0.005_r8 
      m = mm(1)
!
!!$ if ( m > (nsf-1)/2 ) then
!!$  write(itty,*)' INPUT ERROR: mm(1) > (nsf-1)/2, nsf=',nsf
!!$  stop
!!$ end if
!!$!
      lmax(1:nfe) = lmax1
      lmin(1:nfe) = -lmax1
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
!

      return
 9000 format ( 1x, 10a8, / 1x, 4a5 / )
 9003 format( " lmax, lmin=",2i4," m, mdiv=",2i4,3x," n=",e12.4/)
!
!
 9004 format( " toroidal mode number       n = ",f5.2/   &
      " max fourier mode        lmax = ",i5/   &
      " min fourier mode        lmin = ",i5/   &
      " jmax2                        = ",i5/)
 9005 format(/" *** error: increase jmax2")
!
!
 9006 format(/" *** warning: lmax  /=  abs(lmin) can cause" ,&
        " inaccuracies in the inversion of the surface green" ,&
        " function")
!
 9100 format ( 10a8 )
    end subroutine


