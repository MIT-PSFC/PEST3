      SUBROUTINE scrdescur (rmb,zmb,nmb,mpolin,icon,ioflag, &
     &  fsq,lasym0,lzsmlonly,nmbout,rmbout,zmbout,ier)
      USE scrunch_inc1
      USE scrunch_inc2
      IMPLICIT NONE
       INTEGER 	m ! <descur.f>
       INTEGER 	ntoff ! <descur.f>
       INTEGER 	ntype ! <descur.f>
       REAL*8 	time_ou ! <descur.f>
       INTEGER 	modeno ! <descur.f>
       REAL*8 	gout ! <descur.f>22
       REAL*8 	g11 ! <descur.f>
       REAL*8 	gmin ! <descur.f>
       INTEGER 	iter ! <descur.f>
       INTEGER 	imodes ! <descur.f>
       REAL*8 	r10 ! <descur.f>
       INTEGER 	n ! <descur.f>
       INTEGER 	nplane ! <descur.f>
       REAL*8 	time_in ! <descur.f>
       INTEGER 	ier ! <descur.f>
       INTEGER 	i ! <descur.f>
       INTEGER 	irst ! <descur.f>
       INTEGER 	nmbout ! <descur.f>
       INTEGER 	nmb ! <descur.f>
 
       real*8 :: gtest  ! dmc APR 2010

!       THIS IS PROGRAM DESCUR - WHICH USES A STEEPEST DESCENT
!       ALGORITHM TO FIND A LEAST SQUARES APPROXIMATION TO AN
!       ARBITRARY 3-D SPACE CURVE. ANGLE CONSTRAINTS, BASED ON A
!       MINIMUM SPECTRAL WIDTH CRITERION, ARE APPLIED.
!       THE CONSTRAINT IS SATISFIED BY TANGENTIAL VARIATIONS
!       ALONG THE CURVE.
!
!       Written by Steve Hirshman ORNL
!       Modifications by Dick Wieland PPPL
!
!       Reference:   "Optimized Fourier Representations for
!                     three-dimensional magnetic surfaces", S.P. Hirshman
!                     and H.K. Meier, Phys. Fluids 28,1387 (1985)
!***********************************************************************
!       REVISION 1:  January 26, 1989
!                    Compute Fit to individual toroidal planes separately
!                    and then Fourier decompose in phi, using optimized
!                    angle representation
!       REVISION 1.1 September 18, 1990 (Wieland)
!                    Modularize It.
!       REVISION 1.2 October 6, 1992 (Wieland)
!                    Ioflag to optionally turn off IO to screen.
!                    Ier return to indicate error
!                    Ier=1;  Ok
!                    Ier=11; The magnetic axis is not enclosed by boundary
!                    Ier=12; Incorrect initial angle assignment
!                    Ier=13; Time step reduced 100 times without convergence
!                    Ier=99: some other error; message displayed.
!       REVISION 1.3 March 31, 1993 (Wieland)
!                    So far unsuccessful attempt to do a hot restart
!                    C.f. lines with NINIT
!       REVISION 1.4 May 11, 1993 (Wieland)
!                    turn on the Rmns and Zmnc terms
!       REVISION 1.5 Oct 21, 1993 (Wieland)
!                    modify the constraints per SPH
!                    ZC(m=1) = RS(m=1)
!       REVISION 1.6 Feb 7, 1994 (Wieland) final SPH asym update
!                    If asym terms (RC,ZS) are on, lasym=.true.
!       REVISION 1.7 Dec 20, 1994 (Wieland) use LZAKHAROV0 to
!                    prototype UDA with z=z0+z*sin(theta) only
!       REVISION 1.8 Jan 28, 1997 (Terpstra) Change dimensions in
!                    Subroutine Funct1des to * from 1 because this
!                    is flagged as an error on the HP's at GA
!       REVISION 1.9 Apr 25, 1997 (Wieland) clean up interface and
!                    incorporate changes made by JPJ at JET
!
!***********************************************************************
 
!**   Input:  rmb(nmb)        REAL*8 array of boundary R's (any units)
!**           zmb(nmb)        REAL*8 array of boundary Z's (any units)
!**           nmb             integer ... number of boundary points
!**           convention for rmb,zmb is that rmb(1)=rmb(nmb), zmb(1)=zmb(nmb)
!**           mpolin          number of moments to compute (6 means do 0...5)
!**           icon            integer ... dummy, not used
!                (rev dmc -- if icon > 0 icon = LUN for "scrdescur_debug.dat"
!                a debug binary file)
!**           ioflag          logical ... TRUE for printout
!**   Output: fsq             REAL*8 ... resulting "force" residue
!**   Input:  lasym0          logical ... TRUE to calculate Rmc,Rms,Zmc,Zms
!**                                       FALSE to calc only Rmc,Zms terms
!**           nmbout          1st dim of rmbout,zmbout
!**   Output: rmbout(nmbout,2)  rmbout(*,1) = Rmc terms
!**                             rmbout(*,2) = Rms terms
!**           zmbout(nmbout,2)  zmbout(*,1) = Zmc terms
!**                             zmbout(*,2) = Zms terms
!**           ier               integer ... error return (1 = OK)
!**
!**   Fitting the boundary points to :
!**     R = sum(m=0,mpolin-1){Rmc*cos(m*theta) + Rms*sin(m*theta)}
!**     Z = sum(m=0,mpolin-1){Zmc*cos(m*theta) + Zms*sin(m*theta)}
 
!**           * * * * * * * * * * * * * * * * * * * * *
 
        REAL*8 rmb(nmb),zmb(nmb),fsq
        REAL*8 rmbout(0:nmbout,2),zmbout(0:nmbout,2)
        integer mpolin,icon
	logical ioflag,lasym0,lzsmlonly
 
!**********************************************************************
!       CONTROL DATA THAT CAN BE SPECIFIED BY USER
!**********************************************************************
      integer, save :: niter=3000, nstep=100, ninit=0
      REAL*8, save :: ftol=1.D-5
!**********************************************************************
!       DATA USED INTERNALLY (NOT TO BE CHANGED BY CASUAL USER)
!**********************************************************************
      REAL*8, save :: tot_time=0.0D0
!**********************************************************************
!  dmc: debug
        if(icon.gt.0) then
           open(unit=icon,file='scrdescur_debug.dat',status='unknown', &
                form='unformatted')
           call scrpreset_wr(icon)
           write(icon) nmb,mpolin
           write(icon) rmb(1:nmb)
           write(icon) zmb(1:nmb)
           write(icon) lasym0
           write(icon) nmbout
           close(unit=icon)
        endif
!**********************************************************************
        ier=1  ! assume success initially
        fsq=-99.0D0
 
        nfp=1
        deltf=1.0D0
        irst =1
        pexp=4.0D0
        qexp=1.0D0
!**********************************************************************
        lasym = lasym0
        lzakharov = lzsmlonly
	ninit = 0         ! to turn off non-working hot restart
	nresets = 0
	ioflagc = ioflag
	nthetax = nmb
        mpolx = min(17,mpolin)  ! dmc: >17 seems to be unsafe (Feb. 2005)
        if (nmb.gt.ntheta) then
            write(lunmsg,1) nmb,ntheta
 1          format(' Scrunch%% error: # boundary pts (',i4, &
     &      ') exceeds limit (',i4,')')
            CALL scrabwarn(lunmsg,'Scrunch%% error: # boundary pts > limit')
            ier=99
            return
        end if
        if (mpolin.gt.mpol) then
            write(lunmsg,2) mpolin,mpol
 2          format(' Scrunch%% error: # moments )',i4, &
     &      ') exceeds limit (',i4,')')
            CALL scrabwarn(lunmsg,' Scrunch%% error: # moments > limit')
            ier=99
            return
        end if
 
!**********************************************************************
!       MPOL = NUMBER OF POLOIDAL MODES USED IN CURVE-FIT
!
!       MPNT = NUMBER OF R AND Z MODES APPEARING IN FIT
!       NTHETA=NUMBER OF THETA POINTS
!       N2   = TOTAL NUMBER OF RESIDUALS IN FSQ
!
!       STACKING OF XVEC ARRAY (FOR FIXED TOROIDAL PLANE)
!       XVEC(1,mpol):Rmcos            XVEC(1+mpol,2*mpol):Rmsin
!       XVEC(1+2*mpol,3*mpol):Zmcos   XVEC(1+3*mpol,4*mpol):Zmsin
!       XVEC(4*mpol,n2): Theta angle
!**********************************************************************
        twopi=8.D0*atan(1.D0)
!**********************************************************************
!       INITIALIZE FIXED m,n ARRAYS
!**********************************************************************
        if (ninit.eq.0) CALL scrfixaraydes
!**********************************************************************
!       READ IN POINTS ON CURVE TO BE FIT
!**********************************************************************
!       call inptsdes
        do 7 i=1,nmb
            rin3d(i,1) = rmb(i)
            zin3d(i,1) = zmb(i)
 7      continue
!**********************************************************************
!       COMPUTE INITIAL GUESSES (MUST ORDER ANGLES MONOTONICALLY FOR
!       NON-STARLIKE DOMAINS)
!**********************************************************************
        if (ninit.eq.0) CALL scringuessdes (ier)
!**	if (ier.ne.1) return
	if (ier.ne.1) Then
           write(lunmsg,*) ' At error return of descur, ier, fsq =', ier, fsq
           Return
        Endif
!**********************************************************************
!       BEGIN MAIN INTEGRATION LOOP
!**********************************************************************
        CALL cptimr8 (time_in)
        if (ioflagc) write(lunmsg,10)
 10     format(/' ITERATIONS    RMS ERROR    FORCE GRADIENT    <M>', &
     &  '    m #''s <= ON')
 
        do 1000 nplane = 1,nphi2
        if (ioflagc) write(lunmsg,200) nplane
 200    format(/'                  Fitting toroidal plane # ',i3)
!**********************************************************************
!       INITIALIZE M=0 and M=1 MODE AMPLITUDES
!**********************************************************************
	if (ninit.eq.0) then
        do n = 1,n2
        xstore(n)= 0.D0
        xdot(n)  = 0.D0
        xvec(n)  = 0.D0
	end do
      CALL scramplituddes(r0n(nplane),z0n(nplane),angle(1,nplane), &
     &  xvec,xvec(1+mpol),xvec(1+mpol2),xvec(1+mpol3), &
     &  xvec(1+mpol4),rin3d(1,nplane),zin3d(1,nplane),ier)
	end if
!**	if (ier.ne.1) return
	if (ier.ne.1) Then
           write(lunmsg,*) ' At error return of descur, ier, fsq =', ier, fsq
           Return
        Endif
        r10 = .5D0*(abs(xvec(2      )) + abs(xvec(2+ mpol)) &
     &      +     abs(xvec(2+mpol2)) + abs(xvec(2+mpol3)))
        imodes = mpolx
        delt=deltf
	ninit = 1
 
        do 30 iter=1,niter
      CALL scrfunct1des(xvec,xvec(1+mpol),xvec(1+mpol2),xvec(1+mpol3), &
     &  xvec(1+mpol4),gvec,gvec(1+mpol),gvec(1+mpol2), &
     &  gvec(1+mpol3),gvec(1+mpol4),fsq,r10,imodes)
        go to (50,40)iter
        gmin = MIN(gmin,gnorm)
      CALL screvolvedes(g11)
!**********************************************************************
!       RUDIMENTARY TIME STEP CONTROL
!**********************************************************************
        if(gnorm/gmin.gt.1.D6)irst = 2
        if( (irst.eq.2) .or. (gmin.eq.gnorm) ) &
     &     CALL scrrestartdes(irst,ier)
 
        !  DMC, experimenting with remedy for intermittent instability
        !  at higher moment counts

        if(imodes.lt.10) then
           gtest=1.D-4
        else
           gtest=1.D-6
        endif

        if((gnorm.lt.gtest).and.(imodes.lt.mpolx))then

        imodes = imodes + 1
      CALL scrrestartdes(irst,ier)
!**	if (ier.ne.1) return
	if (ier.ne.1) Then
           write(lunmsg,*) ' At error return of descur, ier, fsq =', ier, fsq
           Return
        Endif
        irst = 2
        delt = delt/.95D0
        endif
        if(mod(iter,nstep).eq.0.or.(gnorm.lt.ftol**2))go to 60
        go to 30
 40     g11=gnorm
        go to 30
 50     gmin = gnorm
        imodes = 3
 60     gout = sqrt(gnorm)
        modeno = imodes - 1
        if(iter.eq.1)modeno = mpolx-1
      CALL cptimr8 (time_ou)
        tot_time = tot_time + (time_ou - time_in)
        if (ioflagc) write(lunmsg,110)iter,fsq,gout,specw,modeno
      CALL cptimr8 (time_in)
        if((gnorm.lt.ftol**2.and.imodes.eq.mpolx).or.fsq.lt.ftol) &
     &  go to 80
 30     continue
 110    format(i8,1x,2e16.3,f10.2,i8)
!**********************************************************************
!       STORE FINAL ANSWER FOR POST-PROCESSING
!**********************************************************************
 80     do 120 ntype = 1,4
        ntoff = mpol*(ntype-1)
        do 120 m = 1,mpolx
 120    result(nplane,m,ntype) = xvec(m+ntoff)
 1000   continue
!**********************************************************************
!       FINAL OUTPUT LOOP
!**********************************************************************
        if (ioflagc) write(lunmsg,90)tot_time
 90     format(/,' COMPUTATIONAL TIME = ',1pe12.3,' SECONDS'/)
 320    if (ioflagc) write(lunmsg,330)pexp,qexp
 330    format(' ANGLE CONSTRAINTS WERE APPLIED ',/, &
     &  ' BASED ON RM**2 + ZM**2 SPECTRUM WITH P = ', &
     &  f8.2,' AND Q = ',f8.2,/)
!        call print_it(result(1,1,1),result(1,1,2),
!     >  result(1,1,3),result(1,1,4))
      CALL scrprintdes
        do i=1,mpolx
            rmbout(i-1,1) = rmomb(i,1)
            rmbout(i-1,2) = rmomb(i,2)
            zmbout(i-1,1) = zmomb(i,1)
            zmbout(i-1,2) = zmomb(i,2)
        end do
!       call plotter
 
!!      write(lunmsg,*) ' %nscrunch:  At the end of descur, fsq =', fsq
!!        call scrdesdbg('scrdescur',rmb,zmb,nmb,nmbout,rmbout,zmbout)
 
      End
 
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
