     SUBROUTINE pstmemory(k)
!
! for memory allocation ... aplet  15.11.93
!
 USE pstcom
 USE combla
 USE comfft
 USE l22com
 USE mtrik1
 USE r33com
 USE comggf
 USE newmet
 USE temps
 USE combla
 USE l34com
      
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
 INTEGER :: iok, k
!
 SELECT CASE (k)
!
!23456789012345678901234567890123456789012345678901234567890123456789012

  CASE(1)
!
 ALLOCATE( dpsi(nfe), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate dpsi(nfe)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( dtent(nfe), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate dtent(nfe)'
    stop 'ERROR: cannot allocate in memory' 
 endif
!
 ALLOCATE( psinod(nfe), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate psinod(nfe)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( msub(nfe), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate msub(nfe)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( jsub(nfe), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate jsub(nfe)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( jtot(nfe+2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate jtot(nfe+2)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( lmax(nfe), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate lmax(nfe)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( lmin(nfe), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate lmin(nfe)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
 ALLOCATE( fa(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate fa(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( fpa(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate fpa(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ga(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ga(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( gpa(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate gpa(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( pa(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate pa(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ppa(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ppa(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( psia(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate psia(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( psibig(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate psibig(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( qa(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate qa(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( qpa(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate qpa(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( qppa(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate qppa(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!!$ ALLOCATE( rhoa(nsf0), stat=iok )
!!$ if (iok /= 0 ) then
!!$    print *,'ERROR: cannot allocate rhoa(nsf0)'
!!$    stop 'ERROR: cannot allocate in memory'
!!$ endif
!!$!
!!$ ALLOCATE( alpha(nsf0), stat=iok )
!!$ if (iok /= 0 ) then
!!$    print *,'ERROR: cannot allocate alpha(nsf0)'
!!$    stop 'ERROR: cannot allocate in memory'
!!$ endif
!
 ALLOCATE( di(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate di(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( dr(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate dr(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rgradpsi(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rgradpsi(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ajphi(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ajphi(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( betapol(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate betapol(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( aiaspect(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate aiaspect(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( alqolp(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate alqolp(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( invfft(nths), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate invfft(nths)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( sfft(nths), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate sfft(nths)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( omfft(nths), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate omfft(nths)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( onfft(nths), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate onfft(nths)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
!
 ALLOCATE( scqdp1(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate scqdp1(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( scqdp2(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate scqdp2(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xinf(nths), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xinf(nths)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( zinf(nths), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate zinf(nths)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xwal(nths), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xwal(nths)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( zwal(nths), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate zwal(nths)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( temp1(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate temp1(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( temp2(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate temp2(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( temp3(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate temp3(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( temp4(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate temp4(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( psinew(nsf0+1), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate psinew(nsf0+1)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( psisin(nsg2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate psisin(nsg2)'
    stop 'ERROR: cannot allocate in memory'
 endif

!aplet  17.8.95
 ALLOCATE( sigavr(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate sigavr(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( sigavp(nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate sigavp(nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
 ALLOCATE( qdelp(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate qdelp(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( delta(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate delta(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xjprym(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xjprym(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xdth(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xdth(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( zdth(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate zdth(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xdps(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xdps(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( zdps(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate zdps(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xjacob(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xjacob(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
 ALLOCATE( grpssq(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate grpssq(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( grpsth(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate grpsth(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xsqdps(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xsqdps(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xsq(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xsq(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
! aplet 6/8/96
 ALLOCATE( xa(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xa(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( za(nths,nsf0), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate za(nths,nsf0)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
!
      return
!
 CASE(2)
!
 ALLOCATE( ut(nmtg2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ut(nmtg2)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( vt(nmtg2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate vt(nmtg2)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( u(nfmgb2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate u(nfmgb2)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( va(nfmgb2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate va(nfmgb2)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
      return

CASE(3)

 ALLOCATE( xix0(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xix0(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( chi0(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate chi0(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xilog(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xilog(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( chlog(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate chlog(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xicon(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xicon(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( chcon(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate chcon(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xix1(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xix1(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( chi1(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate chi1(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xix2(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xix2(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( chi2(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate chi2(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xsml0(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xsml0(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( csml0(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate csml0(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xsml1(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xsml1(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( csml1(nfn,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate csml1(nfn,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
 ALLOCATE( gm1(nfn,nfn), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate gm1(nfn,nfn)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ttq(nfn,nfn), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ttq(nfn,nfn)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ttr(nfn,nfn), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ttr(nfn,nfn)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( geem(nfn,nfn), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate geem(nfn,nfn)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( vacmat(nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate vacmat(nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( vacmti(nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate vacmti(nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( alxi1e(nsf0,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate alxi1e(nsf0,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( alxi1o(nsf0,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate alxi1o(nsf0,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( qdchie(nsf0,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate qdchie(nsf0,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( qdchio(nsf0,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate qdchio(nsf0,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ddchie(nsf0,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ddchie(nsf0,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ddchio(nsf0,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ddchio(nsf0,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( gm1dco(nsf0,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate gm1dco(nsf0,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( gm1dce(nsf0,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate gm1dce(nsf0,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( wpot((nfe+1)*nfm*(2*nfm+1)), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate wpot'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( akin((nfe+1)*nfm*(2*nfm+1)), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate akin'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
 ALLOCATE( ry8(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ry8(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( x1frbo(nsf0,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate x1frbo(nsf0,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( x1frbe(nsf0,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate x1frbe(nsf0,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( c1frbo(nsf0,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate c1frbo(nsf0,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( c1frbe(nsf0,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate c1frbe(nsf0,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( x1fpro(nsf0,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate x1fpro(nsf0,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( x1fpre(nsf0,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate x1fpre(nsf0,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( c1fpro(nsf0,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate c1fpro(nsf0,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( c1fpre(nsf0,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate c1fpre(nsf0,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( w0l1ou(nfe,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate w0l1ou(nfe,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( w0l1eu(nfe,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate w0l1eu(nfe,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xisolo(nfe,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xisolo(nfe,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( xisole(nfe,nfm,nsg), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xisole(nfe,nfm,nsg)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rw1(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rw1(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
 ALLOCATE( ry2(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ry2(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rz3(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rz3(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rw4(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rw4(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ry5(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ry5(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rz6(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rz6(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rw7(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rw7(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ry9(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ry9(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rz10(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rz10(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rk1(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rk1(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rk4(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rk4(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( rk7(3,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate rk7(3,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
 ALLOCATE( tw(9,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate tw(9,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( ty(9,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ty(9,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( tz(9,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate tz(9,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( td(9,nfm,nfm), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate td(9,nfm,nfm)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
 ALLOCATE( wreg(nkb,nkc), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate wreg(nkb,nkc)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( wkin(nkb,nkc), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate wkin(nkb,nkc)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( spot(nkc2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate spot(nkc2)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( skin(nkc2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate skin(nkc2)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
 ALLOCATE( amat(nkc,nkc), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate amat(nkc,nkc)'
    stop 'ERROR: cannot allocate in memory'
 endif
!
!
 END SELECT
end SUBROUTINE pstmemory

