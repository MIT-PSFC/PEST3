     SUBROUTINE pstfree_mem(k)
!
! the complementary routine to memory
! a pletzer sep 15 1999
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
!
  CASE(1)
!
 DEALLOCATE( dpsi, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate dpsi'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( dtent, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate dtent'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( psinod, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate psinod'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( msub, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate msub'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( jsub, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate jsub'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( jtot, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate jtot'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( lmax, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate lmax'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( lmin, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate lmin'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( fa, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate fa'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( fpa, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate fpa'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ga, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ga'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( gpa, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate gpa'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( pa, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate pa'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ppa, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ppa'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( psia, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate psia'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( psibig, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate psibig'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( qa, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate qa'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( qpa, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate qpa'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( qppa, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate qppa'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!!$ DEALLOCATE( rhoa, stat=iok )
!!$ if (iok /= 0 ) then
!!$    print *,'ERROR: cannot deallocate rhoa'
!!$    stop 'ERROR: cannot deallocate in FREE_MEM'
!!$ endif
!!$!
!!$ DEALLOCATE( alpha, stat=iok )
!!$ if (iok /= 0 ) then
!!$    print *,'ERROR: cannot deallocate alpha'
!!$    stop 'ERROR: cannot deallocate in FREE_MEM'
!!$ endif
!
 DEALLOCATE( di, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate di'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( dr, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate pbig'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rgradpsi, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rgradpsi'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ajphi, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ajphi'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( betapol, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate betapol'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( aiaspect, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate aiaspect'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( alqolp, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate alqolp'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( invfft, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate invfft'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( sfft, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate sfft'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( omfft, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate omfft'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( onfft, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate onfft'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
!
 DEALLOCATE( scqdp1, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate scqdp1'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( scqdp2, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate scqdp2'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xinf, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xinf'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( zinf, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate zinf'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xwal, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xwal'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( zwal, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate zwal'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( temp1, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate temp1'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( temp2, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate temp2'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( temp3, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate temp3'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( temp4, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate temp4'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( psinew, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate psinew'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( psisin, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate psisin'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!aplet  17.8_r8 .95
 DEALLOCATE( sigavr, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate sigavr'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( sigavp, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate sigavp'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
!
 DEALLOCATE( wreg, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate wreg'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( wkin, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate wkin'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( alxi1e, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate alxi1e'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( alxi1o, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate alxi1o'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( qdchie, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate qdchie'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( qdchio, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate qdchio'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ddchie, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ddchie'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ddchio, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ddchio'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( gm1dco, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate gm1dco'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( gm1dce, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate gm1dce'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( spot, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate spot'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( skin, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate skin'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
 DEALLOCATE( wpot, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate wpot'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( akin, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate akin'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( amat, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate amat'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( xix0, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xix0'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( chi0, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate chi0'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xilog, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xilog'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( chlog, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate chlog'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xicon, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xicon'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( chcon, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate chcon'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xix1, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xix1'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( chi1, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate chi1'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xix2, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xix2'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( chi2, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate chi2'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( qdelp, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate qdelp'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( delta, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate delta'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xjprym, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xjprym'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xdth, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xdth'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( zdth, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate zdth'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xdps, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xdps'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( zdps, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate zdps'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xjacob, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xjacob'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( ry8, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ry8'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( grpssq, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate grpssq'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( grpsth, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate grpsth'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xsml0, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xsml0'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( csml0, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate csml0'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xsml1, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xsml1'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( csml1, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate csml1'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( xsqdps, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xsqdps'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xsq, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xsq'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( gm1, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate gm1'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ttq, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ttq'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ttr, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ttr'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( geem, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate geem'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( vacmat, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate vacmat'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( vacmti, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate vacmti'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
! aplet 6/8/96
 DEALLOCATE( xa, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xa'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( za, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate za'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( x1frbo, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate x1frbo'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( x1frbe, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate x1frbe'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( c1frbo, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate c1frbo'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( c1frbe, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate c1frbe'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( x1fpro, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate x1fpro'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( x1fpre, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate x1fpre'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( c1fpro, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate c1fpro'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( c1fpre, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate c1fpre'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( w0l1ou, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate w0l1ou'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( w0l1eu, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate w0l1eu'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xisolo, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xisolo'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( xisole, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate xisole'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rw1, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rw1'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( ry2, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ry2'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rz3, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rz3'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rw4, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rw4'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ry5, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ry5'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rz6, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rz6'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rw7, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rw7'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ry9, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ry9'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rz10, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rz10'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rk1, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rk1'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rk4, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rk4'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( rk7, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate rk7'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
 DEALLOCATE( tw, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate tw'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( ty, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ty'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( tz, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate tz'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( td, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate td'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
      return
!
 CASE(2)
!
 DEALLOCATE( ut, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate ut'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( vt, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate vt'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( u, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate u'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
 DEALLOCATE( va, stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot deallocate va'
    stop 'ERROR: cannot deallocate in FREE_MEM'
 endif
!
!
      return
!
 END SELECT
end SUBROUTINE pstfree_mem
