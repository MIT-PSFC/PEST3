subroutine pstMetRealloc(nt5, ns, ier)
  !
  ! re-allocate memory for metric quantities
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
 implicit none
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
  integer, intent(in) :: nt5, ns
  integer, intent(out) :: ier

  integer iok

  ier = 0
!
!!$ IF(ALLOCATED(dpsi)) DEALLOCATE(dpsi, STAT=iok)
!!$ ALLOCATE(dpsi(nfe), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate dpsi'
!!$
!!$ IF(ALLOCATED(dtent)) DEALLOCATE(dtent, STAT=iok)
!!$ ALLOCATE(dtent(nfe), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate dtent'
!!$
!!$ IF(ALLOCATED(psinod)) DEALLOCATE(psinod, STAT=iok)
!!$ ALLOCATE(psinod(nfe), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate psinod'
!!$
!!$ IF(ALLOCATED(msub)) DEALLOCATE(msub, STAT=iok)
!!$ ALLOCATE(msub(nfe), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate msub'
!!$
!!$ IF(ALLOCATED(jsub)) DEALLOCATE(jsub, STAT=iok)
!!$ ALLOCATE(jsub(nfe), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate jsub'
!!$
!!$ IF(ALLOCATED(jtot)) DEALLOCATE(jtot, STAT=iok)
!!$ ALLOCATE(jtot(nfe), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate jtot'
!!$
!!$ IF(ALLOCATED(lmax)) DEALLOCATE(lmax, STAT=iok)
!!$ ALLOCATE(lmax(nfe), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate lmax'
!!$
!!$ IF(ALLOCATED(lmin)) DEALLOCATE(lmin, STAT=iok)
!!$ ALLOCATE(lmin(nfe), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate lmin'

 IF(ALLOCATED(fa)) DEALLOCATE(fa, STAT=iok)
 ALLOCATE(fa(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate fa'

 IF(ALLOCATED(fpa)) DEALLOCATE(fpa, STAT=iok)
 ALLOCATE(fpa(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate fpa'

 IF(ALLOCATED(ga)) DEALLOCATE(ga, STAT=iok)
 ALLOCATE(ga(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate ga'

 IF(ALLOCATED(gpa)) DEALLOCATE(gpa, STAT=iok)
 ALLOCATE(gpa(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate gpa'

 IF(ALLOCATED(pa)) DEALLOCATE(pa, STAT=iok)
 ALLOCATE(pa(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate pa'

 IF(ALLOCATED(ppa)) DEALLOCATE(ppa, STAT=iok)
 ALLOCATE(ppa(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate ppa'

!!$ IF(ALLOCATED(psia)) DEALLOCATE(psia, STAT=iok)
!!$ ALLOCATE(psia(ns), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate psia'
!!$
!!$ IF(ALLOCATED(psibig)) DEALLOCATE(psibig, STAT=iok)
!!$ ALLOCATE(psibig(ns), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate psibig'
!!$
 IF(ALLOCATED(qa)) DEALLOCATE(qa, STAT=iok)
 ALLOCATE(qa(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate qa'

 IF(ALLOCATED(qpa)) DEALLOCATE(qpa, STAT=iok)
 ALLOCATE(qpa(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate qpa'

 IF(ALLOCATED(qppa)) DEALLOCATE(qppa, STAT=iok)
 ALLOCATE(qppa(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate qppa'

 IF(ALLOCATED(di)) DEALLOCATE(di, STAT=iok)
 ALLOCATE(di(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate di'

 IF(ALLOCATED(dr)) DEALLOCATE(dr, STAT=iok)
 ALLOCATE(dr(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate dr'

 IF(ALLOCATED(rgradpsi)) DEALLOCATE(rgradpsi, STAT=iok)
 ALLOCATE(rgradpsi(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate rgradpsi'

 IF(ALLOCATED(ajphi)) DEALLOCATE(ajphi, STAT=iok)
 ALLOCATE(ajphi(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate ajphi'

 IF(ALLOCATED(betapol)) DEALLOCATE(betapol, STAT=iok)
 ALLOCATE(betapol(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate betapol'

 IF(ALLOCATED(aiaspect)) DEALLOCATE(aiaspect, STAT=iok)
 ALLOCATE(aiaspect(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate aiaspect'

 IF(ALLOCATED(alqolp)) DEALLOCATE(alqolp, STAT=iok)
 ALLOCATE(alqolp(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate alqolp'

!!$ IF(ALLOCATED(invfft)) DEALLOCATE(invfft, STAT=iok)
!!$ ALLOCATE(invfft(nt5), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate invfft'
!!$
!!$ IF(ALLOCATED(sfft)) DEALLOCATE(sfft, STAT=iok)
!!$ ALLOCATE(sfft(nt5), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate sfft'
!!$
!!$ IF(ALLOCATED(omfft)) DEALLOCATE(omfft, STAT=iok)
!!$ ALLOCATE(omfft(nt5), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate omfft'
!!$
!!$ IF(ALLOCATED(onfft)) DEALLOCATE(onfft, STAT=iok)
!!$ ALLOCATE(onfft(nt5), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate onfft'
!!$
!!$ IF(ALLOCATED(scqdp1)) DEALLOCATE(scqdp1, STAT=iok)
!!$ ALLOCATE(scqdp1(ns), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate scqdp1'
!!$
!!$ IF(ALLOCATED(scqdp2)) DEALLOCATE(scqdp2, STAT=iok)
!!$ ALLOCATE(scqdp2(ns), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate scqdp2'
!!$
 IF(ALLOCATED(xinf)) DEALLOCATE(xinf, STAT=iok)
 ALLOCATE(xinf(nt5), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate xinf'

 IF(ALLOCATED(zinf)) DEALLOCATE(zinf, STAT=iok)
 ALLOCATE(zinf(nt5), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate zinf'

 IF(ALLOCATED(xwal)) DEALLOCATE(xwal, STAT=iok)
 ALLOCATE(xwal(nt5), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate xwal'

 IF(ALLOCATED(zwal)) DEALLOCATE(zwal, STAT=iok)
 ALLOCATE(zwal(nt5), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate zwal'

!!$ IF(ALLOCATED(temp1)) DEALLOCATE(temp1, STAT=iok)
!!$ ALLOCATE(temp1(ns), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate temp1'
!!$
!!$ IF(ALLOCATED(temp2)) DEALLOCATE(temp2, STAT=iok)
!!$ ALLOCATE(temp2(ns), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate temp2'
!!$
!!$ IF(ALLOCATED(temp3)) DEALLOCATE(temp3, STAT=iok)
!!$ ALLOCATE(temp3(ns), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate temp3'
!!$
!!$ IF(ALLOCATED(temp4)) DEALLOCATE(temp4, STAT=iok)
!!$ ALLOCATE(temp4(ns), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate temp4'

!!$ IF(ALLOCATED(psinew)) DEALLOCATE(psinew, STAT=iok)
!!$ ALLOCATE(psinew(ns), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate psinew'
!!$
!!$ IF(ALLOCATED(psisin)) DEALLOCATE(psisin, STAT=iok)
!!$ ALLOCATE(psisin(nsg2), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate psisin'

 IF(ALLOCATED(sigavr)) DEALLOCATE(sigavr, STAT=iok)
 ALLOCATE(sigavr(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate sigavr'

 IF(ALLOCATED(sigavp)) DEALLOCATE(sigavp, STAT=iok)
 ALLOCATE(sigavp(ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate sigavp'

 IF(ALLOCATED(qdelp)) DEALLOCATE(qdelp, STAT=iok)
 ALLOCATE(qdelp(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate qdelp'

 IF(ALLOCATED(delta)) DEALLOCATE(delta, STAT=iok)
 ALLOCATE(delta(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate delta'

 IF(ALLOCATED(xjprym)) DEALLOCATE(xjprym, STAT=iok)
 ALLOCATE(xjprym(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate xjprym'

 IF(ALLOCATED(xdth)) DEALLOCATE(xdth, STAT=iok)
 ALLOCATE(xdth(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate xdth'

 IF(ALLOCATED(zdth)) DEALLOCATE(zdth, STAT=iok)
 ALLOCATE(zdth(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate zdth'

 IF(ALLOCATED(xdps)) DEALLOCATE(xdps, STAT=iok)
 ALLOCATE(xdps(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate xdps'

 IF(ALLOCATED(zdps)) DEALLOCATE(zdps, STAT=iok)
 ALLOCATE(zdps(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate zdps'

 IF(ALLOCATED(xjacob)) DEALLOCATE(xjacob, STAT=iok)
 ALLOCATE(xjacob(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate xjacob'

 IF(ALLOCATED(grpssq)) DEALLOCATE(grpssq, STAT=iok)
 ALLOCATE(grpssq(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate grpssq'

 IF(ALLOCATED(grpsth)) DEALLOCATE(grpsth, STAT=iok)
 ALLOCATE(grpsth(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate grpsth'

 IF(ALLOCATED(xsqdps)) DEALLOCATE(xsqdps, STAT=iok)
 ALLOCATE(xsqdps(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate xsqdps'

 IF(ALLOCATED(xsq)) DEALLOCATE(xsq, STAT=iok)
 ALLOCATE(xsq(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate xsq'

 IF(ALLOCATED(xa)) DEALLOCATE(xa, STAT=iok)
 ALLOCATE(xa(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate xa'

 IF(ALLOCATED(za)) DEALLOCATE(za, STAT=iok)
 ALLOCATE(za(nt5, ns), STAT=iok)
 IF(iok/=0) PRINT*,'***ERROR** could not allocate za'

!!$ IF(ALLOCATED(ut)) DEALLOCATE(ut, STAT=iok)
!!$ ALLOCATE(ut(nmtg2), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate ut'
!!$
!!$ IF(ALLOCATED(vt)) DEALLOCATE(vt, STAT=iok)
!!$ ALLOCATE(vt(nmtg2), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate vt'
!!$
!!$ IF(ALLOCATED(u)) DEALLOCATE(u, STAT=iok)
!!$ ALLOCATE(u(nfmgb2), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate u'
!!$
!!$ IF(ALLOCATED(va)) DEALLOCATE(va, STAT=iok)
!!$ ALLOCATE(va(nfmgb2), STAT=iok)
!!$ IF(iok/=0) PRINT*,'***ERROR** could not allocate va'
  
 ier = iok

end subroutine pstMetRealloc
