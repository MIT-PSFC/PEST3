!....................................................................
      SUBROUTINE pstremap2
!....................................................................
!
!      this routine interpolates equilibrium data from the psi mesh
!      given by psia(i),i = 1,nosurf  to the psi mesh
!      given by psinew(i),i=ist, nusurf.
!
 USE pstcom

 USE combla

 USE mtrik1

 USE newmet

 USE l22com

 USE temps

 use i2mex_mod

!     include 'l33com.inc'  
 USE r33com

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER I
      INTEGER NRAT
      INTEGER J
      INTEGER LGIVUP
      real*8 G
      real*8 G2
      real*8 PP
      real*8 ZY
      real*8 TOPS
      real*8 TOTCUR
      real*8 ZHL
      real*8 ZHR
      real*8 ZDPSI
      real*8 B2SUM
      real*8 VTSUM
      real*8 SSUM
      real*8 SJPHI
      real*8 DSUBR
      real*8 DSUBEE
      real*8 BP2SUM
      real*8 RGRPS
      real*8 XMAX
      real*8 XMIN
      real*8 ZMIN
      real*8 ZMAX
      real*8 ZJACOB
      real*8 ZBTHETA
      real*8 ZMAJR
      real*8 ZMINR
      real*8 B0SC
      INTEGER NR
      real*8 ZQ0S
      real*8 ZLS
      real*8 ZRO
      real*8 ZETA
      real*8 ZS1
      real*8 ZETA3
      real*8 PSISN
      real*8 ZDRMIN
      real*8 ZDELC
      real*8 ZQR
      INTEGER IER
      real*8 ZG34
      real*8 GAMMAR
      real*8 ZG14
      real*8 ZPI8
      real*8 ZCONST
      real*8 ZQC

      real*8 ze(nusurf), zf(nusurf), zh(nusurf)


real(r8), dimension(:), allocatable :: the

!
!     parameter(nsf2 = 2*nsf)
!
!
!
!      interpolate all radial profiles...
!
!     CALL pstmemory(2)
!
 LOGICAL	 	 	 :: ldelta
 LOGICAL	 	 	 :: lmapck
 LOGICAL	 	 	 :: lvacdsk
 REAL*8	 	 	 :: psis
 REAL*8	 	 	 :: vp
 REAL*8	 	 	 :: qs
 REAL*8	 	 	 :: qps
 REAL*8	 	 	 :: pps
 REAL*8	 	 	 :: e
 REAL*8	 	 	 :: f
 REAL*8	 	 	 :: h
 REAL*8	 	 	 :: em
 REAL*8	 	 	 :: ak
 REAL*8	 	 	 :: gee
 REAL*8	 	 	 :: di1
 REAL*8	 	 	 ::  dr1
 REAL*8	 	 	 ::  feq
 REAL*8	 	 	 ::  zor2
 REAL*8	 	 	 ::  zvp
 REAL*8	 	 	 ::  minorRad4pt
 INTEGER	 	 	 :: iok
 INTEGER, DIMENSION(2) 	 :: iop
 COMMON /c2f90ldelta/ ldelta 
 COMMON /c2f90lmapck/ lmapck 
 COMMON /c2f90lvacdsk/ lvacdsk 
 COMMON /c2f90psis/ psis 
 COMMON /c2f90vp/ vp 
 COMMON /c2f90qs/ qs 
 COMMON /c2f90qps/ qps 
 COMMON /c2f90pps/ pps 
 COMMON /c2f90e/ e 
 COMMON /c2f90f/ f 
 COMMON /c2f90h/ h 
 COMMON /c2f90em/ em 
 COMMON /c2f90ak/ ak 
 COMMON /c2f90gee/ gee 
!
      temp1(1) =  0.0_r8 
      temp1(nosurf) =  0.0_r8 
      n1surf = nosurf - 1
!

  nsf0 = nusurf
  !print*,'*** remap *** nusurf=',nusurf
  nsf   = nsf0
  nsf2 = 2*nsf
  !print*,'*** remap *** nsf2=',nsf2
  nfe = 1 + nsf/2
  nfe1=nfe+1
  nfmg   = nfm
  nfmga  = ipolp2 * nfmg+1
  !print*,'ipolp2=',ipolp2
  nfmgb  = ipolp2 * 2 * nfmga
  nfmgb2 = nfmgb
! re-allocate

call pstMetRealloc(nths, nsf, ier)
if(ier/=0) then
   print*,'***ERROR** failed to reallocate memory for metric quantities'
endif


! interpolate equilbrium onto new radial psinew mesh

allocate(the(nths0+1), stat=iok)
call i2mex_getOriT(nths0+1, the, iok)
call i2mex_error(iok)

!!$the = meq_o%x_spl%x1
call pstCompMetric(nths0+1, nusurf, the, psinew, ier)
  if(ier/=0) then
     print*,'***ERROR** occurred in pstRemap2 after pstCompMetric'
  endif
deallocate(the, stat=iok)
!
!
!
      call i2mex_getGlasserGreenJohnsonEFH(nusurf, psinew, ze, zf, zh, ier)
      dr = ze + zf + zh**2
      di = ze + zf + zh - 0.25_r8

      call pstToAxis1(psinew, di)
      call pstToAxis1(psinew, dr)

      if (lextra  .AND.  (imode >  1)) go to 270

      !
      ! print out remapped quantities...
      !

tops =  0.0_r8 
totcur =  0.0_r8 
do surf = 2,nusurf ! changed from 1..nosurf 10/10/99
      zhl =  0._r8 
      if( surf >   1    )  &
           zhl = (psinew(surf) - psinew(surf-1))/ 2.0_r8 
      zhr =  0._r8 
      if( surf < nosurf)  &
           zhr = (psinew(surf+1) - psinew(surf))/ 2.0_r8 
      zdpsi = zhl + zhr
	CALL pstdsubi(surf,b2sum,vtsum,ssum,sjphi,dsubr,dsubee,bp2sum,rgrps)
!!$	di(surf) = dsubee
!!$	dr(surf) = dsubr
	ajphi(surf) = -sjphi
	xmax = maxval( xa(1:mth,surf) )
	xmin = minval( xa(1:mth,surf) )
        zmin = minval( za(1:mth,surf) )
        zmax = maxval( za(1:mth,surf) )
        zjacob = sum( xjacob(1:mth,surf) )/mth
        zbtheta = sqrt( (xmax-xmin)*(zmax-zmin) )/(2 * zjacob)
	aiaspect(surf) = (xmax-xmin)/(xmax+xmin)
	alqolp(surf) =  0.0_r8 
if( abs(pa(surf)*qpa(surf)) >  1.E-10_r8  ) then
	alqolp(surf) = - ppa(surf)*qa(surf)/(pa(surf)*qpa(surf))
end if
        rgradpsi(surf) = rgrps

        tops = tops + ssum * zdpsi
        totcur = totcur - twopi * sjphi * zdpsi
!
! change def of poloidal beta_p from integral over p' to local def
!
!       betapol(surf) =  4._r8  * twopi * tops * pa(surf) / totcur**2
        betapol(surf) =  2._r8 * pa(surf) /zbtheta**2

end do

call pstToAxis1(psinew, ajphi)
call pstToAxis1(psinew, aiaspect)
call pstToAxis1(psinew, alqolp)
call pstToAxis1(psinew, rgradpsi)
call pstToAxis1(psinew, betapol)

do i = 1, nosing
   !
   nrat = 1
   do j = 1, i
      nrat = nrat + (msing(j+1)-msing(j))*msub(j)
   end do
   nr = nrat + msub(i)/2

   !
   CALL pstggjcoef(nr,e,f,h,di1,dr1,ak,gee,em,feq,zor2,zvp,zq0s,zls)

   ! get the four point minor radius which takes the elongation into account
   xmax = maxval( xa(1:mth,nusurf) )
   xmin = minval( xa(1:mth,nusurf) )
   zmin = minval( za(1:mth,nusurf) )
   zmax = maxval( za(1:mth,nusurf) )
   minorRad4pt=(xmax-xmin+zmax-zmin)/4
  
   ! now get the same (and the btheta) at the rational surface
   aklay(i) = ak
   gelay(i) = gee
   qslay(i) = qa(nr)
   qpslay(i) = qpa(nr)
   xmax = maxval( xa(1:mth,nr) )
   xmin = minval( xa(1:mth,nr) )
   zmin = minval( za(1:mth,nr) )
   zmax = maxval( za(1:mth,nr) )
   zjacob = sum( xjacob(1:mth,nr) )/mth
   zbtheta = sqrt( (xmax-xmin)*(zmax-zmin) )/(2 * zjacob)

   !write(99,*) 'Qo and  Xo do not have the eta^1/3 factor multiplied in'
   !write(99,*) 'H     Qo     Xo'
   write(outilc,'(18e18.10)') majorRadius,sqrt(b0SquareCentre),zbtheta, &
                   minorRad4pt,(xmax-xmin+zmax-zmin)/4, xmax-majorRadius, &
                   0.0_r8, 0.0_r8, e, f, gee, h, ak, feq, zor2, zvp, zq0s, zls
enddo
CALL pstempty(outilc)
!
  270 continue
!
!       save metric for use again later.
!
      return
!
 7000 CALL psterrmes(outpst,'remap',lgivup)
      return
      end



