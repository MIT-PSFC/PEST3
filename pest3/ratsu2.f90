!......................................................................
      SUBROUTINE pstratsu2
!
! modified on  9.5.94 to treat non-monotonic q profiles (aplet)
!......................................................................
 USE pstcom
 USE comggf
 IMPLICIT NONE
  integer, parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 ZQMIN
      real*8 ZQMAX
      INTEGER I, jmin, jmax, jpoints
      real*8 ZPSI
      real*8 ZPSL
      real*8 ZPSR
      INTEGER JS, jsl, jsr, jcount, j, jleft, jrigh
      real*8, dimension(:), allocatable :: zs, ze, zf, zh

!
      
 real*8, allocatable :: zpsia(:), zqa(:)
 integer jer, jn, jns

call i2mex_getOriNs(jns, jer)
allocate(zpsia(jns), zqa(jns))

call i2mex_getOriPsi(jns, zpsia, jer)
call i2mex_getQ(jns, zpsia, zqa, jer)

zqmin = minval( zqa(1:jns) )
zqmax = maxval( zqa(1:jns) )
jmax = min( floor(n* zqmax), lmax(1))
jmin = max( ceiling(n* zqmin), lmin(1))

! find all rational surfaces with resonant poloidal mode  jmin <= ms <= jmax

zpsl = zpsia(1)
zpsr = zpsia(jns)
psisin = 0._r8
psisin(1)  = zpsia(1)
jn = int(n)
nosing = 0
jcount = 0

! Look for rational surfaces by scanning q from the axis to the 
! edge. An acceptable rational surface ms/n is such that ms is 
! part of the Fourier spectrum. In the end, whether a rational 
! is actually singular will depend on msin. Msin should be /= 0
! to be selected as a singular surface.

do j = 1, jns-1
   
   if(nosing == nsg2-2) exit 

   jsl = floor(n*zqa(j  ))
   jsr = floor(n*zqa(j+1))
   if(jsr > jmax .or. jsl < jmin-1) cycle
   if(jsl /= jsr) then
      ! Check to see that two rational surfaces arent within one grid step. If
      ! so then decrement the js
      if (abs(jsr-jsl)>1) then
        if (jsr>jsl) then
          jsr=jsl+1
        else 
          jsl=jsr+1 
        endif
      endif
      js = max(jsl, jsr)
      ! found a rational surface!
      jcount = jcount + 1
      if(jcount > nsg) then
         if (nosing > 0) then
            ! ok got least one singular surface
            print *,'WARNING: too many rational surfaces in ratsu2. Try increasing nsg'
            exit
         else
            ! hopeless
            stop 'ERROR: too many rational surfaces in ratsu2. Try increasing nsg'
         endif
      endif
      ! do we really want it?
      if(msin(jcount)==1 .or. (msin(1)==js .and. nosing==0)) then
         ! go locate, exactly
         jleft = max(j-1,1)
         jrigh = min(j+4, jns)
         call i2mex_findRationalSurface(js, jn, &
              & zpsia(jleft), zpsia(jrigh), zpsi, jer)
         call i2mex_error(jer)
         if(jer==0) then
            ! yes
            nosing = nosing + 1
            psisin(nosing+1) = zpsi 
         else
            stop 'error'
         endif
      endif
   endif

enddo

nsing1 = nosing + 1
nsing2 = nsing1 + 1
psisin(nsing2) = zpsia(jns)
jpoints = nosing + 4

allocate(zs(jpoints), ze(jpoints), zf(jpoints), zh(jpoints))

zs(1:3) = psisin(1) + (psisin(2)-psisin(1))* &
     & (/ (real(i-1,r8)/real(4-1,r8), i=1,3) /)
zs(4:jpoints) = psisin(2:nsing2)

call i2mex_getGlasserGreenJohnsonEFH(jpoints, zs, ze, zf, zh, jer)
call i2mex_error(jer)

eflay(1:nosing) = zf(4:3+nosing)
 hlay(1:nosing) = zh(4:3+nosing)
drlay(1:nosing) = ze(4:3+nosing) + zf(4:3+nosing) + zh(4:3+nosing)**2
  xmu           = sqrt(abs(-drlay + (hlay-0.5_r8)**2))

deallocate(zs, ze, zf, zh)
deallocate(zpsia, zqa)

end subroutine pstratsu2



