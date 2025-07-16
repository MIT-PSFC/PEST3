!................................
subroutine pstshape(arat,xrad,ellip,pdelta &
     ,zzmin,zzmax,zxmin,zxmax)
  !......................
  !
  ! find min and max
  !
  USE pstcom
  USE mtrik1
  IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
  real*8 ARAT
  real*8 XRAD
  real*8 ELLIP
  real*8 PDELTA
  real*8 ZZMIN
  real*8 ZZMAX
  real*8 ZXMIN
  real*8 ZXMAX
  real*8 X0
  real*8 ZX1
  real*8 ZX2

  !
  INTEGER, dimension(1)	 :: imin, imax

  zzmin = minval( zinf(1:mth) )
  zzmax = maxval( zinf(1:mth) )
  zxmin = minval( xinf(1:mth) )
  zxmax = maxval( xinf(1:mth) )

  xrad = (zxmax - zxmin)/ 2._r8 
  x0 = (zxmax+zxmin)/ 2._r8 
  arat = xrad / x0
  ellip = (zzmax-zzmin)/(zxmax-zxmin)
  imin = minloc(zinf(1:mth))
  imax = maxloc(zinf(1:mth))
  zx1 = xinf( imax(1) )
  zx2 = xinf( imin(1) )

  pdelta = (zxmax+zxmin-zx1-zx2)/(zxmax-zxmin)
  !
end subroutine pstshape


