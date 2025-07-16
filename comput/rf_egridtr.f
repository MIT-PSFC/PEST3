c
      subroutine rf_egridtr(eminev,zdemin,egrmax,maxen,enerev)
c
c  REAL*8 interface
c
      implicit NONE
      real*8, intent(in) :: eminev      ! min energy (eV) (usually zero)
      real*8, intent(in) :: zdemin      ! energy step at eminev
      real*8, intent(in) :: egrmax      ! max energy (eV)
      integer, intent(in) :: maxen      ! no. of grid zones
      real*8, intent(out) :: enerev(0:maxen) ! energy zone *bdys*
c
c  note 0:maxen indexing of enerev!
c
c------------------------------------------
c  for TRANSP (dmc Feb 96) control the max energy of the non-uniform
c  grid -- mod dmc Jan 1997)
c
c  the energy grid is linearly spaced in steps of (zdemin) until the
c  first step of a logarithmic spacing from the top of the linearly
c  spaced region is larger than the linear step (zdemin) would be.
c
c------------------------------------------
c
      real*8 zEbase,zErat,deove,zlebase,zlemax,zle
      integer i,i0,ileft
c
      real*8, parameter :: ZERO=0.0d0
      real*8, parameter :: ONE=1.0d0
c
c------------------------------------------
c
      enerev(0)=max(ZERO,eminev)
      if(enerev(0).eq.ZERO) then
         enerev(1)=zdemin
         i0=1
      else
         i0=0
      endif

      !  linear part-- until log spacing yields bigger step...

      do
         zEbase=enerev(i0)
         zErat=egrmax/zEbase
         ileft=maxen-i0
         deove=exp(log(zErat)/ileft)
         if(zEbase*(deove-ONE).lt.zdemin) then
            enerev(i0+1)=enerev(i0)+zdemin
            i0=i0+1
         else
            exit
         endif
      enddo

      zlebase = log(zEbase)
      zlemax = log(egrmax)
      do i=i0+1,maxen
         zle=(zlebase*(maxen-i)+zlemax*(i-i0))/(maxen-i0)
         enerev(i)=exp(zle)
      enddo
c
      return
      end
c-----------------------------
      subroutine rf_egridtr_real(eminev,zdemin,egrmax,maxen,enerev)
c
c  REAL interface
c
      implicit NONE
      real, intent(in) :: eminev        ! min energy (eV) (usually zero)
      real, intent(in) :: zdemin        ! energy step at eminev
      real, intent(in) :: egrmax        ! max energy (eV)
      integer, intent(in) :: maxen      ! no. of grid zones
      real, intent(out) :: enerev(0:maxen) ! energy zone *bdys*
c
      real*8 zeminev,zegrmax,zzdemin
      real*8 zenerev(0:maxen)
c
      zeminev=eminev
      zzdemin=zdemin
      zegrmax=egrmax
c
      call rf_egridtr(zeminev,zzdemin,zegrmax,maxen,zenerev)
c
      enerev=zenerev
c
      return
      end
