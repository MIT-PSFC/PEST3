subroutine rfx_antcen(nant,ra,za,rcen,zcen,hh)

  ! compute antenna center and half-height from array info
  ! hh = (zmax-zmin)/2; rcen = r value of location in ra(:),za(:) closest to
  !                            center of structure

  ! mod DMC May 2008: hh is now half the distance btw end pts; 
  !   (rcen,zcen) is the center -- allow off-midplane antenna.

  implicit NONE

  !---------
  integer, intent(in) :: nant
  real*8, intent(in) :: ra(nant),za(nant)
  real*8, intent(out) :: rcen,zcen
  real*8, intent(out) :: hh

  !---------
  real*8 :: zcmin,zctest,zrr,zzz
  integer :: k,icmin

  !---------

  zcmin=(ra(nant)-ra(1))**2 + (za(nant)-za(1))**2
  hh = zcmin/2

  icmin=1

  do k=2,nant
     zctest = (ra(1)+ra(nant)-2*ra(k))**2 + (za(1)+za(nant)-2*za(k))**2
     if(zctest.lt.zcmin) then
        zcmin=zctest
        icmin=k
     endif
  enddo

  rcen=ra(icmin)
  zcen=za(icmin)

end subroutine rfx_antcen
