subroutine pstSetBetmsh(p)

  ! Set global mesh packing exponent. For instance a value
  ! > 1 increases the density of nodes near the edge.

  use comggf
  implicit none
  real*8, intent(in) :: p

  betmsh = p
end subroutine pstSetBetmsh

subroutine pstSetAlfmsh(pdim,p)

  ! Set global mesh packing exponent. For instance a value
  ! > 1 increases the density of nodes.

  use comggf
  implicit none
  integer, intent(in) :: pdim ! no of resonant surfaces + 2 (axis,edge)
  real*8, intent(in) :: p(*)
 
  alfmsh(1) = p(1)
  alfmsh(pdim+2) = p(2) ! here the alfmsh of the edge is moved to its index.
end subroutine pstSetAlfmsh

subroutine pstSetWidmsh(pdim,p)

  ! Set global mesh packing zone width.  The width of the packed
  ! zones in normalized flux space.  Should be <1 for each.  ~0.1

  use comggf
  implicit none
  integer, intent(in) :: pdim ! no of resonant surfaces + 2 (axis,edge)
  real*8, intent(in) :: p(*)

  widmsh(1:pdim) = p(1:pdim)

end subroutine pstSetWidmsh

subroutine pstSetPackpts(pdim,p)

  ! Set the numer of points packed on either side of each surface within dlay.
  ! Should be >=1 for each. 1 leaves as only the 2 points at dlay.

  use comggf
  implicit none
  integer, intent(in) :: pdim ! no of resonant surfaces + 2 (axis,edge)
  integer, intent(in) :: p(*)

  packpts(1:pdim) = p(1:pdim)

end subroutine pstSetPackpts

subroutine pstSetPackmth(p)

  ! Set global mesh packing method.  0 packs at surfaces and axis,
  ! 1 includes the edge.  See pest3.hh

  use comggf
  implicit none
  integer, intent(in) :: p

  packmth = p
end subroutine pstSetPackmth


subroutine pstSetLsymhi(k)

  ! Set switch for large solution support width. A value of 1 
  ! forces the support to be symmetric about singular surfaces.
  ! Default is asymmetric and determined by distance between the
  ! singular surfaces.

  use comggf
  implicit none
  integer, intent(in) :: k

  lsymhi = .false.
  if(k /= 0) lsymhi = .true.
end subroutine pstSetLsymhi

subroutine pstSetDlayb(p)

  ! Set the normalized large support width. Should be < 1.0, that
  ! is the the large soutions should not overlap neighboring singular
  ! surfaces.

  use comggf
  implicit none
  real*8, intent(in) :: p

  dlayb = p
end subroutine pstSetDlayb

subroutine pstSetDlay(p)

  ! Set the normalized large support width. Should be < 1.0, that
  ! is the the large soutions should not overlap neighboring singular
  ! surfaces.

  use pstcom
  implicit none
  real*8, intent(in) :: p

  dlay = p
end subroutine pstSetDlay

subroutine pstSetUniformMesh(k)

  ! Force the mesh to be uniform. Default is packed mesh
  ! near the axis and about singular surfaces.

  use comggf
  integer, parameter :: r8=selected_real_kind(12,100)
  integer, intent(in) :: k

  LPACK=.TRUE.
  if(k==1) then
     LPACK=.FALSE.
     alfmsh = 1.0_r8
  endif
end subroutine pstSetUniformMesh


subroutine pstSetIsolver(k)

  ! Select solver: 0 for built-in, > 0 for
  ! LAPACK.

  use pstcom
  implicit none
  integer, intent(in) :: k

  isolver = k
end subroutine pstSetIsolver

subroutine pstSetB(pb)

  ! Set the wall distance normalized to the minor radius from 
  ! the plasma edge: PB <= 0 for wall on plasma, 1.0 > PB <=10.0 
  ! for conformal wall, and PB > 10.0 for wall at infinity.

  use pstcom
  use l22com
  implicit none
  real*8, intent(in) :: pb

  wall = .false.
  infwal = .false.
  b = pb
  if(b <= 0.0) wall = .true.
  if(b > 10.0) infwal = .true.
end subroutine pstSetB

subroutine pstSetN(j)
  use pstcom
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)

  ! Set the toroidal mode number.

  integer, intent(in) ::  j

  n = real(j, r8)
  return
end subroutine pstSetN

subroutine pstSetMM(jdim, j)
  use pstcom
  implicit none

  ! Set the number or radial nodes. For a convergence study take
  ! J=[N1, N2, ... 0] (last element must be zero). 

  integer, intent(in) :: jdim ! no of runs
  integer, intent(in) :: j(*) ! no of radial elements in run

  integer jsize

  mm = 0
  if(jdim > nrun) print*,'**WARNING** ajusted mm dimension in pstSetMM'
  jsize = min(nrun, jdim)

  mm(1:jsize) = j(1:jsize)
  m = mm(1)
end subroutine pstSetMM


subroutine pstSetMSIN(jdim, j)
  use pstcom
  implicit none

  ! Select the singular surfaces. J is typically an array of 0's and
  ! 1's: 0's for rational surfaces that are skipped and 1's for 
  ! rational surfaces that are singular, counting from the magnetic
  ! axis outwards. For instance, to select the 3/2 and 4/2 surfaces 
  ! for a monotonically increasing q profile from 0.9 to 6.4, J should
  ! be set to [0, 1, 1, 0,...]. 
  

  integer, intent(in) :: jdim ! no of rational surfaces
  integer, intent(in) :: j(*) ! activation flag for each rational surface.

  msin = 0
  if(jdim > nsg) then
     print*,'pstSetMSIN: WARNING jdim > nsg will truncate'
     msin(1:nsg) = j(1:nsg)
  else
     msin(1:jdim) = j(1:jdim)
  endif

end subroutine pstSetMSIN

subroutine pstSetWALL(j)

  ! Set the conducting wall on the plasma edge. This option overrides
  ! pstSetB: 1=wall on plasma.

  use pstcom
  implicit none
  integer, intent(in) ::  j

  wall=.TRUE.
  if(j==0) wall = .FALSE.

  return
end subroutine pstSetWALL

subroutine pstSetINFWALL(j)
  use pstcom
  implicit none

  ! Set wall at inifinity. This option overrides
  ! pstSetB: 1=wall at infinity.

  integer, intent(in) :: j

  infwal=.TRUE.
  if(j==0) infwal = .FALSE.

  return
end subroutine pstSetINFWALL
