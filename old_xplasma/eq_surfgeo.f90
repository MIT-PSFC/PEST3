!
! compute various flux surface parameters derived from the position of 
! the flux surface on the R,Z grid
!
! the midplane intercepts are defined as the major radiaus at the height of the
! centroid of the flux surface.
!
! the rho value will be forced to be larger then  1.e-5*(rho_bdy-rho_axis)
!
! notes from gelong:
!c
!c> Date: Wed, 22 Sep 1999 11:59:40 -0400
!c> From: Glenn Bateman <bateman@fusion.physics.lehigh.edu>
!c> 	Consider a toroidal magnetic surface (also called a flux surface).
!c> Use the following variables
!c>
!c> R_out   = major radius to the outboard edge of the flux surface
!c>           that is, the largest value of major radius anywhere on the surface
!c> R_in    = major radius to the inboard edge (ie, smallest value of major
!c> radius)
!c> R_top   = major radius to the top edge of the flux surface
!c>           that is, to the highest point on the flux surface cross section
!c> R_bottom = major radius to the lowest point on the flux surface
!c> height  = distance between the elevations of the highest point
!c>           and the lowest point
!c> R_dent_in  = major radius to the most indented inboard part of the surface
!c>           in cases where the flux surface is bean shaped
!c>           with a dent on the inboard edge
!c>           Note: R_dent_in = R_in whenever the flux surface is not indented
!c> R_dent_out = major radius to the most indented outboard part of the surface
!c>
!c> r_minor = half width = (R_out - R_in ) / 2     (usually called minor radius)
!c> r_major = ( R_out + R_in ) / 2                 (usually called major radius)
!c>
!c> elongation:	kappa = height / ( R_out - R_in )
!c>
!c> tiangularity:	delta = ( r_major - min ( R_top, R_bottom ) ) / r_minor
!c>
!c> indentation:	indent = max ( ( R_dent_in - R_in ) / rminor ,
!c> 			       ( R_dent_out - R_out ) / rminor )
!c>
!c> 	Note: for normal convex surfaces, indent = 0.0
!c> For bean-shaped surfaces with indentation on the inboard edge, indent > 0
!c> For bean-shaped surfaces with indentation on the outboard edge, indent < 0
!c>
!c> 	For surfaces with the point of the triangle facing outward,
!c> delta > 0.0, while if the triangle points inward, delta < 0.0
!c>
!c> 	For vertically elongated surfaces, kappa > 1.0
!c> For horizontally elongated surfaces, kappa < 1.0.
!c>
!c> squareness based on Holcomb Phys. Plasmas 16, 056116 (2009), see fig. 1
!c>
!c> Z_rmax = elevation where R==R_out
!c>
!c> Mikkelsen squareness attempts to match moment representation
!c>
!c>    R = Rmajor + Rminor*cos(theta + asin(triang)sin(theta))
!c>    Z = Rmajor + elong*Rminor*sin(theta + square*sin(2*theta))
!c>
!c> by finding the point Rmk = R(pi/4.), interpolating from the boundary
!c> data to Zmk and then using the moment representation at Z(pi/4.) to
!c> derive the squareness.
!c>
!c>
!c>
!
subroutine eq_surfgeo(ivec, rho, phi, ntheta, &
     elong, triang, indent, zmidp, rin_midp, rout_midp, limits, ier)

  ! mod DMC Apr 2011: method now available natively in xplasma
  !   see xplasma2/xplasma_rzgeo.f90

  use xplasma_obj_instance
  use eq_module

  implicit none

  integer, intent(in) :: ivec             ! number of surfaces
  real*8,  intent(in) :: rho(ivec)        ! radial coordinate of flux surface
  real*8,  intent(in) :: phi(ivec)        ! toroidal coordinate of flux surface
  integer, intent(in) :: ntheta           ! number of poloidal points to use around the flux surface
                                          ! or 0 for the default of 400

  real*8,  intent(out) :: elong(ivec)    ! surface elongation
  real*8,  intent(out) :: triang(ivec)   ! triangularity
  real*8,  intent(out) :: indent(ivec)   ! indentation

  real*8,  intent(out) :: zmidp (ivec)   ! height of flux surface centroid
  real*8,  intent(out) :: rin_midp(ivec)  ! inner midplane intercept at centroid height
  real*8,  intent(out) :: rout_midp(ivec)  ! outer midplane intercept at centroid height

  real*8,  intent(out) :: limits(20,ivec) ! surface limits used for elong,triang,indent
                                          ! 1  -> R_in
                                          ! 2  -> R_out
                                          ! 3  -> R_top
                                          ! 4  -> Z_top
                                          ! 5  -> R_bottom
                                          ! 6  -> Z_bottom
                                          ! 7  -> R_dent_in
                                          ! 8  -> R_dent_out
                                          ! 9  -> R_centroid
                                          ! 10 -> Z_centroid
                                          ! 11 -> lower Holcomb squareness
                                          ! 12 -> upper Holcomb squareness
                                          ! 13 -> Z_rmax
                                          ! 14 -> Rmk_low
                                          ! 15 -> Zmk_low
                                          ! 16 -> lower Mikkelsen squareness
                                          ! 17 -> Rmk_upp
                                          ! 18 -> Zmk_upp
                                          ! 19 -> upper Mikkelsen squareness

  integer, intent(out) :: ier             ! nonzero on error
    
  !--------------------------------------------------------------

  call xplasma_surfgeo(s, rho, &
       elong, triang, indent, &
       zmidp, rin_midp, rout_midp, &
       ier, &
       ntheta=ntheta, &
       phival=phi, &
       auxdata=limits)

  if(ier.ne.0) then
     call eq_errmsg('non-zero status code in old_xplasma:eq_surfgeo.f90')
  endif

end subroutine eq_surfgeo
