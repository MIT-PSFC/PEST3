subroutine pstGetQBounds(q0, qa, qmin, qmax)

  ! Return the safety factor q on axis (q0), at the edge (qa) 
  ! the minum values of q (qmin) and the maximum value of q
  ! (qmax).

  use i2mex_mod
  implicit none
  real*8, intent(out) :: q0, qa, qmin, qmax

  integer ns, ier
  real*8, dimension(:), allocatable :: psi, q

  call i2mex_getOriNs(ns, ier)
  call i2mex_error(ier)

  allocate(psi(ns), q(ns))
  call i2mex_getOriPsi(ns, psi, ier)
  call i2mex_error(ier)
  
  call i2mex_getQ(ns, psi, q, ier)
  call i2mex_error(ier)

  q0 = q(1)
  qa = q(ns)
  qmin = minval(q)
  qmax = maxval(q)
  deallocate(psi, q)

end subroutine pstGetQBounds

subroutine pstGetGSAverageError(res)
  
  ! Return average Grad-Shafranov Error normalized 
  ! to R J_phi on axis.

  use i2mex_mod
  implicit none
  real*8, intent(out) :: res

  integer nt1, ns, ier, i
  real*8, dimension(:), allocatable :: psi, the
  real*8, dimension(:,:), allocatable :: gserror

  call i2mex_getOriNt1(nt1, ier)
  call i2mex_error(ier)
  call i2mex_getOriNs(ns, ier)
  call i2mex_error(ier)
  
  allocate(the(nt1))
  allocate(psi(ns))
  allocate(gserror(nt1, ns))

  call i2mex_getOriT(nt1, the, ier)
  call i2mex_error(ier)

  call i2mex_getOriPsi(ns, psi, ier)
  call i2mex_error(ier)

  call i2mex_getGSError(nt1, ns, the, psi, gserror, ier)
  call i2mex_error(ier)
  res = SUM(SUM( gserror(1:nt1-1,2:ns), 1))/real((nt1-1)*ns,i2mex_r8)

  call i2mex_2Dx(nt1, ns, the, psi, gserror, 'gserror_pst.dx', ier)
  call i2mex_error(ier)

  deallocate(the, psi, gserror)

end subroutine pstGetGSAverageError
  

subroutine pstGetWlambda(p)

  ! Return the lowest eignevalue of delta-W, th epotential
  ! energy (determines ideal stability).

  use pstcom
  implicit none
  real*8, intent(out) :: p
  p = wlambda
end subroutine 

subroutine pstGetrMagnetic(p)

  ! Return the magnetic axis position [m].

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = rMagnetic
end subroutine 

subroutine pstGetb0SquareCentre(p)

  ! Return the magnetic field strength at the 
  ! geometric centre [T].

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = b0SquareCentre
end subroutine 

subroutine pstGettotalToroidalCurrent(p)

  ! Return the total plasma current [A/mu0].

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = totalToroidalCurrent
end subroutine pstGettotalToroidalCurrent

subroutine pstGetbetaPoloidal(p)

  ! Return the poloidal beta. 

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = betaPoloidal
end subroutine 

subroutine pstGetbetaToroidal(p)

  ! Return the toroidal beta.

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = betaToroidal
end subroutine 

subroutine pstGetTroyonG(p)

  ! Return the Troyon factor, also known 
  ! as beta-Normalized.

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = TroyonG
end subroutine 

subroutine pstGetfourPiInductance(p)

  ! Return the internal inductance li

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = fourPiInductance
end subroutine 

subroutine pstGetMajorRadius(p)

  ! Return the major radius (R at the geometric
  ! centre) [m].

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = majorRadius
end subroutine pstGetMajorRadius

subroutine pstGetMinorRadius(p)

  ! Return the minor radius [m].

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = minorRadius
end subroutine pstGetMinorRadius

subroutine pstGetElongation(p)

  ! Return the cross-section elongation.

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = elongation
end subroutine pstGetElongation


subroutine pstGetTriangularity(p)

  ! Return the cross-section triangularity

  use pstcom
  implicit none
  real*8, intent(out) :: p

  p = triangularity
end subroutine pstGetTriangularity


subroutine pstGetNosurf(j)

  ! Return the number of coarse grain radial surfaces.

  use pstcom
  implicit none
  integer, intent(out) :: j

  j = nosurf
end subroutine pstGetNosurf

subroutine pstGetMth(j)

  ! Return the number of poloidal sections.

  use pstcom
  implicit none
  integer, intent(out) :: j

  j = mth
end subroutine pstGetMth

subroutine pstGetIneg(kneg, nodes, nodes_size)

  ! Return the number 'kneg' of zero-crossings occurring at the node 
  ! values 'nodes'. A value of kneg /=0 indicates ideal instability
  ! by virtue of Newcomb's criterion.

  use pstcom
  implicit none
  integer, intent(out) :: kneg, nodes(*)
  integer, intent(in) :: nodes_size ! size of nodes

  integer i
  kneg = nIdealInstabilities
  do i = 1, minval( (/ nIdealInstabilities, nodes_size, 100 /) )
     nodes(i) = idealInstabilityNode(i)
  enddo
end subroutine pstGetIneg

subroutine pstGetNosing(josing)

  ! Return the number of singular surfaces.

  use pstcom
  implicit none
  integer, intent(out) :: josing

  josing = nosing
end subroutine pstGetNosing

subroutine pstGetDelpr(kdim, ksize, re_delpr, im_delpr)

  ! Return the Delta' matching matrix in real/imag parts. 
  ! kdim: leading dimension of re_delpr and im_delpr
  ! ksize: actual size (in general the no of singular surfaces nosing).
  ! *NOTE* To get the matching data in [1/psi] multiply these data
  !        by CMATCH.
  
  use pstcom
  implicit none
  integer, intent(in) :: kdim, ksize
  real*8, intent(out) :: re_delpr(kdim,*), im_delpr(kdim,*)

  integer i, j
  do j=1, ksize
     do i = 1, ksize
        re_delpr(i, j) =  real(delpr(i, j))
        im_delpr(i, j) = aimag(delpr(i, j))
     enddo
  enddo
end subroutine pstGetDelpr

subroutine pstGetGampr(kdim, ksize, re_gampr, im_gampr)

  ! Return the Gama' matching matrix in real/imag parts. 
  ! kdim: leading dimension of re_gampr and im_gampr
  ! ksize: actual size (in general the no of singular surfaces nosing).
  ! *NOTE* 
  ! 1. To get the matching data in [1/psi] multiply these data
  !    by CMATCH.
  ! 2. These are the raw data after renormalization, which differ from
  !    the small to large solution coefficients by a term XSMNUS that 
  !    has been subtracted off to regularize the Gamma's.
  
  use pstcom
  implicit none
  integer, intent(in) :: kdim, ksize
  real*8, intent(out) :: re_gampr(kdim, *), im_gampr(kdim, *)

  integer i, j
  do j=1, ksize
     do i = 1, ksize
        re_gampr(i, j) =  real(gampr(i, j))
        im_gampr(i, j) = aimag(gampr(i, j))
     enddo
  enddo
end subroutine pstGetGampr

subroutine pstGetError_Delpr(kdim, ksize, err_delpr)

  ! Return Delta' error estimate from linear least square fit in the
  ! no of radial surfaces.

  use pstcom
  implicit none
  integer, intent(in) :: kdim, ksize
  real*8, intent(out) :: err_delpr(kdim, *)

  integer i, j
  do j=1, ksize
     do i = 1, ksize
        err_delpr(i, j) =  error_delpr(i, j)
     enddo
  enddo
end subroutine pstGetError_Delpr

subroutine pstGetError_Gampr(kdim, ksize, err_gampr)

  ! Return Gamma' error estimate from linear least square fit in the
  ! no of radial surfaces.

  use pstcom
  implicit none
  integer, intent(in) :: kdim, ksize
  real*8, intent(out) :: err_gampr(kdim, *)

  integer i, j
  do j=1, ksize
     do i = 1, ksize
        err_gampr(i, j) =  error_gampr(i, j)
     enddo
  enddo
end subroutine pstGetError_Gampr

subroutine pstGetRsnorm(ksize, pnorm)

  ! Return rs normalization factor. Multiply the PEST3 
  ! matching data by this factor to get an approximate 
  ! normalized value in the minor radius.

  use pstcom
  use comggf
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: pnorm(*)

  pnorm(1:ksize) = 0.0
  if(ksize <= nosing) pnorm(1:ksize) = rsnorm(1:ksize)
end subroutine pstGetRsnorm

subroutine pstGetXsmnus(ksize, psmnus)
  
  ! Return the subtraction coefficient of small solution of
  ! odd parity (Eq. 61 in PBD). This is the large solution 
  ! coefficient of order one that is used to renormalized the large 
  ! solution to give a well-behaved Frobenius expansion in the zero
  ! beta limit. 
  ! 

  use pstcom
  use comggf
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: psmnus(*)

  psmnus(1:ksize) = 0.0
  if(ksize <= nosing) psmnus(1:ksize) = xsmnus(1:ksize)
end subroutine pstGetXsmnus

subroutine pstGetCmatch(ksize, pmatch)

  ! Return 2 mu f normalization factor in Delta'. Use this factor
  ! to obtain the matching data in [1/psi].

  use pstcom
  use comggf
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: pmatch(*)

  pmatch(1:ksize) = 0.0
  if(ksize <= nosing) pmatch(1:ksize) = cmatch(1:ksize)
end subroutine pstGetCmatch
  
subroutine pstGetPsisin(ksize, ppsi0, ppsis, ppsia)

  ! Return the values of the poloidal flux on axis (ppsi0), at the 
  ! singular surfaces (ppsis), and at the edge (ppsia). 

  use pstcom
  use comggf
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: ppsi0, ppsis(*), ppsia

  ppsi0 = psisin(1)
  ppsis(1:ksize) = 0.0
  if(ksize <= nosing) ppsis(1:ksize) = psisin(2:ksize+1)
  ppsia = psisin(nosing+2)
end subroutine pstGetPsisin

subroutine pstGetXmu(ksize, pmu)

  ! Return mu = sqrt(-D_I) at the singular surfaces, where
  ! D_I are the ideal Mercier indices.
  
  use pstcom
  use comggf
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: pmu(*)

  pmu(1:ksize) = 0.0
  if(ksize <= nosing) pmu(1:ksize) = xmu(1:ksize)
end subroutine pstGetXmu
  
subroutine pstGetDrlay(ksize, p)

  ! Return the resistive Mercier indices at the singular
  ! surfaces.

  use pstcom
  use comggf
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: p(*)

  p(1:ksize) = 0.0
  if(ksize <= nosing) p(1:ksize)= drlay(1:ksize)
end subroutine pstGetDrlay

subroutine pstGetQslay(ksize, p)

  ! Return the safety factor values at the singular
  ! surfaces.

  use pstcom
  use comggf
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: p(*)

  p(1:ksize) = 0.0
  if(ksize <= nosing) p(1:ksize)= qslay(1:ksize)
end subroutine pstGetQslay

subroutine pstGetHlay(ksize, p)

  ! Return the Glasser-Greene-Johnson coefficents H at the singular
  ! surfaces.

  use pstcom
  use comggf
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: p(*)

  p(1:ksize) = 0.0
  if(ksize <= nosing) p(1:ksize)= hlay(1:ksize)
end subroutine pstGetHlay

subroutine pstGetEflay(ksize, p)

  ! Return the Glasser-Greene-Johnson coefficents F at the singular
  ! surfaces.

  use pstcom
  use comggf
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: p(*)

  p(1:ksize) = 0.0
  if(ksize <= nosing) p(1:ksize)= eflay(1:ksize)
end subroutine pstGetEflay

subroutine pstGetXiLogTerm(ksize, re_p, im_p)

  ! Return the zero-beta log term (re and imag parts)

  use pstcom
  use comggf  
  implicit none
  integer, intent(in) :: ksize
  real*8, intent(out) :: re_p(*), im_p(*)

  re_p(1:ksize) = 0.0 ; im_p(1:ksize) = 0.0
  if(ksize <= nosing) then
     re_p(1:ksize)= real(XiLogTerm(1:ksize))
     im_p(1:ksize)= aimag(XiLogTerm(1:ksize))
  endif
end subroutine pstGetXiLogTerm
