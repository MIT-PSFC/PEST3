!-------------------------------------------------------
!  this source file contains the following routines:
!
!    scrunch_rzmom:  input {R,Z} data, return modified {R,Z} data and
!                    corresponding Fourier moments.
!
!    scrunch_rzmom77:  input {R,Z} data, return modified {R,Z} data and
!                    corresponding Fourier moments; no. of moments 
!                    kept can be different than moments array dimension
!
!    scrunch_rz:     input {R,Z} data, return modified {R,Z} data;
!                    scrunch_rzmom is used; Fourier moments not saved.
!
!    scrunch_moms:   input Fourier moments, return modified Fourier 
!                    moments; scrunch_rzmom is used, {R,Z} data is
!                    generated temporarily but not saved.
!
!-------------------------------------------------------
subroutine scrunch_moms(nth, xrho,ns, rmc1d,rms1d,zmc1d,zms1d, imom,iopt, ier)
!
  implicit NONE
!
  integer, intent(in) :: nth    ! size of theta grid to use for {R,Z} expansion
                                ! inth=max(51,nth) actually used...
  integer, intent(in) :: ns     ! no. of surfaces
  real*8, intent(in) :: xrho(ns) !  = sqrt(phi/philim) (sqrt. tor. flux)
 
  integer, intent(in) :: imom   ! no. of moments not counting 0'th moment
 
  real*8, intent(inout) :: rmc1d(0:imom,ns) !    rmc1d  R cos moments
  real*8, intent(inout) :: rms1d(0:imom,ns) !    rms1d  R sin moments
  real*8, intent(inout) :: zmc1d(0:imom,ns) !    zmc1d  Z cos moments
  real*8, intent(inout) :: zms1d(0:imom,ns) !    zms1d  Z sin moments

  integer, intent(in) :: iopt   ! option (0=standard, 1=equal arc theta)
  !   iopt in {1,3,5,7} => equal arc; otherwise use standard scrunch theta.
 
  integer, intent(out) :: ier   ! completion code: 0 = OK
!
!------------------ 
  integer inth,is,ith,im,jopt
  real*8, dimension(:,:), allocatable :: r,z
  real*8, parameter :: twopi = 6.2831853071795862D+00
  real*8 ztheta
  real*8 snthtk(imom),csthtk(imom)   ! sin(im*theta),cos(im*theta)
!------------------
!  allocate {R,Z}
!
  inth=max(51,nth)
!
  allocate(r(inth,ns),z(inth,ns))
!
! use expansion formula to fill in the values

  jopt=iopt
  if(jopt.gt.1) jopt=jopt-2

  do is=1,ns
     do ith=1,inth
        ztheta=(ith-1)*twopi/(nth-1)
        call scrSINCOS(ztheta,imom,snthtk,csthtk)
        r(ith,is)=rmc1d(0,is)
        z(ith,is)=zmc1d(0,is)
        do im=1,imom
           r(ith,is)=r(ith,is)+rmc1d(im,is)*csthtk(im)+rms1d(im,is)*snthtk(im)
           z(ith,is)=z(ith,is)+zmc1d(im,is)*csthtk(im)+zms1d(im,is)*snthtk(im)
        enddo
     enddo
  enddo

  call scrunch_rzmom(r,z,inth, xrho,ns, rmc1d,rms1d,zmc1d,zms1d, imom,jopt,ier)

  deallocate(r,z)

end subroutine scrunch_moms
!-------------------------------------------------------
subroutine scrunch_rz(r,z,nth,xrho,ns,imom,iopt,ier)
!
!  transform r(1:nth,1:ns),z(1:nth,1:ns) to a representation
!  parametrized by the scruncher-optimized theta parameter
!
!  the 1st dimensions of {r,z} contain the present theta variation;
!  r(1,ns)=r(nth,ns) and z(1,ns)=z(nth,ns) is expected to hold
!
!  r(1:nth,1) and z(1:nth,1) contain the magnetic axis (replicated nth times)
!  which should not vary.
!
!  the plasma boundary is in r(1:nth,ns), z(1:nth,ns).
!
!  in the transformation, use truncated Fourier series, with the 0'th
!  moment + `imom' higher moments.
!
!  updown asymmetric surfaces are allowed.
!
!    --dmc 5 Feb 2001--
!
!
  implicit NONE
 
  integer, intent(in) :: nth    ! no. of theta pts
  integer, intent(in) :: ns     ! no. of surfaces
 
  real*8, intent(inout) :: r(nth,ns),z(nth,ns)
  real*8, intent(in) :: xrho(ns) !  = sqrt(phi/philim) (sqrt. tor. flux)
 
  integer, intent(in) :: imom   ! no. of moments not counting 0'th moment

  integer, intent(in) :: iopt   ! option: 0=standard, 1=equal arc: th=0:2pi
                                ! option: 2=standard, 3=equal arc: th=-pi:pi
                                !     use odd value of nth to avoid FFT re-eval
                                !     even nth forces FFT re-eval as final step
                                ! option: 5=equal arc: th=0:2pi; re-eval from
                                !         Fourier coefficients as final step
                                ! option: 7=equal arc: th=-pi:pi; re-eval from
                                !         Fourier coefficients as final step
 
  integer, intent(out) :: ier   ! completion code: 0 = OK
!                             ier=489 means axis pts don't all match
!                             ier=490 means periodicity test failed
!
!
!---------------------------
! local:
!
! Fourier moments description of surfaces
!
  real*8, dimension(:,:), allocatable :: rmc1d,rms1d,zmc1d,zms1d
!
!---------------------------
!
 
  allocate(rmc1d(0:imom,ns),rms1d(0:imom,ns))
  allocate(zmc1d(0:imom,ns),zms1d(0:imom,ns))
 
  call scrunch_rzmom(r,z,nth,xrho,ns, rmc1d,rms1d,zmc1d,zms1d, imom,iopt,ier)
 
  deallocate(rmc1d,rms1d,zmc1d,zms1d)
 
  return
end subroutine scrunch_rz
 
subroutine scrunch_rzmom77(r,z,nth,xrho,ns, rmc1dd,rms1dd,zmc1dd,zms1dd, &
     imomd,imom,iopt,ier)

  ! this routine is like scrunch_rzmom (below) except that the output
  ! moments array dimension and the actual no. of moments can be different.

  ! this routine has the same arguments as scrunch_rz (above) *plus*
  ! the moments sets:
  !
  !    rmc1dd  R cos moments
  !    rms1dd  R sin moments
  !    zmc1dd  Z cos moments
  !    zms1dd  Z sin moments
 
  implicit NONE
 
  integer, intent(in) :: nth    ! no. of theta pts
  integer, intent(in) :: ns     ! no. of surfaces
 
  real*8, intent(inout) :: r(nth,ns),z(nth,ns)
  real*8, intent(in) :: xrho(ns) !  = sqrt(phi/philim) (sqrt. tor. flux)
 
  integer, intent(in) :: imomd  ! no. of moments **array dimension**
  integer, intent(in) :: imom   ! no. of moments **used** not counting 0'th moment
 
  real*8, intent(out) :: rmc1dd(0:imomd,ns)
  real*8, intent(out) :: rms1dd(0:imomd,ns)
  real*8, intent(out) :: zmc1dd(0:imomd,ns)
  real*8, intent(out) :: zms1dd(0:imomd,ns)

  integer, intent(in) :: iopt   ! option: 0=standard, 1=equal arc: th=0:2pi
                                ! option: 2=standard, 3=equal arc: th=-pi:pi
                                !     use odd value of nth to avoid FFT re-eval
                                !     even nth forces FFT re-eval as final step
                                ! option: 5=equal arc: th=0:2pi; re-eval from
                                !         Fourier coefficients as final step
                                ! option: 7=equal arc: th=-pi:pi; re-eval from
                                !         Fourier coefficients as final step

  integer, intent(out) :: ier   ! completion code: 0 = OK
 
  !                             ier=489 means axis pts don't all match
  !                             ier=490 means periodicity test failed
 
  !---------------------------------------------
  real*8, dimension(:,:), allocatable :: rmc1d,rms1d,zmc1d,zms1d
  integer :: imomi
  !---------------------------------------------

  imomi=min(imom,imomd)

  allocate(rmc1d(0:imomi,ns),rms1d(0:imomi,ns))
  allocate(zmc1d(0:imomi,ns),zms1d(0:imomi,ns))

  rmc1dd=0; rms1dd=0; zmc1dd=0; zms1dd=0

  call scrunch_rzmom(r,z,nth,xrho,ns, rmc1d,rms1d,zmc1d,zms1d, imomi,iopt, ier)
  if(ier.eq.0) then
     rmc1dd(0:imomi,ns)=rmc1d(0:imomi,ns)
     rms1dd(0:imomi,ns)=rms1d(0:imomi,ns)
     zmc1dd(0:imomi,ns)=zmc1d(0:imomi,ns)
     zms1dd(0:imomi,ns)=zms1d(0:imomi,ns)
  endif
  deallocate(rmc1d,rms1d,zmc1d,zms1d)

end subroutine scrunch_rzmom77
 
subroutine scrunch_rzmom(r,z,nth,xrho,ns, rmc1d,rms1d,zmc1d,zms1d, imom, iopt,&
     ier)
 
  ! this routine has the same arguments as scrunch_rz (above) *plus*
  ! the moments sets:
  !
  !    rmc1d  R cos moments
  !    rms1d  R sin moments
  !    zmc1d  Z cos moments
  !    zms1d  Z sin moments
 
  implicit NONE
 
  integer, intent(in) :: nth    ! no. of theta pts
  integer, intent(in) :: ns     ! no. of surfaces
 
  ! numerical surfaces, modified on output:
  real*8, intent(inout) :: r(nth,ns),z(nth,ns)

  real*8, intent(in) :: xrho(ns) !  = sqrt(phi/philim) (sqrt. tor. flux)
 
  integer, intent(in) :: imom   ! no. of moments not counting 0'th moment
 
  real*8, intent(out) :: rmc1d(0:imom,ns)
  real*8, intent(out) :: rms1d(0:imom,ns)
  real*8, intent(out) :: zmc1d(0:imom,ns)
  real*8, intent(out) :: zms1d(0:imom,ns)

  integer, intent(in) :: iopt   ! option: 0=standard, 1=equal arc: th=0:2pi
                                ! option: 2=standard, 3=equal arc: th=-pi:pi
                                !     use odd value of nth to avoid FFT re-eval
                                !     even nth forces FFT re-eval as final step
                                ! option: 5=equal arc: th=0:2pi; re-eval from
                                !         Fourier coefficients as final step
                                ! option: 7=equal arc: th=-pi:pi; re-eval from
                                !         Fourier coefficients as final step
                                ! (all other values equivalent to iopt=0).

  integer, intent(out) :: ier   ! completion code: 0 = OK
 
  !                             ier=489 means axis pts don't all match
  !                             ier=490 means periodicity test failed
 
  !---------------------------------------------
  ! local:
  ! indices, etc.
 
  integer i,j,k,im
  real*8, parameter :: twopi = 6.2831853071795862D+00
  real*8, parameter :: onepi = 3.1415926535897931D+00
  real*8, parameter :: zero = 0.0D+00
  real*8 ztheta
  real*8 snthtk(imom),csthtk(imom)   ! sin(im*theta),cos(im*theta)
 
  logical :: ifft_eval
  logical :: ipi_shift
  integer :: itest,icen

  real*8, dimension(:), allocatable :: zwk

  !-------------------------------------------------
  !  executable code...
 
  ier=0
  do j=2,nth
     if((r(j,1).ne.r(1,1)).or.(z(j,1).ne.z(1,1))) then
        ier=489
     endif
  enddo
  if(ier.ne.0) return
 
  do i=2,ns
     if((r(1,i).ne.r(nth,i)).or.(z(1,i).ne.z(nth,i))) then
        ier=490
     endif
  enddo
  if(ier.ne.0) return
 
  if((iopt.eq.1).or.(iopt.eq.3).or.(iopt.eq.5).or.(iopt.eq.7)) then

     call eq_arc_scrunch(r,z,nth,ns, rmc1d,rms1d,zmc1d,zms1d, imom,ier)
     ! r & z defined over (0:2pi) now

     if(iopt.eq.1) then
        ifft_eval=.FALSE.
        ipi_shift=.FALSE.

     else if(iopt.eq.3) then
        ipi_shift=.TRUE.
        itest = nth/2
        itest = 2*itest
        if(nth.eq.itest) then
           ! even number
           ifft_eval=.TRUE.
        else
           ! odd number; FFT re-eval avoided
           ifft_eval=.FALSE.
        endif
     else if(iopt.eq.5) then
        ifft_eval=.TRUE.
        ipi_shift=.FALSE.
     else if(iopt.eq.7) then
        ifft_eval=.TRUE.
        ipi_shift=.TRUE.
     endif
        
  else
     call scriscrunchx(xrho,r,z,nth,ns, rmc1d,rms1d,zmc1d,zms1d, 2,ns, &
          imom,ier)

     ifft_eval=.TRUE.
     ipi_shift= (iopt.eq.2)

  endif

  if(ier.ne.0) then
     ! ZERO all outputs
     rmc1d = zero
     rms1d = zero
     zmc1d = zero
     zms1d = zero
     r = zero
     z = zero

  else
     rms1d(0,:)=zero    ! make sure these unused array elements are defined
     zms1d(0,:)=zero

     if(ifft_eval) then

        ! regenerate {R,Z} from scrunched moments

        ! start at 2nd surface (axis is not moved).

        do i=2,ns
           do j=1,nth-1
              if(ipi_shift) then
                 ztheta= -onepi + (j-1)*twopi/(nth-1)
              else
                 ztheta=(j-1)*twopi/(nth-1)
              endif
              call scrSINCOS(ztheta,imom,snthtk,csthtk)
              r(j,i)=rmc1d(0,i)
              z(j,i)=zmc1d(0,i)
              do im=1,imom
                 r(j,i)=r(j,i)+rmc1d(im,i)*csthtk(im)+rms1d(im,i)*snthtk(im)
                 z(j,i)=z(j,i)+zmc1d(im,i)*csthtk(im)+zms1d(im,i)*snthtk(im)
              enddo
           enddo
           r(nth,i)=r(1,i)
           z(nth,i)=z(1,i)
        enddo

     else if(ipi_shift) then
        !  th=0:2pi -> th=-pi:pi w/o FFT re-eval.  nth must be odd.

        allocate(zwk(nth))
        icen = nth/2 + 1
        do i=2,ns
           zwk(1)=r(icen,i)
           zwk(nth)=zwk(1)
           k=icen
           do j=2,nth-1
              k=k+1
              if(k.gt.nth) k=2
              zwk(j)=r(k,i)
           enddo
           r(1:nth,i)=zwk

           zwk(1)=z(icen,i)
           zwk(nth)=zwk(1)
           k=icen
           do j=2,nth-1
              k=k+1
              if(k.gt.nth) k=2
              zwk(j)=z(k,i)
           enddo
           z(1:nth,i)=zwk
        enddo
        deallocate(zwk)

     endif
  endif
 
end subroutine scrunch_rzmom
