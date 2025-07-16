subroutine eq_getjac(ivec,zrho,zchi,zphi,irzmode,jtens,detj,ierr)

  use xplasma_obj_instance
  use eq_module

  !  get metric tensor and determinant, given mag. coordinate

  implicit NONE

!  input:
  integer ivec                      ! vector dimensioning
  real*8 zrho(abs(ivec)),zchi(abs(ivec)),zphi(abs(ivec)) ! (rho,chi,phi)
  integer irzmode                   ! =1: d[R,Z,Lphi]/d[rho,chi,phi]

!  output:
  real*8 jtens(3,3,abs(ivec))            ! d[x,y,z]/d[rho,chi,phi]
  real*8 detj(abs(ivec))                 ! determinant of jtens

  integer ierr                      ! error code, 0=OK

  ! irzmode=0:
  !  jtens(1:3,1,j) = (dx/drho,dx/dchi,dx/dphi) at j'th vector element
  !  jtens(1:3,2,j) = (dy/drho,dy/dchi,dy/dphi) at j'th vector element
  !  jtens(1:3,3,j) = (dz/drho,dz/dchi,dz/dphi) at j'th vector element

  ! irzmode=1:
  !  jtens(1:3,1,j) = (dR/drho,dR/dchi,dR/dphi) at j'th vector element
  !  jtens(1:3,2,j) = (dZ/drho,dZ/dchi,dZ/dphi) at j'th vector element
  !  jtens(1:3,3,j) = R*(dphi/drho,dphi/dchi,1) at j'th vector element
  !               for axisymmetric cases jtens(1:3,3,j) = (0,0,R)
  !
  ! for compatibility with existing f77 xplasma codes:
  ! if irzmode=2, jtens(1:3,3) = (Z,0,R) is returned.

  !---------------------------

  logical :: ccwflag,axisymm
  integer :: jvec,j
  real*8 :: zr(abs(ivec)),zz(abs(ivec)),zdetj(abs(ivec)),zcos,zsin
  real*8 :: drdrho(abs(ivec)),dzdrho(abs(ivec))
  real*8 :: drdchi(abs(ivec)),dzdchi(abs(ivec))
 
  !--------------------------------
 
  jvec=abs(ivec)
  ccwflag = (ivec.gt.0)

  call xplasma_global_info(s,ierr, axisymm=axisymm)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?error detected in eq_getjac:'
     call xplasma_error(s,ierr,lunerr)
  else if(.not.axisymm) then
     ierr=1
     write(lunerr,*) ' ?eq_getjac: only available for axisymmetric geometry!'
  endif
  if(ierr.eq.0) then

     if(irzmode.eq.2) then
        call xplasma_rzjac(s,zrho,zchi,ierr, &
             ccwflag=ccwflag, &
             z=zz,r=zr,rzdetj=zdetj, &
             drdrho=drdrho,dzdrho=dzdrho, &
             drdtheta=drdchi,dzdtheta=dzdchi)
     else
        call xplasma_rzjac(s,zrho,zchi,ierr, &
             ccwflag=ccwflag, &
             r=zr,rzdetj=zdetj, &
             drdrho=drdrho,dzdrho=dzdrho, &
             drdtheta=drdchi,dzdtheta=dzdchi)
     endif

     if(ierr.eq.0) then

        detj=zr*zdetj

        if(irzmode.eq.0) then
           do j=1,jvec
              zcos=cos(zphi(j))
              zsin=sin(zphi(j))

              jtens(1,1,j)=drdrho(j)*zcos
              jtens(2,1,j)=drdchi(j)*zcos
              jtens(3,1,j)=-zr(j)*zsin

              jtens(1,2,j)=drdrho(j)*zsin
              jtens(2,2,j)=drdchi(j)*zsin
              jtens(3,2,j)=zr(j)*zcos

              jtens(1,3,j)=dzdrho(j)
              jtens(2,3,j)=dzdchi(j)
              jtens(3,3,j)=0
           enddo
           
        else
           do j=1,jvec
              jtens(1,1,j)=drdrho(j)
              jtens(2,1,j)=drdchi(j)
              jtens(3,1,j)=0

              jtens(1,2,j)=dzdrho(j)
              jtens(2,2,j)=dzdchi(j)
              jtens(3,2,j)=0

              if(irzmode.eq.2) then
                 jtens(1,3,j)=zz(j)
              else
                 jtens(1,3,j)=0
              endif
              jtens(2,3,j)=0
              jtens(3,3,j)=zr(j)
           enddo

        endif

     else
        write(lunerr,*) ' ?error detected in eq_getjac xplasm_rzjac call:'
        call xplasma_error(s,ierr,lunerr)
     endif

  endif

  if(ierr.ne.0) then
     jtens = 0.0d0
     detj  = 0.0d0
  endif

end subroutine eq_getjac
