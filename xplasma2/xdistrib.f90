subroutine xdistrib(iorder, nx,x,fprime, ibcxmin,bcxmin,ibcxmax,bcxmax,  &
     &  ntarg,ztarg,xtarg, fnorm, tol, ier)
!
!  find a set of roots F(xtarg)=ztarg  (to rel. tolerance tol), where
!  F is constructed from fprime as follows:
!
!     f'(x) = interpolating function of order "iorder" (see below) of
!              {[x(1),fprime(1)],...,[x(nx),fprime(nx)]}
!             **assumed positive definite** fprime(1:nx) are checked
!              and satisfying the given boundary conditions
!
!     (nx,x(1:nx),fprime(1:nx), iorder & BCs given, x strict ascending).
!
!     Then,
!               integral[x(1) to x](f'(x')dx')
!     F(x) =  ---------------------------------  (ranging from 0.0 to 1.0)
!             inregral[x(1) to x(nx)](f'(x')dx')
!
!     is constructed from f', and, the vectorized solver "zridderx" is
!     used to find the roots xtarg(1:ntarg) satisfying
!
!         (F(xtarg(i)) = ztarg(i), i=1,ntarg)
!
!     to x accuracy tol*(x(nx)-x(1)) and F accuracy tol.
!
!     (ztarg(1:ntarg) are given, represent the normalized integral values,
!     and need not be ordered, though 0.0<=ztarg(i)<=1.0 is expected)
!
!     the normalizing integral[x(1) to x(nx)](f'(x')dx') is returned as
!     `fnorm'.
!
!---------------------------------------------
!
!     outputs: xtarg(1:ntarg), ier
!
!     no error messages are written here; call xdistrib_error(lun,ier)
!     to write a message on unit lun
!
!---------------------------------------------
!
!  the interpolating function -- constructed with pspline library routines
!
!    iorder=-1 -- step function:  btw x(j) & x(j+1) f = fprime(j)
!      (boundary condition arguments ignored)
!
!    iorder=0  -- piecewise linear interpolant
!      (boundary condition arguments ignored)
!
!    iorder=1  -- Akima Hermite interpolation
!      (boundary condition codes -1, 0, 1 are OK, but ibcxmin=ibcxmax
!       is expected!)
!
!    iorder=2 or 3 -- spline
!      (boundary condition codes -1, 0, 1, 2, ... 7 all OK)
!
!    ibcxmin -- boundary condition code for left endpoint x(1)
!     bcxmin -- associated parameter (if 1<=ibcxmin<=4)
!    ibcxmax -- boundary condition code for right endpoint x(nx)
!     bcxmax -- associated parameter (if 1<=ibcxmax<=4)
!
!    boundary condition codes are as described in "cspline" in the
!    pspline library.  see:  http://w3.pppl.gov/~pshare/help/pspline_hlp.html
!
!    The most common BC choices are:
!
!     ibcxmin=ibcxmax=-1 -- periodic boundary condition
!        f''(x(1))=f''(x(nx))  f'''(x(1))=f'''(x(nx))
!
!     ibcxmin(max)=0     -- default boundary condition
!        "not a knot" for Spline, divided difference for Hermite
!
!     ibcxmin(max)=1, bcxmin(max)= f'' at indicated boundary
!        specify the slope of f' (i.e. f'') in the argument.
!
!     except for periodic, the boundary conditions are independent
!     of eachother.
!
!---------------------------------------------
!
  use xdistrib_mod
!
  implicit NONE
!
!  ***arguments***
!
  integer, intent(in) :: iorder       ! interpolating function order
  integer, intent(in) :: nx           ! dimension of x & fprime arrays
  real*8, intent(in) :: x(nx)         ! x coordinate -- STRICT ASCENDING
  real*8, intent(in) :: fprime(nx)    ! fprime -- all values .GE. ZERO
!
!  BC data -- ignored for some values of iorder
!
  integer, intent(in) :: ibcxmin,ibcxmax  ! boundary condition codes
  real*8, intent(in) :: bcxmin,bcxmax ! boundary condition data
!
!  target data
!
  integer, intent(in) :: ntarg        ! no. of target points
  real*8, intent(in) :: ztarg(ntarg)  ! F targets all values BTW 0.0 and 1.0
!
!  roots found
!
  real*8, intent(out) :: xtarg(ntarg)
!
!  normalizing integral
!
  real*8, intent(out) :: fnorm
!
  real*8, intent(in) :: tol           ! rel. accuracy tolerance
!
  integer, intent(out) :: ier         ! completion code, 0=OK
!
!    error codes -- see xdistrib_error
!
!-----------------------------------------------------------------
!  local...
!
  integer i,j,ilinx,ipx,iwarn,ifail
  integer, parameter, dimension(3) :: ict = (/ 1,0,0 /)
  integer, parameter :: ivec=10
  real*8 zero,one
  real*8 fmin,fmax,xmin,xmax
!
  real*8 feqi_1int
!
  real*8 wk(max(ivec,nx))
!
  real*8 fspv(ivec)  ! for test evaluations of interpolating function
!
  logical lmask(ntarg)
  real*8 zero_vec(ntarg),one_vec(ntarg),dum_vec(ntarg)
!
!-----------------------------------------------------------------
 
  external xdistrib_sub            ! for zridderx
 
!-----------------------------------------------------------------
!
  ier = 0
  zero = ier
  one = 1
  xtarg = zero
!
!  ...simple error checks...
!
  if(ntarg.le.0) ier=98
  if(nx.le.1) ier=99
  if(ier.ne.0) return
!
  if((iorder.lt.-1).or.(iorder.gt.3)) then
     ier=100
     return
  endif
!
  do i=1,nx-1
     if(x(i).ge.x(i+1)) then
        ier=101
        return
     endif
     if((fprime(i).eq.zero).and.(fprime(i+1).eq.zero)) ier=106
  enddo
  if(ier.ne.0) return
!
  xmin=x(1)
  xmax=x(nx)
!
  fmin=fprime(1)
  fmax=fprime(1)
!
  do i=2,nx
     fmin=min(fmin,fprime(i))
     fmax=max(fmax,fprime(i))
  enddo
  if(fmin.lt.zero) ier=102
  if((iorder.eq.-1).and.(fmin.eq.zero)) ier=106
!
  if(iorder.eq.1) then
     if((ibcxmin.lt.-1).or.(ibcxmin.gt.1)) ier=103
     if((ibcxmax.lt.-1).or.(ibcxmax.gt.1)) ier=104
     if(ibcxmin.ne.ibcxmax) ier=109
  else if(iorder.gt.1) then
     if((ibcxmin.lt.-1).or.(ibcxmin.gt.7)) ier=103
     if((ibcxmax.lt.-1).or.(ibcxmax.gt.7)) ier=104
  endif
  if(ier.ne.0) return
!
  do i=1,ntarg
     if((ztarg(i).lt.zero).or.(ztarg(i).gt.one)) then
        ier=105
        return
     endif
  enddo
!
  if(ier.ne.0) return
!
!---------------------------------------------------------------
! set up interpolating function
!
  allocate(splcoef(ncdim(iorder),nx))
  if(iorder.eq.-1) then
     splcoef(1,1:nx-1)=fprime(1:nx-1)
     splcoef(1,nx)=splcoef(1,nx-1)
  else
     splcoef(1,1:nx)=fprime(1:nx)
  endif
!
! compute hermite/spline coefficients
!
  if(iorder.eq.3) then
     call r8cspline(x,nx,splcoef,ibcxmin,bcxmin,ibcxmax,bcxmax, &
          wk,nx,ilinx,ier)
     if(ier.ne.0) ier=111
!
  else if(iorder.eq.2) then
     call r8mkspline(x,nx,splcoef,ibcxmin,bcxmin,ibcxmax,bcxmax, &
          ilinx,ier)
     if(ier.ne.0) ier=111
!
  else if(iorder.eq.1) then
     ipx=0
     if(ibcxmin.eq.-1) ipx=1
     if(ibcxmin.eq.1) then
        ipx=2
        splcoef(2,1)=bcxmin
        splcoef(2,nx)=bcxmax
     endif
     call r8akherm1p(x,nx,splcoef,ilinx,ipx,ier)
     if(ier.ne.0) ier=111
  endif
!
  if(ier.ne.0) go to 1000
!
  nx_spl=nx
  iorder_spl=iorder
!
!  setup interpolation zone lookup data
!
  if(ibcxmin.eq.-1) then
     ipx=1
  else
     ipx=0
  endif
!
  allocate(xpkg(nx,4))
!
  call r8genxpkg(nx,x,xpkg,ipx,0,1,tol,-3,ier)
  if(ier.ne.0) then
     ier=112
     go to 1000
  endif
!
!  negative evaluations test
!  evaluate ivec pts btw each successive node pair and check if the
!  min. is negative.
!
  if(iorder.ge.1) then
     do i=1,nx-1
 
        wk(1:ivec)=x(i) + (x(i+1)-x(I))*(/ (j/(one*(ivec+1)), j=1,ivec) /)
 
        if(iorder.eq.3) then
           call r8spvec(ict,ivec,wk,ivec,fspv,nx,xpkg,splcoef,iwarn,ier)
        else if(iorder.eq.2) then
           call r8vecspline(ict,ivec,wk,ivec,fspv,nx,xpkg,splcoef,iwarn,ier)
        else if(iorder.eq.1) then
           call r8vecherm1(ict,ivec,wk,ivec,fspv,nx,xpkg,splcoef,iwarn,ier)
        endif
 
        if(ier.ne.0) then
           ier=120
        else
           if(minval(fspv).lt.zero) ier=121
        endif
     enddo
  endif
!---------------------------------------------------------------
! integrate fprime*dx through each zone
!
  allocate(ftabl(nx))
!
  ftabl(1)=zero
  do i=1,nx-1
     ftabl(i+1)=ftabl(i)+feqi_1int(iorder,x,splcoef,nx,i,x(i+1)-x(i))
  enddo
  fnorm=ftabl(nx)
!
!---------------------------------------------------------------
! find the roots...
!
  lmask=.false.                   ! mask array
  do i=1,ntarg
     if(ztarg(i).eq.zero) then
        lmask(i)=.true.           ! left end point
        xtarg(i)=x(1)
     else if(ztarg(i).eq.one) then
        lmask(i)=.true.           ! right end point
        xtarg(i)=x(nx)
     endif
     zero_vec(i)=x(1)
     one_vec(i)=x(nx)
  enddo
!
  call zridderx(ntarg,lmask,zero_vec,one_vec,tol*(x(nx)-x(1)),tol, &
       xdistrib_sub, xtarg, ifail, ntarg, ztarg, 1, dum_vec, 1)
!
  if(ifail.ne.0) ier=120
!
!---------------------------------------------------------------
!  exit...
 
1000 continue
 
  if(allocated(xpkg)) deallocate(xpkg)
  if(allocated(splcoef)) deallocate(splcoef)
  if(allocated(ftabl)) deallocate(ftabl)
!
  return
!
end subroutine xdistrib
 
 
subroutine xdistrib_sub(ntarg,lmask,xtarg,atarg, &
     & idum0,ztarg,idum1,zdum,idum2)
  !
  ! evaluation subroutine for zridderx
  !
  !    looking for zeroes of F(x(k)) =
  !                    int(x0 to x(k))[f'*dx']/int(x0 to xn)[f'*dx'] - z(k)
  !
 
  use xdistrib_mod
  implicit NONE
 
  ! ------------------------ arguments --------------------------------
 
  integer, intent(in) :: ntarg         ! input vector dimension
  logical, intent(in) :: lmask(ntarg)  ! masking vector
  real*8, intent(in) :: xtarg(ntarg)   ! input x values to evaluate
  real*8, intent(out) :: atarg(ntarg)  ! output function eval. results.
  real*8, intent(in) :: ztarg(ntarg)   ! data for function eval.
 
  integer, intent(in) :: idum0,idum1,idum2  ! placeholders
  real*8, intent(in) :: zdum(ntarg)    ! placeholder
  ! ------------------------ local stuff ------------------------------
 
  integer nvec,jvec,i,imode,iwarn
  real*8 xvec(ntarg)
  integer iv(ntarg)
  real*8 dxn(ntarg),zdum2(ntarg),zdum3(ntarg)
 
  real*8 feqi_1int                     ! integration subfunction
 
  ! ------------------------ executable code --------------------------
 
  ! --> find the active points
 
  nvec=0
  do i=1,ntarg
     if(.not.lmask(i)) then
        nvec=nvec+1
        xvec(nvec)=xtarg(i)
     endif
  enddo
 
  ! --> table lookup on these
 
  imode=1
 
  iv=0
  call r8xlookup(nvec,xvec,nx_spl,xpkg,imode,iv,dxn,zdum2,zdum3,iwarn)
 
  jvec=0
  do i=1,ntarg
     if(.not.lmask(i)) then
        jvec=jvec+1
        atarg(i)=ftabl(iv(jvec)) + &
             & feqi_1int(iorder_spl,xpkg,splcoef,nx_spl,iv(jvec),dxn(jvec))
        atarg(i)=atarg(i)/ftabl(nx_spl) - ztarg(i)
     endif
  enddo
 
  return
 
end subroutine xdistrib_sub
 
 
subroutine xdistrib_error(lun,ier)
!
!  report errors
!
  select case(ier)
 
     case(98)
        write(lun,*) '?xdistrib: at least one target element required.'
 
     case(99)
        write(lun,*) '?xdistrib: at least two x array elements required.'
 
     case(100)
        write(lun,*) '?xdistrib: "iorder" argument out of range.'
 
     case(101)
        write(lun,*) '?xdistrib:  x array not strict ascending.'
 
     case(102)
        write(lun,*) '?xdistrib:  fprime array has negative element.'
 
     case(103)
        write(lun,*) '?xdistrib:  invalid bdy cond: left end point.'
 
     case(104)
        write(lun,*) '?xdistrib:  invalid bdy cond: right end point.'
 
     case(105)
        write(lun,*) '?xdistrib:  target element not btw zero and one.'
 
     case(106)
        write(lun,*) '?xdistrib:  fprime array has consecutive zero elements.'
 
     case(109)
        write(lun,*) '?xdistrib:  when using Hermite, iorder=1,'
        write(lun,*) '            ibcxmin .eq. ibcxmax is required.'
 
     case(111)
        write(lun,*) '?xdistrib:  spline setup routine returned an error.'
 
     case(112)
        write(lun,*) '?xdistrib:  zone lookup setup routine returned an error.'
 
     case(120)
        write(lun,*) '?xdistrib:  error code during interpolation test.'
 
     case(121)
        write(lun,*) '?xdistrib:  interpolation test yields negative value.'
 
     case default
 
  end select
 
  return
 
end subroutine xdistrib_error
 
