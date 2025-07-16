subroutine eqi_irhofun(id_axis,zlbl,inprof,zprof,iflag,zsm,id,ierr)

  use xplasma_definitions
  use eqi_rzbox_module

  !  make XPLASMA profile -- integrated quantity; smooth by 1/2 zone width
  !  to insure smooth derivative from spline; preserve integrated result.

  !  the input quantity is I(x) = integral[0 to x]{f(x)*(dV/dx)*dx}
  !  where V is volume (iflag=1), area (iflag=2), or a suitable user defined 
  !  quantity (iflag=<xplasma profile id>)

  !  the output quanity is a spline of a smoothed version of I(x) which
  !  preserves the integral at x=1.

  !  mod DMC Apr 2009 -- recoded (too many artifacts in the old method)
  !  synopsis:
  !    infer f(x), smooth this, recompute I(x), apply renormalization

  !  details: we actually treat a positive definite function
  !    g(x) = f(x) - min(f), min(f)=<minimum value of f(x) over x in [0:1]>
  !    g(x) is smoothed, its integral is preserved, I(min(f)) is added back in.

  implicit NONE
  integer, intent(in) :: id_axis    ! axis id -- must be rho or akin to rho
  character*(*), intent(in) :: zlbl ! name for profile function to be created
  integer, intent(in) :: inprof     ! size of the integrated data profile
  real*8, intent(in) :: zprof(inprof) ! the integrated data provided
                         ! if inprof = size(x_axis) zprof(1)=0 is expected
                         ! if inprof = size(x_axis)-1 the axial data point
                         !    is presumed to be omitted.
  integer, intent(in) :: iflag      ! =1 -- volume normalization;
                         ! derivative evaluations -> W/m^3, #/sec/m^3, etc.
                                    ! =2 -- area normalization;
                         ! derivative evaluations -> A/m^2 (current density).
                                    ! = <id> -- user provided normalization.

  real*8, intent(in) :: zsm         ! smoothing width
  !  a multiplier of max(1,zsm/(x(2)-x(1))) is applied to the smoothing
  !  width (= the local grid spacing); if zsm=0.0 use the standard 1/2 zone 
  !  width smooth, i.e. smoothing delta = local grid spacing.

  integer, intent(out) :: id        ! id of stored profile (if successful)
  integer, intent(out) :: ierr      ! completion code, 0=OK.

  !  if an error occurs and ierr is set, id=0 will be returned.

  !-------------------------
  !  local variables:
  integer isize                ! orig. gridsize, including mag. axis
  integer isizp                ! augmented gridsize = isize + 1

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  real*8, parameter :: ZERO = 0.0_R8
  real*8, parameter :: HALF = 0.5_R8
  real*8, parameter :: ONE = 1.0_R8
  real*8, parameter :: TWO = 2.0_R8
  real*8, parameter :: THREE = 3.0_R8
  real*8, parameter :: QUARTER = 0.25_R8
  real*8, parameter :: zltol = 1.0e-9_R8

  real*8 :: zsmul

  !  the x axis grid and corresponding volumes (or analogous V(x)).
  real*8, dimension(:), allocatable :: zxgrid0,zv0

  real*8 :: zfmin,zfmax,zdv,zvol,zgnorm

  !  the derived function on a shifted grid, and its smoothed version:
  real*8, dimension(:), allocatable :: zxgrid,zg,zgsm

  !  the original & reconstructed profiles I(x); latter is splined
  real*8, dimension(:), allocatable :: zi0,zi1
  
  !  delta & epsln profiles for smoothing...
  real*8, dimension(:), allocatable :: zdelta,zepsln
  real*8 :: zdum1,zx
  integer :: ix,i,j,idum2,idum3

  integer :: idum,ier0
  character*120 zmsg

  !-------------------------

#ifdef __DEBUG

  integer :: inbdbg
  common/nbdbg_com1/ inbdbg   ! zeroed in nbstart; set .gt.0 for plots

  integer :: inzons,inxx,inext,iprev
  real*8, dimension(:), allocatable :: zxx,zf1,zf2

#endif
  
  id = 0
  ierr = 0

  call xplasma_grid_size(sp,id_axis,isize,ierr)
  if((ierr.ne.0).or.(isize.eq.0)) then
     zmsg=' '
     write(zmsg,*) ' ?eqi_irhofun: invalid id_axis: ',id_axis
     call xplasma_errmsg_append(sp,trim(zmsg))
     ierr=501
     return
  endif

  if(inprof.eq.isize) then
     if(zprof(1).ne.ZERO) then
        zmsg=' '
        write(zmsg,*) ' ?eqi_irhofun: integrated profile expected, value'
        call xplasma_errmsg_append(sp,trim(zmsg))
        zmsg=' '
        write(zmsg,*) '  on axis should be ZERO, instead got: ',zprof(1)
        call xplasma_errmsg_append(sp,trim(zmsg))
        ierr=9999
        return
     endif
  else if(inprof.ne.(isize-1)) then
     zmsg=' '
     write(zmsg,*) ' ?eqi_irhofun: # of points in profile (inprof): ',inprof
     call xplasma_errmsg_append(sp,trim(zmsg))
     zmsg=' '
     write(zmsg,*) '  inconsistent with gridsize (axis id, size): ', &
          id_axis,isize
     call xplasma_errmsg_append(sp,trim(zmsg))
     ierr=9999
     return
  endif

  ! OK........

  allocate(zxgrid0(isize),zv0(isize),zi0(isize),zi1(isize))
  call xplasma_grid(sp,id_axis,zxgrid0,ier0)

  zi0(1)=ZERO
  if(inprof.eq.(isize-1)) then
     zi0(2:isize)=zprof
  else
     zi0 = zprof
  endif

  if(max(abs(zxgrid0(1)),abs(zxgrid0(isize)-ONE)).gt.zltol) then
     zmsg=' '
     write(zmsg,*) ' ?eqi_irhofun: base grid does not cover range [0,1]: ', &
          zxgrid0(1),zxgrid0(isize)
     call xplasma_errmsg_append(sp,trim(zmsg))
     ierr=9999
     deallocate(zv0,zi0,zi1,zxgrid0)
     return
  endif

  ! smoothing width factor: based on central zone spacing, zone spacing
  ! could vary radially

  ! minimum smoothing based on grid spacing; maximum is (1/4) the profile
  ! width, with the multiplier corresponding to this estimated at rho=0:

  zsmul = max(ONE,min( zsm/(zxgrid0(2)-zxgrid0(1)), &
       QUARTER*(zxgrid0(isize)-zxgrid0(1))/(zxgrid0(2)-zxgrid0(1)) ))

  ! V(x)

  if(iflag.eq.1) then
     call xplasma_volume(sp,zxgrid0,zv0,ier0)
  else if(iflag.eq.2) then
     call xplasma_area(sp,zxgrid0,zv0,ier0)
  else
     call xplasma_eval_prof(sp,iflag,zxgrid0,zv0,ier0)
  endif

  zvol = zv0(isize)

  do ix=2,isize
     zdv = zv0(ix)-zv0(ix-1)
     if(zdv.le.ZERO) then
        zmsg=' '
        write(zmsg,*) &
             ' ?eqi_irhofun: volume weighting not monotonic increasing:'
        call xplasma_errmsg_append(sp,trim(zmsg))
        zmsg=' '
        write(zmsg,'(" V(",i4,")=",1pe11.4," V(",i4,")=",1pe11.4)') &
             ix-1,zv0(ix-1),ix,zv0(ix)
        call xplasma_errmsg_append(sp,trim(zmsg))
        ierr=9999
        deallocate(zv0,zi0,zi1,zxgrid0)
        return
     endif
  enddo

  ! form working grid: zone centers + extra point at bdys

  isizp=isize+1

  allocate(zxgrid(isizp),zg(isizp),zgsm(isizp))
  allocate(zepsln(isizp),zdelta(isizp))
  zxgrid(1)=ZERO
  zxgrid(isizp)=ONE
  do ix=2,isize
     zxgrid(ix)=HALF*(zxgrid0(ix-1)+zxgrid0(ix))
     zdelta(ix)=zsmul*(zxgrid0(ix)-zxgrid0(ix-1))
     zdv = zv0(ix)-zv0(ix-1)
     zg(ix) = (zi0(ix)-zi0(ix-1))/zdv
     if(ix.eq.2) then
        zfmin = zg(ix)
        zfmax = zg(ix)
     else
        zfmin = min(zfmin,zg(ix))
        zfmax = max(zfmax,zg(ix))
     endif
  enddo

  if(zfmin.eq.zfmax) then
     ! no radial variation
     zi1 = zi0

  else
     zg(1)=zg(2)
     zg(isizp)=zg(isize)

     zdelta(1)=zdelta(2)
     zdelta(isizp)=zdelta(isize)

     do ix=1,isizp
        zepsln(ix)=zfmax-zfmin
        zg(ix)=zg(ix)-zfmin
     enddo

     zgnorm = zi0(isize) - zfmin*zvol

     call r8filtr6(zxgrid,zg,zgsm,isizp,zepsln,isizp,zepsln,0,zdelta, &
       0,ZERO,0,ZERO,ONE,ONE,zdum1,idum2,idum3)

     ! reconstruct integrated profile

     zi1(1)=ZERO
     do ix=2,isize
        zdv = zv0(ix)-zv0(ix-1)
        zi1(ix)=zi1(ix-1) + zdv*zgsm(ix)
     enddo

     ! check normalization (zg,zgsm were constructed positive definite)

     zgnorm = zgnorm/zi1(isize)

     ! apply
     do ix=2,isize
        zx=zxgrid0(ix)
        zi1(ix)=zi1(ix)*zgnorm
     enddo

     ! add in constant part
     zvol=ZERO
     do ix=2,isize
        zdv=zv0(ix)-zv0(ix-1)
        zvol=zvol+zdv
        zi1(ix)=zi1(ix)+zfmin*zvol
     enddo
  endif

  !  now create the function as a spline; f'(x)=0 on axis, not a knot @bdy

  call xplasma_create_1dprof(sp,zlbl,id_axis,zi1,id,ierr, &
       ispline=2,ibca=1,zbca=ZERO)

#ifdef __DEBUG
  if((ierr.eq.0).and.(id.gt.0).and.(inbdbg.gt.0)) then
     inzons=isize-1
     allocate(zxx(10*inzons+1),zf1(10*inzons+1),zf2(10*inzons+1))

     zxx(1)=zxgrid0(1)
     zf1(1)=zi0(1)
     inext=1
     do i=1,inzons
        iprev=inext
        inext=iprev+10
        zxx(inext)=zxgrid0(1+i)
        zf1(inext)=zi0(1+i)
        do j=1,9
           zxx(iprev+j)=zxx(iprev)+j*(zxgrid0(1+i)-zxgrid0(i))/10
           zf1(iprev+j)=zi0(i)+j*(zi0(i+1)-zi0(i))/10
        enddo
     enddo

     call xplasma_eval_prof(sp,id,zxx,zf2,ierr)

     call r8_grafx2(zxx,zf1,zf2,10*inzons+1,'x',' ', &
          'eqi_irhofun debug plot',' ',zlbl)

     deallocate(zxx,zf1,zf2)
  endif
#endif

  deallocate(zxgrid0,zv0,zi0,zi1)
  deallocate(zxgrid,zg,zgsm,zepsln,zdelta)

end subroutine eqi_irhofun
