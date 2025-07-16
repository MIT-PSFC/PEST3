subroutine pstInit(kmax1, kmeshMax, kinputFormat, kinputPath, kin_len, ptime, &
     & newQ, iakima, x2axis, x2edge, np1, nr, separatrix, ier)

  ! Initialize PEST3. This will read the equilibrium, store it as cubic splines,
  ! rescale it if necessary and allocate the memory for PEST3.

  USE pstcom
  use l22com
  use i2mex_mod
  IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
  integer, intent(inout) :: kmax1 ! Fourier modes ..-kmax1...+kmax1. This is an 
                                  ! input variable. However, if the equilibrium
                                  ! resolution is judged to be insufficient 
                                  ! then corrective action is taken to reduce 
                                  ! kmax1.
  integer, intent(in) :: kmeshMax ! max finite element mesh size
  integer, intent(in) :: kinputFormat ! equilibrium input format 
  character*(*), intent(in) :: kinputPath ! file or MDSPlus path
  integer, intent(in) :: kin_len ! length of above char string
  real(r8), intent(in) :: ptime ! MDSPlus equilibrium time in [sec]
  real(r8), intent(in) :: newQ ! new q value at edge, or if negative set q to -newQ on axis.
  integer, intent(in) :: iakima ! set to 1 if Jacobian should be Akima splined
  real(r8), intent(in) :: x2axis ! to axis extrapolation fraction
  real(r8), intent(in) :: x2edge ! to axis extrapolation fraction
  integer, intent(in) :: np1 ! No of poloidal rays for contouring
  integer, intent(in) :: nr ! No of radial surface for contouring
  real(r8), intent(in) :: separatrix ! How close to the separatrix
  
                               ! a zero value leaves the q and g profiles unchanged.
  integer, intent(out) :: ier ! error flag (0=ok)

  integer :: nfm2, nfm41, mpol, nt1

  ier = 0

  call pstsetlnk (1)
  !
  CALL pstdefglo
  !
  !      read equilibrium
  !
  inputFormat = kinputFormat ! 
  ! -1 -> chease INP1
  !  0 -> TRANSP UFILES
  ! +1 -> chease inp1.cdf
  ! +2 -> jsolver eqdsk.cdf 
  ! +3 -> G-eqdsk file from efit
  ! +4 -> G-eqdsk file from efit thru ESC-Q
  ! +5 -> Menard input format in psipqgRZ MKS format
  ! +6 -> Belova input format in freeqbe netCDF format
  ! +7 -> read ceo.cdf as produced by tesc for instance
  ! etc...
  
  i2mex_to_axis_fraction = x2axis
  i2mex_to_edge_fraction = x2edge
  i2mex_nt1 = np1
  i2mex_ns = nr
  i2mex_LAST_NORM_SURFACE_IS = separatrix
  i2mex_xzakima = iakima

  if(inputFormat==0 .and. len(kinputPath)>=4) then

     if(kinputPath(1:4)=='MDS+') then

        ! MDSPlus

        write(*,*)'Access MDSPlus data ', TRIM(kinputPath(1:kin_len))
        write(*,'(a,f10.4,a)') &
             & ' at time ', ptime,' [sec]'
        call i2mex_fromMDSPlus(TRIM(kinputPath(1:kin_len)), &
             & ptime, 65, 302, i2mex_counterclockwise, ier)  
        call i2mex_error(ier)
        ! refine equilibrium
        mpol = 4
        call i2mex_refineESCPAndQ(mpol, ier)
        call i2mex_error(ier)

     endif

  else 


     select case(inputFormat)

     case(:-1)
        call i2mex_readEquilibrium(inputFormat, kinputPath(1:kin_len),  &
             & i2mex_clockwise, ier)
        call i2mex_error(ier)
     case(0)

        write(*,*)'Access UFILE data ', TRIM(kinputPath(1:kin_len))
        write(*,'(a,f10.4,a)') &
             & ' at time ', ptime,' [sec]'
        call i2mex_fromMDSPlus(TRIM(kinputPath(1:kin_len)), &
             & ptime, 65, 302, i2mex_counterclockwise, ier)  
        call i2mex_error(ier)
        ! refine equilibrium
        mpol = 4
        call i2mex_refineESCPAndQ(mpol, ier)
        call i2mex_error(ier)

     case(1:7)
        call i2mex_readEquilibrium(inputFormat, kinputPath(1:kin_len),  &
             & i2mex_clockwise, ier)
        call i2mex_error(ier)
     case(8)
        call i2mex_readEquilibrium(inputFormat, kinputPath(1:kin_len),  &
             & i2mex_clockwise, ier)
        call i2mex_error(ier)
        ! refine equilibrium
 !       mpol = 4
 !       call i2mex_refineESCPAndQ(mpol, ier)
 !       call i2mex_error(ier)
     end select

  endif

  if(ier/=0) then
     print*,'***ERROR** occurred when reading input file ', &
          & TRIM(kinputPath(1:kin_len))
     return 
  endif

!!$     print*,'refine equilibrium'
!!$         ! refine equilibrium
!!$         mpol =4 
!!$         call i2mex_refineESCPAndQ(mpol, ier)
!!$         call i2mex_error(ier)

  if(newQ /= 0.0_r8) then
     if(newQ > 0.0_r8) then
        ! adjust q at edge
        call i2mex_scaleQ(i2mex_o%psi_edge, newQ, ier)
        call i2mex_error(ier)
     else 
        ! adjust q on axis
        call i2mex_scaleQ(0.0_r8, abs(newQ), ier)
        call i2mex_error(ier)
     endif
  endif

!!$  call pstCheckMercier(ier)
!!$  if(ier/=0) then
!!$     print*,'***ERROR** occurred: Plasma is Mercier unstable '
!!$     return 
!!$  endif
  

  call i2mex_getOriNt1(nt1, ier)
  if(2*kmax1+1 > (nt1-1)/2 - 1) then
     print *,'--WARNING-- Too many Fourier modes ', 2*kmax1+1 
     print *,'            given the available equilibrium resolution of ', nt1-1
     kmax1 = ((nt1-1)/2 - 2)/2
     print *,'            Number of Fourier modes now reduced to ', 2*kmax1+1
  endif

  CALL pstinglo(kmeshMax)
  !
  CALL pstdefmod

  !
  !      input data from cards
  !
!!$      CALL pstcardmo

  lmax1 = kmax1
  lmax(1:nfe) =  lmax1
  lmin(1:nfe) = -lmax1
  jmax2 = 2*lmax1 + 1
  ! set dims here

  nfn = jmax2
  nfm = nfn
  nsg3=ipolp*2*(nfm+2)
  nthsfm=4 * nthss * nfm
  nkd=4*nfm
  nke=2*nfm*nfe
  nfv0=nfm * nfm
  nmt=ipolp* nfm * ( nfe-2 ) + nfm
  nmt2=nmt*nmt
  nmtg = nfm*nfe
  nfmgbx = 2*(nfm+1)
  nfm21=2*nfm + 1
  nfm2=nfm21*nfm21
  nfm41=2*nfm21
  nfm8=nfm41*nfm41
  nkb=ipolp*nfm + 1
  nkc = ipolp * 2 * nfm + 2
  nad = (nkc)*(nkc+1)/2
  nkc2 = nkc*nkc

  CALL pstmemory(3) ! allocate memory for fourier modes


  return
end subroutine pstInit

