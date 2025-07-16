      SUBROUTINE pstinglo(kmeshMax)
!......................................................................
!      this SUBROUTINE pstreads control logicals for the package
!      and opens necessary disk files..
!
!......................................................................
!
 USE pstcom
 use i2mex_mod
 use newmet
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
  integer, intent(in) :: kmeshMax ! max finite lement mesh size
      INTEGER NMAT
      INTEGER LNSF0
      INTEGER LNTHS0
      INTEGER LNFM
      INTEGER LNFN
      INTEGER ID1
      INTEGER IS
      INTEGER ID3

      integer cdfid, ier, iok,i
!
 CHARACTER*6                       :: mpout1, mapdsk

 real(r8), dimension(:), allocatable :: the, psi

  nths0 = i2mex_o%nt1 - 1
  mth = nths0
  mth1 = mth + 1
  mth2 = mth + 2
  nths = nths0 + 5
  dth = twopi / mth  

    
  nsf0  = i2mex_o%ns
  if( nsf0 < (kmeshMax-1)*2 + 3) then
     
     ! increase equilibrium mesh resolution
     nsf0 = (kmeshMax-1)*2 + 4
  else
     ! equilibrium mesh is sufficiently fine
     nsf0  = i2mex_o%ns
  endif

  
  allocate(psi(nsf0    ), stat=iok)

  if(nsf0 > i2mex_o%ns) then
    ! uniform mesh in psi
    psi = i2mex_o%psi_edge* (/ (real(i-1,r8)/real(nsf0-1,r8), i=1,nsf0) /)
  else
    ! take original mesh
    call i2mex_getOriPsi(i2mex_o%ns, psi, ier)
    call i2mex_error(ier)
  endif

  nsf   = nsf0
  nsf2 = 2*nsf
  nths = nths0 + 5
  nthsh = nths/2
  nthss= nths -3

  nfe = 1 + nsf/2
  nfe1=nfe+1
  nths2=2*nths
  nthsq=4*nthss*nthss+1  

  ! allocate memory for metric quantities
  CALL pstmemory(1)

  allocate(the(nths0+1), stat=iok)
  
  call i2mex_getOriT(nths0+1, the, ier)
  call i2mex_error(ier)

  ! back off near axis
  psi(1) = psi(1) + (psi(2)-psi(1))*0.01_r8

  psia = psi
  psibig = twopi*(psi - psi(nsf))
  nosurf = nsf0
  dth = twopi / mth
  r   = 1.0_r8
  psimin = minval(psibig)
  psilim = maxval(psibig)
  dth = twopi / mth
  r2     = r * r
  r4     = r2 * r2
  r6     = r4 * r2
  psidif = ( psilim-psimin ) / twopi


  call pstCompMetric(nths0+1, nsf, the, psi, ier)
  if(ier/=0) then
     print*,'***ERROR** occurred in pstInglo after pstCompMetric'
  endif

  deallocate(the, psi)
  xma = SUM(xa(1:nths0,1))/real(nths0,r8)

!!$      nsf0 = nsf
!!$      nths0= nths - 5

!!$      if(nsf0 < 2*maxval(mm) ) then
!!$         print *,'INGLO: nsf0 =',nsf0,' < 2*maxval(mm): must decrease mm = ', mm
!!$         stop 'FATAL ERROR in INGLO'
!!$      endif
      if(nths0 < 4*jmax2) then
         print *,'INGLO: NTHS0 < 4*jmax2: must decrease number of Fourier modes = ',jmax2
         stop 'FATAL ERROR in INGLO'
      endif
!!$      nfn = nths0/2 + 1
!!$      nfm = nfn
!
!!$      nsf2 = 2*nsf
!!$      nths = nths0 + 5
!!$      nfe = 1 + nsf/2
!!$      nthsh = nths/2
!!$      nthss= nths -3
!!$      nfe1=nfe+1
!!$      nths2=2*nths
!!$      nthsq=4*nthss*nthss+1
!!$!
!!$      lnsf0 = nsf0
!!$      lnths0 = nths0
!!$      lnfm = nfm
!!$      lnfn = nfn
!!$!
!!$!
!!$      CALL pstmemory(1)
!!$      CALL pstmemory(3)
!
     dpsi(1)= one / mdiv / m
!
!      assign disk files...
!
!      channels for standard i/o....
!
      call pstsetlnk ( 2 )
!
!
 1000 format(' **error** possible mismatch in parameter values',/,   'file size expected=',i6,'disk file size is=',i6)
      return
 7000 CALL psterrmes(outpst,'inglo',is)
      end


