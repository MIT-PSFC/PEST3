subroutine pstlaeigen(lambda)
  !
  ! compute eigenvalue using LAPACK routine
  ! A. Pletzer March 17 2000
  use  pstcom
  use l34com
  implicit none
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

  real*8, intent(out) :: lambda

  complex*16, dimension(:,:), allocatable :: lapack_AB, lapack_Q, lapack_Z
  complex*16, dimension(:), allocatable ::  lapack_WORK
  real*8, dimension(:), allocatable :: lapack_W, lapack_RWORK
  real*8 lapack_VL, lapack_VU, lapack_ABSTOL
  integer, dimension(:), allocatable :: lapack_IWORK, lapack_IFAIL
  integer lapack_N, lapack_KD, lapack_LDAB, lapack_LDQ, &
       & lapack_IL, lapack_IU, lapack_M, &
       & lapack_LDZ, lapack_INFO
  character*1 lapack_JOBZ, lapack_RANGE, lapack_UPLO

  integer mf
  integer j1, j2, nadres

  integer jm, nlong, jc, jf1, jf2, jl1, jl2

  mf = jsub(1)


  lapack_JOBZ='N' ! 'N' eigenvalues only
  lapack_RANGE='I' ! 'I' eigenvalues in (VL, VU] are sought
  lapack_UPLO='U' ! 'U' upper triangle is stored
  lapack_N = mp*mf       ! matrix order
  lapack_KD= 2*mf-1      ! no of super diagonals
  lapack_LDAB = lapack_KD + 1
  lapack_LDQ=1 !lapack_N  !JOBZ = 'N', the array Q is not referenced.
  lapack_VL = -0.01_r8  ! lower eigenvalue bound
  lapack_VU = 0.01_r8  
  lapack_IL = 1 ! no of eigenvalues
  lapack_IU = 1 
  lapack_ABSTOL = 1.e-10_r8
  lapack_LDZ = 1
  allocate(lapack_AB(lapack_LDAB, lapack_N))
  !allocate(lapack_Q(lapack_LDQ,lapack_N))
  allocate(lapack_Q(1,1))
  allocate(lapack_W(lapack_N))
  allocate(lapack_Z(lapack_LDZ,lapack_IL))
  allocate(lapack_WORK(lapack_N))
  allocate(lapack_RWORK(7*lapack_N))
  allocate(lapack_IWORK(5*lapack_N))
  allocate(lapack_IFAIL(lapack_N))

! do jc=1,(nfe+1)*nfm*(2*nfm+1)
!   write(99,*) wpot(jc)
! enddo

 lapack_AB=(0.0_r8, 0.0_r8)
 nlong = mf*(2*mf+1) ! triangle length
 
 do jm = 1, mp-1
    jc = (jm-1)*nlong + 1
    do jf1 = 1, mf
       do jf2 = 1, 2*mf-jf1+1
          jl1 = 2*mf + 1 - jf2
          jl2 = (jm-1)*mf + jf1 + jf2 - 1
          lapack_AB( jl1, jl2 ) = wpot(jc)
          jc = jc + 1
       enddo
    enddo
 enddo
 ! last triangle is special
 jm = mp
 jc = (jm-1)*nlong + 1 - (mf*(mf+1))/2 ! adjust address back to addvac mod.
 do jf1 = 1, mf
    do jf2 = 1, mf-jf1+1
       jl1 = 2*mf + 1 - jf2
       jl2 = (jm-1)*mf + jf1 + jf2 - 1
       lapack_AB( jl1, jl2 ) = wpot(jc)
       jc = jc + 1
    enddo
 enddo

  call ZHBEVX( lapack_JOBZ, lapack_RANGE, lapack_UPLO, lapack_N, &
       & lapack_KD, lapack_AB, lapack_LDAB, lapack_Q, lapack_LDQ, &
       & lapack_VL, lapack_VU, lapack_IL, &
       &         lapack_IU, lapack_ABSTOL, lapack_M, lapack_W, &
       & lapack_Z, lapack_LDZ, lapack_WORK, lapack_RWORK, &
       lapack_IWORK, lapack_IFAIL, &
       lapack_INFO )

  if(lapack_INFO /=0) then
     print*,'***ZHBEVX Failure after attempting to pstsolve eigenvalue problem'
     print*,'INFO=', lapack_INFO         
     print*,'IFAIL=', lapack_IFAIL
     stop 'ERROR: failure in solving eigenproblem in LAEIGEN'
  endif

  lambda = lapack_W(1)

  deallocate(lapack_AB)
  deallocate(lapack_Q)
  deallocate(lapack_W)
  deallocate(lapack_Z)
  deallocate(lapack_WORK)
  deallocate(lapack_RWORK)
  deallocate(lapack_IWORK)
  deallocate(lapack_IFAIL)

end subroutine pstlaeigen
