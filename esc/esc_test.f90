program test1

  implicit none

  integer, parameter :: r8=selected_real_kind(12,100)

  integer, parameter :: nt1=17, ns=21
  character*32 label
  real(r8) s_in(ns) ! sqrt(norm. toroidal flux)
  real(r8) pin(ns) ! pressure in mu0 * Pascal
  real(r8) piter(ns) ! pressure in mu0 * Pascal
  real(r8) jparallel(ns) ! parallel current <J.B>/r0<B. grad phi> in [mu0 A/m^2]
  real(r8) qin(ns) ! safety factor
  real(r8) rbound(nt1), zbound(nt1) ! the boundary points [m^2]
  real(r8) rbtor ! the covariant toroidal magnetic field at the plasma edge [T m]
  real(r8) r0 ! the major radius
  integer ifail ! failure flag (0 = ok)
  integer ijob  ! ijob=0: fresh start
  ! ijob=1: save equilibrium state in file <label>.es
  ! ijob=2: load  equilibrium state from file <label>.es
  ! ijob=3: load from and save into file <label>.es
  integer mpol  ! number of poloidal modes (typically 3-4, should be <=8)
  real(r8) dummy  ! not presently used

  real(r8), parameter :: pi=3.14159265358979312_r8
  real(r8) t(nt1), titer
  real(r8) z0, a, kappa, delta, q0, qa, b0, p0
  real(r8), dimension(:), allocatable :: psibar, phibar, g, iota, q, p, s
  real(r8), dimension(:,:), allocatable :: Rcos, Rsin, Zcos, Zsin
  real(r8), dimension(:,:), allocatable :: R, Z
  integer i, j, np1, na1, niter, nf
  integer, parameter :: clockwise = 1, counterclockwise = 0

  ! Profile and geometry specification

  r0=0.7_r8; z0=0.0_r8; a=0.6_r8; kappa=2.0_r8; delta=0.7_r8

  t = 2.0_r8*pi*(/ (real(i-1,r8)/real(nt1-1,r8), i=1, nt1) /)
  rbound = r0 +       a*cos( t + delta*sin(t) )
  zbound = z0 + kappa*a*sin( t                )

  q0 = 1.0_r8; qa = 10.0_r8; b0 = 1.0_r8; p0=0.04_r8;

  rbtor = r0*b0
  s_in  = (/ (real(i-1,r8)/real(ns-1,r8), i=1, ns) /)
  pin= p0*(1.0_r8 - s_in**2)
  qin = q0*(1.0_r8 + (qa/q0 - 1.0_r8)*s_in**4)
  jparallel = 2.0*b0*(1.0_r8 - s_in**2)/(r0*q0)


  write(*,*)' '
  write(*,*)'============================================================='
  write(*,*)'Test 1: Input profiles are p and j//. First start with zero'
  write(*,*)'pressure, then increase to prescribed value.'
  write(*,*)'============================================================='
  write(*,*)' '

  ! Initialization

  label='nstx.12345Z00' ! <MACHINE>.<RUNID>
  call esc0(label)

  niter = 3
  do i = 0, niter-1
     piter = real(i,r8)*pin/real(niter-1,r8)

     write(*,'(i4,a,f10.6)') i,' 2 p0/B0^2->',2._r8*piter(1)/b0**2 

     call esc_p_jparallel(piter, jparallel, rbound, zbound, nt1, &
     & rbtor, r0)

     ! Execute ESC

     ijob = 0 
     mpol = 4
     call esc(ifail, ijob, mpol, dummy)

     ! Test against failure

     if(ifail /= 0) then
        print *, 'ESC-failure ', ifail
        print *, 'possible causes are:'
        print *, 'ifail=2 problems near the boundary'
        print *, 'ifail=4 problems near the axis'
        print *, 'ifail=8 problems in the middle'
        print *, 'odd ifail numbers are usually tolerable'
     endif

  enddo

  ! Extract data. Many routines of the type escGetXXXX are available
  ! to extract the data; check file get.c for a complete list. All
  ! quantities are re-interpolated on the uniform sqrt(toroidal
  ! field) mesh.

  np1 = 129
  na1 = 11
  allocate(psibar(na1), phibar(na1), g(na1), iota(na1), q(na1), p(na1), s(na1))
  nf = 8 ! number of Fourier modes
  allocate(Rcos(nf, na1), Rsin(nf, na1))
  allocate(Zcos(nf, na1), Zsin(nf, na1))
  allocate(R(np1, na1), Z(np1, na1))
  call escGetPsibar(psibar, na1) ! poloidal flux/2*pi
  call escGetPhibar(phibar, na1) ! toroidal flux/2*pi
  call escGetG(g, na1) ! covariant toroidal B-function
  call escGetIota(iota, na1) ! 1/q 
  call escGetQ(q, na1) ! q
  call escGetP(p, na1) ! pressure
  call escGetFourierR(Rcos, Rsin, nf, na1) ! the grid
  call escGetFourierZ(Zcos, Zsin, nf, na1)  
  call escGetR(R, np1, na1, counterclockwise)  
  call escGetZ(Z, np1, na1, counterclockwise) 
  ! ... 

  ! Finally, clean-up

  call esc1
  
  write(*,*)' '
  write(*,*) '  Phibar    Psibar       g         q       iota'
  do i = 1, na1
     write(*,'(5f10.6)') phibar(i), psibar(i), g(i), q(i), iota(i)
  enddo

  write(*,*)''
  write(*,*) 'R-COS '
  write(*,'(8i8)') (i, i=0,nf-1)
  do i = 1, na1
     write(*,'(8f8.5)') (Rcos(j,i), j=1,nf)
  enddo

  write(*,*)''
  write(*,*) 'R-SIN '
  write(*,'(8i8)') (i, i=0,nf-1)
  do i = 1, na1
     write(*,'(8f8.5)') (Rsin(j,i), j=1,nf)
  enddo

  write(*,*)''
  write(*,*) 'Z-COS '
  write(*,'(8i8)') (i, i=0,nf-1)
  do i = 1, na1
     write(*,'(8f8.5)') (Zcos(j,i), j=1,nf)
  enddo

  write(*,*)''
  write(*,*) 'Z-SIN '
  write(*,'(8i8)') (i, i=0,nf-1)
  do i = 1, na1
     write(*,'(8f8.5)') (Zsin(j,i), j=1,nf)
  enddo

  s = sqrt(phibar/phibar(na1))
  call  dumpprof(na1, s, p, 0, 'p & j//: p(s)')
  call  dumpprof(na1, s, q, 1, 'p & j//: q(s)')
  call  dumpprof(na1, s, g, 1, 'p & j//: g(s)')
  call  dumpprof(na1, s, psibar, 1, 'p & j//: psibar(s)')
  call  dumpprof(na1, s, phibar, 1, 'p & j//: phibar(s)')

  call  dumpmesh(np1,na1, R, Z, 0, 'p & j//')

  write(*,*)' '
  write(*,*)'============================================================='
  write(*,*)'Test 2: Input profiles are p and q. This is achieved by first'
  write(*,*)'prescribing an initial j// guess before switching to the q '
  write(*,*)'mode. '
  write(*,*)'============================================================='
  write(*,*)' '

  ! cold start

  write(*,*)''
  write(*,*) 'Run in j// mode...'
  piter = 0.0_r8
  niter = 2
  call esc0(label)
  do i = 1, niter
     piter = real(i-1,r8)*pin/real(niter-1,r8)
     call esc_p_jparallel(piter, jparallel, rbound, zbound, nt1, &
          & rbtor, r0)
     call esc(ifail, 0, mpol, dummy)
     if(ifail /= 0) then
        print *, 'ESC-failure ', ifail
        print *, 'possible causes are:'
        print *, 'ifail=2 problem near the boundary'
        print *, 'ifail=4 problem near the axis'
        print *, 'ifail=8 problem in the middle'
        print *, 'odd ifail numbers are usually tolerable'
     endif
  enddo

  write(*,*)''
  write(*,*) 'Run in q mode...'
  do i = 1, niter
     piter = real(i-1,r8)*pin/real(niter-1,r8)
     call esc_p_q(piter, qin, rbound, zbound, nt1, &
     & rbtor, r0)
     call esc(ifail, 1, mpol, dummy)
     if(ifail /= 0) then
        print *, 'ESC-failure ', ifail
        print *, 'possible causes are:'
        print *, 'ifail=2 problem near the boundary'
        print *, 'ifail=4 problem near the axis'
        print *, 'ifail=8 problem in the middle'
        print *, 'odd ifail numbers are usually tolerable'
     endif
  enddo

  write(*,*)''
  write(*,*) 'Done! Access the data and and dump results in *.mtv files'

  call escGetPsibar(psibar, na1) ! poloidal flux/2*pi
  call escGetPhibar(phibar, na1) ! toroidal flux/2*pi
  call escGetG(g, na1) ! covariant toroidal B-function
  call escGetIota(iota, na1) ! 1/q 
  call escGetQ(q, na1) ! q
  call escGetP(p, na1) ! pressure
  call escGetR(R, np1, na1, counterclockwise)  
  call escGetZ(Z, np1, na1, counterclockwise) 
  ! ... 

  ! Finally, clean-up

  call esc1

  call  dumpmesh(np1,na1, R, Z, 1, 'p & q')
  call  dumpprof(na1, s, p, 1, 'p & q: p(s)')
  call  dumpprof(na1, s, q, 1, 'p & q: q(s)')
  call  dumpprof(na1, s, g, 1, 'p & q: g(s)')
  call  dumpprof(na1, s, psibar, 1, 'p & q: psibar(s)')
  call  dumpprof(na1, s, phibar, 1, 'p & q: phibar(s)')

  write(*,*)' '
  write(*,*)'============================================================='
  write(*,*)'Test 3: Find the beta limit, store the sequence of '
  write(*,*)'equilibria.'
  write(*,*)'============================================================='
  write(*,*)' '

  piter = 0.0_r8
  write(*,*) 'Run in j// mode...'
  call esc0(label)
  niter = 3
  do i=1, niter
     piter = real(i,r8)*pin/real(niter-1,r8)
     write(*,'(i4,a,f10.6)') i,' 2 p0/B0^2->',2._r8*piter(1)/b0**2
     call esc_p_jparallel(piter, jparallel, rbound, zbound, nt1, &
          & rbtor, r0)
     call esc(ifail, 0, mpol, dummy)
     if(ifail /= 0) then
        print *, 'ESC-failure ', ifail
        print *, 'possible causes are:'
        print *, 'ifail=2 problems near the boundary'
        print *, 'ifail=4 problems near the axis'
        print *, 'ifail=8 problems in the middle'
        print *, 'odd ifail numbers are usually tolerable'
     endif
  enddo
  
  ifail = 0
  niter = 10
  i = 1
  write(*,*)''
  write(*,*) 'Run in q mode...'

  do while (ifail == 0 .and.  i <= 1000)
     piter = real(i-1,r8)*pin/real(niter-1,r8)
     write(*,'(i4,a,f10.6)') i,' 2 p0/B0^2->',2._r8*piter(1)/b0**2
     call esc_p_q(piter, qin, rbound, zbound, nt1, &
     & rbtor, r0)
     titer = real(i,r8)
     call esc(ifail, 1, mpol, titer)
     if(ifail /= 0) then
        print *, 'ESC-failure ', ifail
     endif
     i = i + 1
  enddo

!!$  write(*,*) 'Max beta limit reached!'
  
  call esc1 ! clean-up

  write(*,*)' '
  write(*,*)'============================================================='
  write(*,*)'Test 4: Load the last "good" equilibrium from *.es file and'
  write(*,*)'rerun ESC. '
  write(*,*)'============================================================='
  write(*,*)' '

  call esc0(label)

  piter = real(i-2,r8)*pin/real(niter-1,r8)
  call esc_p_q(piter, qin, rbound, zbound, nt1, &
     & rbtor, r0)

  call esc(ifail, 2, mpol, titer-1._r8)
  
  write(*,*)''
  write(*,*) 'Done! Access the data and and dump results in *.mtv files'

  call escGetPsibar(psibar, na1) ! poloidal flux/2*pi
  call escGetPhibar(phibar, na1) ! toroidal flux/2*pi
  call escGetG(g, na1) ! covariant toroidal B-function
  call escGetIota(iota, na1) ! 1/q 
  call escGetQ(q, na1) ! q
  call escGetP(p, na1) ! pressure
  call escGetR(R, np1, na1, counterclockwise)  
  call escGetZ(Z, np1, na1, counterclockwise) 

  call  dumpmesh(np1,na1, R, Z, 1, 'p & q beta-max')
  call  dumpprof(na1, s, p, 1, 'p & q beta-max: p(s)')
  call  dumpprof(na1, s, q, 1, 'p & q beta-max: q(s)')
  call  dumpprof(na1, s, g, 1, 'p & q beta-max: g(s)')
  call  dumpprof(na1, s, psibar, 1, 'p & q beta-max: psibar(s)')
  call  dumpprof(na1, s, phibar, 1, 'p & q beta-max: phibar(s)')

  call esc1 ! clean-up


  

  deallocate(psibar, phibar, g, iota, q, p, s)
  deallocate(Rcos, Rsin)
  deallocate(Zcos, Zsin)
  deallocate(R, Z)

  write(*,*) 'Successful end of ESC test program'
  write(*,*) 'Results can now be visualized if you have acess to '
  write(*,*) 'the plotting program "plotmtv" by typing'
  write(*,*) 'plotmtv prof.mtv'
  write(*,*) 'plotmtv mesh.mtv'

end program test1


!*******************************************************************************
! the following routines dump arrays into ascii files for 
! graphics postprocessing, in a format required by the 
! plotmtv package. If you have plotmtv installed on your
! computer, type "plotmtv <file>.mtv" to visualize the 
! results.

      SUBROUTINE dumpmesh(nt1, ns, x, z, nappend, title)
      IMPLICIT NONE
      integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)    
      REAL(r8) x(*), z(*)
      INTEGER nt1,ns,i,iu
      integer, intent(in) :: nappend ! set to 1 if in append mode
      CHARACTER*(*), intent(in) :: title
      iu = 90      
      if(nappend==0) then
         OPEN(iu,FILE='mesh.mtv',FORM='formatted')
      else
         OPEN(iu,FILE='mesh.mtv',FORM='formatted', POSITION='append')
      endif

      WRITE(iu,*)'$DATA=COLUMN'
      WRITE(iu,*)'%toplabel = "'//title//'"'
      WRITE(iu,*)'%xlabel = "R"'
      WRITE(iu,*)'%ylabel = "Z"'
      WRITE(iu,*)'%linelabel = ""'
      WRITE(iu,*)'%linetype = 1'
      WRITE(iu,*)'%linecolor = 1'
      WRITE(iu,*)'%markertype = 0'
      WRITE(iu,*)'%markercolor = -1'
      WRITE(iu,*)'x y'
      DO i=1,nt1*ns
         WRITE(iu,*)x(i), z(i)
      END DO
      CLOSE(iu)
      RETURN
      END SUBROUTINE dumpmesh
      
      SUBROUTINE dumpprof(ns, x, y, nappend, title)
      IMPLICIT NONE
      integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)    
      INTEGER, intent(in) :: ns
      REAL(r8), intent(in) :: x(*), y(*)
      integer, intent(in) :: nappend ! set to 1 if in append mode
      CHARACTER*(*), intent(in) :: title

      integer i, iu

      iu = 90

      if(nappend==0) then
         OPEN(iu,FILE='prof.mtv',FORM='formatted')
      else
         OPEN(iu,FILE='prof.mtv',FORM='formatted', POSITION='append')
      endif

      WRITE(iu,*)'$DATA=COLUMN'
      WRITE(iu,*)'%toplabel="',title,'"'
      WRITE(iu,*)'%xlabel = ""'
      WRITE(iu,*)'%ylabel = ""'
      WRITE(iu,*)'%linelabel = ""'
      WRITE(iu,*)'%linetype = 1'
      WRITE(iu,*)'%linecolor = 1'
      WRITE(iu,*)'%markertype = 12'
      WRITE(iu,*)'%markercolor = 0'
      WRITE(iu,*)'x y'
      DO i=1,ns
         WRITE(iu,*)x(i), y(i)
      END DO
      WRITE(iu,*)' '
      CLOSE(iu)
      RETURN
      END SUBROUTINE dumpprof
