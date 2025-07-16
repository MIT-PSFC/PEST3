!...................................................................
      SUBROUTINE pstvacds2
!...................................................................
!     fetches the matric vacmat from disk. allows us to use chance's
!     stand-alone vacuum calculation and all it's features without
!     having to re-write the pest code each time. march 1993
 USE pstcom

 USE l22com

 USE mtrik1
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER I
      INTEGER J


!aplet  17.5_r8 .94      common /ivac/ vacmti(nfm,nfm)
 LOGICAL	 	 	 :: ldelta
 LOGICAL	 	 	 :: lmapck
 LOGICAL	 	 	 :: lvacdsk
!!$ REAL*8, DIMENSION(nths) 	 :: xwal
!!$ REAL*8, DIMENSION(nths) 	 :: zwal
 REAL*8, DIMENSION(nfv0) 	 :: vacpst
 COMMON /c2f90ldelta/ ldelta 
 COMMON /c2f90lmapck/ lmapck 
 COMMON /c2f90lvacdsk/ lvacdsk 
      write(6,200)
 200  format('skip vacuum calculation, enter diag vacmat')
      jmax1 = lmax(1) - lmin(1) + 1
      do i = 1, jmax1
      do j = 1, jmax1
      vacmat(i,j) =  0.0_r8 
      vacmti(i,j) =  0.0_r8 
      end do
      end do
!
      do i = 1, jmax1
      write(6,*) 'i = ', i
      read(5,*) vacmat(i,i)
      end do
!
!     print the matrix
!     
      CALL pstmatwrt(vacmat,nfm,nfm,jmax1,jmax1,"vacmat         " )
      CALL pstmatwrt(vacmti,nfm,nfm,jmax1,jmax1,"vacmti         " )
!
      return
      end

