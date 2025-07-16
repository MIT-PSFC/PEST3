!      6.1_r8    output
!.......................................................................
      SUBROUTINE pstoutput(nwz)
!
! corrected for up-down asym. case ap  17.10_r8 .95
! eigvag should be real? => no correction
!
!.......................................................................
 USE pstcom

 USE combla

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER NWZ
      INTEGER MAT2
      INTEGER I
      INTEGER NADRE2
      INTEGER IS


!................
!
!      now write to disk for post processing
!
 REAL*8, DIMENSION(10) 	 :: eigvag
 COMMON /c2f90eigvag/ eigvag 
      if(nwz  ==  0)go to 700
!
!***
      mat2 = mat*2
!***
!      write(outev) real(mat), real(nwz),( eigvag(i), i=1,nwz )
!
      do i = 1, nwz
      nadre2 = ( i-1 ) * mat2 + 1
      CALL pstzrd ( outvec,xt(1),mat2,nadre2,is, 7000 )
!      write(outev) (xt(j), j= 1, mat)
      end do
!
  700 continue
!
!!$      call timer(outpst,'end of 23')
      CALL pstempty (outpst)
!
!
 9110 format(1x," code bloc 23", 1x," cpu:i/o:sys: =",a10,2x, &
 " total microseconds=",i10/)
      return
 7000 CALL psterrmes(outpst,'output',is)
      end


