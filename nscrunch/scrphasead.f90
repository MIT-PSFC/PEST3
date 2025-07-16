      SUBROUTINE scrphasead(rm,ym,immax,imom,alpha)
!
!  adjust moments set by subtracting a phase shift alpha from the
!  theta parametrization
!
      IMPLICIT NONE
       REAL*8 	yjn2 ! <phasead.f>
       REAL*8 	yjn1 ! <phasead.f>
       REAL*8 	rjn2 ! <phasead.f>
       REAL*8 	rjn1 ! <phasead.f>
       INTEGER 	im ! <phasead.f>
       INTEGER 	imom ! <phasead.f>
       REAL*8 	alpha ! <phasead.f>
       INTEGER 	immax ! <phasead.f>
      REAL*8 rm(0:immax,2),ym(0:immax,2)
!
      REAL*8 zsn(100),zcs(100)
!
!------------------------------------
!
      CALL scrsincos(alpha,imom,zsn,zcs)
!
!  0'th moments are unaffected
!
      do im=1,imom
         rjn1= rm(im,1)*zcs(im) +rm(im,2)*zsn(im)
         rjn2=-rm(im,1)*zsn(im) +rm(im,2)*zcs(im)
         yjn1= ym(im,1)*zcs(im) +ym(im,2)*zsn(im)
         yjn2=-ym(im,1)*zsn(im) +ym(im,2)*zcs(im)
         rm(im,1)=rjn1
         rm(im,2)=rjn2
         ym(im,1)=yjn1
         ym(im,2)=yjn2
      enddo
      return
      end
