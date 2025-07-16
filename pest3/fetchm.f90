!**********************************************************************
      SUBROUTINE pstfetchm(is)
!**********************************************************************
!      purpose: to read in the last big block from disk and load into 2
!              big matrices amat,bmat for manipulation.
!
 USE pstcom

 USE l22com

 USE l21com

 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IS
      INTEGER NCOLMN
      INTEGER LENGTH
      INTEGER LGIVUP
      INTEGER NADRES
      INTEGER INRD
      INTEGER IKG
      INTEGER I
      INTEGER J


!
!aplet  18.5_r8 .94
!aplet  18.5_r8 .94
!
!.......................................................................
!
!      set address
!
 REAL*8, DIMENSION(nkc,nkc) 	 :: zrmat
 REAL*8, DIMENSION(nkc,nkc) 	 :: zimat
      ncolmn = 2*ipol * jsub(m)
      length = (ncolmn*(ncolmn+1))/2
!
      lgivup = 1
      nadres = m * length + 1
!
!     now read into spot, holds both pot and kin, in turn
!
      inrd = inpot
      if(is /= 1)inrd = inkin
!***      CALL pstzrd(inrd,spot(1),length,nadres,lgivup, 7000)
 if( inrd == inpot ) then
 spot(1:length) = wpot(nadres:nadres+length-1)
 else
 spot(1:length) = akin(nadres:nadres+length-1)
 end if
!
! load upper triangle...
!
      ikg = 0
      do 10 i = 1, ncolmn
      do 10 j = i, ncolmn
      ikg = ikg + 1
 10   amat(i,j) = spot(ikg)
!
!       symmetrise amat
!
      do 20 i = 2, ncolmn
      do 20 j = 1, i-1
!***   20 amat(i,j) = amat(j,i)
 20   amat(i,j) = conjg( amat(j,i) )
!
      do i = 1, ncolmn
      do j = 1, ncolmn
      zrmat(i,j) =  real( amat(i,j) )
      zimat(i,j) = aimag( amat(i,j) )
      end do
      end do
!
      if ( checkd ) &
  CALL pstmatwrt(zrmat,nkc,nkc,ncolmn,ncolmn,"re amat with no vac " )
!
      if ( checkd ) &
  CALL pstmatwrt(zimat,nkc,nkc,ncolmn,ncolmn,"im amat with no vac " )
!
!
!     check for i/o completion
!
      return
!      error jump
 7000 CALL psterrmes(outpst,'fetchm',lgivup)
      end

