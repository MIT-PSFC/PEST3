!***********************************************************************
      SUBROUTINE pstwritem(is)
!******************************************************************
!
!      this is the inverse of fetchm
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
      INTEGER I
      INTEGER J
      INTEGER LGIVUP
      INTEGER NADRES
      INTEGER IKG
      INTEGER INRD
      INTEGER JTRANG


!aplet  18.5_r8 .94
!aplet  18.5_r8 .94
!......................................................................
!
 REAL*8, DIMENSION(nkc,nkc) 	 :: zrmat
 REAL*8, DIMENSION(nkc,nkc) 	 :: zimat
      ncolmn = 2*ipol * jsub(m)
!
      length = (ncolmn*(ncolmn+1))/2
!
!      set address for writes.
!
      do i = 1, ncolmn
      do j = 1, ncolmn
      zrmat(i,j) =  real( amat(i,j) )
      zimat(i,j) = aimag( amat(i,j) )
      end do
      end do
!
!
      if ( checkd ) &
  CALL pstmatwrt(zrmat,nkc,nkc,ncolmn,ncolmn,"re amat with vacuum " )
!
      if ( checkd ) &
  CALL pstmatwrt(zimat,nkc,nkc,ncolmn,ncolmn,"im amat with vacuum " )
!
      lgivup=1
      nadres = (m-1) * length + 1
!
!      now pick up the potential/kinetic energy block
!
      ikg = 0
      do 10 i = 1, ncolmn
      do 10 j = i, ncolmn
      ikg = ikg + 1
   10 spot(ikg) = amat(i,j)
!
!      write to disk
!
!      first check on ioc. inrd=inpot for pot, and = inkin for kin
!
      inrd = inpot
      if(is  /=  1) inrd = inkin
!***       CALL pstzwr(inrd,spot(1),length,nadres,lgivup, 7000)
if ( inrd == inpot ) then
 wpot(nadres:nadres+length-1) = spot(1:length)
else
 akin(nadres:nadres+length-1) = spot(1:length)
end if
!
!      now we must fix the overlapping triangle of the m-2'th block.
!      set the index parameters. use spot to hold both the pot
!      and kin triangles, check indexing
!
      jtrang = ipol * jsub(m-1)
!
!     now load into array for write
!
      ikg = 0
      do 30 i = 1, jtrang
      do 30 j = i,jtrang
      ikg = ikg + 1
      spot(ikg) = amat(i,j)
   30 continue
!
!      now for the address. use ikg the number of words to be written
!
      nadres = (m-1) * length - ikg + 1
!
!***      CALL pstzwr(inrd,spot(1),ikg,nadres,lgivup, 7000)
if ( inrd == inpot ) then
 wpot(nadres:nadres+ikg-1) = spot(1:ikg)
else
 akin(nadres:nadres+ikg-1) = spot(1:ikg)
end if
!
      return
!      error jump
 7000 CALL psterrmes(outpst,'writem',lgivup)
      end

