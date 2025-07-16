!
!----------------------------------------------------------------------
!      1.5   initialization of indexing.
!.......................................................................
!
      SUBROUTINE pstindint
!.......................................................................
 USE pstcom

 USE l22com

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)


!
!      choice of linear or cubic basis functions...
!
      ibas = 1
!aplet  20.12.95      if ( lcub ) ibas = 2
      ipol = 1
!      1.5.1   set up rank of polarization blocks.
!             this version adds 2 extra rows and columns to each block
!             to allow for the adddition of matrix elements from small
!             singular solutions.......rg april  80.
!
      do 1511 m1 = 1, mp
         jsub(m1)  = ibas * ( lmax(1) - lmin(1) + 1 )
 1511 continue
!
!      1.5.2   set pointers to finite element blocks.
!
!
      jtot(1) = 1
      jtot(2) = jtot(1) + jsub(1)
      do 1521 m1 = 3, mp
      jtot(m1) = jtot(m1-1) + ipol * jsub(m1-1)
 1521 continue
      jtot( mp+1 ) = jtot( mp ) + ipol * jsub( mp )
!
!      1.5.3   find rank of matrix
!
      mat = jtot( mp ) + ipol * jsub( mp ) - 1
      if ( .not. checki )   return
 
      write ( outmod, 9000 )   ( jtot(m1),m1=1,mp)
      write ( outmod, 9001 )   mat
!
 9000 format ( 2x, "check of polarization block pointers",/ &
,(1x,    "jtot =", 20i4 ))
 9001 format ( 1x, "total matrix rank =", i4 )
      return
      end


