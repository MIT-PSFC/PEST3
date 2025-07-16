          SUBROUTINE pstbshift(a)
!     lcm (bshift)
 USE pstcom
 USE comivi
 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IS
      INTEGER NC
      INTEGER N2
      INTEGER N22
      INTEGER N23
      INTEGER N21
      INTEGER NCN21
      INTEGER NN2
      INTEGER NL2
      INTEGER NCNN21
      INTEGER NCNN2
      INTEGER NC1NN2
      INTEGER I
      INTEGER NADRES
      INTEGER J
      INTEGER JJ


!***         dimension a(1)
!
 COMPLEX*16, DIMENSION(*) 	 :: a
      is = 1
         nc=1
         n2=nlong/2
         n22=nc*n2
         n23=(nc+1)*n2
         n21=n2+1
         ncn21=nc*n2+1
         nn2=nlong-n2
         nl2=nc*nlong
         ncnn21=nc*(nlong-n2)+1
         ncnn2=nc*(nlong-n2)
         nc1nn2=(nc+1)*(nlong-n2)
!
!
      do 3 i = 1, ng
      nadres = (i-1) * nlong + 1
!!      CALL pstzrd(nda,a(1),n2,nadres,is, 7000)
!!      CALL pstzrd(ndb,a(n21),n2,nadres,is, 7000)
 a(1:n2) = wpot(nadres:nadres+n2-1)
 a(n21:n21+n2-1) = akin(nadres:nadres+n2-1)
!
      do 1 j = 1, n2
      jj = j + n22
      a(j) = a(j) - al0 * a(jj)
    1 continue
!
!!    CALL pstzwr(nds,a(1),n2,nadres,is, 7000)
 if ( nds == inpot ) then
    wpot(nadres:nadres+n2-1) = a(1:n2)
 else if ( nds == inkin ) then
    akin(nadres:nadres+n2-1) = a(1:n2)
 else
    ashift(nadres:nadres+n2-1) = a(1:n2)
 end if
!
      nadres = nadres + n2
!
!!      CALL pstzrd(nda,a(1),nn2,nadres,is, 7000)
!!      CALL pstzrd(ndb,a(ncnn21),nn2,nadres,is, 7000)
 a(1:nn2) = wpot(nadres:nadres+nn2-1)
 a(ncnn21:ncnn21+nn2-1) = akin(nadres:nadres+nn2-1)
!
      do 2 j = 1, nn2
      jj = j + ncnn2
      a(j) = a(j) - al0 * a(jj)
    2 continue
!
!!    CALL pstzwr(nds,a(1),nn2,nadres,is, 7000)
 if ( nds == inpot ) then
    wpot(nadres:nadres+nn2 -1) = a(1:nn2)
 else if ( nds == inkin ) then
    akin(nadres:nadres+nn2 -1) = a(1:nn2)
 else
    ashift(nadres:nadres+nn2 -1) = a(1:nn2)
 endif
!
    3 continue
      return
!
!       trouble
 7000 CALL psterrmes(ndlt, 'bshift', is)
         end


