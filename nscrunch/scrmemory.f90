      SUBROUTINE scrmemory
 
      USE mintrp
      USE scrunch_inc1
      USE scrunch_inc2
!...begin heap allocation
      IMPLICIT NONE

       INTEGER ialloc

       ialloc = 0
       ALLOCATE( rin(nu,ncmax), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating rin(nu,ncmax)')
       END IF
       ALLOCATE( zin(nu,ncmax), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating zin(nu,ncmax)')
       END IF
       ALLOCATE( xvec(n2), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating xvec(n2)')
       END IF
       ALLOCATE( gvec(n2), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating gvec(n2)')
       END IF
       ALLOCATE( xdot(n2), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating xdot(n2)')
       END IF
       ALLOCATE( xstore(n2), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating xstore(n2)')
       END IF
       ALLOCATE( m1(mnd), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating m1(mnd)')
       END IF
       ALLOCATE( n1(mnd), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating n1(mnd)')
       END IF
       ALLOCATE( rmomb(mpol,2), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating rmomb(mpol,2)')
       END IF
       ALLOCATE( zmomb(mpol,2), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating zmomb(mpol,2)')
       END IF
       ALLOCATE( mm(mpol), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating mm(mpol)')
       END IF
       ALLOCATE( dm1(mpol), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating dm1(mpol)')
       END IF
       ALLOCATE( faccon(mpol), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating faccon(mpol)')
       END IF
       ALLOCATE( xmpq(mpol,4), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating xmpq(mpol,4)')
       END IF
       ALLOCATE( result(nphi2,mpol,4), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating result(nphi2,mpol,4)')
       END IF
       ALLOCATE( rin3d(ntheta,nphi), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating rin3d(ntheta,nphi)')
       END IF
       ALLOCATE( zin3d(ntheta,nphi), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating zin3d(ntheta,nphi)')
       END IF
       ALLOCATE( angle(ntheta,nphi), STAT=ialloc )
       IF(ialloc .NE. 0) THEN
          CALL SCRABORTER('ERROR allocating angle(ntheta,nphi)')
       END IF
!...end heap allocation
       RETURN
       END
