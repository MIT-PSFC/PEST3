!
!...................................................................
      SUBROUTINE pstvacdsk
!...................................................................
!     fetches the matric vacmat from disk. allows us to use chance's
!     stand-alone vacuum calculation and all it's features without
!     having to re-write the pest code each time. march 1993
 USE pstcom

 USE l22com

 USE mtrik1
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER IOVAC
      INTEGER NDSK
      INTEGER LEN
      INTEGER IIFF
      INTEGER LGIVUP
      INTEGER NADRES
      INTEGER LMN
      INTEGER LMX
      INTEGER LENGTH
      INTEGER J12
      INTEGER J2
      INTEGER J1
      INTEGER I
      real*8 OUTPEST


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
      write(12,200)
      write(outmod,200)
 200  format('skip vacuum calculation here, using vacmat from fort.36 ')
      iovac = 36
      ndsk = 1
      CALL pstzop(iovac,"vacmat",len,ndsk,iiff,999)
      lgivup = 1
      nadres = 1
      CALL pstzrd(iovac,lmn,1,nadres,lgivup,999)
      nadres = nadres + 1
      CALL pstzrd(iovac,lmx,1,nadres,lgivup,999)
      if(lmin(1)  /=  lmn  .OR.  lmax(1)  /=  lmx) then
      write(outmod,100)lmn,lmx,lmin(1),lmax(1)
      write(6,100)lmn,lmx,lmin(1),lmax(1)
!
!
  100 format(' error***** mismatch in the l-values', &
 /,'expecting',3x,    2i4, ' vacmat has', 2i4)
      stop 'ERROR in vacuum'
      endif
      nadres = nadres + 1
      CALL pstzrd(iovac,xwal(1),mth2,nadres,lgivup,999)
      nadres = nadres + mth2
      CALL pstzrd(iovac,zwal(1),mth2,nadres,lgivup,999)
      nadres = nadres + mth2
      length = jmax1**2
      CALL pstzrd(iovac,vacpst(1),length,nadres,lgivup,999)
!
      j12 = 1
!
      do 700 j2 = 1, jmax1
      do 700 j1 = 1, jmax1
      vacmat(j1,j2) = vacpst(j12)
      j12 = j12 + 1
 700  continue
      if( .not. lsymz)then
         nadres = nadres + length
         CALL pstzrd(iovac,vacpst(1),length,nadres,lgivup,999)
!     
         j12 = 1
!     
         do 750 j2 = 1, jmax1
            do 750 j1 = 1, jmax1
               vacmti(j1,j2) = vacpst(j12)
               j12 = j12 + 1
 750     continue
      endif
      
      CALL pstzcl( iovac, 999 )
!     
!     print the matrix
!     
      CALL pstmatwrt(vacmat,nfm,nfm,jmax1,jmax1,"vacmat         " )
!
!     draw the wall
!
! aplet  17.5_r8 .94
!     if(xwal(1)  >    0.0_r8 )then
!        xmx =  0.5_r8  * ( xinf(1) + xinf(mth2/2))
!        CALL pstdrawc2(zwal,xwal,zinf,xinf,1,mth1,"z","x",xmx,xma,zma)
!     endif
      if( imode ==  1 ) then
      write(outmod, 1470)
      do i = 1, mth1
      write(outmod, 1480) i, xinf(i), zinf(i), xwal(i), zwal(i)
      end do
      end if
!
!
 1470 format(//' plasma and wall boundaries:' &
 / 1x,' i ',1x,5x,'xinf',  11x,'zinf',11x,'xwal',11x,'zwal')
!
 1480 format(1x,i3,4(1x,e14.7))
!
!aplet  17.5_r8 .94
!
      return
  999 continue
      CALL psterrmes(outpest,'vacdsk')
      stop 'ERROR in vacuum'
      end

