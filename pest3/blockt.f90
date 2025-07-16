!------------------------------------------------------------
      SUBROUTINE pstblockt
!------------------------------------------------------------
!      this subroutine packs delta-w into disk file potdsk...
!
 USE pstcom

 USE l21com

 USE l34com

 USE r33com

 USE l22com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER NCOLMN
      INTEGER LENGTH
      INTEGER LENGT2
      INTEGER NOVLAP
      INTEGER JDIAG
      INTEGER NADRES
      INTEGER NEQJS
      INTEGER ID1
      INTEGER IS
      INTEGER ID3
      INTEGER LGIVUP
      INTEGER KG
      INTEGER I
      INTEGER J
      INTEGER K2
      INTEGER KKK
      INTEGER K1
      INTEGER II


!
!      set constants..
!
      ncolmn = 2 * ipol * jsub(m1)
      length = ( ncolmn*(ncolmn+1) ) / 2
      lengt2 = 2*length
      novlap = ncolmn / 2
      jdiag = length - ( novlap*(novlap+1) ) / 2
!
      if ( m1  /=  1 )  go to 40
!
      nadres = 1
      neqjs = 0
      if ( (imode >  1)  .OR.  liter )  go to 40
!
   10 continue
!!$      CALL pstzop(inpot,'potdsk ',neqjs,id1,is,7000)
!!$      if(.not.lmarg) CALL pstzop(inkin,'kindsk',neqjs,id3,lgivup,7000)
!
   40 continue
!
!      store into spot
!
       if ( m1  ==  1 )  go to 55
!
!      add last triangle of spot
!
      kg = jdiag
      do i = 1, novlap
         do j = i, novlap
            kg = kg + 1
            spot(kg) = wreg(i,j)
            skin(kg) = wkin(i,j)
         end do
      end do
!
!      put spot on disk..
!
      lgivup = 1
!***      nadres = ( m1-2 )*length + 1
!***      CALL pstzwr( inpot,spot(1),length,nadres,lgivup,7000)
!***      if(.not.lmarg) 
!***     . CALL pstzwr(inkin,skin(1),length,nadres,lgivup,7000)
      nadres = ( m1-2 )*length + 1
      wpot(nadres:nadres+length-1) = spot(1:length)
      if(.not.lmarg) &
  akin(nadres:nadres+length-1) = skin(1:length)
!
      if ( .not. check1 )  go to 55
      if ( m1  ==  mbeg )  write(outmod,9001)
 9001 format(2x,' second block at singular surface ' )
      write(outmod,9000) m1
 9000 format(1x,' re potential energy at m1 = ',i4 )
      k2 = 0
      kkk = 0
      do 9003 i = 1, ncolmn
      k1 = k2 + 1
      k2 = k1 + ncolmn - i
      kkk = kkk + 1
      write(outmod,9002) ( zero,ii=1,kkk-1 ),  &
                    ( REAL(spot(kg)), kg=k1,k2 )
 9003 continue
      write(outmod,9004) m1
 9004 format(1x,' im potential energy at m1 = ',i4 )
      k2 = 0
      kkk = 0
      do 9005 i = 1, ncolmn
      k1 = k2 + 1
      k2 = k1 + ncolmn - i
      kkk = kkk + 1
      write(outmod,9002) ( zero,ii=1,kkk-1 ),  &
                    (aimag(spot(kg)), kg=k1,k2 )
 9005 continue
 9002 format(20(1x,10e12.4,/))
!
   55 continue
!
      if ( m1  ==  mp )  go to 29
!
!      fill up first part of spot..
!
      kg = 0
      do  i = 1, novlap
         do  j = i, ncolmn
            kg = kg + 1
            spot(kg) = wreg(i,j)
            skin(kg) = wkin(i,j)
         end do
      end do
!
   29 continue
!
!       repeat last block...
!
!***      nadres = nadres + length
!***      CALL pstzwr( inpot,spot(1),length,nadres,lgivup,7000)
!***      if(.not.lmarg) 
!***     . CALL pstzwr(inkin,skin(1),length,nadres,lgivup,7000)
      nadres = nadres + length
      wpot(nadres:nadres+length-1) = spot(1:length)
      if(.not.lmarg) &
  akin(nadres:nadres+length-1) = skin(1:length)
!
      return
!
 7000 CALL psterrmes(outpst,'blockt',lgivup )
      end


