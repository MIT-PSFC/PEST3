         SUBROUTINE pstconalr(eps)

 USE pstcom

 USE combla

 USE comivi
 
 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      real*8 EPS
      INTEGER IS
      INTEGER NADRES
      INTEGER N1
      INTEGER JN
      INTEGER IKD
      INTEGER IB
      INTEGER IM
      INTEGER IELEM
      INTEGER JM
      INTEGER MLT


!
!
!     read first radial block
!
      is = 1
      nadres = 1
!!      CALL pstzrd(nds,a(1),nlong,nadres,is, 7000)
      if ( nds == inpot ) then
         a(1:nlong) = wpot(nadres:nadres+nlong-1)
      else if ( nds == inkin ) then
         a(1:nlong) = akin(nadres:nadres+nlong-1)
      else
         a(1:nlong) = ashift(nadres:nadres+nlong-1)
      end if
!
         if (mf  ==  0) go to 95
!
!     scan over ng-1 radial intervals
!
         n1=ng-1
         do  70  jn=1,n1
         if (jn  ==  1) go to 45
!
!     set first overlapping triangle
!
         ikd=0
         ib=0
         m=mf+ml+1
         do  40  im=1,mf
         m=m-1
         ielem=ikd
         do  30  jm=im,mf
         ielem=ielem+1
         ib=ib+1
         a(ielem)=c(ib)
  30     continue
         ikd=ikd+m
  40     continue
!
!     decompose radial block
!
  45     continue

         CALL pstcald(a,eps,mf,ml,nsing,neg, jn)
!
         if (nsing  ==  -1 ) return
!
!     write decomposed radial block
!
!!      CALL pstzwr(nds,a(1),nlong,nadres,is, 7000)
  if (nds == inpot ) then
     wpot(nadres:nadres+nlong-1) = a(1:nlong)
  else if( nds == inkin ) then
     akin(nadres:nadres+nlong-1) = a(1:nlong)
  else
     ashift(nadres:nadres+nlong-1) = a(1:nlong)
  end if
!
!     store last overlapping triangle into b
!
         mlt=nlong-(mf*(mf+1))/2
         ib=0
         do  60  im=1,mf
         do  50  jm=im,mf
         ib=ib+1
         mlt=mlt+1
         c(ib)=a(mlt)
  50     continue
  60     continue
!
!     read next radial block
!
      nadres = nadres + nlong
!!      CALL pstzrd(nds,a(1),nlong,nadres,is, 7000)
      if ( nds == inpot ) then
         a(1:nlong) = wpot(nadres:nadres+nlong-1)
      else if ( nds == inkin ) then
         a(1:nlong) = akin(nadres:nadres+nlong-1)
      else
         a(1:nlong) = ashift(nadres:nadres+nlong-1)
      end if
!
  70     continue
!
!     set last overlapping triangle
!
         m=mf+ml+1
         ikd=0
         ib=0
         do  90  im=1,mf
         m=m-1
         ielem=ikd
         do  80  jm=im,mf
         ielem=ielem+1
         ib=ib+1
         a(ielem)=c(ib)
  80     continue
         ikd=ikd+m
  90     continue
!
!     decompose last radial block
!
  95     continue

         CALL pstcald(a,eps,0,mf+ml,nsing,neg, jn)
!
         if (nsing  ==  -1 ) return
!
!     write last decomposed radial block
!
!!      CALL pstzwr(nds,a(1),nlong,nadres,is, 7000)
 if ( nds == inpot ) then
    wpot(nadres:nadres+nlong-1) = a(1:nlong)
 else if ( nds == inkin ) then
    akin(nadres:nadres+nlong-1) = a(1:nlong)
 else
    ashift(nadres:nadres+nlong-1) = a(1:nlong)
 endif
    
 900  format (/,' potential energy',/,(10e12.4))
!
         return
!
!       trouble
 7000 CALL psterrmes(ndlt,'conalr',is)
         end


