      subroutine ck_sepseg(xsegs,nsegs,norder,iorder)
c
      implicit NONE
      integer, intent(in) :: nsegs      ! no. of segments
      real, intent(in) :: xsegs(2,nsegs) ! segments; xsegs(1,i).le.xsegs(2,i)
      integer, intent(out) :: norder    ! ordering flag
                                        ! +1 -- ordering OK
                                        !  0 -- segments not disjoint
                                        ! -1 -- invalid segment
      integer, intent(out) :: iorder(nsegs) ! ordering
c
c  returning ordering of disjoint segments if possible
c  a segment is a pair of numbers x1,x2 satisfying x1.le.x2
c  (x1a,x2a) is disjoint from (x1b,x2b) if x1a.gt.x2b or x2a.lt.x1b
c
c  if all segments are disjoint, their ordering is returned in iorder(...)
c  so that xsegs(1:2,iorder(1)) is the segment with lowest values,
c          xsegs(1:2,iorder(2)) the next lowest, etc.
c
c  caution: N^2 algorithm is used, number of segments is assumed to be small!
c
      integer ii,jj,inside,iord(nsegs)
c--------------------------
c
      norder=-1
      do ii=1,nsegs
         if(xsegs(1,ii).gt.xsegs(2,ii)) return  ! invalid segment
      enddo
c
c  simple answer if only one segment
c
      if(nsegs.eq.1) then
         norder=1
         iorder(norder)=1
         return
      endif
c
c  scan multiple segments -- n**2 algorithm
c
      norder=0
      iord=1
      inside=0
c
      do jj=1,nsegs
         do ii=1,nsegs
            if(ii.eq.jj) cycle
            if((xsegs(2,jj).lt.xsegs(1,ii)).or.
     >         (xsegs(1,jj).gt.xsegs(2,ii))) then
               if(xsegs(1,jj).gt.xsegs(2,ii)) then
                  iord(jj)=iord(jj)+1
               endif
            else
               inside=ii
               exit
            endif
         enddo
         if(inside.ne.0) exit
      enddo
c
      if(inside.ne.0) then
         continue                       ! not disjoint
      else
         norder=nsegs
         iorder=iord
      endif
c
      return
      end
