      integer function iscmp1(zname,ablis,n)
C
C  look for a match of ZNAME in list ABLIS
C  return the index in ablis of a matching name, or 0 if not found
C
      implicit NONE
      integer n
      character*(*) zname
      character*(*) ablis(n)
C
C  local
C
      integer lmax
      parameter (lmax=512)
C
      integer ila,ilz,i,iscmp0,idummy
C
      character*1 dbuff(lmax)
      integer str_length
C
C-----------------------------
C
      iscmp1 = 0
C
      ila=min(lmax,len(ablis(1)))
      ilz=str_length(zname)
C
      if((ilz.eq.0).or.(ilz.gt.ila)) return
C
C  length OK
C
      do i=1,ila
         dbuff(i)=' '
         if(i.le.ilz) then
            dbuff(i)=zname(i:i)
            call uupper(dbuff(i))
         endif
      enddo
C
      iscmp1 = iscmp0(dbuff,ablis,n,ila,idummy)
C
      return
      end
