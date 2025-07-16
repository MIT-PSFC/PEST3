      subroutine mkixrid(arid, zrid, lrid)
C
C 11/04/96 CAL
C
C Make runid for TRRUN.ISAM file
C To take advantage of index key ordering,
C shot numbers must be aligned
C e.g.  12345A01
C      123456A100
C
#ifdef __VMS
      integer len
      character*(*) arid, zrid
      integer       lrid, lrid1, lshot
C
      len = 10
      call checkid(arid,lrid,lshot)
      zrid=' '
      if(lrid.eq.0) return
C
	if(lrid .eq. len) then	
           lrid1 = 1
	else if (lrid-lshot .eq. 4) then
	  lrid1=len-lrid+1	
	else if (lrid .eq. lshot) then
	  lrid1=len-lrid+1	
        else
	  lrid1=len-lrid	
	ENDIF	
	zrid(lrid1:)=arid(1:lrid)	
C
#endif
        return
        end
 
 
 
 
 
 
