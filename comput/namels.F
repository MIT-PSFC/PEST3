#include "fpreproc/f77_dcomplx.h"
#include "fpreproc/fp_const.h"
      subroutine namels(iZ,iA,cname)
c names an element with charge and mass iZ and iA
c  (dmc 26 Apr 2001:  use Rob Andre's impurity labeler.  Note:  do not
c  REAL*8 convert this routine!)
 
      use periodic_table_mod
 
C============
C idecl:  explicitize implicit REAL declarations:
!
      IMPLICIT NONE
      REAL*8, PARAMETER :: ONE = 1.0D0 ! _fgtok_constant_
!
C============
      integer iZ  ! charge of the element
      integer iA  ! mass of the element
      real   r4A  ! atomic mass (real*4)  ! real*4 for periodic table software
      character cname*(*) ! character string name of the element
 
c default:  from periodic table...
      character*12 tmp
 
      integer iAmin,iZnum
      REAL   Astd  ! real*4 for periodic table software
 
c check: standard mass
 
      if(iZ .gt. 2) then
 
         iZnum=iZ
         do
            call standard_amu(iZnum,Astd)
            iAmin=Astd+ONE
            if(iAmin.ge.iA) exit
            iZnum=iZnum+1
         enddo
 
c impurity label
 
         r4A = iA
         tmp = to_periodic_table(iZnum,r4A,iZ,0)
         cname = tmp
 
      endif
 
c identify hydrogen isotopes:
 
      if(iZ .eq. 1)then
          if(iA .eq. 1) cname = 'H'
          if(iA .eq. 2) cname = 'D'
          if(iA .eq. 3) cname = 'T'
      endif
 
      if(iZ .eq. 2)then
          if(iA .eq. 3) cname = 'He3'
          if(iA .eq. 4) cname = 'He4'
      endif
 
      return
      end
