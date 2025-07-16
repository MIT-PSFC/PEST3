#include "fpreproc/library_names.h"
#include "fpreproc/blas_names.h"
      SUBROUTINE scrrestartdes(irst,ier)
      USE scrunch_inc1
      IMPLICIT NONE
       INTEGER 	ier ! <restartdes.f>
       INTEGER 	l ! <restartdes.f>
       INTEGER 	irst ! <restartdes.f>
!**********************************************************************
!       This routine either stores an accepted value of the local solution
!       (irst=1) or reset the present solution to a previous value (irst=2)
!**********************************************************************
        go to (15,25)irst
 
 15     CALL Scopy(n2,xvec,1,xstore,1)
        return
 
 25     do 20 l = 1,n2
        xdot(l) = 0.D0
 20     xvec(l) = xstore(l)
 
        delt = .95D0* delt
        irst = 1
        nresets = nresets + 1
        if(nresets .ge. 100)then
        write(lunmsg,*)' Time step reduced 100 times without convergence'
        ier=13
        endif
        return
        end
 
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
! 29Oct2006 fgtok
