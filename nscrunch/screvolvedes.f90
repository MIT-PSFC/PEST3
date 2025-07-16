      SUBROUTINE screvolvedes(g11)
      USE scrunch_inc1
      IMPLICIT NONE
       INTEGER 	i ! <evolvedes.f>
       REAL*8 	fac ! <evolvedes.f>
       REAL*8 	b1 ! <evolvedes.f>
       REAL*8 	otav ! <evolvedes.f>
       REAL*8 	dtau ! <evolvedes.f>
       REAL*8 	g11 ! <evolvedes.f>
       REAL*8 	ftest ! <evolvedes.f>
       REAL*8 	bmax ! <evolvedes.f>
        data bmax /0.15D0/
 
        ftest=gnorm/g11
        dtau=abs(1.D0-ftest)
        g11=gnorm
        otav=dtau/delt
        dtau=delt*otav+.001D0
        if(dtau.gt.bmax)dtau=bmax
        b1=1.D0-.5D0*dtau
        fac=1.D0/(1.D0+.5D0*dtau)
        do 10 i=1,n2
        xdot(i)=fac*(xdot(i)*b1-delt*gvec(i))
 10     xvec(i)=xvec(i)+xdot(i)*delt
        return
        end
 
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
