      SUBROUTINE scrfixaraydes
      USE scrunch_inc1
      IMPLICIT NONE
       INTEGER 	m ! <fixaraydes.f>
       INTEGER 	n ! <fixaraydes.f>
       INTEGER 	nn0 ! <fixaraydes.f>
       INTEGER 	ntor ! <fixaraydes.f>
       INTEGER 	l ! <fixaraydes.f>
!**********************************************************************
!       This routine store toroidal and poloidal mode number arrays
!**********************************************************************
        l = 0
        ntor = max(1,nphi-1)
        nn0 = 1 - (ntor+1)/2
        do 5 n=1,ntor
 5      nn(n) = (nn0 + (n-1))*nfp
        do 40 m=1,mpolx
        mm(m)=m-1
        do 40 n=1,ntor
        if(mm(m).eq.0 .and. nn(n) .lt. 0) go to 40
        l=l+1
 30     continue
        m1(l)=mm(m)
        n1(l)=nn(n)
 40     continue
        mpnt=l
        dnorm=2.D0/nthetax
        do 60 m=1,mpolx
        dm1(m)=m-1
        faccon(m) = .125D0*dnorm/(1.D0+dm1(m))**pexp
        if(m.eq.1) go to 60
        xmpq(m,1) = dm1(m)**pexp
        xmpq(m,2) = dm1(m)**qexp
        xmpq(m,4) = xmpq(m,1)
        if(m.le.2) xmpq(m,4) = 0.D0
        xmpq(m,3) = sqrt(xmpq(m,4))
 60     continue
        return
        end
 
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
