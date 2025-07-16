      SUBROUTINE scrprintdes
      USE scrunch_inc1
      USE scrunch_inc2
      IMPLICIT NONE
       INTEGER 	mn ! <printdes.f>
        REAL*8 rmc(nphi2,1),rms(nphi2,1),zmc(nphi2,1),zms(nphi2,1)
        REAL*8 integrate(nphi2),argi(nphi2)
        if (ioflagc) write(lunmsg,30)
 30     format( &
     &  /'c  mb  nb    RmncB       RmnsB       ZmncB       ZmnsB')
        do 70 mn=1,mpnt
        if (ioflagc) write(lunmsg,80) mn-1,0,xvec(mn),xvec(mpol+mn), &
     &    xvec(mpol2+mn),xvec(mpol3+mn)
        rmomb(mn,1)=xvec(mn)
        rmomb(mn,2)=xvec(mpol+mn)
        zmomb(mn,1)=xvec(mpol2+mn)
        zmomb(mn,2)=xvec(mpol3+mn)
 70     continue
 80     format(i5,i4,4e12.4)
        return
        end
 
