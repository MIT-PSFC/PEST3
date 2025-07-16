      SUBROUTINE scrorderdes(rval,zval,xaxis,yaxis,ier)
       USE scrunch_inc1
      IMPLICIT NONE
       INTEGER 	ier ! <orderdes.f>
       REAL*8 	dy ! <orderdes.f>
       REAL*8 	dx ! <orderdes.f>
       REAL*8 	yaxis ! <orderdes.f>
       REAL*8 	y ! <orderdes.f>
       REAL*8 	xaxis ! <orderdes.f>
       REAL*8 	x ! <orderdes.f>
       REAL*8 	residue ! <orderdes.f>
       REAL*8 	save_z ! <orderdes.f>
       REAL*8 	save_r ! <orderdes.f>
       INTEGER 	next ! <orderdes.f>
       INTEGER 	j ! <orderdes.f>
       REAL*8 	shortest ! <orderdes.f>
       INTEGER 	i1 ! <orderdes.f>
       INTEGER 	ip1 ! <orderdes.f>
       INTEGER 	i ! <orderdes.f>
       REAL*8 	old_dist ! <orderdes.f>
        REAL*8 rval(*),zval(*),new_dist,tempr(ntheta),tempz(ntheta)
!**********************************************************************
!       Program ORDER : orders points on a magnetic surface at a
!       fixed toroidal plane and assigns right-handed circulation
!       around flux surface. XAXIS, YAXIS:  Polar-type axis (must lie
!       inside curve to check sign of rotation)
!**********************************************************************
        old_dist = 1.D20
        do 10 i = 1,nthetax-1
        ip1 = i + 1
        i1 = i
        shortest = 1.D20
 15              do 20 j = ip1,nthetax
                 if(i1.gt.1)old_dist = (rval(i1-1)-rval(j))**2 &
     &                              + (zval(i1-1)-zval(j))**2
                 new_dist =(rval(i1)-rval(j))**2 + (zval(i1)-zval(j))**2
                 if((new_dist.le.old_dist).and. &
     &           (new_dist.lt.shortest))then
                         next = j
                         shortest = new_dist
                 endif
 20     continue
!**********************************************************************
!       Swap nearest point (next) with current point (ip1)
!**********************************************************************
        if(shortest.ge.1.D10)then
        save_r = rval(i-1)
        rval(i-1) = rval(i)
        rval(i)= save_r
        save_z = zval(i-1)
        zval(i-1) = zval(i)
        zval(i)= save_z
        i1 = i1 - 1
        ip1 = ip1 - 1
        go to 15
        endif
        save_r = rval(ip1)
        rval(ip1) = rval(next)
        rval(next)= save_r
        save_z = zval(ip1)
        zval(ip1) = zval(next)
        zval(next)= save_z
 10     continue
!**********************************************************************
!       Check that xaxis,yaxis is inside surface and
!       ascertain that the angle rotates counterclockwise
!       using Cauchy's theorem in "complex"-plane
!**********************************************************************
        residue = 0.D0
        do 30 i = 1,nthetax-1
        x = .5D0*(rval(i)+rval(i+1)) - xaxis
        y = .5D0*(zval(i)+zval(i+1)) - yaxis
        dx= rval(i+1) - rval(i)
        dy= zval(i+1) - zval(i)
        residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.D-10)
 30     continue
        x = .5D0*(rval(1)+rval(nthetax)) - xaxis
        y = .5D0*(zval(1)+zval(nthetax)) - yaxis
        dx= rval(1) - rval(nthetax)
        dy= zval(1) - zval(nthetax)
        residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.D-10)
 
        if( residue .lt.(-.90D0*twopi))then
        do 40 i = 2,nthetax
        j = nthetax - i + 2
        tempr(i) = rval(j)
        tempz(i) = zval(j)
 40     continue
        do 50 i = 2,nthetax
        rval(i) = tempr(i)
        zval(i) = tempz(i)
 50     continue
        else if( abs(residue) .lt. (.90D0*twopi) )then
        write(lunmsg,*)' The magnetic axis is not enclosed by boundary '
        ier=11
	return
        endif
 
        return
        end
 
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
