      subroutine mkquin2(sub,x,nx,th,nth,fquin)
C
C  create a data set for interpolation, from evaluation of function and
C  derivatives & 2nd derivatives
C
      external sub                      ! passed subroutine (x,th,output)
      real x(nx)                        ! x coordinate array
      real th(nth)                      ! th coordinate array
C
      real fquin(0:5,nx,nth)            ! function data / spline coeff array
C
C
C  sub's interface:  subroutine sub(xi,th,ra(6))
C    ra(1) -- value of fcn f
C    ra(2) -- df/dx
C    ra(3) -- df/dth
C    ra(4) -- d2f/dx2
C    ra(5) -- d2f/dth2
C    ra(6) -- d2f/dx.dth
C
C----------------------------
C
      do ix=1,nx
         do ith=1,nth
            call sub(x(ix),th(ith),fquin(0,ix,ith))
         enddo
      enddo
C
      return
      end
