      SUBROUTINE scrpreset(jmax,ns,imom)
 
      USE mintrp
      USE scrunch_inc1
      IMPLICIT NONE
      INTEGER jmax, ns, imom
      ntheta = jmax
      mpol = imom+1
      mnd=mpol*(nphi+1)
      mpol2=2*mpol
      mpol3=3*mpol
      mpol4=4*mpol
      n2=mpol4+ntheta
      nu      = jmax
      ntheta0 = jmax
      ncmax   = ns
      moption = 3
      numc = ns
      facmid = 0.0D0  ! magic numbers
      facedg = 0.0D0  ! don't change
      facdx0 = 0.0D0  ! (defaults set in scriscrunch).
      RETURN
      END SUBROUTINE scrpreset
 
      subroutine scrpreset_wr(icon)
 
        ! write on LUN icon the info needed to call scrpreset
 
        USE mintrp
        USE scrunch_inc1
      IMPLICIT NONE
        integer, intent(in) :: icon
 
        write(icon) ncmax,nu,mpol
 
      end subroutine scrpreset_wr
! 29Oct2006 fgtok -s r8_precision.sub "r8con.csh conversion"
