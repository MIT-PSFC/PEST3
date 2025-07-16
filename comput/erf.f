      real function erf(realarg)
C
C  single precision error function
C     simply making use of double precision function in r8slatec
C
      external derf             ! link r8slatec if this routine is called!
C
      real realarg
C
      double precision darg,dans
C
      darg=realarg
      dans=derf(darg)
C
      erf=dans
      return
      end
