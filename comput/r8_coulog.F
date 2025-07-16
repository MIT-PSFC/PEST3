#include "fpreproc/fp_const.h"
      subroutine r8_coulog(et,at,zt,denb,tempb,ab,zb,cln,nb,Bmag,
     >                    llstix)
C
C (see the warning below.)
C
C calculates the coulomb logarithm for a test particle of energy et
C (kev), atomic mass at and atomic charge zt, colliding with NB
C maxwellian background species with densities denb (/m**3),
C temperatures tempb (keV), atomic masses ab, atomic charges zb.  The
C electrons should be included as one of the background species, with
C Z=-1 and A=(1./1836.1).  The list of coulomb logarithms is returned
C in cln.  The magnetic field should be input as bmag (Tesla).
C
C Note that the test particle loglambda's produced by this subroutine
c are not symmetric, i.e. it is not true that
c loglam_alpha_beta=loglam_beta_alpha.  If such symmetry is desired,
c one should set Et=(3/2)T_t and call this routine twice with species
c alpha and beta switched, and take the average value of log_lambda
c which is returned.
C
C Written by Greg Hammett, 16-May-1989.
C
C The exact details of the coulomb logarithm are probably not worth
C worrying about too much.  The recipe I use here is similiar to the
C one given in the NRL plasma formulary.  So it includes the quantum
C mechanics correction which is usually important for collisions with
C electrons.  Rather than completely ignoring the Debye shielding by
C particles with vthermal < v, it contains a smooth transition which I
C have rigorously shown from the Balescu-Lenard operator is correct in
C the limits vthermal>>v and vthermal<<v.  However, the Balescu-Lenard
C operator predicts a substantial amount of Cerenkov radiation of
C plasma waves if the particle speed exceeds the electron thermal
C speed.  I am ignoring this Cerenkov radiation (and its associated
C drag).  To look at it properly, Nat Fisch claims that one would have
C to consider the reabsorption of the Cerenkov radiation as well.
C
C Warning:  GWH 3/1/91:  the formulas I am using here provide very
C little Debye shielding if the test particles are much faster than the
c thermal electrons.  In that case, one would need to worry about whether or
C not this is correct.  There appear to be fundamental problems with
C the Balescu-Lenard operator.
C
C The main references for this stuff are D.V. Sivukhin, Reviews of
C Plasma Physics, VOl. 4 (1966) (especially good on quantum corrections,
C physical insight), Trubnikov, Reviews of Plasma Physics, Vol. 1
C (1963) (simple, complete derivation of Fokker-Planck equation,
C complete with justification for the Coulomb Logarithm in terms of
C deriving the collision operator from the dielectric response of the
C plasma), Rob Goldston's Ph.D. thesis (references on strong magnetic
C field corrections, collisions with multi-electron impurities),
C Krommes' class notes on the derivation of the Landau operator from the
C Balescu-Lenard operator, and my notes on Krommes' notes relating to
C the problem when the test particle is much faster than the thermal
C electron velocity.
C
C============
C idecl:  explicitize implicit INTEGER declarations:
!
      IMPLICIT NONE
      REAL*8, PARAMETER :: ZERO = 0.0D0 ! _fgtok_constant_
      REAL*8, PARAMETER :: ONE = 1.0D0 ! _fgtok_constant_
      REAL*8, PARAMETER :: TWO = 2.0D0 ! _fgtok_constant_
      REAL*8, PARAMETER :: C1000 = 1000.0D0 ! _fgtok_constant_
!
      INTEGER nb,i
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 sum,omega2,vrel2,rmax,rmincl,rminqu,rmin
C============
      REAL*8 et,at,zt,Bmag
      REAL*8 denb(nb),tempb(nb),ab(nb),zb(nb),cln(nb)
 
      Logical llstix       ! tbt 7/10/95
 
c
c first calculate the maximum impact parameter, rmax:
c
c rmax=( sum_over_species omega2/vrel2 )**(-1/2)
c
c where,
c omega2 = omega_p**2 + omega_c**2
c vrel2=T/m+2.*E_t/m_t
c
      sum=ZERO
      do 100 i=1,nb
         omega2=1.74D0*zb(i)**2/ab(i)*denb(i)
     1	+9.18D15*zb(i)**2/ab(i)**2*Bmag**2
         vrel2=9.58D10*(tempb(i)/ab(i) + TWO*et/at)
         sum=sum+omega2/vrel2
 100  continue
      rmax=sqrt(ONE/sum)
 
c next calculate rmin, including quantum corrections.  The classical
C rmin is:
c
C rmincl = e_alpha e_beta / (m_ab vrel**2)
C
C where m_ab = m_a m_b / (m_a+m_b) is the reduced mass.
c vrel**2 = 3 T_b/m_b + 2 E_a / m_a
c (Note:  the two different definitions of vrel2 used in this code
c are each correct for their application.)
c
C The quantum rmin is:
C
C rminqu = hbar/( 2 exp(0.5) m_ab vrel)
C
C and the proper rmin is the larger of rmincl and rminqu
C
      do 200 i=1,nb
         vrel2=9.58D10*(3*tempb(i)/ab(i)+2*et/at)
         rmincl=0.13793D0*abs(zb(i)*zt)*(ab(i)+at)/ab(i)/at/vrel2
         rminqu=1.9121D-8*(ab(i)+at)/ab(i)/at/sqrt(vrel2)
         rmin=max(rmincl,rminqu)
         cln(i)=log(rmax/rmin)
         if(cln(i) .lt. ONE) then
            write(6,*) 'warning from COULOG: coulomb logarithms < 1!'
            cln(i)=ONE
         endif
200   continue
 
c debug section:  set all log(lambda)'s to a simple formula for
c log(lambda_e) to benchmark with Stix's solutions:
 
      if(llstix) then
      do 300 i=1,nb
         cln(i)=24.0D0-log(sqrt(denb(1)/1.D6)/
     >      (tempb(1)*C1000))
 300  continue
      endif
 
      return
      end
! 11May2001 fgtok -s r8_precision.sub "r8con.csh conversion"
! 11May2001 fgtok
