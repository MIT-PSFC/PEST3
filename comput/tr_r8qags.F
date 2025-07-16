#include "fpreproc/library_names.h"
c freeware integrator adapted for TRANSP: added auxilliary argument
c for integrand (dmc Apr 2003)
c
      subroutine tr_r8qags(f,aux,a,b,epsabs,epsrel,
     *   result,abserr,neval,ier,
     *   limit,lenw,last,iwork,work)
c***begin prologue  dqags
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             (end-point) singularities, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & prog. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral  i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x,aux). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            aux(*) - auxilliary arguments for f, passed thru
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account. however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table. it is presumed that
c                             the requested tolerance cannot be
c                             achieved, and that the returned result is
c                             the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero.except when limit or lenw is invalid,
c                             iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end
c                    with ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, detemines the
c                    number of significant elements actually in the work
c                    arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers
c                    to the error estimates over the subintervals
c                    such that work(limit*3+iwork(1)),... ,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2),
c                    and k = limit+1-last otherwise
c
c            work  - double precision
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end-points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end-points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last)
c                     contain the error estimates.
c
c***references  (none)
c***routines called  dqagse,xerror
c***end prologue  dqags
c
c
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 a,abserr,b,epsabs,epsrel,f,result,work,aux(*)
      integer ier,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  dqags
      ier = 6
      neval = 0
      last = 0
      result = 0.0_R8
      abserr = 0.0_R8
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for dqagse.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call tr_r8qagse(f,aux,a,b,epsabs,epsrel,limit,result,abserr,neval,
     *  ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
cxx      if(ier.ne.0) call xerror(26habnormal return from dqags,26,ier,lvl)
cxx      return
      end
      subroutine tr_r8qagse(
     *   f,aux,a,b,epsabs,epsrel,limit,result,abserr,neval,
     *   ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqagse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             (end point) singularities, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x,aux). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            aux(*) - auxilliary arguments for f, passed thru
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b)
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved, and that the
c                             returned result is the best which can be
c                             obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             epsabs.le.0 and
c                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
c                             result, abserr, neval, last, rlist(1),
c                             iord(1) and elist(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left end points
c                     of the subintervals in the partition of the
c                     given integration range (a,b)
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            rlist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced in the
c                     subdivision process
c
c***references  (none)
c***routines called  d1mach,dqelg,dqk21,dqpsrt
c***end prologue  dqagse
c
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,correc,abs,defabs,defab1,defab2,d1mach,max,
     *  dres,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,errmax,
     *  error1,error2,erro12,errsum,ertest,f,oflow,resabs,reseps,result,
     *  res3la,rlist,rlist2,small,uflow,aux(*)
      integer id,ier,ierro,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     *  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     * res3la(3),rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine dqelg (rlist2 should be of dimension
c            (limexp+2) at least).
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2 containing
c                       the part of the epsilon table which is still
c                       needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left interval
c           *****2    - variable for the right interval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered up
c                       to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine is
c                       attempting to perform extrapolation i.e. before
c                       subdividing the smallest interval we try to
c                       decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagse
      epmach = d1mach(4)
c
c            test on validity of parameters
c            ------------------------------
      ier = 0
      neval = 0
      last = 0
      result = 0.0_R8
      abserr = 0.0_R8
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0_R8
      elist(1) = 0.0_R8
      if(epsabs.le.0.0_R8.and.
     *   epsrel.lt.max(0.5E+02_R8*epmach,0.5E-28_R8))
     *   ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      ierro = 0
      call tr_r8qk21(f,aux,a,b,result,abserr,defabs,resabs)
c
c           test on accuracy.
c
      dres = abs(result)
      errbnd = max(epsabs,epsrel*dres)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.1.0E+02_R8*epmach*defabs.and.abserr.gt.errbnd)
     *   ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     *  abserr.eq.0.0_R8) go to 140
c
c           initialization
c           --------------
c
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      nres = 0
      numrl2 = 2
      ktmin = 0
      extrap = .false.
      noext = .false.
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ksgn = -1
      if(dres.ge.(0.1E+01_R8-0.5E+02_R8*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest error
c           estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5_R8*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call tr_r8qk21(f,aux,a1,b1,area1,error1,resabs,defab1)
        call tr_r8qk21(f,aux,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 15
        if(abs(rlist(maxerr)-area12).gt.0.1E-04_R8*abs(area12)
     *  .or.erro12.lt.0.99_R8*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max(epsabs,epsrel*abs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of subintervals
c           equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(max(abs(a1),abs(b2)).le.(0.1E+01_R8+0.1E+03_R8*epmach)*
     *  (abs(a2)+0.1E+04_R8*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
   30   call tr_r8qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(errsum.le.errbnd) go to 115
c ***jump out of do-loop
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(abs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over the
c           larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
c ***jump out of do-loop
          if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call tr_r8qelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1E-02_R8*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = max(epsabs,epsrel*abs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 100
c
c           prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5_R8
        erlarg = errsum
        go to 90
   80   small = abs(b-a)*0.375_R8
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if(ier+ierro.eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0_R8.and.area.ne.0.0_R8) go to 105
      if(abserr.gt.errsum) go to 115
      if(area.eq.0.0_R8) go to 130
      go to 110
  105 if(abserr/abs(result).gt.errsum/abs(area)) go to 115
c
c           test on divergence.
c
  110 if(ksgn.eq.(-1).and.max(abs(result),abs(area)).le.
     * defabs*0.1E-01_R8) go to 130
      if(0.1E-01_R8.gt.(result/area).or.(result/area).gt.0.1E+03_R8
     * .or.errsum.gt.abs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0_R8
      do 120 k = 1,last
         result = result+rlist(k)
  120 continue
      abserr = errsum
  130 if(ier.gt.2) ier = ier-1
  140 neval = 42*last-21
  999 return
      end
      subroutine tr_r8qk21(f,aux,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk21
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  21-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x,aux). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c            aux(*) - auxilliary arguments for f, passed thru
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 21-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 10-point gauss rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk21
c
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 a,absc,abserr,b,centr,abs,dhlgth,max,min,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     *  resg,resk,reskh,result,uflow,wg,wgk,xgk,aux(*)
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 21-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 10-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 10-point gauss rule
c
c           wgk    - weights of the 21-point kronrod rule
c
c           wg     - weights of the 10-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.066671344308688137593568809893332_R8/
      data wg  (  2) / 0.149451349150580593145776339657697_R8/
      data wg  (  3) / 0.219086362515982043995534934228163_R8/
      data wg  (  4) / 0.269266719309996355091226921569469_R8/
      data wg  (  5) / 0.295524224714752870173892994651338_R8/
c
      data xgk (  1) / 0.995657163025808080735527280689003_R8/
      data xgk (  2) / 0.973906528517171720077964012084452_R8/
      data xgk (  3) / 0.930157491355708226001207180059508_R8/
      data xgk (  4) / 0.865063366688984510732096688423493_R8/
      data xgk (  5) / 0.780817726586416897063717578345042_R8/
      data xgk (  6) / 0.679409568299024406234327365114874_R8/
      data xgk (  7) / 0.562757134668604683339000099272694_R8/
      data xgk (  8) / 0.433395394129247190799265943165784_R8/
      data xgk (  9) / 0.294392862701460198131126603103866_R8/
      data xgk ( 10) / 0.148874338981631210884826001129720_R8/
      data xgk ( 11) / 0.000000000000000000000000000000000_R8/
c
      data wgk (  1) / 0.011694638867371874278064396062192_R8/
      data wgk (  2) / 0.032558162307964727478818972459390_R8/
      data wgk (  3) / 0.054755896574351996031381300244580_R8/
      data wgk (  4) / 0.075039674810919952767043140916190_R8/
      data wgk (  5) / 0.093125454583697605535065465083366_R8/
      data wgk (  6) / 0.109387158802297641899210590325805_R8/
      data wgk (  7) / 0.123491976262065851077958109831074_R8/
      data wgk (  8) / 0.134709217311473325928054001771707_R8/
      data wgk (  9) / 0.142775938577060080797094273138717_R8/
      data wgk ( 10) / 0.147739104901338491374841515972068_R8/
      data wgk ( 11) / 0.149445554002916905664936468389821_R8/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 10-point gauss formula
c           resk   - result of the 21-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk21
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5_R8*(a+b)
      hlgth = 0.5_R8*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 21-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0_R8
      fc = f(centr,aux)
      resk = wgk(11)*fc
      resabs = abs(resk)
      do 10 j=1,5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,aux)
        fval2 = f(centr+absc,aux)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
      do 15 j = 1,5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,aux)
        fval2 = f(centr+absc,aux)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   15 continue
      reskh = resk*0.5_R8
      resasc = wgk(11)*abs(fc-reskh)
      do 20 j=1,10
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0_R8.and.abserr.ne.0.0_R8)
     *   abserr = resasc*
     *   min(0.1E+01_R8,(0.2E+03_R8*abserr/resasc)**1.5_R8)
      if(resabs.gt.uflow/(0.5E+02_R8*epmach)) abserr = max
     *   ((epmach*0.5E+02_R8)*resabs,abserr)
      return
      end
! 15Oct2002 fgtok -s qag.table "quadpack renames for imbedding in TRANSP system"
! 22Oct2003 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Oct2003 fgtok
