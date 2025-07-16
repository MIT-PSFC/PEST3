subroutine eqi_bbox(zrhoi, zBmin,zchiBmin, zBmax,zchiBmax, ierr)

  use xplasma_definitions
  use eqi_rzbox_module

  !  given a flux surface (zrho) find Bmin and Bmax -- axisymmetry assumed.

  IMPLICIT NONE

  !  input:

  REAL*8 zrhoi                      ! surface to be examined

  !  output:
  REAL*8 zBmin                      ! Bmin on surface
  real*8 zchiBmin                   ! (chi) of minimum
  REAL*8 zBmax                      ! Bmax on surface
  real*8 zchiBmax                   ! (chi) of maximum

  integer ierr                      ! completion code, 0=OK

  external feq_Bmmx

  !---------------------------------------------------

  REAL*8 zchi(nsrch1)             ! search points (poloidal angle coord)
  REAL*8 zB(0:1,nsrch1)           ! Z,dZ/dchi @ each search pt

  REAL*8 zchi0,zchi1,zchimmx

  REAL*8 zeps,zeta,ztmp,zrho

  integer i,nsrch,id_Bmod
  integer iminB,imaxB,iflag

  REAL*8 zrhovec(nsrch1),zbvec(nsrch1,2)

  logical imask,axisymm

  real*8, parameter :: C2PI = 6.2831853071795862D+00
  real*8, parameter :: ZERO = 0.0d0
  !---------------------------------------------------

  zrho=max(1.0d-12,zrhoi)

  imask=.FALSE.
  ierr=0

  call xplasma_global_info(sp,ierr, axisymm=axisymm)
  if(ierr.ne.0) return

  if(.not.axisymm) then
     ierr=107
     call xplasma_errmsg_append(sp, &
          ' ?eqi_bbox:  non-axisymmetry not yet supported!')
  endif
  if(ierr.ne.0) return

  call xplasma_find_item(sp,'Bmod',id_bmod,ierr)
  if(ierr.ne.0) then
     call xplasma_errmsg_append(sp,' ?eqi_bbox:  mod(B) spline not found.')
  endif
  if(ierr.ne.0) return

  zchi0=ZERO
  zchi1=C2PI
  nsrch=nsrch0

  do i=1,nsrch
     zrhovec(i)=zrho
     zchi(i)=zchi0+(i-1)*(zchi1-zchi0)/nsrch
  enddo
  zchi(nsrch1)=zchi1
  zrhovec(nsrch1)=zrhovec(1)

  call eqi_fBmod(nsrch1,zrhovec,zchi,2,nsrch1,zBvec,ierr)
  if(ierr.ne.0) return

  do i=1,nsrch1
     zb(0:1,i)=zbvec(i,1:2)
     if(i.eq.1) then
        zBmin=zB(0,1)
        zBmax=zBmin
        iminB=1
        imaxB=1
        zchiBmin=zchi(1)
        zchiBmax=zchi(1)
     else
        if(zB(0,i).gt.zBmax) then
           imaxB=i
           zBmax=zB(0,i)
           zchiBmax=zchi(i)
        endif
        if(zB(0,i).lt.zBmin) then
           iminB=i
           zBmin=zB(0,i)
           zchiBmin=zchi(i)
        endif
     endif
  enddo

  !  find exact min & max using rootfinder on derivative

  zeps=C2PI*2.0d-15
  zeta=zBmax*2.0d-15

  !  cobble for if extrema is at end point and there is a small numerical
  !  discrepancy (could happen for problems with "symmetry"

  iflag=0
  if(zB(1,1)*zB(1,nsrch1).lt.0.0) then
     if(abs(zB(1,1)-zB(1,nsrch1)).lt.(10*zeta)) then
        zeta=1.1*max(zeta,abs(zB(1,1)),abs(zB(1,nsrch1)))
     else
        call xplasma_errmsg_append(sp, &
             ' %eqi_bbox: dB/dchi inconsistent at theta extrema.')
        iflag=1
     endif
  endif

  do i=1,nsrch
     if(zB(1,i).gt.ZERO .and. zB(1,i+1).lt.ZERO) then
        call zridderx(1,imask,zchi(i),zchi(i+1),zeps,zeta, &
             feq_Bmmx,zchimmx,ierr, 1,zrho,1,ztmp,1)
        if(ierr.ne.0) return
        if(ztmp.gt.zBmax) then
           zBmax=ztmp
           zchiBmax=zchimmx
        endif
     endif
     if(zB(1,i).lt.ZERO .and. zB(1,i+1).gt.ZERO) then
        call zridderx(1,imask,zchi(i),zchi(i+1),zeps,zeta, &
             feq_Bmmx,zchimmx,ierr, 1,zrho,1,ztmp,1)
        if(ierr.ne.0) return
        if(ztmp.lt.zBmin) then
           zBmin=ztmp
           zchiBmin=zchimmx
        endif
     endif
  enddo

  ierr=0
  if(iflag.eq.1) ierr=9999

end subroutine eqi_bbox

!-----------------------------------------------
subroutine feq_Bmmx(ivec,imask,zchi,zanswer, &
     ivecd,zinput,ninput,zoutput,noutput)

  use xplasma_definitions
  use eqi_rzbox_module

  IMPLICIT NONE

  !  function for zridderx -- dR/dchi @ some rho; return mod(B); also
  !  will find extrema by looking for dB/dchi=0

  integer ivec,ivecd
  logical imask(ivec)
  integer ninput,noutput
  real*8 zanswer(ivec)
  REAL*8 zchi(ivec),zinput(ivecd,ninput),zoutput(ivecd,noutput)

  !---------------------------

  REAL*8 zrho_use(ivec),zchi_use(ivec),zvals(ivec,0:1)
  integer ict_vald(10),ierr,i,jvec

  !---------------------------

  jvec=0
  do i=1,ivec
     if(.not.imask(i)) then
        jvec=jvec+1
        zrho_use(jvec)=zinput(i,1)
        zchi_use(jvec)=zchi(i)
     endif
  enddo
  if(jvec.eq.0) return              ! nothing to be done

  call eqi_fBmod(jvec,zrho_use,zchi_use,2,jvec,zvals(1:jvec,0:1),ierr)
  if(ierr.ne.0) return

  jvec=0
  do i=1,ivec
     if(.not.imask(i)) then
        jvec=jvec+1
        zoutput(i,1)=zvals(jvec,0)  ! B field
        zanswer(i)=zvals(jvec,1)    ! dB/dchi
     endif
  enddo

end subroutine feq_Bmmx

!---------------------------------------------------
subroutine eqi_bbvec(zrho,ivec,znormb,zbmin,zbmax,zbvals,zthbr1,zthbr2,ierr)

  !  find vector of B values and corresponding theta locations on
  !  a rho surface (axisymmetry assumed).  Only non-singular flux
  !  surfaces should be specified (i.e. avoid magnetic axis).

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE

  real*8, intent(in) :: zrho  ! flux surface
  integer, intent(in) :: ivec ! size of input/output vectors
  real*8, intent(in) :: znormb(ivec) ! ascending sequence in range [0,1]

  !  znormb(j) indicates j'th B value wanted, per the formula
  !      znormb(j) = (zbvals(j)-Bmin)/(Bmax-Bmin)
  !  where Bmin, Bmax are the min & max mod(B) on the zrho surface.

  !  i.e. znormb(1) = 0 refers to Bmin;
  !       znormb(ivec) = 1 refers to Bmax.

  !  Bmin and Bmax will be found during execution of this routine.

  real*8, intent(out) :: zbmin,zbmax   ! min & max mod(B) on surface
  real*8, intent(out) :: zbvals(ivec)  ! B values corresponding to znormb
  real*8, intent(out) :: zthbr1(ivec)  ! "lower branch" theta sequence
  real*8, intent(out) :: zthbr2(ivec)  ! "upper branch" theta sequence

  integer, intent(out) :: ierr      ! status code (0=OK)

  !  B(thrbr1(j))=B(thrbr2(j))=zbvals(j) will be satisfied.

  !  if theta(Bmax) > theta(Bmin) then
  !      {thbr1} in [theta(Bmax)-2*pi,theta(Bmin)]
  !      {thbr2} in [theta(Bmin),theta(Bmax)]

  !  if theta(Bmax) < theta(Bmin) then
  !      {thbr1} in [theta(Bmax),theta(Bmin)]
  !      {thbr2} in [theta(Bmin),theta(Bmax)+2*pi]

  !-------------------------------------------------------------
  external f_eqi_bfind

  logical imask(ivec)
  integer i,ichk,ibr,j
  real*8 zeps,zeta,ztol,zrhoi,ztest
  real*8 zthbmin,zthbmax(2)

  real*8 zth_srch0(nsrch0),zb_srch0(nsrch0),zrho_srch0(nsrch0)
  real*8 ztha(ivec),zthb(ivec),zthans(ivec),zinput(ivec,2)
  real*8 zdum(ivec,1)

  logical :: axisymm

  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ONE = 1.0d0
  real*8, parameter :: C2PI = 6.2831853071795862D+00
  !-------------------------------------------------------------

  !  perform error checks...

  zbmin=0
  zbmax=0
  zbvals=0
  zthbr1=0
  zthbr2=0

  ierr=0
  call xplasma_global_info(sp,ierr,axisymm=axisymm)
  if(ierr.ne.0) return

  if(.not.axisymm) then
     ierr=107
     call xplasma_errmsg_append(sp, &
          ' ?eqi_bbvec:  non-axisymmetry not yet supported!')
     return
  endif

  if(ivec.le.0) then
     call xplasma_errmsg_append(sp, &
          ' ?eqi_bbvec: "ivec" is not a positive number.')
     ierr=9999
  endif
  if(ierr.gt.0) return

  ztol=2.0d-15

  zrhoi=max(100*ztol,min(ONE,zrho))

  if((maxval(znormb).gt.1.0d0+ztol).or.(minval(znormb).lt.-ztol)) then
     call xplasma_errmsg_append(sp,' ?eqi_bbvec: "znormb" not in range [0,1].')
     ierr=ierr+1
  endif

  do i=1,ivec-1
     if(znormb(i+1).le.znormb(i)) then
        ierr=ierr+1
        call xplasma_errmsg_append(sp, &
             ' ?eqi_bbvec: "znormb" not strict ascending.')
        exit
     endif
  enddo

  if(ierr.ne.0) return
  !---------------------------------------
  !  OK.  get Bmin & Bmax & corresponding theta values

  call eqi_bbox(zrhoi,zbmin,zthbmin,zbmax,zthbmax(1),ierr)
  if(ierr.ne.0) return

  zeps=C2PI*ztol
  zeta=ztol*zBmax

  if(zthbmax(1).gt.zthbmin) then
     zthbmax(2)=zthbmax(1)
     zthbmax(1)=zthbmax(2)-C2PI
  else
     zthbmax(2)=zthbmax(1)+C2PI
  endif

  !  the target B values:

  do i=1,ivec
     zbvals(i) = zBmin + znormb(i)*(zBmax-zBmin)
  enddo

  !  set up for vectorized root finder call to find theta values
  !  corresponding to desired B values

  imask=.FALSE.

  !  check output extrema

  ierr=0
  ichk=0

  if(znormb(1).le.ztol) then
     imask(1)=.TRUE.
     ichk=ichk+1
     zbvals(1)=zBmin
     zthbr1(1)=zthBmin
     zthbr2(1)=zthBmin
  endif

  if(znormb(ivec).ge.1.0d0-ztol) then
     imask(ivec)=.TRUE.
     ichk=ichk+1
     zbvals(ivec)=zBmax
     zthbr1(ivec)=zthbmax(1)
     zthbr2(ivec)=zthbmax(2)
  endif

  if(ichk.ge.ivec) return

  !  set up initial search grid: both sides...

  zrho_srch0=zrhoi
  do ibr=1,2
     do i=1,nsrch0
        zth_srch0(i)=(zthbmin*(nsrch0-i)+zthbmax(ibr)*(i-1))/(nsrch0-1)
     enddo
     call eqi_fBmod(nsrch0,zrho_srch0,zth_srch0,1,nsrch0,zb_srch0,ierr)
     if(ierr.ne.0) exit
     do i=1,ivec
        do j=1,nsrch0-1
           ztest=(zb_srch0(j)-zbvals(i))*(zb_srch0(j+1)-zbvals(i))
           if(ztest.le.ZERO) then
              ztha(i)=min(zth_srch0(j),zth_srch0(j+1))
              zthb(i)=max(zth_srch0(j),zth_srch0(j+1))
              if(abs(zb_srch0(j)-zbvals(i)).lt.zeta) then
                 if(ibr.eq.1) then
                    zthbr1(i)=zth_srch0(j)
                 else
                    zthbr2(i)=zth_srch0(j)
                 endif
                 imask(i)=.TRUE.
              else if(abs(zb_srch0(j+1)-zbvals(i)).lt.zeta) then
                 if(ibr.eq.1) then
                    zthbr1(i)=zth_srch0(j+1)
                 else
                    zthbr2(i)=zth_srch0(j+1)
                 endif
                 imask(i)=.TRUE.
              endif
              exit                  ! take the first valid interval
           endif
        enddo
        if(ibr.eq.1) then
           zthans(i)=zthbr1(i)
        else
           zthans(i)=zthbr2(i)
        endif
        zinput(i,1)=zrhoi
        zinput(i,2)=zbvals(i)
     enddo

     call zridderx(ivec,imask,ztha,zthb,zeps,zeta,f_eqi_bfind, &
          zthans,ierr,ivec,zinput,2,zdum,1)
     if(ierr.ne.0) exit

     if(ibr.eq.1) then
        zthbr1=zthans
     else
        zthbr2=zthans
     endif
  enddo

end subroutine eqi_bbvec

!-----------------------------------------------
subroutine f_eqi_Bfind(ivec,imask,zchi,zanswer, &
     ivecd,zinput,ninput,zoutput,noutput)

  use xplasma_definitions
  use eqi_rzbox_module

  IMPLICIT NONE

  !  function for zridderx -- dR/dchi @ some rho; return mod(B); also
  !  will find extrema by looking for dB/dchi=0

  integer ivec,ivecd
  logical imask(ivec)
  integer ninput,noutput
  real*8 zanswer(ivec)
  REAL*8 zchi(ivec),zinput(ivecd,ninput),zoutput(ivecd,noutput)

  !---------------------------

  REAL*8 zrho_use(ivec),zchi_use(ivec),zvals(ivec,1)
  integer ierr,i,jvec

  !---------------------------

  jvec=0
  do i=1,ivec
     if(.not.imask(i)) then
        jvec=jvec+1
        zrho_use(jvec)=zinput(i,1)
        zchi_use(jvec)=zchi(i)
     endif
  enddo
  if(jvec.eq.0) return              ! nothing to be done

  call eqi_fBmod(jvec,zrho_use,zchi_use,1,jvec,zvals,ierr)
  if(ierr.ne.0) return

  jvec=0
  do i=1,ivec
     if(.not.imask(i)) then
        jvec=jvec+1
        zanswer(i)=zvals(jvec,1)-zinput(i,2)  ! B field - Btarget
     endif
  enddo

end subroutine f_eqi_Bfind
