c-------------------------------------------------------
c ported from f77 xplasma dmc Feb. 2006
c
      subroutine eqi_rzbox(zrhoi,zphi,
     >     iwant_Rmin,zRmin,zth_Rmin,iwant_Rmax,zRmax,zth_Rmax,
     >     iwant_Zmin,zZmin,zth_Zmin,iwant_Zmax,zZmax,zth_Zmax,
     >     ierr)
c
      use xplasma_definitions
      use eqi_rzbox_module
c
c  given a flux surface (zrho) find a minimal box (zRmin,zRmax)x(zZmin,zZmax)
c  (all toroidal angles) which incloses the surface
c
c  this amounts to finding the Rmin,Rmax,Zmin,Zmax of the surface.
c
      IMPLICIT NONE
c
      REAL*8 zrhoi                      ! surface to be boxed
      REAL*8 zphi                       ! toroidal angle
c
c  note: zphi not currently used-- reserved for future non-azisymmetry upgrade
c
      logical :: iwant_Rmin             ! TRUE for Rmin output
      logical :: iwant_Rmax             ! TRUE for Rmax output
      logical :: iwant_Zmin             ! TRUE for Zmin output
      logical :: iwant_Zmax             ! TRUE for Zmax output
c
c  output as requested-- approximate values *always* returned; set iwant
c  switch for exact (to mach precision) value...
c
      REAL*8 zRmin                      ! Rmin of enclosing box in (R,Z)
      REAL*8 zRmax                      ! Rmax of enclosing box in (R,Z)
      REAL*8 zZmin                      ! Zmin of enclosing box in (R,Z)
      REAL*8 zZmax                      ! Zmax of enclosing box in (R,Z)

c  and corresponding theta (or "chi") values, accurate as requested:
      REAL*8 zth_Rmin
      REAL*8 zth_Rmax
      REAL*8 zth_Zmin
      REAL*8 zth_Zmax
C
      integer ierr                      ! completion code, 0=OK
C
      external feq_Rmmx,feq_Zmmx
C
C---------------------------------------------------
c
      REAL*8 zchi(nsrch1)             ! search points (poloidal angle coord)
      REAL*8 zR(0:1,nsrch1)           ! R,dR/dchi @ each search pt
      REAL*8 zZ(0:1,nsrch1)           ! Z,dZ/dchi @ each search pt
      real*8 zRZbuf(nsrch1,4)         ! interp. buffer
      integer iflist(4),idchis(4)
c
      REAL*8 zchi0,zchi1,zchimmx
c
      REAL*8 zeps,zeta,ztmp,zrho
c
      integer i
      integer iminR,iminZ,imaxR,imaxZ
      integer iminR2,iminZ2,imaxR2,imaxZ2
C
      logical imask
C
      REAL*8 zrhovec(nsrch1)
c-----------------------------------
      integer :: lRnum,lZnum,iertmp
c
      real*8 :: bdy_tol
      logical :: axisymm,bflag
      integer :: id_bdy_info
c
      real*8, parameter :: c2pi = 6.2831853071795862D+00
      real*8, parameter :: csmall = 1.0d-14
      real*8, parameter :: czero = 0.0d0
      real*8, parameter :: cone  = 1.0d0
C---------------------------------------------------
C
      ierr=0
      zrmin=0
      zrmax=0
      zzmin=0
      zzmax=0
C
      zrho=max(100*csmall,zrhoi)
C
      call xplasma_global_info(sp,ierr, axisymm=axisymm, bdytol=bdy_tol)
      if(ierr.ne.0) return
C
      id_bdy_info = 0
      bflag=.FALSE.
      if(abs(cone-zrho).lt.bdy_tol) then
         bflag = axisymm
         call xplasma_find_item(sp,'__RZ_BDY_INFO',id_bdy_info,ierr,
     >        nf_noerr=.TRUE.)
      endif
C
      imask=.FALSE.
c
      zchi0=0
      zchi1=c2pi
C
      do i=1,nsrch1
         zrhovec(i)=zrho
         zchi(i)=zchi0+(i-1)*(zchi1-zchi0)/(nsrch0)
      enddo
      zchi(nsrch1)=zchi1
c
      call xplasma_common_ids(sp,ierr, id_R=lRnum,id_Z=lZnum)
      if((lRnum.eq.0).or.(lZnum.eq.0)) then
         ierr=2504
         return
      endif

      iflist(1)=lRnum
      iflist(2)=lZnum
      iflist(3)=lRnum
      iflist(4)=lZnum
      idchis(1:2)=0
      idchis(3:4)=1

      if(bflag.and.(id_bdy_info.gt.0)) then
         call eqi_rzbox_binfo_get(cone,zrzbuf,ierr)
      else
         call xplasma_RZeval_2d(sp,iflist,
     >        xplasma_rho_coord,zrhovec,xplasma_theta_coord,zchi,
     >        zrzbuf,ierr, ideriv2s=idchis)
      endif
      if(ierr.ne.0) return
c
      if(bflag.and.(id_bdy_info.eq.0)) then
         call eqi_rzbox_binfo_set(cone,zrzbuf,id_bdy_info,ierr)
      endif
      if(ierr.ne.0) return
c
      do i=1,nsrch1
         zR(0,i)=zrzbuf(i,1)    ! R
         zZ(0,i)=zrzbuf(i,2)    ! Z
         zR(1,i)=zrzbuf(i,3)   ! dR/dchi
         zZ(1,i)=zrzbuf(i,4)   ! dZ/dchi
         if(i.eq.1) then
            zRmin=zR(0,1)
            zRmax=zRmin
            zZmin=zZ(0,1)
            zZmax=zZmin
            iminR=1
            imaxR=1
            iminZ=1
            imaxZ=1
         else
            if(zR(0,i).gt.zRmax) then
               imaxR=i
               zRmax=zR(0,i)
            endif
            if(zR(0,i).lt.zRmin) then
               iminR=i
               zRmin=zR(0,i)
            endif
            if(zZ(0,i).gt.zZmax) then
               imaxZ=i
               zZmax=zZ(0,i)
            endif
            if(zZ(0,i).lt.zZmin) then
               iminZ=i
               zZmin=zZ(0,i)
            endif
         endif
      enddo
C
      zth_Rmin=zchi(iminR)
      zth_Rmax=zchi(imaxR)
      zth_Zmin=zchi(iminZ)
      zth_Zmax=zchi(imaxZ)
C
C  check the branch cut
C  the derivatives should match there, but may fail to do so exactly
C  for numerical reasons.  If the derivatives change sign, it is because
C  there is an extrema there and the derivatives are extremely close to
C  zero.
C
      if(zr(1,1)*zr(1,nsrch1).le.czero) then
         zr(1,1)=czero
         zr(1,nsrch1)=zr(1,1)
      else
         zr(1,1)=(zr(1,1)+zr(1,nsrch1))/2
         zr(1,nsrch1)=zr(1,1)
      endif
C
      if(zz(1,1)*zz(1,nsrch1).le.czero) then
         zz(1,1)=czero
         zz(1,nsrch1)=zz(1,1)
      else
         zz(1,1)=(zz(1,1)+zz(1,nsrch1))/2
         zz(1,nsrch1)=zz(1,1)
      endif
C
C  find exact min & max using rootfinder on derivative
C
      zeps=c2pi*csmall
      zeta=csmall*max(zRmax,abs(zZmin),abs(zZmax))/c2pi
C
      do i=1,nsrch0
         if(iwant_Rmax .and. 
     >        zR(1,i).gt.czero .and. zR(1,i+1).lt.czero) then
            call zridderx(1,imask,zchi(i),zchi(i+1),zeps,zeta,
     >         feq_Rmmx,zchimmx,ierr,
     >         1,zrho,1,ztmp,1)
            if(ierr.ne.0) then
               ierr=9999         ! not expected
               return
            endif
            if(ztmp.gt.zRmax) then
               zRmax=ztmp
               zth_Rmax=zchimmx
            endif
         endif
         if(iwant_Rmin .and.
     >        zR(1,i).lt.czero .and. zR(1,i+1).gt.czero) then
            call zridderx(1,imask,zchi(i),zchi(i+1),zeps,zeta,
     >         feq_Rmmx,zchimmx,ierr,
     >         1,zrho,1,ztmp,1)
            if(ierr.ne.0) then
               ierr=9999         ! not expected
               return
            endif
            if(ztmp.lt.zRmin) then
               zRmin=ztmp
               zth_Rmin=zchimmx
            endif
         endif
         if(iwant_Zmax .and.
     >        zZ(1,i).gt.czero .and. zZ(1,i+1).lt.czero) then
            call zridderx(1,imask,zchi(i),zchi(i+1),zeps,zeta,
     >         feq_Zmmx,zchimmx,ierr,
     >         1,zrho,1,ztmp,1)
            if(ierr.ne.0) then
               ierr=9999         ! not expected
               return
            endif
            if(ztmp.gt.zZmax) then
               zZmax=ztmp
               zth_Zmax=zchimmx
            endif
         endif
         if(iwant_Zmin .and. 
     >        zZ(1,i).lt.czero .and. zZ(1,i+1).gt.czero) then
            call zridderx(1,imask,zchi(i),zchi(i+1),zeps,zeta,
     >         feq_Zmmx,zchimmx,ierr,
     >         1,zrho,1,ztmp,1)
            if(ierr.ne.0) then
               ierr=9999         ! not expected
               return
            endif
            if(ztmp.lt.zZmin) then
               zZmin=ztmp
               zth_Zmin=zchimmx
            endif
         endif
      enddo
C
      ierr=0
C
      return
      end
C---------------------------------------------------------------------
      subroutine feq_Rmmx(ivec,imask,zchi,zanswer,
     >   ivecd,zinput,ninput,zoutput,noutput)
c
      use xplasma_definitions
      use eqi_rzbox_module
C
      IMPLICIT NONE
C
C  function for zridderx -- dR/dchi @ some rho; return R also
C  will find extrema by looking for dR/dchi=0
C
      integer ivec,ivecd
      logical imask(ivec)
      integer ninput,noutput
c
c  note do not want "intent(out)"-- not all elements will be set...
c
      REAL*8 zchi(ivec),zanswer(ivec),
     >   zinput(ivecd,ninput),zoutput(ivecd,noutput)
c
c---------------------------
c
      REAL*8 zrho_use(ivec),zchi_use(ivec),zvals(ivec,0:1)
      integer i,jvec,ierr,lRnum
c
      integer :: iflist(2),iderivs(2)
c---------------------------
c
      jvec=0
      do i=1,ivec
         if(.not.imask(i)) then
            jvec=jvec+1
            zrho_use(jvec)=zinput(i,1)
            zchi_use(jvec)=zchi(i)
         endif
      enddo
      if(jvec.eq.0) return              ! nothing to be done
c
      call xplasma_common_ids(sp,ierr, id_R=lRnum)

      iflist=lRnum
      iderivs(1)=0
      iderivs(2)=1

      call xplasma_RZeval_2d(sp,iflist,
     >     xplasma_rho_coord,zrho_use(1:jvec),
     >     xplasma_theta_coord,zchi_use(1:jvec),
     >     zvals(1:jvec,0:1),ierr, ideriv2s=iderivs)
      if(ierr.ne.0) then
         ierr=9999   ! not expected here
         return
      endif

      jvec=0
      do i=1,ivec
         if(.not.imask(i)) then
            jvec=jvec+1
            zoutput(i,1)=zvals(jvec,0)  ! R
            zanswer(i)=zvals(jvec,1)    ! dR/dchi
         endif
      enddo
c
      return
      end
C-----------------------------------------------
      subroutine feq_Zmmx(ivec,imask,zchi,zanswer,
     >   ivecd,zinput,ninput,zoutput,noutput)
c
      use xplasma_definitions
      use eqi_rzbox_module
C
      IMPLICIT NONE
C
C  function for zridderx -- dZ/dchi @ some rho; return Z also
C  will find extrema by looking for dZ/dchi=0
C
      integer ivec,ivecd
      logical imask(ivec)
      integer ninput,noutput
c
c  note do not want "intent(out)"-- not all elements will be set...
c
      REAL*8 zchi(ivec),zanswer(ivec),
     >   zinput(ivecd,ninput),zoutput(ivecd,noutput)
c
c---------------------------
c
      REAL*8 zrho_use(ivec),zchi_use(ivec),zvals(ivec,0:1)
      integer i,jvec,ierr,lZnum
c
      integer :: iflist(2),iderivs(2)
c---------------------------
c
      jvec=0
      do i=1,ivec
         if(.not.imask(i)) then
            jvec=jvec+1
            zrho_use(jvec)=zinput(i,1)
            zchi_use(jvec)=zchi(i)
         endif
      enddo
      if(jvec.eq.0) return              ! nothing to be done
c
      call xplasma_common_ids(sp,ierr, id_Z=lZnum)

      iflist=lZnum
      iderivs(1)=0
      iderivs(2)=1

      call xplasma_RZeval_2d(sp,iflist,
     >     xplasma_rho_coord,zrho_use(1:jvec),
     >     xplasma_theta_coord,zchi_use(1:jvec),
     >     zvals(1:jvec,0:1),ierr, ideriv2s=iderivs)
      if(ierr.ne.0) then
         ierr=9999   ! not expected here
         return
      endif

      jvec=0
      do i=1,ivec
         if(.not.imask(i)) then
            jvec=jvec+1
            zoutput(i,1)=zvals(jvec,0)  ! Z
            zanswer(i)=zvals(jvec,1)    ! dZ/dchi
         endif
      enddo
c
      return
      end
C-----------------------------------------------
