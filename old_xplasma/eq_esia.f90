subroutine eq_esia(lun_esia,ins,int1,ztflux,id_pmhd,id_nen) 

  implicit NONE

  ! This routine writes the MHD information in the ``compact'' form
  ! required by Leonid's ESI architecture. (One code where
  ! the equilibrium file is used is in Roscoe's new  
  ! (non-Boozer) version of orbit.)
  !  Kumar -- 04202004

  ! integer, parameter :: ins = 21 ! Fixed by ESC
  ! integer, parameter :: int1=65
  integer, parameter :: r8=selected_real_kind(12,100)
  integer ierr,luntrm,lnout
  integer lun_esia       ! logical unit where to write file

  real*8 :: r8time = 1.0d34

  integer :: ik,jk,imid,id_tflux=0,id_pflux=0,id_pmhd,id_q=0,id_g=0,id_nen
  integer :: ins,int1,id_ne=0,ifail,iflag,impol
  integer :: ifcns(3)
  real*8  :: ztime  ! ESC related; not presently used
  real*8 :: zdrho,zpi,zpi2,zxi_axis,zxi_bdy,zraxis,zdum1,zdum2,ztflux
  real*8 :: zrho(ins),zpmhd(ins),zq(ins) 
  real*8 :: zg(ins),zg_r(ins),zg_rr(ins)
  real*8 :: zne(ins),zne_r(ins)
  real*8 :: ztf_r(ins),ztf_rr(ins),zpf_r(ins),zpf_rr(ins)
  real*8 :: zpr_r(ins),zpr_rr(ins),zpf_rb(ins)
  real*8 :: zp(ins),zp_r(ins),zt(ins),zt_r(ins)  
  real*8, dimension(:,:), allocatable :: zrhon,zchin 
  real*8, dimension(:,:), allocatable :: hfun
  real*8, dimension(:,:,:), allocatable :: zfun
  integer :: iregion(int1)
  real*8 :: zrhob(int1),zphi(int1),zchi(int1),zR(int1),zZ(int1)

  luntrm = lnout(0)

! Prepare Rho and Chi coordinates
!  (same as a and gq in Leonid's represenatation)

     zrho(1)=0.0_r8
     zdrho=1.0_r8/(float(ins-1))
     do ik=2,ins
       zrho(ik)=zrho(ik-1)+zdrho
     enddo
     zrho(ins)=1.0_r8

     zpi = 3.14159265358979312_r8 
     zpi2=2.0_r8*zpi
     do ik=1,int1
      zchi(ik)=zpi2*(ik-1)/(int1-1)
     enddo

! Prepare g and g' 
! (same as F and F' in ESI representation]

     ! g=R*Bt (T*m)

     call eq_gfnum('G', id_g)

     call eq_rgetf(ins,zrho,id_g,0,zg,ierr)
     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("G",...)'
     endif

     call eq_rgetf(ins,zrho,id_g,1,zg_r,ierr)
     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("G_r",...)'
     endif

     call eq_rgetf(ins,zrho,id_g,2,zg_rr,ierr)
     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("G_rr",...)'
     endif

! Prepare gF'/a, gF'', gY'/a, (gY'/a)' (ESC) 

!    First get derivatives of poloidal flux

     call eq_gfnum('PSI',id_pflux)
     call eq_rgetf(ins,zrho,id_pflux,1,zpf_r,ierr)
     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("PF_R",...)'
     endif

!    Store a back-up of unaltered derivative of poloidal flux

     zpf_rb(1:ins)=zpf_r(1:ins)

     call eq_gfnum('PSI',id_pflux)
     call eq_rgetf(ins,zrho,id_pflux,2,zpf_rr,ierr)
     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("PF_RR",...)'
     endif

!    Calculate derivative of (psi'/a) using
!    quotient rule
 
     do ik=2,ins
      zpf_rr(ik)=-(zrho(ik)*zpf_rr(ik)-zpf_r(ik))/(zrho(ik)*zrho(ik))
     enddo

! Following based on observations of Leonid's data
!  (Zero second derivative is enforced at the center)

     zpf_rr(1)=0.0_r8

     do ik=2,ins
      zpf_r(ik)=-zpf_r(ik)/(zrho(ik))
     enddo
     zpf_r(1)=2.0_r8*zpf_r(2)-zpf_r(3)

! Now derivative of toroidal flux

!    (1/X)d(Phi)/dX=2*tflux from the
!    definition of X=sqrt(Phi/tflux)

     ztf_r=2.0_r8*ztflux/zpi2

     ztf_rr=0.0_r8

! Prepare T,Ta,P, an Pa (ESC) 

     do ik=2,ins
      zt(ik)=-zg(ik)*zg_r(ik)/zpf_rb(ik)
     enddo
     zt(1)=2.0_r8*zt(2)-zt(3)

     do ik=2,ins
      zt_r(ik)= -(zpf_rb(ik)*(zg(ik)*zg_rr(ik)+zg_r(ik)*zg_r(ik))             &
                -zg(ik)*zg_r(ik)*zpf_rr(ik))/(zpf_rb(ik)*zpf_rb(ik))
     enddo
!    zt_r(1)=2.0_r8*zt_r(2)-zt_r(3)

! Following based on observations of Leonid's data
!  (Zero derivative is enforced at the center)

     zt_r(1)=0.0_r8


     ! Pressure in Pascals 

     call eq_rgetf(ins,zrho,id_pmhd,1,zpr_r,ierr)
     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("zpr_r",...)'
     endif

     call eq_rgetf(ins,zrho,id_pmhd,2,zpr_rr,ierr)
     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("zpr_rr",...)'
     endif


     ! convert pressure to mu0 Pascal
     ! (4*pi*1.0E-7*1.0E6)
     ! 1.0E6 is for converting 1/cm^3 to 1/m^3 
     ! zpr_r=0.4_r8*zpi*zpr_r
     ! zpr_rr=0.4_r8*zpi*zpr_rr

    zpr_r=0.4_r8*zpi*(1.0e-6*zpr_r)
    zpr_rr=0.4_r8*zpi*(1.0e-6*zpr_rr)

     do ik=2,ins
      zp(ik)=-zpr_r(ik)/zpf_rb(ik)
     enddo
     zp(1)=2.0_r8*zp(2)-zp(3)

     ! Second derivative using quotient rule

     do ik=2,ins
      zp_r(ik)=-(zpf_rb(ik)*zpr_rr(ik)-zpf_rr(ik)*zpr_r(ik))/(zpf_rb(ik)*zpf_rb(ik))
     enddo

!    zp_r(1)=2.0_r8*zp_r(2)-zp_r(3)

! Following based on observations of Leonid's data
!  (Zero derivative is enforced at the center)

     zp_r(1)=0.0_r8

! Prepare r, z and their derivatives

     call eq_gfnum('R',ifcns(1))
     call eq_gfnum('Z',ifcns(2))
     call eq_gfnum('BMOD',ifcns(3))

     allocate(zrhon(int1, ins))
     allocate(zchin(int1, ins))
     allocate(hfun(int1, ins))
     allocate(zfun(int1,ins,12))

     zrhon=spread(zrho,dim=1,ncopies=int1)
     zchin=spread(zchi,dim=2,ncopies=ins) 

     call eq_hrhochi(int1*ins,zrhon,zchin,3,ifcns,                    &
         &                (/1,1,1,0,0,1/),int1*ins,zfun,ierr)
     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("R,Z, |B| etc.",...)'
     endif

     hfun(:,:)=0.0_r8

     zne(:)=0.0_r8
     zne_r(:)=0.0_r8

! Prepare Ne and d(Ne)/da

     call eq_rgetf(ins,zrho,id_nen,0,zne,ierr)

     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("NE",...)'
     endif

     call eq_rgetf(ins,zrho,id_nen,1,zne_r,ierr)

     if(ierr.ne.0) then
       write(luntrm,*)  &
         ' ?eq_esia:  error in eq_rgetf("NE_r",...)'
     endif

     zne(1:ins)=zne(1:ins)*1E-20
     zne_r(1:ins)=zne_r(1:ins)*1E-20

     write (lun_esia,'(a)') '!!! Do not edit this file'
     write (lun_esia,'(i3,a,i3,a)') int1,' x ',ins,' - numbers of poloidal x radial data points'
     write (lun_esia,'(a)') 'gq:'

     write (lun_esia,'(4(e24.16))') (zchi(ik),ik=1,int1)

     write (lun_esia,'(TR23,a,TR23,a,TR22,a)') 'a','F','Fa'
     write (lun_esia,'(3(e24.16))') (zrho(ik),zg(ik),zg_r(ik),ik=1,ins)

     write (lun_esia,'(TR19,a,TR20,a,TR19,a,TR16,a)') "gF'/a","gF''","gY'/a","(gY'/a)'"
     write (lun_esia,'(4(e24.16))') (ztf_r(ik),ztf_rr(ik),zpf_r(ik),zpf_rr(ik),ik=1,ins)

     write (lun_esia,'(TR23,a,TR20,a,TR23,a,TR22,a)') "T","Ta''","P","Pa"
     write (lun_esia,'(4(e24.16))') (zt(ik),zt_r(ik),zp(ik),zp_r(ik),ik=1,ins)

     write (lun_esia,'(TR23,a,TR20,a,TR19,a,TR14,a)') "r","r'_a","r'_gq","r''_{a,gq}"
     write (lun_esia,'(4(e24.16))') ((zfun(ik,jk,1),zfun(ik,jk,2),zfun(ik,jk,3),        &
                                 zfun(ik,jk,4),ik=1,int1),jk=1,ins)

     write (lun_esia,'(TR23,a,TR20,a,TR19,a,TR14,a)') "z","z'_a","z'_gq","z''_{a,gq}"
     write (lun_esia,'(4(e24.16))') ((-zfun(ik,jk,5),-zfun(ik,jk,6),-zfun(ik,jk,7),        &
                                 -zfun(ik,jk,8),ik=1,int1),jk=1,ins)

     write (lun_esia,'(TR23,a,TR20,a,TR19,a,TR14,a)') "B","B'_a","B'_gq","B''_{a,gq}"
     write (lun_esia,'(4(e24.16))') ((zfun(ik,jk,9),zfun(ik,jk,10),zfun(ik,jk,11),        &
                                 zfun(ik,jk,12),ik=1,int1),jk=1,ins)

     write (lun_esia,'(TR18,a,TR13,a,TR12,a,TR9,a)') "gh'_gq",         &
                                         "gh''_{a,gq}","gh''_{gq,gq}","gh'''_{a,gq,gq}"
     write (lun_esia,'(4(e24.16))') ((hfun(ik,jk),hfun(ik,jk),hfun(ik,jk),                &
                                 hfun(ik,jk),ik=1,int1),jk=1,ins)

     write (lun_esia,'(a,TR10,a,TR6,a)') "Te=Ti=15keV N=3.26e+20 a",   &
                                                  "ne [10^20/m^3]","dne/da [10^20/m^3]"     
     write (lun_esia,'(3(e24.16))') (zrho(ik),zne(ik),zne_r(ik),ik=1,ins)

     
     deallocate(zrhon)
     deallocate(zchin)
     deallocate(hfun)
     deallocate(zfun)

end subroutine eq_esia 
