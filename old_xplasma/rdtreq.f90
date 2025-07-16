subroutine rdtreq_hdr(lun,runid,time,tavg,pcur,tflux,ntheta,nsurf)
!
!  read ascii file containing MHD equilibrium data from TRANSP.
!  header information -- runid, scalars & no. of theta pts & surfaces.
!
  implicit NONE
!
  integer, intent(in) :: lun           ! lun of file open for reading
!
  character*(*) runid                  ! runid
!
  real, intent(out) :: time            ! data time
  real, intent(out) :: tavg            ! +/- averaging time
  real, intent(out) :: pcur            ! plasma current, Amps
  real, intent(out) :: tflux           ! enclosed toroidal flux, Webers
!
  integer, intent(out) :: ntheta       ! no. of theta grid points
  integer, intent(out) :: nsurf        ! no. of flux surfaces
!
!--------------------------------------------
!
  character(80) :: line
  integer i,i1,i2,str_length
!
!--------------------------------------------
! format stmts
!
1001 format(1x,a)
1002 format(1x,1pe13.6)
1003 format(1x,i6)
!
!--------------------------------------------
!
  do
     read(lun,1001) line
     if(line(1:12).eq.'TRANSP runid') then
 
        read(lun,1001) line
 
        ! use trailing non-blank word as runid label
 
        i2=str_length(line)
        i1=1
        do i=i2,1,-1
           if(line(i:i).eq.' ') then
              i1=i+1
              exit
           endif
        enddo
 
        runid=line(i1:i2)
 
     else if(line(1:5).eq.'time:') then
        read(lun,1002) time
     else if(line(1:5).eq.'tavg:') then
        read(lun,1002) tavg
     else if(line(1:5).eq.'pcur:') then
        read(lun,1002) pcur
     else if(line(1:6).eq.'tflux:') then
        read(lun,1002) tflux
     else if(line(1:7).eq.'ntheta:') then
        read(lun,1003) ntheta
     else if(line(1:6).eq.'nsurf:') then
        read(lun,1003) nsurf
        exit
     else
        call errmsg_exit('?rdtreq_hdr:  unexpected line:  '//line)
     endif
  enddo
  return
end subroutine rdtreq_hdr
!============================================================================
subroutine rdtreq_data(lun,ntheta,nsurf,theta,rho,psi,pmhd,qmhd,gmhd,R,Z)
!
!  read ascii file containing MHD equilibrium data from TRANSP.
!  profile information -- dimensions from header
!
  implicit NONE
!
  integer, intent(in) :: lun           ! lun of file open for reading
  integer, intent(in) :: ntheta        ! theta dimension
  integer, intent(in) :: nsurf         ! flux surface dimension
!
  real, intent(out) :: theta(ntheta)   ! theta grid
  real, intent(out) :: rho(nsurf)      ! rho grid
  real, intent(out) :: psi(nsurf)      ! psi (Webers/rad)
  real, intent(out) :: pmhd(nsurf)     ! pressure (Pascals)
  real, intent(out) :: qmhd(nsurf)     ! q profile (dimensionless)
  real, intent(out) :: gmhd(nsurf)     ! g profile (R*B, m*Tesla)
!
  real, intent(out) :: R(ntheta,nsurf) ! R(theta,psi)
  real, intent(out) :: Z(ntheta,nsurf) ! Z(theta,psi)
!
!--------------------------------------------
!
  integer i,j
  character(80) :: line
!
!--------------------------------------------
! format stmts
!
1001 format(1x,a)
1004 format(5(1x,1pe13.6))
!
!--------------------------------------------
!
  do
     read(lun,1001) line
     if(line(1:6).eq.'theta:') then
        read(lun,1004) (theta(i),i=1,ntheta)
     else if(line(1:4).eq.'rho:') then
        read(lun,1004) (rho(i),i=1,nsurf)
     else if(line(1:4).eq.'psi:') then
        read(lun,1004) (psi(i),i=1,nsurf)
     else if(line(1:5).eq.'pmhd:') then
        read(lun,1004) (pmhd(i),i=1,nsurf)
     else if(line(1:5).eq.'qmhd:') then
        read(lun,1004) (qmhd(i),i=1,nsurf)
     else if(line(1:5).eq.'gmhd:') then
        read(lun,1004) (gmhd(i),i=1,nsurf)
     else if(line(1:2).eq.'R:') then
        read(lun,1004) ((R(i,j),i=1,ntheta),j=1,nsurf)
     else if(line(1:2).eq.'Z:') then
        read(lun,1004) ((Z(i,j),i=1,ntheta),j=1,nsurf)
        exit
     else
        call errmsg_exit('?rdtreq_data:  unexpected line:  '//line)
     endif
  enddo
  return
end subroutine rdtreq_data
 
