      subroutine getmpa(zrmc,zymc,mimom,imom,zrin,zrout,zymp)
C
C  DMC 18 March 1994:
C  given the updown asymmetric moments of a surface, define an
C  approximate midplane to the surface.  Return the inner and
C  outer R intercepts of this midplane with the surface, and
C  return the height of the surface above y=0.
C
C  rev dmc 25 Mar 1994:  use
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER imom,mimom,inth,imm,ith
!============
! idecl:  explicitize implicit REAL declarations:
      REAL zth,zr,zy
!============
      REAL zrmc(0:mimom,2)  ! R moments
      REAL zymc(0:mimom,2)  ! Y moments
C
C  mimom gives the moments 1st array dimension; imom gives the number
C   of *active* moments; imom.le.mimom ...
C
      REAL zrin            ! inner intercept R (output)
      REAL zrout           ! outer intercept R (output)
      REAL zymp            ! midplane vertical displacement (output)
C
      parameter (inth=400)
      REAL zrcon(inth),zycon(inth)
      REAL zcos(0:imom),zsin(0:imom)
C---------------------------------------------------------------------
C
      zcos(0)=1
      zsin(0)=0
      do ith=1,inth-1
        zth=(6.2831853071795862*(ith-1))/(inth-1)
        zr=0.0
        zy=0.0
        call sincos(zth,imom,zsin(1:imom),zcos(1:imom))
        do imm=0,imom
          zr=zr+zrmc(imm,1)*zcos(imm)+zrmc(imm,2)*zsin(imm)
          zy=zy+zymc(imm,1)*zcos(imm)+zymc(imm,2)*zsin(imm)
        enddo
        zrcon(ith)=zr
        zycon(ith)=zy
      enddo
C
      zrcon(inth)=zrcon(1)
      zycon(inth)=zycon(1)
C
      call getmpa_ry(zrcon,zycon,inth,zrin,zrout,zymp)
C
      return
      end
C---------------------------------------------------------------------
C
      subroutine getmpa_ry(zrcon,zycon,inth,zrin,zrout,zymp)
C
C  find midplane elevation and R intercepts based on centroid, starting
C  from a closed contour of (R,Y) pairs
C
      implicit NONE

      integer, intent(in) :: inth
      real, intent(in) :: zrcon(inth),zycon(inth)  ! closed contour

      real, intent(out) :: zrin,zrout    ! midplane intercepts (approx.)
      real, intent(out) :: zymp          ! midplane elevation (centroid)
C
C  local:
C
      real :: zrmp,zr,zy,zrp,zyp,zytest,zrans
      integer :: ith
C
C  find centroid
C
      call plcentr(zrcon,zycon,inth,zrmp,zymp)
C
C  find approximate midplane intercept locations
C
      zr=zrcon(1)
      zy=zycon(1)
C
      do ith=2,inth
c
        zrp=zr
        zyp=zy
        zr=zrcon(ith)
        zy=zycon(ith)
c
        zytest=(zy-zymp)*(zymp-zyp)
        if(zytest.ge.0.0) then
          zrans=(zr*(zymp-zyp)+zrp*(zy-zymp))/(zy-zyp)
          if(zrans.gt.zrmp) then
            zrout=zrans
          else
            zrin=zrans
          endif
        endif
      enddo
C
      return
      end
