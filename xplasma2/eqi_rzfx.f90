subroutine eqi_rzfx(r,nr,z,nz,inside,fcn)

  !  extrapolate based on a nearest neighbour method weighted average
  !  method the fcn to beyond the limiter regions

  implicit NONE
  integer nr,nz                     ! grid dimensions
  real*8 r(nr),z(nz)                ! grid
  logical :: inside(nr,nz)          ! T means inside limiter; F: outside
  real*8 fcn(nr,nz)                 ! fcn needing extrapolation

  !  where inside(i,j)=.TRUE., the fcn is already defined

  !  ad hoc method:
  !  form list of values that border the extrapolation region; average
  !  over these values with inverse square of distance weighting to derive
  !  the extrapolation

  !-----------------------------------------------------------
  !  arrays amply big enough for "edge value list"...

  real*8 rlist(nr*nz),zlist(nr*nz),fval(nr*nz)
  real*8 zwt,zr2sum,zfr2sum

  integer i,j,k,ip,im,jp,jm,ilist
  !-----------------------------------------------------------

  ilist=0

  do j=1,nZ
     jm=max(1,j-1)
     jp=min(nZ,j+1)

     do i=1,nR
        im=max(1,i-1)
        ip=min(nR,i+1)

        if(inside(i,j)) then
           if((.not.inside(i,jm)).or.(.not.inside(i,jp)).or. &
                (.not.inside(im,j)).or.(.not.inside(ip,j))) then

              !  add edge point to list

              ilist=ilist+1
              rlist(ilist)=r(i)
              zlist(ilist)=z(j)
              fval(ilist)=fcn(i,j)

           endif
        endif

     enddo
  enddo

  !  ok use list for extrapolation

  do j=1,nZ
     do i=1,nR

        if(.not.inside(i,j)) then

  !  extrapolate here

           zfr2sum=0.0d0
           zr2sum=0.0d0
           do k=1,ilist
              zwt=1.0d0/((r(i)-rlist(k))**2+(z(j)-zlist(k))**2)
              zr2sum=zr2sum+zwt
              zfr2sum=zfr2sum+fval(k)*zwt
           enddo

           fcn(i,j)=zfr2sum/zr2sum

        endif
     enddo
  enddo

end subroutine eqi_rzfx
