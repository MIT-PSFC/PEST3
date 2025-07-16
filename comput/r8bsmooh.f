      subroutine r8bsmooh(zprof,znorm)
 
      real*8 zprof(*),znorm(*)
 
      call r8bsmoop(zprof,znorm)   ! patch holes
 
      call r8bsmoo(zprof,znorm)    ! smooth result
 
      return
      end
C------------------------------------------------------------
      subroutine r8bsmoop(zprof,znorm)
C
C  patch "holes" (zero values) in profile,
C  then smooth.
C
C  done for MC radial profiles that are formed by summing of the
C  profile itself and its normalization
C
C    zprof(j) = [MC sum]f*dw
C    zwsum(j) = [MC sum]dw
C
C  and at normalization time
C
C    if (zwsum(j).gt.0.0) then
C         zprof(j)=zprof(j)/zwsum(j)
C    else
C         zprof(j)=0.0
C    endif
C
C  which can leave "holes" in the profile in case of lousy statistics.
C
C  *** SO ***
C  fill these holes with flat extrapolation & linear interpolation
C  before smoothing with "r8bsmoo"
C
C-----------------------------
C
      use r8bsmoo_mod
      implicit NONE
C
      real*8 zprof(MJ)
      real*8 znorm(MJ)
C
      real*8 z1(mj),zx1,zf1,zx2,zf2,zx,zf
C
      integer j,jnz,jj,isrch,jsave1,jsave2
C
C-----------------------------
C
C  0.  normalize
C
      do j=lcentr,ledge
         zprof(j)=zprof(j)/znorm(j)
         z1(j)=1.0
      enddo
C
C  1.  fill in from left
C
      do j=lcentr,ledge
         if(zprof(j).ne.0.0) go to 10
      enddo
C
      go to 1000                        ! all zero:  exit now...
C
 10   continue
      jnz=j
      do jj=lcentr,jnz-1
         zprof(jj)=zprof(jnz)
      enddo
C
C  2.  fill in from right
C
      do j=ledge,lcentr,-1
         if(zprof(j).ne.0.0) go to 20
      enddo
C
 20   continue
      jnz=j
      do jj=ledge,jnz+1,-1
         zprof(jj)=zprof(jnz)
      enddo
C
C  3.  fill in remaining gaps with linear interpolation
C
      isrch=1
      do j=lcentr+1,ledge-1
         if(isrch.eq.1) then
C  find start of zeroed section
            if(zprof(j).eq.0.0) then
               jsave1=j-1
               zx1=xi(jsave1,2)
               zf1=zprof(jsave1)
               isrch=2
            endif
         else if(isrch.eq.2) then
C  find end of zeroed section
            if(zprof(j).ne.0.0) then
               jsave2=j
               zx2=xi(jsave2,2)
               zf2=zprof(jsave2)
               isrch=1
C  apply patch
               do jj=jsave1+1,jsave2-1
                  zx=xi(jj,2)
                  zf=zf1+(zf2-zf1)*(zx-zx1)/(zx2-zx1)
                  zprof(jj)=zf
               enddo
            endif
         endif                          ! isrch mode
      enddo                             ! j loop
C
C  un-normalize
C
      do j=lcentr,ledge
         zprof(j)=zprof(j)*znorm(j)
      enddo
C
 1000 continue
      return
      end
