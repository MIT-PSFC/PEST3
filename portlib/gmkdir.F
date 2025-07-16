#define BUFSIZ 300
!
      subroutine gmkdir(root_in,rpath_in,ier)
!
!  make a multi-level subdirectory starting from root
!
!  unix example:  call gmkdir('HOME','.cache/sub1/sub2/sub3',ier)
!    => $HOME/.cache/sub1/sub2/sub3 directory is created (if it does
!       not already exist).  Each higher level directory is also created
!       as needed.
!
!  vms example:  call gmkdir('SYS$LOGIN','[zcache.sub1.sub2.sub3]',ier)
!    => subdirectory chain [zcache.sub1.sub2.sub3] is created under the
!       user's login directory (SYS$LOGIN) -- each subdirectory level is
!       created as needed.
!
      implicit NONE
!
!  input:
!
      character*(*) root_in             ! root -- start point (must exist)
      character*(*) rpath_in            ! subdirectory chain
!
!  output:
      integer ier                       ! error code, 0 = OK
!
!-------------------------------------------------------------
!
      character*BUFSIZ start
      character*BUFSIZ next
!
      character*250 root,rpath
!
!  DMC -- procedure to check directories /u & /p can be unreliable
!  and is not needed.  If ' ' is passed for root, start check a 
!  little further down the path than "/".
!
      integer, parameter :: iroot_param_srch = 6
      integer :: ii,istat
!
      character*512 zcmd
!-------------------------------------------------------------
!
      if(root_in.ne.' ') then
         root = root_in
         rpath = rpath_in
      else
         if(rpath_in(1:1).eq.'/') then
            do ii=iroot_param_srch,1,-1
               if(rpath_in(ii:ii).eq.'/') then
                  root=rpath_in(1:ii)
                  rpath=rpath_in(ii+1:)
                  exit
               endif
            enddo
         else if(rpath_in(1:1).eq.'~') then
            if(rpath_in(2:2).eq.'/') then
               root='HOME'
               rpath=rpath_in(3:)
            else
               ! (I don't think this can work: ~<username>)
               root=' '
               rpath = rpath_in
            endif
         else
            call sget_cwd(root)
            rpath = rpath_in    ! path relative to $cwd
         endif
      endif

!!      write(0,*) ' ** root_in = "',trim(root_in),'"'
!!      write(0,*) ' ** rpath_in = "',trim(rpath_in),'"'
!!      write(0,*) ' ** root = "',trim(root),'"'
!!      write(0,*) ' ** rpath = "',trim(rpath),'"'

      ier=0
      call gmkdir_init(root,start,rpath,next)

!!      write(0,*) ' gmkdir_init: '
!!      write(0,*) '   root = "',trim(root),'"'
!!      write(0,*) '   start = "',trim(start),'"'
!!      write(0,*) '   rpath = "',trim(rpath),'"'
!!      write(0,*) '   next = "',trim(next),'"'

 10   continue
      if((ier.eq.0).and.(next.ne.' ')) then
         call gmkdir_next(start,next,ier)

!!!         write(0,*) ' gmkdir_next: '
!!!         write(0,*) '   start = "',trim(start),'"'
!!!         write(0,*) '   next = "',trim(next),'"'
!!!         write(0,*) '   ier = ',ier

         go to 10
      endif
!
      return
      end
!
!-------------------------------------------------------------
!
      subroutine gmkdir_init(root,start,rpath,next)
!
      implicit NONE
!
      character*(*) root                ! root (as passed)
      character*(*) start               ! starting directory (standardized)
!
      character*(*) rpath               ! relative path (as passed)
      character*(*) next                ! starting rel. path (standardized)
!-------------
      integer str_length
      integer :: iclen,ic,ic1,ic2,istat,iplen,io
!
      character*BUFSIZ ztest
      character*15 cpid
      character*30 tmpfile
!---------------------------------
!
      start=' '
      next=' '
!
#ifdef __UNIX
!
!  root directory...
!
      call sget_pid_str(cpid,iplen)
      tmpfile = 'gmkdir_'//cpid(1:iplen)//'.tmp'
!
      if(root.ne.' ') call ufilnam(root,' ',start)
!
!  relative path from there...
!
      iclen=str_length(rpath)
!
!  does requested directory already exist?
!
      call find_io_unit(io)
!
      ztest = trim(start)//'/'//trim(rpath)//'/'//trim(tmpfile)
      open(unit=io,file=trim(ztest),status='unknown',iostat=istat)
      if(istat.eq.0) then
                                !  open successful; directory must exist.
         close(unit=io,status='delete',iostat=istat)
         start = trim(start)//'/'//trim(rpath)
         next = ' '
         return
      endif
!
!  work up through rpath parents...
!
      ic=iclen
      do
         ic=ic-1
         if(ic.le.0) exit
         if(rpath(ic:ic).ne.'/') cycle

         ztest = trim(start)//'/'//rpath(1:ic)//'/'//trim(tmpfile)
         open(unit=io,file=trim(ztest),status='unknown',iostat=istat)
         if(istat.eq.0) then
                                !  open successful; directory must exist.
            close(unit=io,status='delete',iostat=istat)
            start = trim(start)//'/'//rpath(1:ic)
            next = rpath(ic+1:)
            return
         endif
      enddo
!
!  try start directory
!
      ztest = trim(start)//'/'//trim(tmpfile)
      open(unit=io,file=trim(ztest),status='unknown',iostat=istat)
      if(istat.eq.0) then
                                !  open successful; directory must exist.
         close(unit=io,status='delete',iostat=istat)
         next=rpath
         return
      endif
!
!  this code reached if directory path (start) is not fully in place.
!  (this would be unusual).
!
      ic=0
 10   continue
      ic=ic+1
      if(ic.gt.iclen) go to 1000        ! reached the end
      if((ichar(rpath(ic:ic)).ne.0).and.(rpath(ic:ic).ne.' ').and.
     >   (rpath(ic:ic).ne.'~').and.(rpath(ic:ic).ne.'/')) then
         ic1=ic
      else
         go to 10
      endif
!
      ic=iclen+1
 20   continue
      ic=ic-1
      if(rpath(ic:ic).ne.'/') then
         ic2=ic
      else
         go to 20
      endif
!
      next=rpath(ic1:ic2)
      go to 1000
#endif
!-------------------
#ifdef __VMS
!
!  root directory...
!
      if(root.ne.' ') then
         itest=index(root,'[')
         if(itest.eq.0) itest=index(root,']')
         if(itest.eq.0) itest=index(root,'<')
         if(itest.eq.0) itest=index(root,'>')
         if(itest.ne.0) then
            start=root                  ! use "as is"
            go to 500
         endif
!
!  may be a logical name that has to be translated (e.g. SYS$LOGIN)
!
         iclen=str_length(root)
         ic=0
 110     continue
         ic=ic+1
         if((ichar(root(ic:ic)).ne.0).and.(root(ic:ic).ne.' ')) then
            ic1=ic
         else
            go to 110
         endif
!
         ic2=iclen
         if(root(ic2:ic2).eq.':') ic2=ic2-1
         call mpi_sget_env(root(ic1:ic2),ztest,iertmp)
         if(ztest.eq.' ') then
            start=root                  ! translation failed, use "as is"
         else
            ilz=str_length(ztest)
            if(ilz.gt.2) then
               if(ztest(ilz-1:ilz).eq.'.]') then
                  start=ztest(1:ilz-2)//']' ! rooted name, omit trailing .
               else if(ztest(ilz-1:ilz).eq.'.>') then
                  start=ztest(1:ilz-2)//'>' ! rooted name, omit trailing .
               else
                  start=ztest(1:ilz)    ! use translated value verbatim
               endif
            else
               start=ztest(1:ilz)       ! use translated value verbatim
            endif
         endif
      endif
!
!  relative path subdirectory...
!
 500  continue
      iclen=str_length(rpath)
      ic=0
!
 510  continue
      ic=ic+1
      if(ic.gt.iclen) go to 1000        ! reached the end
      if((ichar(rpath(ic:ic)).ne.0).and.(rpath(ic:ic).ne.' ').and.
     >   (rpath(ic:ic).ne.'[').and.(rpath(ic:ic).ne.'<').and.
     >   (rpath(ic:ic).ne.'.')) then
         ic1=ic
      else
         go to 510
      endif
!
      ic=iclen+1
 520  continue
      ic=ic-1
      if((rpath(ic:ic).ne.']').and.(rpath(ic:ic).ne.'>').and.
     >   (rpath(ic:ic).ne.'.')) then
         ic2=ic
      else
         go to 520
      endif
!
      next=rpath(ic1:ic2)
#endif
 1000 continue
      return
      end
!----------------------------------------------------------------------
!
      subroutine gmkdir_next(start,next,ier)
!
! check the next directory
!
      implicit NONE
!
      character*BUFSIZ start            ! current root
      character*BUFSIZ next             ! current relative path
!
      integer ier
!
!-----------------
#include "fpreproc/byte_declare.h"
!
      integer jsystem,str_length,cmkdir
      integer :: ils,ilx,inext,ild,iln,istat
!
      character*1 zdelim,zterm
      character*BUFSIZ dirfile,newdir
!
      character*512 zcmd
      BYTE_DECLARE, dimension(:), allocatable :: cpath
!-----------------
!
#ifdef __UNIX
      zdelim='/'
#else
      zdelim='.'
#endif
!
      ils=str_length(start)
      ilx=str_length(next)
!
      inext=index(next,zdelim)-1
      if(inext.le.0) then
         inext=ilx
      endif
!
      if(start.ne.' ') then
#ifdef __UNIX
         dirfile=start(1:ils)//next(1:inext)
#else
         dirfile=start(1:ils)//next(1:inext)//'.DIR'
         zterm=start(ils:ils)
         newdir=start(1:ils-1)//zdelim//next(1:inext)//zterm
#endif
      else
!
!  blank root
!
#ifdef __UNIX
         dirfile=next(1:inext)
#else
         dirfile=next(1:inext)//'.DIR'
         zterm=']'
         newdir='['//zdelim//next(1:inext)//zterm
#endif
      endif
      ild=str_length(dirfile)
#ifdef __UNIX
      newdir=dirfile(1:ild)//'/'
      ild=ild+1
#endif
      iln=str_length(newdir)
!-------------------------------------------
! make the shell command
!   ...mod DMC: call cmkdir(<path>)
#ifdef __UNIX
      if(allocated(cpath)) deallocate(cpath)
      allocate(cpath(ild+1))
      call cstring(newdir(1:ild),cpath,'2C')
      istat = cmkdir(cpath)
#else
      zcmd='if f$search("'//dirfile(1:ild)//'").eqs."" then $ '//
     >   'create/directory '//newdir(1:ild)
      istat = jsystem(zcmd)   ! VMS only
#endif
!
#ifdef __UNIX
      ier=istat
#else
      if(mod(istat,2).ne.1) then
         ier=istat
      else
         ier=0
      endif
#endif
!--------------------------------------------
!  new root and rel. path
!
      start=newdir
      if(inext.lt.ilx) then
         next=next(inext+2:ilx)
      else
         next=' '
      endif
!
      return
      end
 
