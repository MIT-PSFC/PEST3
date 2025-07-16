!#define _TEST
   module mtrfile
!. mtrfile                           File handling operations
!
!
!  Author                            jim.conboy@ccfe.ac.uk
!  Version                           1.00,  06Jan2012
!  Modifications
!
!  1.00  06Jan2012                   jim.conboy@ccfe.ac.uk
!                                    Written
!----_^---------------------------------=========================================|
!
!
   use mtrcons, only  : cdlm, cdkdlm, cfpath, cspace, csp8, char0
   use logmod
   implicit  none
   private   cfnapp
!
!..GLOBAL
!
   integer, parameter                :: ncfnm = 1024    ! Max chars in file path
   character(len=8)                  :: cbl8
   character(len=ncfnm)              :: cerr          & ! error msg buffer      
                                       ,cbuf            ! working space
!
   contains

  subroutine filename( cfnm, istat, c0, c1, c2, c3, c4 )
  !.. filename                         construct a filename from strings
  !                                    codesys/source/portlib/mtrfile.filename.
  !
  !
  !  Author                            jim.conboy@ccfe.ac.uk
  !  Version                           1.00,  06Jan2012
  !  Modifications
  !
  !  1.00  06Jan2012                   jim.conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  character(len=*), intent(inout)    :: cfnm           ! file name
  integer, intent(out)               :: istat          ! return status, 0=OK
  character(len=*), intent(in)       :: c0, c1         ! path elements
  character(len=*), intent(in),    &
                    optional         :: c2,c3,c4
  !
  !..local
  integer                            :: lfnm         & ! length cfnm
                                       ,l1             ! index, last char in cfnm
  !
  character(len=*),parameter         :: cSr='portlib/mtrfile.filename'
  !----------------------------------::================!.........................|
  !
  lfnm = len( cfnm )
  !
  if( c0 == '+' )                               then
     l1 = len_trim(cfnm) 
     if( cfnm(l1:l1) == cdlm )  l1 = l1-1
                                                else
     l1   = 0
     cfnm = repeat(csp8, len(cfnm)/8)
     call cfnapp( cfnm, lfnm, l1, c0, istat )
     if( istat /= 0 )                         then
         print *,cSr,' Arg 1 istat = ',istat 
                       return ;               endif
                                                endif
  !
  call cfnapp( cfnm, lfnm, l1, c1, istat )
  if( istat /= 0 )                              then 
      print *,cSr,' Arg 2 istat = ',istat 
                       return ;                 endif
  !
  !--       3rd field
  if( .not. present( c2 ))                      then
                                                return 
                                                else
     call cfnapp(  cfnm, lfnm, l1, c2, istat )
     if( istat /= 0 )                        then 
      print *,cSr,' Arg 3 istat = ',istat 
                       return ;              endif
                                                endif
  !
  !--       4th field
  if( .not. present( c3 ))                      then
                                                return 
                                                else
     call cfnapp(  cfnm, lfnm, l1, c3, istat )
     if( istat /= 0 )                        then
      print *,cSr,' Arg 1 istat = ',istat 
                       return ;              endif 
                                                endif
  !
  !--       5th field
  if( .not. present( c4 ))                      then
                                                return 
                                                else
     call cfnapp(  cfnm, lfnm, l1, c4, istat )
     if( istat /= 0 )                        then 
                       return ;              endif
                                                endif
  !
  end subroutine filename

  subroutine cfnapp( cpath, lpath, ix, cdir, istat )
  !.. cfnapp                           append single path element to file path
  !                                    codesys/source/portlib/mtrfile.cfnapp.
  !
  !
  !  Author                            <author>@<place>
  !  Version                           1.00,  00mmm2012
  !  Modifications
  !
  !  1.00  00mmm2012                   <author>@<place>
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  character(len=1024),intent(inout) :: cpath
  integer, intent(in)                :: lpath          ! length of the path string
  integer, intent(inout)             :: ix             ! to last used character
  character(len=*), intent(in)       :: cdir           ! path element to append
  integer, intent(out)               :: istat          ! return status, 0=OK
                                                       ! > lpath :  # chars needed
                                                       ! < lpath :  ptr to invalid character
  !
  !..local
  integer                            :: ldir           ! length of element to add
  !
  character(len=*),parameter         :: cSr='portlib/mtrfile.cfnapp'
  !----------------------------------::================!.........................|
  !
  ldir = len_trim(cdir)
  if( ldir == 0 )                                       return  ! ??
  if( cdir(ldir:ldir) == char(0))   ldir = ldir - 1
  !d print *,cSr,ldir,cdir(:ldir)
  !d call cdump( cSr//cdir, cdir, istat )
  !
  if( ix+ldir+1 > lpath )                               then
!    call errorlog( cSr,'filename buffer overflow' )
    istat = ix+ldir+1 ;                                 return
                                                        endif
  !
  if( ix /= 0 .and. cdir(1:1) /= cdlm )                 then
     ix = ix + 1
     cpath(ix:ix) = cdlm
                                                        endif
  !
  istat = verify( cdir(:ldir), cfpath )
  if( istat /= 0 )                                      then
     print *,cSr,' -E- ',istat
     write(cerr,'(a,z2,2a)')  'Invalid char  = (Z'      &
                            , ichar(cdir(istat:istat))  &
                            ,') follows ',cdir(1:min(1000,istat-1))
!     call warnlog( cSr, cerr )
     istat = istat + ix                                 ! address of char in path
                                                        endif
  !
  cpath(ix+1:ix+ldir)  = cdir(:ldir)
  ix = ix + ldir
  !
  if( cpath(ix:ix) == cdlm )                            then
     cpath(ix:ix)  = cspace
     ix = ix - 1
                                                        endif
  !
  end subroutine cfnapp

  function get_cwd( copt )
  !.. get_cwd                             Return cwd 
  !                                       codesys/source/portlib/mtrfile.get_cwd.
  !
  !
  !     Author                            <author>@<place>
  !     Version                           1.00,  17Jan2012
  !     Modifications
  !
  !     1.00  17Jan2012                   <author>@<place>
  !                                       Written
  !----_^---------------------------------=========================================|
  !
      implicit  none
  !
  !..GLOBAL
  !--   variable                       :: name           ! description
  !
  !..Arguments
  !
      character(len=ncfnm)           :: get_cwd        !
      character(len=*), optional     :: copt           ! 'full|last' 
  !
  !..local
      integer                        :: ildisk        &
                                       ,ilcwd         &
                                       ,ixdlm         &
                                       ,ixcwd
  !
      character(len=*),parameter     :: cFn='portlib/mtrfile.get_cwd'
  !----_^------------------------------::================!.........................|
  !
      cbuf = repeat( '        ', len(cbuf)/8 )
  !
  !   get current working directory (VMS or UNIX)
  !
      ixcwd = 33
      call showdefl( cbuf(1:ixcwd-1), ildisk     &
                    ,cbuf(ixcwd:   ), ilcwd   )
      if( cbuf(1:1) /= cspace )                          then
         ildisk = len_trim( cbuf(:ixcwd-1) )
         if( cbuf(ildisk:ildisk) /= cdkdlm )           then
            ildisk = ildisk + 1
            cbuf(ildisk:ildisk) = cdkdlm
                                                       endif
            cbuf(:ixcwd-1) =  adjustr( cbuf(:ixcwd-1))
            ixcwd = ixcwd - ildisk
                                                          endif
   !
      if( present( copt ))                                then
         if( copt == 'last' )                           then
            ixdlm = index( cbuf, cdlm, .true. )
            if( ixdlm /= 0 )     ixcwd = ixdlm 
                                                        endif
                                                          endif
      get_cwd = trim( cbuf(ixcwd:))

      end function get_cwd

  subroutine ioerr( lunerr, ioerrno, file, iolun, status  )
  !.. ioerr                            Diagnostics for IO error
  !                                    codesys/source/portlib/ioerr.
  !
  !
  !  Author                            Jim.Conboy@ccfe.ac.uk
  !  Version                           1.20  12Jan2012
  !  Modifications
  !
  !  1.20  12Jan2012                   Jim.Conboy@ccfe.ac.uk
  !                                    included in mtrfile
  !                                    add cStatus arg
  !
  !  1.10  11Nov2011                   Jim.Conboy@ccfe.ac.uk
  !                                    use logmod routines, if lunerr < 0
  !
  !  1.00  31Oct2011                   Jim.Conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  use logmod
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  integer, intent(in)                :: lunerr        &! lun to write error info
                                       ,ioerrno        ! io error code
  character(len=*), intent(in), &
                         optional    :: file           ! file attempted
  integer,   intent(in), optional    :: iolun          ! Lun for which error occured
  character(len=*), intent(in)     &
                  , optional         :: status         ! requested status
  !
  !..local
  character(len=1024)                :: cbuf           ! io buffer
  !
  character(len=*),parameter         :: cSr='portlib/ioerr'
  !----------------------------------::================!.........................|
  !
  if( ioerrno == 0 )                                     return
  !
#ifdef _LAHEY
  call iostat_msg( ioerrno, cbuf )
#else
  cbuf = "No error description available " 
#endif

  if( lunerr > -1 )                                      then 
      write(lunerr,'(a,i8,2a)') &
        ' IO error ',ioerrno,': ',trim(cbuf)
                                                         else
!      call warnlog( trim(cbuf), ioerrno )
                                                         endif
!
  if( present(file) .and. file(1:1) /= '/' )             then
!         relative file name - where are we ?? 
     call sget_cwd( cbuf )
     if( lunerr > -1 )                                then
      write(lunerr,'(2a)') ' cwd : ',trim(cbuf)  
                                                      else
!      call infolog( 'cwd :', trim(cbuf) )
                                                      endif
                                                         endif
!    
  if( present( file ) .and. present( iolun ))            then
    write( cbuf, '(3a,i4)' ) &
          '    for ',trim(file),' on unit ',iolun
  elseif( present( iolun ) )                             then
    write( cbuf, '( a,i4)' ) '    on unit ',iolun
  elseif( present( file ) )                              then
    write( cbuf, '(2a)' ) '    for ',trim(file)
  else
                                                         return
                                                         endif
!
  if( present(status ) )  &
    cbuf = cbuf // ' Status '// trim( status )
  if( lunerr > -1 )                                      then
      write(lunerr,'(2a)') '          ',trim(cbuf)  
                                                         else
      call infolog( trim(cbuf) )
                                                         endif
!
  end subroutine ioerr


  subroutine mkdir( cdir, istat, cdir_full, set )
  !.. mkdir                            Create a directory ( tree); definition may
  !                                    contain environment variables
  !                                    codesys/source/portlib/mkdir.f90
  !
  !                                    fpopen is used to expand the directory name using 
  !                                    echo, so reproducing the (Bourne) shell syntax 
  !                                    
  !  Author                            jim.conboy@ccfe.ac.uk
  !  Version                           1.10  04Jan2012
  !  Modifications
  !
  !  1.10  04Jan2012                   jim.conboy@ccfe.ac.uk
  !                                    Protect against trailing file delimiter
  !
  !  1.00  20May2011                   jim.conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  ! character(len=1024)              :: cbuf
  !
  !..Arguments
  !
  character(len=*), intent(in)       :: cdir           ! the directory name
  integer, intent(out)               :: istat          ! return status, 0=OK
  character(len=*), intent(out),    &
                    optional         :: cdir_full      ! full dir nm
  logical, intent(in), optional      :: set            ! .t to change directory
  !
  !..local
  integer                            :: i
  integer                            :: lx             ! len ( cdir_ex )
  logical                            :: lexist
  character(len=1)                   :: char
  character(len=256)                 :: cdir_ex &
                                       ,cbl256
  !
  character(len=1),parameter         :: cd = '/'
  character(len=*),parameter         :: cSr='portlib/mtrfile.mkdir'
  !----------------------------------::================!.........................|
  !
  print *,cSr,' -I- ',trim(cdir)
  cbl256  = repeat( cbl8, len(cbl256)/8)
  cdir_ex = repeat( cbl8, len(cdir_ex)/8)
  !
  cdir_ex = trim(cdir)
  call parse( cdir_ex, istat )
  ! 
  lx = len_trim(cdir_ex)
  if( cdir_ex(lx:lx) == cd )  lx = lx - 1
  !d print *,cSr,' -I- ',cdir_ex(:lx)
  !
  call fpopen( 'mkdir -vp '//cdir_ex(:lx)//char(0), cbuf )  ! & create directory
  !
  !   Test for success - a test that the file is in fact a directory
  !                      would be better
  !
  inquire( file=cdir_ex(:lx), exist=lexist )                ! .. did we ?
  if( lexist )                                      then
    istat = 0
    if( present( cdir_full )) &
        cdir_full(:len(cdir_full))  = cdir_ex(:lx)//cbl256
    if( present( set ))                         then
        !d print *,cSr,' len ',len_trim(cdir_ex)
        if( set ) call set_cwd( trim(cdir_ex), istat )
                                                endif
                                                    else
    print *,cSr,' -E- could not create ',trim(cdir)
    istat = 1
                                                    endif
  !
  return
  end subroutine mkdir

  subroutine mklink( cfnm, clnm, istat )
  !.. mklink                           Create a ( unix file system ) link
  !                                    codesys/source/portlib/mtrfile.mklink.
  !
  !
  !  Author                            Jim.Conboy@ccfe.ac.uk
  !  Version                           1.00,  13Dec2011
  !  Modifications
  !
  !  1.00  13Dec2011                   Jim.Conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  character(len=*), intent(in)       :: cfnm           ! the target file
  character(len=*), intent(inout)    :: clnm           ! the link name
  integer, intent(out)               :: istat          ! return status, 0=OK
  !
  !..local
  integer                            :: l              ! buffer index 
  !
  character(len=512)                 :: cbuf         & ! buffer ( for  Unix command )
                                       ,cout           ! output ( from Unix command )
  character(len=*),parameter         :: cSr='portlib/mtrfile.mklink'
  !----------------------------------::================!.........................|
  !
  cbuf = 'ln -snf '                  ! -n no deref; -f delete existing =>
                                     ! any existing link will be ( silently ) replaced
  l = len_trim(cbuf)
#ifndef __MPI
  cbuf(l+2:) = trim(cfnm) // cspace // trim(clnm )
#else
  cbuf(l+2:) = trim(cfnm)
  call parse( cbuf(l+1:), istat )
  l = len_trim(cbuf)
  cbuf(l+1:) = cspace//trim(clnm)
  call parse( cbuf(l+2:), istat )
#endif
  !
!  call infolog( cSr, trim(cbuf) )
  cout = repeat( csp8, len(cout)/8 )
  call fpopen( trim(cbuf)//' 2>&1'//char(0), cout )
  l = len_trim(cout )
  if(      len_trim(cout) > 2    &
     .and. cout(1:2) /= cspace//char0  )   then   ! Is zero terminated.. 
!     call errorlog( csr, trim(cout) )
     call cdump( 'mklink.cout', cout(:l), l )
     istat = 1
                                           else
     istat = 0
                                           endif
  !
  end subroutine mklink


  subroutine qset_cwd( creq, istat, ccwd )
  !.. qset_cwd                         Change cwd, unless there already..
  !                                    codesys/source/portlib/mtrfile.qset_cwd.
  !
  !
  !  Author                            <author>@<place>
  !  Version                           1.00,  19Jan2012
  !  Modifications
  !
  !  1.00  19Jan2012                   <author>@<place>
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  character(len=*), intent(in)       :: creq           ! requested directory
  integer, intent(out)               :: istat          ! return status, 0=OK
  character(len=*), intent(out),    &
                    optional         :: ccwd           ! where we ended up..
  !
  !..local
  integer                            :: ixreq        & !
                                       ,lenreq       & ! length creq ( w/out trailing / )
                                       ,lencwd         !        cwd
  !
  character(len=ncfnm)               :: creqp          ! expanded directory name
  !
  character(len=*),parameter         :: cSr='portlib/mtrfile.qset_cwd'
  !----------------------------------::================!.........................|
  !
  cbuf = get_cwd( 'full' )
  lencwd = len_trim(cbuf) 
  if( cbuf(lencwd:lencwd) == cdlm ) lencwd = lencwd - 1
  !
  creqp  = creq
  call parse( creqp, istat )
  lenreq = len_trim(creqp)
  if( creqp(lenreq:lenreq) == cdlm )  lenreq = lenreq - 1
  ! 
  ixreq = index( cbuf(:lencwd), creqp(:lenreq), .true. )
  !d print *,cSr,lencwd,lenreq,ixreq
  !d print *,cbuf(:lencwd)
  !d print *,creqp(:lenreq)
  if( ixreq + lenreq -1 == lencwd )                     then
  !        Already there
!     call warnlog( cSr, ' already in '//trim(creqp) )
     istat = 1
                                                        else
     call set_cwd( creqp, istat )
     if( istat /= 0 )                               then
       write(cerr,'(a,i6,a)') &
          'set_cwd error ',istat,': '
!       call warnlog( csr, trim(cerr)//trim(creqp) )
                                                    endif
                                                        endif
  !
  cbuf = get_cwd( 'full' )
  lencwd = len_trim(cbuf) 
  if( cbuf(lencwd:lencwd) == cdlm ) lencwd = lencwd - 1  
  ! 
  ixreq = index( cbuf(:lencwd), creqp(:lenreq), .true. )

  !d print *,lenreq,cbuf(:lenreq)
  !d print *,lencwd,cbuf(:lencwd),ixreq

  if( ixreq + lenreq -1 /= lencwd )                        then
!     call errorlog( cSr, ' in '//trim(cbuf) )
     istat = 8
                                                        else
     call infolog( cSr, ' in '//trim(cbuf) )
     istat = 0
                                                        endif
  !
  if( present( ccwd ))                                  then
     ccwd = trim( cbuf )
                                                        endif

  end subroutine qset_cwd



      subroutine set_cwd( dir, ierr )
      character(len=*)           :: dir                 ! directory to cd to...
      integer                    :: ierr                ! error status code, 0 = normal
 
#include "fpreproc/byte_declare.h"
      BYTE_DECLARE dirbuf(512)
      integer, parameter         :: ildmax = 512 
      integer                    :: chdir           &
                                   ,linux_chdir
      integer str_length,ild
      character(len=12)          :: clen
      character(len=*), parameter   :: cSr = 'portlib/mtrfile.set_cwd'
!
      ild=max(1,len_trim(dir))
      !d print *,cSr, ' length ',ild
#ifdef __VMS
      write(6,*) ' ?sset_cwd:  VMS:  cannot change directory.'
      ierr = 1                          ! VMS:  no can do
#endif
#ifdef __UNIX
 
#ifdef __CRAY
      call pxfchdir(dir,0,ierr)
#elif __HP || __IBM || __RS6000
      ierr = chdir(dir(1:ild)//char(0))
#else
      if( ild .lt. ildmax ) then
         call cstring(dir(1:ild),dirbuf,'2C')
         ierr = linux_chdir(dirbuf)
      else
         write(clen,'(I8)')  ild
!         call errorlog( 'portlib/mtrfile.set_ccwd',  ' path too long: '//clen )
         write(6,*) 'set_cwd : too long ',dir(1:ildmax)
         ierr=1
      endif
#endif
 
#endif  /*  __UNIX  */
#ifdef _DEBUG
      if( ierr == 0 )  write(6,*) 'portlib/mtrfile.set_cwd : ',trim(dir)
#endif
 
      return
      end  subroutine set_cwd
 

  subroutine tf_inquire( ilun, cfnm_a, istat_a )
  !.. tf_inquire                       File status enquiry
  !                                    codesys/source/portlib/mtrfile.tf_inquire.
  !
  !
  !  Author                            Jim.Conboy@ccfe.ac.uk
  !  Version                           1.00,  04Jan2012
  !  Modifications
  !
  !  1.00  04Jan2012                   Jim.Conboy@ccfe.ac.uk
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  use logmod
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  integer, intent(in), optional      :: ilun           ! unit for enquire
  character(len=*), intent(in)     &
                  , optional         :: cfnm_a         ! File name 
  integer, intent(out), optional     :: istat_a        ! return status, 0=OK
  !
  !..local
  integer, parameter                 :: ncfnm = 512  & ! len file name buffer
                                       ,lunio = 6      ! for report
  !
  integer                            :: istat        & ! return status
                                       ,iflb         & ! file length ( bytes )
                                       ,lbuf         & ! ptr into cbuf
                                       ,lfnm           ! len(cfnm)
  !
  logical                            :: lex          & ! .T if file exists
                                       ,lnm          & !            has a name
                                       ,lopen          !            is opened
  !
  character(len=20)                  :: caccs        & ! access method
                                       ,cact         & ! action { R | W | RW }
                                       ,cread        & ! read status  {yes|no|unknown}
                                       ,cwrite         ! write status   "

  character(len=ncfnm)               :: cfnm         & ! file name
                                       ,cbuf           ! IO buffer
  !
  character(len=*),parameter         :: cSr='portlib/mtrfile.tf_inquire'
  !----------------------------------::================!.........................|
  !
  if( present(ilun))                                    then
     inquire(        unit     = ilun   &
                    ,iostat   = istat  &
                    ,access   = caccs  &
                    ,action   = cact   &
                    ,exist    = lex    &
!                    ,flen     = iflb   &
                    ,name     = cfnm   & 
                    ,named    = lnm    &
                    ,opened   = lopen  &
                    ,read     = cread  &
                    ,write    = cwrite &
                                         )
     !
     if( istat /= 0 )                                then
        write(cbuf,*) ' iostat = ',istat  &
                     ,' for unit ',ilun,' :'
                                                     else
        if( .not. lnm )  cfnm = '<unnamed file>'
        lfnm  = min( len_trim(cfnm), ncfnm-12-9-6 ) 
        write(cbuf,'(3a,i4)') ' Status for ',cfnm(:lfnm)  &
                     ,' on unit ',ilun
                                                     endif
  !
  else if( present( cfnm_a ))                           then
     inquire(        file     = cfnm_a &
                    ,iostat   = istat  &
                    ,access   = caccs  &
                    ,action   = cact   &
                    ,exist    = lex    &
!                    ,flen     = iflb   &
                    ,name     = cfnm   & 
                    ,named    = lnm    &
                    ,opened   = lopen  &
                    ,read     = cread  &
                    ,write    = cwrite &
                                         )
     !
     if( istat /= 0 )                                then
        write(cbuf,*) ' iostat = ',istat  &
                     ,' for file ',trim(cfnm),' :'
                                                     else
        write(cbuf,*) ' Status for ',trim(cfnm)        
                                                     endif
                                                        else
     if( present( istat_a))    istat_a = -1
                                                        return
                                                        endif
  !
  if( present( istat_a))    istat_a = istat
  !
  if( istat /= 0 )                                      then
        lbuf = len_trim(cbuf)
#ifdef _LAHEY
        call iostat_msg( istat, cbuf(lbuf+1:) )
#else
        cbuf(lbuf+1:) = "No error description available " 
#endif
!        call errorlog( cSr, trim(cbuf))
                                                        return
                                                        endif
  !
  write(lunio,*) trim(cbuf)
  if( .not. lex ) write(lunio,*) '    File does not exist ! '
  write(lunio,100 ) trim(caccs), trim(cact), iflb, lopen &
                   ,trim( cread), trim(cwrite)
  !
  100 format( 4x,'Access ',t24,a  &
             /4x,'Action ',t24,a  &
             /4x,'Length ',t24,i9,' bytes' &
             /4x,'Open ? ',t24,l  &
             /4x,'Read ? ',t24,a, t36,'Write ? ',t48,a &
            )
  !
  end subroutine tf_inquire

   end module mtrfile

#ifdef _TEST

  program test_trfile
  !.. test_trfile                      As on the tin
  !                                    codesys/source/portlib/mtrfile.test_trfile.
  !
  !  Usage
  ! >cpp --traditional =1/portlib/mtrfile.f90 ./mtrfile.f90
  ! >lf95 --nfix mtrfile.f90 -o ttrfile -I $CS/mod -l:$CS/lib/portlib.a
  !
  !  Test VMS function ( on Unix )
  ! >cpp --traditional -D_VMS_TEST =1/portlib/showdefl.for ./showdefl.for
  ! >lf95 -c showdefl.for
  ! >lf95 --nfix mtrfile.f90 showdefl.o -o ttrfile -I $CS/mod -l:$CS/lib/portlib.a
  !
  !
  !  Author                            <author>@<place>
  !  Version                           1.00,  17Jan2012
  !  Modifications
  !
  !  1.00  17Jan2012                   <author>@<place>
  !                                    Written
  !----_^---------------------------------=========================================|
  !
  use mtrfile
  use logmod
  !
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..local
  integer                            :: i            & ! loop index 
                                       ,istat
  character(len=256)                 :: ccwd
  !
  character(len=*),parameter         :: cPr='portlib/mtrfile.test_trfile'
  !----------------------------------::================!.........................|
  !
  call openlog( 'testtrfile.log' )
  call setloglevel( 'inf' )
  !
  ccwd = get_cwd( 'full' )
  print *,cPr,' -I- "',trim(ccwd),'"'
  !
  ccwd = get_cwd( 'last' )
  print *,cPr,' -I- "',trim(ccwd),'"'
  !
  ccwd = 'TranspTemp'
  call mklink( '/common/transp_shared/Code/transp/JET_56/tmp', ccwd, istat )
  !
  call mkdir( 'Test', istat )
  !
  call qset_cwd( 'Test', istat, ccwd  )
  print *,cPr,' in ',trim(ccwd )
  !
  call qset_cwd( 'Test', istat, ccwd  )
  print *,cPr,' in ',trim(ccwd )
  !
  end program test_trfile
#endif
