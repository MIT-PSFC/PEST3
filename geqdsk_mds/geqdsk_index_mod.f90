module geqdsk_index_mod

  implicit NONE
  save

  !  compound data type contains all information in an EFIT_INDEX file
  !  sample contents:
  !    "ntimes=4"
  !    "path='/pub/outgoing/pshare'"
  !    "time=0.15   filename=g008500.00150"
  !    "time=0.20   filename=g008500.00200"
  !    "time=0.25   filename=g008500.00250"
  !    "time=0.294  filename=g008500.00294"

  !  file contains list of EFIT G-eqdsk files, each associated with a time
  !  value.

  type geqdsk_index
     character*128 :: path      ! path to index file itself
     character*128 :: ifname    ! index filename (not including directory path)

     integer :: ntimes          ! no. of time points "ntimes"
     integer :: nguard          ! initialization state (-999 = OK)
     integer :: ntime2          ! no. of time points (file records)

     real*8, dimension(:), allocatable :: times         ! time vector
     character*20, dimension(:), allocatable :: ctimes  ! time vector (char)
     character*128, dimension(:), allocatable :: efiles ! associated filenames

  end type geqdsk_index

  !---------------------------------------------------------
  ! a single instance is provided here...

  logical :: geqdsk_index1_init = .FALSE.
  type (geqdsk_index) :: geqdsk_index1

  !---------------------------------------------------------

  private :: fn_split,parse,rshift,g1val  ! private methods...

  contains

    subroutine geqdsk_index_preset(geqi)

      !  put geqdsk_index object into defined but empty state

      type(geqdsk_index), intent(inout) :: geqi ! geqdsk_index object

      geqi%path=' '
      geqi%ifname=' '
      geqi%ntimes=0
      geqi%ntime2=0
      geqi%nguard=-999
      if(allocated(geqi%times)) deallocate(geqi%times)
      if(allocated(geqi%ctimes)) deallocate(geqi%ctimes)
      if(allocated(geqi%efiles)) deallocate(geqi%efiles)

    end subroutine geqdsk_index_preset

    subroutine geqdsk_index_clear(geqi)

      !  clear contents of geqdsk_index object

      type(geqdsk_index), intent(inout) :: geqi ! geqdsk_index object

      if(geqi%nguard.eq.-999) then
         geqi%path=' '
         geqi%ifname=' '
         geqi%ntimes=0
         geqi%ntime2=0
         if(allocated(geqi%times)) deallocate(geqi%times)
         if(allocated(geqi%ctimes)) deallocate(geqi%ctimes)
         if(allocated(geqi%efiles)) deallocate(geqi%efiles)
      else
         call geqdsk_index_preset(geqi)
      endif
    end subroutine geqdsk_index_clear

    subroutine geqdsk_index_msg(geqi,ierr,ilunmsg)

      type(geqdsk_index), intent(inout) :: geqi ! geqdsk_index object
      integer, intent(in) :: ierr     ! error code
      integer, intent(in) :: ilunmsg  ! logical unit for message

      !  CONTENTS of geqi are cleared at end of this routine ***

      if(ierr.eq.1) then
         write(ilunmsg,*) ' ?geqdsk_index_mod: blank filename.'
      else if(ierr.eq.2) then
         write(ilunmsg,*) ' ?geqdsk_index_mod: file open failure.'
      else if(ierr.eq.3) then
         write(ilunmsg,*) ' ?geqdsk_index_mod: file record parse failure.'
      else if(ierr.eq.4) then
         write(ilunmsg,*) ' ?geqdsk_index_mod: time value before ntimes defined.'
      else if(ierr.eq.5) then
         write(ilunmsg,*) ' ?geqdsk_index_mod: time records out of order.'
      else if(ierr.eq.6) then
         write(ilunmsg,*) ' ?geqdsk_index_mod: no. of times does not match ntimes.'
         write(ilunmsg,*) '  "ntimes" in file: ',geqi%ntimes
         write(ilunmsg,*) '  From counting EFITS: ',geqi%ntime2
      else
         write(ilunmsg,*) ' ?geqdsk_index_mod: error (unspecified).'
      endif

      call geqdsk_index_clear(geqi)

    end subroutine geqdsk_index_msg

    subroutine geqdsk_index_init(ilun,filename,geqi,ierr)

      !  initialize a geqdsk_index object from an ascii file
      !  if there is a leading string "EFIT_INDEX:" ignore it.

      integer, intent(in) :: ilun               ! fortran i/o unit no. to use.
      character*(*), intent(in) :: filename     ! file containing data for geqi
      type(geqdsk_index), intent(inout) :: geqi ! geqdsk_index object
      integer, intent(out) :: ierr              ! completion code, 0=OK

      !-----------------------
      character*128 ifpath
      character*128 ifname
      character*128 zbuf,zval1,zval2
      character*10 :: cival
      character*20 :: cgval
      integer i,j,ilen,istat,itimes,itype
      !-----------------------
      !  get filename...

      call fn_split(filename,ifname,ifpath,ierr)  ! split into path & filename
      if(ierr.ne.0) return

      ! does it match object current object contents?
      if((ifpath.eq.geqi%path).and.(ifname.eq.geqi%ifname)) then
         return  ! object already initialized from file
      endif

      ! no: initialize object from file.  First, open the file...
      call geqdsk_index_clear(geqi)  ! and clear old object contents

      open(unit=ilun,file=trim(ifpath)//trim(ifname),status='old',iostat=istat)
      if(istat.ne.0) then
         close(unit=ilun,iostat=istat)
         ierr=2
         return  ! open failure.
      endif

      ! opened OK; read contents

      itimes=-1

      do
         read(ilun,'(A)',end=888) zbuf
         call parse(zbuf,itype,zval1,zval2)

         if(itype.eq.0) then
            ierr=3
            exit

         else if(itype.eq.1) then
            ! ntimes=<no. of times> record
            if(zval2.ne.' ') then
               ierr=3
               exit
            endif
            call rshift(zval1,cival)
            read(cival,'(i10)',iostat=istat) geqi%ntimes
            if(istat.ne.0) then
               ierr=3
               exit
            endif
            itimes=0  ! OK can read times now...
            allocate(geqi%times(geqi%ntimes))
            allocate(geqi%ctimes(geqi%ntimes))
            allocate(geqi%efiles(geqi%ntimes))

         else if(itype.eq.2) then
            ! path=<file-path> record  -- ignore
            continue

         else if(itype.eq.3) then
            ! time=<time-val> filename=<filename>
            if(itimes.lt.0) then
               ierr=4
               exit
            endif
            itimes=itimes+1

            if(itimes.gt.geqi%ntimes) cycle

            geqi%ctimes(itimes)=zval1
            geqi%efiles(itimes)=zval2
            call rshift(zval1,cgval)
            read(cgval,'(g20.0)',iostat=istat) geqi%times(itimes)
            if(istat.ne.0) then
               ierr=3
               exit
            endif
            if(itimes.gt.1) then
               if(geqi%times(itimes-1).ge.geqi%times(itimes)) then
                  ierr=5
                  exit
               endif
            endif
         endif
      enddo
888   continue
      close(unit=ilun)
      geqi%ntime2 = itimes
      if(ierr.eq.0) then
         if(itimes.ne.geqi%ntimes) then
            ierr=6
         endif
      endif

      geqi%path=ifpath
      geqi%ifname=ifname
    end subroutine geqdsk_index_init

    subroutine fn_split(filename,ifname,ifpath,ierr)

      ! split into path & filename

      character*(*), intent(in) :: filename
      character*(*), intent(out) :: ifname
      character*(*), intent(out) :: ifpath
      integer, intent(out) :: ierr

      integer str_length,i,ilen,istart
      character*10 ztest

      ierr=0
      ifname=' '
      ifpath=' '
      ilen=str_length(filename)
      if(ilen.eq.0) then
         ierr=1  ! non-blank input expected
         return
      endif

      ! exclude leading "EFIT_INDEX:" if present...

      istart=1
      i=index(filename,':')
      if(i.gt.0) then
         ztest=filename(1:10)
         call uupper(ztest)
         if(ztest.eq.'EFIT_INDEX') then
            istart=i+1
         endif
      endif

      if(index(filename,'/').eq.0) then
         ifname=filename(istart:)
         ifpath=' '
      else
         do i=ilen,istart,-1
            if(filename(i:i).eq.'/') then
               if(i.eq.ilen) then
                  ierr=1  ! trailing "/" not expected...
               else
                  ifname=filename(i+1:ilen)
                  ifpath=filename(istart:i)
               endif
               exit
            endif
         enddo
      endif

    end subroutine fn_split

    subroutine rshift(zinval,zrsval)
      character*(*), intent(in) :: zinval  ! input value
      character*(*), intent(out) :: zrsval ! value right shifted

      integer :: str_length,irlen,iilen

      zrsval=' '
      iilen=str_length(zinval)
      if(iilen.eq.0) return

      irlen=len(zrsval)
      if(iilen.gt.irlen) return

      zrsval(irlen-iilen+1:irlen)=zinval(1:iilen)
    end subroutine rshift

    subroutine parse(zbuf,itype,zval1,zval2)
      character*(*), intent(in) :: zbuf    ! line read from file
      integer, intent(out) :: itype
      character*(*), intent(out) :: zval1,zval2

      !  parse line:
      !   type 1:   ntimes=<no. of times>
      !   type 2:   path=<path-string>
      !   type 3:   time=<time-value>  filename=<filename-string>

      character*10 :: ztest
      character*128 :: part2
      integer :: ieqs,i,ii,ifin

      zval1=' '
      zval2=' '
      itype=0

      ieqs=index(zbuf,'=')
      if(ieqs.eq.0) return

      ztest=' '
      ii=0
      do i=1,ieqs-1
         if(zbuf(i:i).ne.' ') then
            ii=ii+1
            if(ii.gt.len(ztest)) exit
            ztest(ii:ii)=zbuf(i:i)
         endif
      enddo

      call uupper(ztest)
      if(ztest.eq.'NTIMES') then
         itype=1
         call g1val(zbuf,ieqs+1,zval1,ifin)
      else if(ztest.eq.'PATH') then
         itype=2
         call g1val(zbuf,ieqs+1,zval1,ifin)
      else if(ztest.eq.'TIME') then
         itype=3
         call g1val(zbuf,ieqs+1,zval1,ifin)
         ieqs=index(zbuf(ifin:),'=')
         if(ieqs.eq.0) then
            itype=0
            zval1=' '
         else
            ieqs=ieqs+ifin-1
            ii=0
            do i=ifin,ieqs-1
               if(zbuf(i:i).ne.' ') then
                  ii=ii+1
                  if(ii.gt.len(ztest)) exit
                  ztest(ii:ii)=zbuf(i:i)
               endif
            enddo
            call uupper(ztest)
            if(ztest.ne.'FILENAME') then
               itype=0
               zval1=' '
            else
               call g1val(zbuf,ieqs+1,zval2,ifin)
            endif
         endif
      else
         itype=0
      endif

    end subroutine parse

    subroutine g1val(zbuf,istart,zval,ifin)
      character*(*), intent(in) :: zbuf  ! input buffer
      integer, intent(in) :: istart      ! scan start
      character*(*), intent(out) :: zval ! value extracted
      integer, intent(out) :: ifin       ! position beyond end of value

      integer :: ii,jj,str_length,ilen

      ilen=str_length(zbuf)

      zval=' '
      ii=istart-1
      jj=0
      do 
         ii=ii+1
         if(ii.gt.ilen) exit

         if((zbuf(ii:ii).ne.' ').and.(zbuf(ii:ii).ne."'").and. &
              (zbuf(ii:ii).ne.'"')) then
            jj=jj+1
            zval(jj:jj)=zbuf(ii:ii)
         else
            if((zbuf(ii:ii).eq.' ').and.(jj.gt.0)) exit
         endif
      enddo
      ifin=ii

    end subroutine g1val

end module geqdsk_index_mod
