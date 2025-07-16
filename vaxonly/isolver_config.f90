!
! Based on a tokamak name and possibly some other parameter, return
! the name of the appropriate Isolver state file.  Usually this state file
! will be found in the $ISOLVERDIR directory.  This subroutine is
! placed here to avoid a dependency on Isolver.  It must be hand maintained.
!
! The config argument can have different formats
!   ''                ->  iso_<tok>.nc, this is the default tokamak state file
!   'name: <name>'    ->  save file will have name iso_<tok>-<name>.nc
!   'shot: <n>'       ->  return a save file appropriate for this shot number, <n> can be empty to use shot=999999999
!   'date: <date>'    ->  return a save file appropriate for the tokamak at this date
!                         <date> -> YY,        two digit year
!                                -> YYYY,      four digit year
!                                -> MmmYYYY,   month and year as in Apr2010
!                                -> ddMmmYYYY, day, month and year  01Ap2010
!   'file: <file>'    -> save file will have the file name <file>, ignore the tokamak
!   '<name>'          -> same as 'name: <name>', iso_<tok>-<name>.nc
!
! The tokamak configurations are defined in $JCODESYSDIR/src/transp/main/jython/isolver/toks
! See also ../isoTokamakConfig.py and the config member in isoTokamakConfig and its subclasses
!
subroutine isolver_config(tok, config, savename, ier)
  implicit none

  character*(*), intent(in)  :: tok       ! name of the tokamak
  character*(*), intent(in)  :: config    ! possibly empty configuration
  character*(*), intent(out) :: savename  ! save file name or a short error message if there was an error
  integer,       intent(out) :: ier       ! nonzero on error

  character(len=len(tok))    :: ltok  ! lower case tokamak
  character(len=len(config)) :: fig   ! left adjusted config
  character*10               :: tag   ! stuff up to the colon in lower case
  character*9                :: cshot ! shot number as character

  integer :: ix       ! index of :
  integer :: iy       ! index after : which is not beyond the string
  integer :: ic       ! tag, 0->name:, 1->shot:, 2->date:
  integer :: ishot    ! shot number
  integer :: dd       ! day
  integer :: mm       ! month
  integer :: year     ! four digit year
  integer :: k        ! temp

  ier = 0

  ltok = adjustl(tok)
  call ulower(ltok)

  fig = adjustl(config)

  if (len_trim(fig)==0) then                
     goto 100                 ! no configuration, use default 
  end if
   
  !
  ! --- generic extraction of config info ---
  !
  ic = 0
  
  dd=1
  mm=1
  year=1

  ix = index(fig,':')

  if (ix>0) then
     ! -- fig has a colon tag --
     tag = fig(1:ix)
     
     call ulower(tag)

     iy = min(len(fig), ix+1)

     if (tag == "name:") then
        ic  = 0
        fig = adjustl(fig(iy:))

     else if (tag == "shot:") then
        ic    = 1
        cshot = adjustl(fig(iy:))

        if (len_trim(cshot)==0) then
           ishot = 999999999          !  'shot: ' is equivalent to last shot
        else
           read(cshot,'(i9)',err=200) ishot
        end if

     else if (tag == "date:") then
        ic  = 2
        fig = adjustl(fig(iy:))

        k = len_trim(fig)
        
        if (k<2) then          ! use default
           goto 100      
  
        else if (k==2) then    ! YY
           mm = 1
           dd = 1
           read(fig(1:2),'(i2)',err=200) year

           if (year<70) then
              year = 2000+year
           else
              year = 1900+year
           end if

        else if (k==4) then    ! YYYY
           mm = 1
           dd = 1
           read(fig(1:4),'(i4)',err=200) year

        else if (k==7) then    ! MmmYYYY
           dd = 1

           mm = getMonth(fig(1:3))
           if (mm==0) goto 200

           read(fig(4:7),'(i4)',err=200) year

        else if (k==9) then    ! ddMmmYYY
           read(fig(1:2),'(i2)',err=200) dd

           mm = getMonth(fig(3:5))
           if (mm==0) goto 200

           read(fig(6:9),'(i4)',err=200) year

        else
           savename = "bad date"
           ier=1 ; return
        end if

     else if (tag == "file:") then
        savename = adjustl(fig(iy:))
        return

     else
        savename = "unknown config tag"
        ier=1 ; return
     end if
  else
     ! -- fig is the name for the name: tag --
     ic = 0
  end if

  ! -- name: <fig> ---
  if (ic==0) then
     call ulower(fig)  ! always lower case, see isoTokamakConfig.py

     savename = "iso_"//trim(ltok)//"-"//trim(fig)//".nc"
     return
  end if

  ! --- tokamak specific ---
  if (ltok=="rga" .and. (ic==1 .or. ic==2)) then  
     !
     ! dummy tokamak for testing and as an example for future additions
     ! go backward in time
     !
     if ((ic==1 .and. ishot>=212005) .or. (ic==2 .and. isAfter(7,4,2010))) then
        savename = "iso_rga.nc"             ! tokamak configuration after last upgrade 07Apr2010 and current configuration
        
     else if ((ic==1 .and. ishot>=127042) .or. (ic==2 .and. isAfter(24,10,2000))) then
        savename = "iso_rga-upgrade1.nc"    ! tokamak configuration after first upgrade 24Oct2000
        
     else
        savename = "iso_rga-orig.nc"        ! tokamak as originally built
     end if
     return
  else if (ltok=="nstx" .and. (ic==1 .or. ic==2)) then  
     !
     ! NSTX
     ! go backward in time
     !
     if ((ic==1 .and. ishot>=142633) .or. (ic==2 .and. isAfter(1,1,2013))) then
        savename = "iso_nstx-upgrade.nc"  ! tokamak configuration after upgrade 31Dec2013 
     else
        savename = "iso_nstx-orig.nc"     ! tokamak as originally built
     end if
     return
  end if
  
  goto 100  ! use default


  ! -- default --
100 continue
  if (ltok=="nstx") then
     savename = "iso_"//trim(ltok)//"-orig.nc"     ! NSTX default
  else if (ltok=="fnsf") then
     savename = "iso_"//trim(ltok)//"-r16v18.nc"   ! FNSF default
  else
     savename = "iso_"//trim(ltok)//".nc"
  end if
  return

  ! -- error during read --
200 continue
  ier = 1
  savename = "config read error"
  return
  
contains
  !
  ! ----- getMonth -----
  ! convert the three letter month in the format Ccc such as Apr to the integer month
  !
  integer function getMonth(ccc)
    character*3, intent(in) :: ccc  ! three letter month
    
    character*3  :: c
    character*36 :: mmm="janfebmaraprmayjunjulaugsepoctnovdec"
    
    integer :: i

    c = ccc
    call ulower(c)

    i = index(mmm,c)
    if (i>0) then
       getMonth = (i-1)/3 + 1
    else
       getMonth = 0
    end if
  end function getMonth
  
  !
  ! ----- isAfter -----
  ! Return true if the configure date dd,mm,year is after or equal to the
  ! argument epoch edd, emm, eyear
  !
  logical function isAfter(edd, emm, eyear)
     integer, intent(in) :: edd   ! day of month at epoch
     integer, intent(in) :: emm   ! month at epoch
     integer, intent(in) :: eyear ! yeat at epoch
     
     if (year>eyear)  then
        isAfter = .true.  ; return
     else if (year<eyear) then
        isAfter = .false. ; return
     end if
     
     if (mm>emm) then
        isAfter = .true.  ; return
     else if (mm<emm) then
        isAfter = .false. ; return
     end if
     
     if (dd>=edd) then
        isAfter = .true.
     else
        isAfter = .false.
     end if
  end function isAfter
end subroutine isolver_config
