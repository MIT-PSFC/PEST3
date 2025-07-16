 subroutine fnmdf( cfnm, cext )
!
!  complete filename in cfnm with cext, if not already given -
!
!  eg
!  cfnm = 'work/54321C23'
!  call fnmdf( cfnm, 'TR.DAT' )
!
!  would return cfnm = 'work/54321C23/54321C23TR.DAT'
!
!  Mods
!    29Mar2009   jim.conboy@jet.uk
!                New at build date 'Jan 24 04:30:36', JET vn 44
!
!--------------------------------------------------------------------
!
!..Args :-
!
 character(len=*), intent(inout)  ::  cfnm       ! (part of ? ) filename
 character(len=*), intent(in)     ::  cext       ! requested extension
!
 integer                          ::  lfnm &     ! length of cfnm
                                     ,lfnt &     ! ..  w/out trailing spaces
                                     ,lastdot  & ! ix of last dot
                                     ,lastslash &! ix of last slash
                                     ,ncdir    & ! nchar in directory name
                                     ,ncap       ! nchar to append
!
 character(len=1)                 ::  cdot='.' &
                                    , cbl =' ' &
                                    , cslash='/' !  filename furniture

!--------------------------------------------------------------------
!
 lfnm      =  len(cfnm)
 lfnt      =  len_trim(cfnm)
 lastdot   =  index( cfnm, cdot, .true. )
 if( cfnm(lfnt:lfnt) == cslash )   then
     cfnm(lfnt:lfnt) = cbl
     lfnt = lfnt - 1
 endif
 lastslash =  index( cfnm, cslash, .true. )
!
!d-  print *,'fnmdf :',lastslash,lastdot,lfnt,lfnm
!
 if ( lastdot > 1 .and. lastdot > lastslash )  return
!
 ncdir = lfnt  - lastslash
 ncap  = ncdir + 1 + len_trim(cext)
 if( lfnt + ncap > lfnm )                       then
    print *,'portlib%fnmdf -E- buffer too short ',lfnm
    ix = min(lfnm,lfnt+1)
    cfnm(ix:ix) = '*'
                                                return 
                                                endif
!
 cfnm(lfnt+1:lfnt+1)  = cslash
 cfnm(lfnt+2:)        = cfnm(lastslash+1:lfnt)
 cfnm(lfnt+2+ncdir:)  = trim(cext)
!
!d-  print *,'fnmdf : "',trim(cfnm),'"'
 return
 end
