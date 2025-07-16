      subroutine mdserr(lun, zstr, stat)
 
      implicit none
 
      character*(*) zstr        ! context message string
      integer stat              ! mdsplus status code
      integer lun               ! lun
 
      integer MDS_ERRSTR, str_length
C
      integer il, istat
      character*128 mds_msg
C
      write(lun,*) zstr,stat,':'
      istat=mds_errstr(stat,mds_msg)
      il=str_length(mds_msg)
      write(lun,*) mds_msg(1:il)
C
      return
      end
