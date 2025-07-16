#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_remove_tag(tag)
      implicit none
      character*(*) tag
      mds_remove_tag = 0
      return
      end
#else
c-----------------------------------------------------------------------
c mds_remove_tag
c remove mdsplus tag
c
c input
c   tag   -- character
c returns
c   1          -- on success
c   odd number -- on error
c
c-----------------------------------------------------------------------
      integer function mds_remove_tag(tag)
      implicit none
      character*(*) tag
 
      integer     status, cmd_status, cmd_stat2
      integer     dsc_cmd_status, dsc_cmd_stat2
      integer     retlen
 
      integer istat, icmdstat
      character*4097 strarg
 
      integer descr, MdsValue, mds__long
      integer mdsdclfdo
      character*1    bcks
 
      bcks=char(92)                 ! x5c = '\'
 
      status = 0
 
      dsc_cmd_status = descr(mds__long(), cmd_status, 0)
C Check if open for edit
      strarg = 'GETDBI("OPEN_FOR_EDIT")' // CHAR(0)
      status = MdsValue(strarg, dsc_cmd_status, 0, retlen)
      istat    = mod(status, 2)
      icmdstat = mod(cmd_status, 2)
 
      if (istat.eq.1) then
C If failure, open for edit:
         if (cmd_status.ne.1) then
             dsc_cmd_stat2 = descr(mds__long(), cmd_stat2, 0)
#ifdef __VMS
            strarg = 'TREESHR->TREE$OPEN_TREE_EDIT($EXPT,$SHOT)'
     >      //CHAR(0)
#else
C           strarg = 'tcl("edit $EXPT/shot=$SHOT")'//CHAR(0)
            print *,
     >      ' ?mds_remove_tag - Error: tree not open for edit'
#endif
            status = MdsValue(strarg, dsc_cmd_stat2, 0, retlen)
            istat    = mod(status, 2)
            icmdstat = mod(cmd_stat2, 2)
         endif
         if ((istat.eq.1) .and. (icmdstat.eq.1)) then
               strarg = 'tcl("remove tag ' //
     >                  tag // ' ")' //CHAR(0)
C     print *,strarg(1:80)
            status = MdsValue(strarg, dsc_cmd_status, 0, retlen)
            istat    = mod(status, 2)
            icmdstat = mod(cmd_status, 2)
            if (istat.eq.1) then
               status = cmd_status
            endif
         endif
      endif
      go to 100
 
 100  continue
 
      mds_remove_tag = status
 
      return
      end
 
#endif
