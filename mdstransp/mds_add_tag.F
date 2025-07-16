#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_add_tag(node,tag)
      implicit none
      character*(*) node, tag
      mds_add_tag = 0
      return
      end
#else
c-----------------------------------------------------------------------
c mds_add_tag
c add mdsplus tag
c
c input
c   node   -- character
c   tag   -- character
c returns
c   1          -- on success
c   odd number -- on error
c
c-----------------------------------------------------------------------
      integer function mds_add_tag(node, tag)
      implicit none
      character*(*) node, tag
 
      integer     status, cmd_status, cmd_stat2
      integer     dsc_cmd_status, dsc_cmd_stat2
      integer     retlen
 
      integer istat, icmdstat
      character*4097 strarg
 
      integer mds__leniftrim
      integer descr, MdsValue, mds__long, mds__cstring
      integer mdsdclfdo
      integer ndsc, cdsc, nid, lt
      character*1    bcks
 
      bcks=char(92)                 ! x5c = '\'
 
      status = 0
 
C Check if tag already exists
      nid=0
      ndsc = descr(mds__long(), nid, 0)
      lt=max(1, mds__leniftrim(tag))+1
      call mds__copywnull(strarg, bcks//tag)
      cdsc = descr(mds__cstring(), strarg, 0, lt)
      status = MdsValue("getnci($,'NID_NUMBER')",cdsc,ndsc,0,retlen)
C      print*,"getnci",status,nid,retlen
 
      if (mod(status,2) .eq. 1) then
         mds_add_tag=1
         return
      endif
 
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
     >      ' ?mds_add_tag - Error: tree not open for edit'
#endif
            status = MdsValue(strarg, dsc_cmd_stat2, 0, retlen)
            istat    = mod(status, 2)
            icmdstat = mod(cmd_stat2, 2)
         endif
         if ((istat.eq.1) .and. (icmdstat.eq.1)) then
               strarg = 'tcl("add tag ' //node // ' ' //
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
 
      mds_add_tag = status
 
      return
      end
 
#endif
