#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_add_node(node, usage, stat)
      implicit none
      character*(*) node, usage
      integer stat
      mds_add_node = -1
      return
      end
#else
      integer function mds_add_node(node, usage, stat)
C Create a MDSplus Node
C Check first if node already exists
      implicit none
      character*(*) node
      character*(*) usage
      integer       stat
 
      integer   mds__long, mds__cstring
 
      integer   IDTYPE_CSTRING
      parameter(IDTYPE_CSTRING=14)
 
      integer  descr, mdsvalue
      integer  nid
      integer  cmd_sts, cmd_sts2, retl
      integer  i,       add
      integer dsc, dsc2, ndsc, cdsc, nl, ulen
      character*4097 strarg
      character*1    bcks
      volatile nid, cmd_sts, cmd_sts2
C---------------------------------------------------------------------------
      bcks=char(92)                 ! x5c = '\'
 
C  /* Check if node already exists */
 
      nid=0
      ndsc = descr(mds__long(), nid, 0)
      nl=index(node,' ')-1
      if (nl .le. 0) nl=len(node)
      cdsc=descr(IDTYPE_CSTRING,node,0,nl)
      stat = MdsValue("getnci($,'NID_NUMBER')",cdsc,ndsc,0,retl)
 
      if (mod(stat,2) .eq. 1) then
         mds_add_node=nid
         return
      endif
 
      add=1
      i=0
      do while ( i .le. 1 .and. add .eq. 1)
         i=i+1
C Add Node
         dsc = descr(mds__long(), cmd_sts, 0)
         if (node(1:1).eq.bcks) then
            strarg = 'tcl("add node ' // bcks // node(1:nl) //
     >           '/usage=' // usage // '")' //CHAR(0)
         else
            strarg = 'tcl("add node ' // node(1:nl) //
     >           '/usage=' // usage // '")' //CHAR(0)
         endif
         stat = MdsValue(strarg, dsc, 0, retl)
         if (mod(stat,2) .ne. 1) then
            mds_add_node = -1
            print*,'?MdsValue(tcl("add node error',node(1:nl),stat
            return
         else
C Check status of tcl add
            if (mod(cmd_sts,2) .eq. 1) then
               mds_add_node = 1
               add = 0
C           not open for edit: open
            else if (cmd_sts .eq.  265388274 ) then
               print*,"  -- NOT OPEN FOR EDIT"
               cmd_sts2 = 0
               dsc2 = descr(mds__long(), cmd_sts2, 0)
               stat = MdsValue("TreeOpenEdit($EXPT,$SHOT)",
     >              dsc2,0,retl)
C               print*,"TreeOpenEdit cmd_sts.. ",cmd_sts2,retl,stat
C check status of MdsValue
               if (mod(stat,2) .ne. 1) then
                  print*,"?mds_add_node -- TreeOpenEdit Error ",stat
                  mds_add_node = -1
                  return
C check status of TreeOpenEdit
               else if ( mod(cmd_sts2,2) .ne. 1) then
                  print*,"?mds_add_node -- TreeOpenEdit cmd_sts.. ",
     >                   cmd_sts2
C     no Privilege
                  if (cmd_sts2 .eq. 265392034) then
                     print*,"?mds_add_node -- No Priv to edit",cmd_sts2
                  else
                     print*,"?mds_add_node --TreeOpenEdit Error ",
     >                      cmd_sts2
                  endif
                  mds_add_node = -1
                  return
               endif ! end TreeOpenEdit
            else if (cmd_sts .eq. 265388168) then
               print*,'?mds_add_node: already there ',node(1:nl)
               mds_add_node=1
               return
            else
               print*,'?mds_add_node error ',stat,cmd_sts
               mds_add_node = 0
               return
            endif ! end not open for edit
         endif
      enddo
      return
      end
#endif
