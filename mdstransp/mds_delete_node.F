#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_delete_node(node)
      implicit none
      character*(*) node
      mds_delete_node = 0
      return
      end
#else
c-----------------------------------------------------------------------
c mds_delete_node
c
c input
c   node       -- character
c returns
c   1          -- on success
c   odd number -- on error
c
c-----------------------------------------------------------------------
      integer function mds_delete_node(node)
      implicit none
      character*(*) node
 
      integer       stat
 
      integer   mds__long, mds__cstring
 
      integer   IDTYPE_CSTRING
      parameter(IDTYPE_CSTRING=14)
 
      integer  descr, mdsvalue
      integer  nid,   cmd_sts
      integer  ndsc, cdsc, nl, retl
      character*4097 strarg

      volatile nid, cmd_sts

! Check if node exists 
 
      nid=0
      ndsc = descr(mds__long(), nid, 0)
      nl=len_trim(node)
      cdsc=descr(IDTYPE_CSTRING,node,0,nl)
      stat = MdsValue("getnci($,'NID_NUMBER')",cdsc,ndsc,0,retl)
! node exists: delete 
      if (mod(stat,2) .eq. 1) then
         ndsc = descr(mds__long(), cmd_sts, 0)
         strarg = 'tcl("delete node ' // node(1:nl) //
     >           '/NOCONFIRM")' //CHAR(0)
         stat = MdsValue(strarg, ndsc, 0, retl)         
         if (mod(stat,2) .ne. 1) then
            mds_delete_node = -1
            print*,'?mds_delete_node -E- MdsValue staus= ',stat
            print*,'     command=',strarg(1:len_trim(strarg)-1)
            return
         else
! Check status of tcl delete
            if (mod(cmd_sts,2) .eq. 1) then
               mds_delete_node = 1
            else
               mds_delete_node = cmd_sts
               print*,'?mds_delete_node -E- failed to delete'//
     &             node(1:nl) // ' status = ',stat
            endif
         endif         
      endif
      return
      end
 
#endif
