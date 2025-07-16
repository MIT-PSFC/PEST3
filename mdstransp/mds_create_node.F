#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_create_node(node, usage)
      implicit none
      character*(*) node, usage
      mds_create_node = 0
      return
      end
#else
c-----------------------------------------------------------------------
c mds_create_node
c add mdsplus node
c
c input
c   node   -- character
c   usage  -- integer
c returns
c   1          -- on success
c   odd number -- on error
c
c               TREE$K_USAGE_ANY       0
c               TREE$K_USAGE_STRUCTURE 1
c               TREE$K_USAGE_ACTION    2
c               TREE$K_USAGE_DEVICE    3
c               TREE$K_USAGE_DISPATCH  4
c               TREE$K_USAGE_NUMERIC   5
c               TREE$K_USAGE_SIGNAL    6
c               TREE$K_USAGE_TASK      7
c               TREE$K_USAGE_TEXT      8
c               TREE$K_USAGE_WINDOW    9
c               TREE$K_USAGE_AXIS     10
c               TREE$K_USAGE_SUBTREE  11
c               TREE$K_USAGE_COMPOUND_DATA 12
c
c 03/24/1999 ler converted to Unix
c 05/23/1999 ler *** ERROR IN CODE
C 12/14/1999 cal make it work
c-----------------------------------------------------------------------
      integer function mds_create_node(node, usage)
      implicit none
      character*(*) node, usage
 
      integer     status, cmd_status, cmd_stat2
      integer     dsc_cmd_status, dsc_cmd_stat2
      integer     retlen
 
      integer istat, icmdstat
      character*4097 strarg
 
      integer descr, MdsValue, mds__long
      integer mdsdclfdo
      integer nl, nu, nn
      character*1    bcks
      volatile      cmd_status,     cmd_stat2
C---------------------------------------------------------------------------
      bcks=char(92)                 ! x5c = '\'
 
c      write(cusage,'(i2)',err=99) usage
      status = 0
      cmd_status = 0
      dsc_cmd_status = descr(mds__long(), cmd_status, 0)
C Check if open for edit
      strarg = 'GETDBI("OPEN_FOR_EDIT")' // CHAR(0)
cxx      print*,'mds_create_node: ',strarg(1:28)
      status = MdsValue(strarg, dsc_cmd_status, 0, retlen)
cxx      print *,'status =',status,',cmd_status =',cmd_status
      istat    = mod(status, 2)
      icmdstat = mod(cmd_status, 2)
 
      if (istat .ne. 1 .or. icmdstat .ne. 1) then
         print *,'?mds_create_node: GETDBI sts = ',status
         print *,'?mds_create_node: GETDBI cmd_sts = ',
     >            cmd_status
      endif
      if (istat.eq.1 .or. status.eq.265388200) then
C If failure, open for edit:
         cmd_stat2 = 0
         if (icmdstat.ne.1) then
            dsc_cmd_stat2 = descr(mds__long(), cmd_stat2, 0)
#ifdef __VMS
            strarg = 'TREESHR->TREE$OPEN_TREE_EDIT($EXPT,$SHOT)'
     >               //CHAR(0)
#else
            strarg = 'tcl("EDIT "//$EXPT//"/shot="//$SHOTNAME)'
     >               // CHAR(0)
#endif
            print*,'mds_create_node: ',strarg(1:28)
            status = MdsValue(strarg, dsc_cmd_stat2, 0, retlen)
            print *,'status =',status,',cmd_stat2 =',cmd_stat2
            istat    = mod(status, 2)
            icmdstat = mod(cmd_stat2, 2)
         endif
         nl=len_trim(node)
         nu=len_trim(usage)


         if ((istat.eq.1) .and. (icmdstat.eq.1)) then
            if (node(1:1).eq.bcks) then
               strarg = 'tcl("add node ' // bcks // node //
     >              '/usage=' // usage // '")' //CHAR(0)
            else
               strarg = 'tcl("add node ' // node(1:nl) //
     >              '/usage=' // usage(1:nu) // '")' //CHAR(0)
            endif
            nn=len_trim(strarg)
            status = MdsValue(strarg(1:nn), dsc_cmd_status, 0, retlen)
cxx            print *,'status =',status,', cmd_status =',cmd_status
            istat    = mod(status, 2)
            icmdstat = mod(cmd_status, 2)
            if (istat.eq.1) then
               status = cmd_status
            endif
         else
            print *,'?mds_create_node: failed on: ',strarg(1:60)
            print *,'?mds_create_node: EDIT status = ',status
            print *,'?mds_create_node: EDIT cmd_stat = ',cmd_stat2
         endif
      endif
      go to 100
C 99   continue
C      write(6,*) '?mds_create_node:  illegal usage argument:  ',usage
 100  continue
 
      mds_create_node = status
 
      return
      end
 
#endif
