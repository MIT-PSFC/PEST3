#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_setserver(server)
      implicit none
      character*(*) server
      mds_setserver = 0
      return
      end
#else
c-----------------------------------------------------------------------
c mds_setserver
c select remote MDSplus server
c
c Input:
c   server -- MDSplus server name or "local"
c Returns:
c          -- mdsplus status codes
c             odd_number success
c             even_number failure
c
c *** WARNING ***
c o Saved current server name only 256 characters long
c o Assume consistent case in server names
c o Assume "local" in lower case
c
c
c 03/03/1999 ler Converting to unix usable
c 04/12/1999 ler Now: MdsConnect returns socket value
c                     MdsDisconnect exists
c-----------------------------------------------------------------------
      integer function mds_setserver(server)
      implicit none
 
      character*(*) server
      integer     status
 
      character*256 current_server
      data          current_server/' '/
      integer       socket
      data          socket/0/
 
      character*4097 str_server
      integer l
 
      integer  mdsfconnect
      integer mds__leniftrim
 
      save socket
C-----------------------------------------------------------------------
C You might not wish to disconnect
C      if (socket.ne.0) then
C        call mdsfdisconnect()
C      endif
 
      l = max(1, mds__leniftrim(server))
      if (server(1:l).eq.'local' .or.
     >    server(1:l).eq.'LOCAL') then
        socket = 0
        current_server = server(1:l)
        status = 1
      else
        call mds__copywnull(str_server, server)
        socket = mdsfconnect(str_server)
C Error = -1
        if (socket.le.0) then
          current_server = ' '
          status = 0
        else
          current_server = server(1:l)
          status = 1
        end if
      end if
 
      mds_setserver = status
      return
      end
#endif
 
