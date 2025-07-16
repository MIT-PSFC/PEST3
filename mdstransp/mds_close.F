#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_close(tree, shot)
      implicit none
      character*(*) tree
      integer       shot
      mds_close = 0
      end
#else
c-----------------------------------------------------------------------
c mds_close
c close a tree
c
c Returns:
c   1          on succcess
c   odd number on failure
c
c 03/17/1999 ler Converting to unix usable
c-----------------------------------------------------------------------
      integer function mds_close(tree, shot)
      implicit none
 
      character*(*) tree
      integer       shot
 
      integer status
      character*4097 strarg
 
      integer MdsClose

      call mds__copywnull(strarg, tree)
      status = MdsClose(strarg, shot)
      mds_close = status
 
      return
      end
#endif
 
