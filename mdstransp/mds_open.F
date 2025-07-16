#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
        integer function mds_open(tree, shot)
        implicit none
        character*(*) tree
        integer*4     shot
        character(len=*),parameter  :: cr='mds_open'

        print *,cr,' -D__NOMDSPLUS > dummy mdstransp/mds_open '
        mds_open = 0

        end
#else
c-----------------------------------------------------------------------
c mds_open
c open a tree
c
c 03/03/1999 ler Converting to unix usable
c-----------------------------------------------------------------------
      integer function mds_open(tree, shot)
      implicit none
 
      character*(*) tree
      integer       shot
 
      integer l, status
      character*4097 strarg
 
      integer MdsOpen
 
      call mds__copywnull(strarg, tree)
      status = MdsOpen(strarg, shot)
 
      mds_open = status
 
      return
      end
#endif
