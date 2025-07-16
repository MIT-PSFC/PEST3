      subroutine trmds_parse(zbuf,mds_server,mds_tree,mds_id1,mds_id2,
     >   mds_tstr,ierr)
c
c  parse an MDS+ tree id string
c  of form
c
c     <server_name>:<tree_name>(<id1>[,<id2>][;t=<time_val>])
c
c  where
c
c      server_name -- is an MDSplus server path with possible port no.
c          examples:
c                     local
c                     birch
c                     birch.pppl.gov
c                     europa.pppl.gov:8501
c                     cmoda.psfc.mit.edu
c
c      tree_name -- an MDSplus tree name
c
c      id1,id2 -- if id2 is absent, id1 = MDSplus shot no.
c                 if id2 is present, id1 = (TRANSP) tok.yy code
c                                    id2 = (TRANSP) runid
c
c      time_val -- if ;t=<time_val> qualifier is present, time_val
c                                    is the time of interest, seconds.
c
      implicit NONE
c
c  input:
c
      character*(*) zbuf                ! string to parse
c
c  output:
c
      character*(*) mds_server          ! server
      character*(*) mds_tree            ! tree name
      character*(*) mds_id1             ! id1 -- shot# or tok.yy if id2 present
      character*(*) mds_id2             ! id2 -- runid (if present) or blank
      character*(*) mds_tstr            ! time_val-- time, seconds, if nonblank
c
      integer ierr                      ! completion code, 0=OK, else parse err
c
c-----------------------------------------------------------
      integer ilz,icln,ilparen,icomma,isemi,ieqs,i
      integer is1,is2,it1,it2,ia1,ia2
      integer str_length
      character*10 zwd
c-----------------------------------------------------------
c
      mds_server=' '
      mds_tree=' '
      mds_id1=' '
      mds_id2=' '
      mds_tstr=' '
c
c  find length, left parenthesis (must be present)
c
      ilz=str_length(zbuf)
      ilparen=index(zbuf,'(')
      if(ilparen.eq.0) then
         ierr=91
         return
      endif
c
c  find terminating right parenthesis (must be present)
c
      if(zbuf(ilz:ilz).ne.')') then
         ierr=99
         return
      endif
c
c  find right terminus of server string
c
      do i=ilparen,2,-1
         if(zbuf(i:i).eq.':') then
            icln=i
            go to 100
         endif
      enddo
c
c  no server string found
c
      ierr=92
      return
c
 100  continue
c
c  server name...
      call trmds_pshrink(zbuf,1,icln-1,is1,is2,ierr)
      if(ierr.ne.0) return
      mds_server=zbuf(is1:is2)
      call uupper(mds_server)
c
c  tree name...
      call trmds_pshrink(zbuf,icln+1,ilparen-1,it1,it2,ierr)
      if(ierr.ne.0) return
      mds_tree=zbuf(it1:it2)
      call uupper(mds_tree)
c
c  optional time-of-interest
c
      isemi=index(zbuf(ilparen:ilz),';')
      if(isemi.ne.0) then
         isemi=isemi+ilparen-1
         call trmds_pshrink(zbuf,isemi+1,ilz-1,it1,it2,ierr)
         if(ierr.ne.0) return
c
c  zbuf(it1:it2) contains non-blank string between ; and closing paren
c    this string must have form "t = <value>"
c
         ieqs=index(zbuf(it1:it2),'=')
         if(ieqs.eq.0) then
            ierr=93
            return
         endif
         ieqs=it1+ieqs-1
c
c  LHS
         call trmds_pshrink(zbuf,it1,ieqs-1,ia1,ia2,ierr)
         if(ierr.ne.0) return
         zwd=zbuf(ia1:ia2)
         call uupper(zwd)
         if((zwd.ne.'T').and.(zwd.ne.'TIME')) then
            ierr=94
            return
         endif
c
c  RHS
c
         call trmds_pshrink(zbuf,ieqs+1,it2,ia1,ia2,ierr)
         if(ierr.ne.0) return
c
         mds_tstr=zbuf(ia1:ia2)
         call uupper(mds_tstr)
c
         ilz=isemi
      endif
c
c  tree id args
c
      icomma=index(zbuf(ilparen:ilz),',')
      if(icomma.eq.0) then
         call trmds_pshrink(zbuf,ilparen+1,ilz-1,ia1,ia2,ierr)
         if(ierr.ne.0) return
         mds_id1=zbuf(ia1:ia2)
         call uupper(mds_id1)
      else
         icomma=ilparen+icomma-1
         call trmds_pshrink(zbuf,ilparen+1,icomma-1,ia1,ia2,ierr)
         if(ierr.ne.0) return
         mds_id1=zbuf(ia1:ia2)
         call uupper(mds_id1)
         call trmds_pshrink(zbuf,icomma+1,ilz-1,ia1,ia2,ierr)
         if(ierr.ne.0) return
         mds_id2=zbuf(ia1:ia2)
         call uupper(mds_id2)
      endif
c
      return
      end
c-----------------------------------------------
      subroutine trmds_pshrink(zbuf,ilim1,ilim2,iact1,iact2,ierr)
c
c  parsing primitive
c  ...find first and last non-blank in zbuf(ilim1:ilim2)
c  ...return ierr.eq.0 if successful
c  ...return ierr.eq.1 if ilim1.gt.ilim2 or ilim1.lt.1 or
c     ilim2.gt.len(zbuf)
c
c  ...non-blank range returned in iact1,iact2
c
      implicit NONE
c
      character*(*) zbuf
      integer ilim1,ilim2
      integer iact1,iact2
      integer ierr
c
c--------------------
c
      ierr=1
      if(ilim1.gt.ilim2) return
      if(ilim1.lt.1) return
      if(ilim2.gt.len(zbuf)) return
c
      if(zbuf(ilim1:ilim2).eq.' ') return
c
c  ok, valid range, non-blanks are present
c
      ierr=0
      do iact1=ilim1,ilim2
         if(zbuf(iact1:iact1).ne.' ') go to 10
      enddo
 10   continue
      do iact2=ilim2,ilim1,-1
         if(zbuf(iact2:iact2).ne.' ') go to 20
      enddo
 20   continue
c
      return
      end
