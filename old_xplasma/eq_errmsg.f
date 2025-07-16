      subroutine eq_errmsg(msg)
c
      use xplasma_obj_instance
      use eq_module
c
      IMPLICIT NONE
c
      character*(*) msg                 ! error message to write
c
c  write error message
c
      integer str_length,ilmsg,ier
      character*80 zmodel
c---------------------------------------------------------------------
c
      write(lunerr,*) ' '
      if(eq_premsg_text.ne.' ') then
         ilmsg=str_length(eq_premsg_text)
         write(lunerr,*) eq_premsg_text(1:ilmsg)
      else
         call xplasma_global_info(s,ier, initLabel=zmodel)
         if(ier.ne.0) zmodel = '?? -- xplasma not initialized.'
         ilmsg=str_length(zmodel)
         write(lunerr,*)
     >      ' (model:  ',zmodel(1:ilmsg),') eq_errmsg report:'
      endif
c
      ilmsg=max(1,str_length(msg))
      write(lunerr,*) msg(1:ilmsg)
c
c  print any saved messages in the xplasma object
c
      call xplasma_error(s,9999,lunerr)
c
      return
      end
