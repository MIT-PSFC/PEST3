      subroutine bad_exit
C
#ifdef __UNIX
#ifdef __MPI
      write(6,*) 
     >' %bad_exit:  generic f77 error exit call (errset_mpi status 999)'
#else
      write(6,*) 
     >' %bad_exit:  generic f77 error exit call (errset status 999)'

#endif
#else
      write(6,*) ' %bad_exit:  generic f77 error exit call.'
#endif
C
      call err_end
      call fort_flush(6)
C
#ifdef __VMS
      call lib$stop( %val(2) )
#endif
#ifdef __UNIX
      call errset_mpi(-1,999)
#endif
C
      STOP                              ! this line not reached.
      END
