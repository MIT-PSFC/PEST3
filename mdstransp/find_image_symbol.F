C Find_image_symbol - routine to call LIB$FIND_IMAGE_SYMBOL with the
C  sig_to_ret condtion handler.  This is done because lib$find_image_symbol
C  signals errors instead of returning them.
C
C  Josh Stillerman 12/5/94
C
#ifdef __VMS
        options /extend
#endif
 
      Integer*4 function Find_image_symbol(image, routine, addr)
 
#ifdef __VMS
      Implicit NONE
      Character*(*) image
      Character*(*) routine
      Integer*4 addr
#endif
 
      Integer*4 STATUS
 
#ifdef __VMS
      External LIB$SIG_TO_RET
      Integer*4 LIB$FIND_IMAGE_SYMBOL
#endif
 
      status = 0
 
#ifdef __VMS
      call LIB$ESTABLISH(LIB$SIG_TO_RET)
      STATUS = LIB$FIND_IMAGE_SYMBOL(image, routine, addr)
#endif
 
      Find_image_symbol = status
      return
 
      end
