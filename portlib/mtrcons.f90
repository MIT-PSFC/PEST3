   module mtrcons
!. mtrcons                           Various constants collected together, to enable
!                                    a degree of consistency between routines
!
!
!  Author                            jim.conboy@ccfe.ac.uk
!  Version                           1.00,  06Jan2012
!  Modifications
!
!  1.00  06Jan2012                   jim.conboy@ccfe.ac.uk
!                                    Written
!----_^---------------------------------=========================================|
!
!
   implicit  none
!
!..GLOBAL
!
   integer, parameter           :: R8      = selected_real_kind(12,100)           !  64 bit
   !
   character(len=*), parameter  :: calph_L = 'abcdefghijklmnopqrstuvwxyz'  &
                                  ,calph_U = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'  &
                                  ,cnum    = '0123456789'                  &
                                  ,cschar  = '`¬!"£$%^&*()_:@~;''#?,./'    &
                                  ,coper   = '+-*/='                       &
                                  ,cbrkt   = '(){}[]<>'
   !
   character(len=1), parameter  :: cdlm    = '/'   & ! file path delimiter
                                  ,cdkdlm  = ':'   & ! disk name delimiter ( VMS )
                                  ,cdol    = '$'   & ! dollar
                                  ,cdot    = '.'   & ! point
                                  ,cus     = '_'   & ! U/score  
                                  ,cq      = ''''  & ! single quote
                                  ,cqq     = '"'   & ! double quote
                                  ,coq     = '`'   & ! open quote
                                  ,cspace  = ' '   &
                                  ,char0   = char(0)
   !
   character(len=8), parameter  :: csp8    = '        '
   !
   character(len=*), parameter  :: cfpath  = calph_L // calph_U // cnum // &
                                             cdlm // cdot // cus // cdol
   !
   character(len=*), parameter              :: cM = 'portlib/mtrcons'
!  contains

   end module mtrcons
