module geqdsk_aux
  !
  ! support module for geqdsk_mod  ---  static items, constants
  ! see geqdsk_mod.f90
  !   (dmc Jan 2001)
  !
  implicit NONE
  SAVE

  integer :: neqdsk = 77          ! LUN for G-EQdisk file i/o
 
  integer :: lunmsg = 6           ! LUN for messages
 
  integer, parameter :: auxr8 = selected_real_kind(12,100)
 
  integer, parameter :: file_access = 1
  integer, parameter :: mds_full_access = 2
  integer, parameter :: mds_slice_access = 3
 
  integer, parameter :: ctable_max = 20
 
  ! units conversion tables to support G-EQdisk MDSplus i/o
 
  type ctable
 
     integer :: n
     character*16 :: id
     character*32 :: units_string(ctable_max)
     real(auxr8) :: conversion_factor(ctable_max)
 
  end type ctable
 
  type(ctable) :: time_units,length_units,bfield_units,current_units
  type(ctable) :: psi_units,fpol_units,pressure_units,energy_units
  type(ctable) :: pprime_units,ffprime_units,no_units

  logical :: ctable_init_flag = .FALSE.
 
contains
 
  subroutine ctable_init
 
    implicit NONE
 
    !  **hand maintained** tables --
    !  **goal is to support conversions for G-EQDSK quantities from all
    !    EFIT MDSplus trees at all experimenal sites
    !       (dmc 1 Feb 2001)
    !
    !    first entry for each table should be the "standard" MKS label
 
    if(ctable_init_flag) return
 
    call ctable_clear(time_units,'time')
    call ctable_add(time_units, 's', 1.0_auxr8)
    call ctable_add(time_units, 'sec', 1.0_auxr8)
    call ctable_add(time_units, 'seconds', 1.0_auxr8)
    call ctable_add(time_units, 'ms', 1.0e-3_auxr8)
    call ctable_add(time_units, 'msec', 1.0e-3_auxr8)
 
    call ctable_clear(length_units,'length')
    call ctable_add(length_units, 'm', 1.0_auxr8)
    call ctable_add(length_units, 'meters', 1.0_auxr8)
    call ctable_add(length_units, 'cm', 1.0e-2_auxr8)
 
    call ctable_clear(bfield_units,'B_field')
    call ctable_add(bfield_units, 't', 1.0_auxr8)
    call ctable_add(bfield_units, 'gauss', 1.0e-4_auxr8)
    call ctable_add(bfield_units, 'kgauss', 1.0e-1_auxr8)
 
    call ctable_clear(current_units,'current')
    call ctable_add(current_units, 'a', 1.0_auxr8)
    call ctable_add(current_units, 'ka', 1.0e3_auxr8)
 
    call ctable_clear(psi_units,'pol_flux_stream_fcn')
    call ctable_add(psi_units, 'Wb/rad', 1.0_auxr8)
    call ctable_add(psi_units, 'v s / rad', 1.0_auxr8)
 
    call ctable_clear(fpol_units,'fpol')
    call ctable_add(fpol_units, 'mT', 1.0_auxr8)
    call ctable_add(fpol_units, 'Tm', 1.0_auxr8)
    call ctable_add(fpol_units, 'm*T', 1.0_auxr8)
    call ctable_add(fpol_units, 'T*m', 1.0_auxr8)
    call ctable_add(fpol_units, 'cm*gauss', 1.0e-6_auxr8)
 
    call ctable_clear(pressure_units,'Pressure')
    call ctable_add(pressure_units, 'Pa', 1.0_auxr8)
    call ctable_add(pressure_units, 'N/m^2', 1.0_auxr8)
 
    call ctable_clear(energy_units,'Energy')
    call ctable_add(pressure_units, 'J', 1.0_auxr8)
    call ctable_add(pressure_units, 'Jles', 1.0_auxr8)
 
    call ctable_clear(pprime_units,'PPrime')
    call ctable_add(pprime_units, 'Pa*rad/Wb', 1.0_auxr8)
    call ctable_add(pprime_units, 'Pa/Wb/rad', 1.0_auxr8)
    call ctable_add(pprime_units, '(N/m^2)/(V s / rad)', 1.0_auxr8)
    call ctable_add(pprime_units, '(N/m^2)/(Wb/rad)', 1.0_auxr8)
 
    call ctable_clear(ffprime_units,'ffprime')
    call ctable_add(ffprime_units, '(mT)^2*rad/Wb', 1.0_auxr8)
    call ctable_add(ffprime_units, '(Tm)^2*rad/Wb', 1.0_auxr8)
    call ctable_add(ffprime_units, '(mT)^2/(Wb/rad)', 1.0_auxr8)
    call ctable_add(ffprime_units, '(Tm)^2/(Wb/rad)', 1.0_auxr8)
    call ctable_add(ffprime_units, '(Tm)^2/(V s / rad)', 1.0_auxr8)
    call ctable_add(ffprime_units, 'Rad/m', 1.0_auxr8)
 
    call ctable_clear(no_units,'no_units')
    call ctable_add(no_units, ' ', 1.0_auxr8)
    call ctable_add(no_units, '-', 1.0_auxr8)
 
    ctable_init_flag = .TRUE.

  end subroutine ctable_init
 
!----------------------------------------------
 
  subroutine ctable_clear(ct,lbl)
 
    implicit NONE
 
    ! clear a conversion table
 
    type(ctable), intent(out) :: ct
    character*(*), intent(in) :: lbl
 
    ct%id = lbl
    ct%n = 0
    ct%units_string = ' '
    ct%conversion_factor = 0.0_auxr8
 
  end subroutine ctable_clear
 
!----------------------------------------------
 
  subroutine ctable_add(ct, string, factor)
 
    implicit NONE
 
    !
    ! add an entry (string, conversion factor) to a units conversion table
    !
 
    type(ctable), intent(inout) :: ct    ! conversion table object
    character*(*), intent(in) :: string  ! new units label
    real(auxr8), intent(in) :: factor    ! conversion factor
 
    integer i,j,str_length
    character*32 cbuf
 
    if( ct%n .ge. ctable_max ) then
       write(lunmsg,*) ' warning on conversion table for "',ct%id,'".'
       write(lunmsg,*) '%geqdsk_aux:  table full, cannot add string: ', &
            & string,' factor: ',factor
       write(lunmsg,*) ' ctable_max = ',ctable_max, &
            & ' -> increase ctable_max parameter in geqdsk_aux module.'
 
    else if(ct%n.lt.0) then
       write(lunmsg,*) ' warning on conversion table for "',ct%id,'".'
       write(lunmsg,*) '?geqdsk_aux:  program error, call ctable_clear first!'
 
    else
 
       ! strip out white space; convert to uppercase
 
       j=0
       cbuf=' '
       do i=1,str_length(string)
          if((string(i:i).ne.' ').and.(string(i:i).ne.char(0)).and. &
               & string(i:i).ne.char(9)) then
             j=j+1
             cbuf(j:j)=string(i:i)
          endif
       enddo
       call uupper(cbuf)
 
       ct%n = ct%n + 1
       ct%units_string(ct%n) = cbuf
       ct%conversion_factor(ct%n) = factor
 
    endif
 
  end subroutine ctable_add
 
!----------------------------------------------
 
  subroutine ctable_check(ct, string, factor)
 
    !
    ! check a units string against contents of conversion table;
    ! output corresponding conversion factor
    ! ...if no match found, output warning and assume factor=1.0
    !
    ! caution:  side effect:  string is converted to uppercase
    !
 
    implicit NONE
 
    type(ctable), intent(in) :: ct       ! conversion table object
    character*(*), intent(inout) :: string  ! units label to be tested
    real(auxr8), intent(out) :: factor   ! conversion factor
 
    integer imatch
 
    if((ct%n.le.0).or.(ct%n.gt.ctable_max)) then
       write(lunmsg,*) ' warning on conversion table for "',ct%id,'".'
       write(lunmsg,*) '?geqdsk_aux:  program error, call ctable_init first.'
       write(lunmsg,*) ' conversion factor = 1.0 assumed.'
       factor = 1.0_auxr8
 
    else
 
       call ctable_match(ct, string, imatch)
 
       if(imatch.gt.0) then
          factor = ct%conversion_factor(imatch)
 
       else
          write(lunmsg,*) '%geqdsk_aux:  no match for units label:  "', &
               & string,'"'
          write(lunmsg,*) ' warning from conversion table for "',ct%id,'".'
          write(lunmsg,*) ' conversion factor = 1.0 assumed.'
          factor = 1.0_auxr8
       endif
 
    endif
 
  end subroutine ctable_check
 
!----------------------------------------------
 
  subroutine ctable_match(ct, string, imatch)
 
    !
    ! check a units string against contents of conversion table;
    ! output index to matching table entry
    ! or imatch=0 if there is no match
    ! or imatch=-1 if there is an error
    !
    ! caution:  side effect:  string is converted to uppercase
    !
    ! this routine writes no warning messages
    !
 
    implicit NONE
 
    type(ctable), intent(in) :: ct       ! conversion table object
    character*(*), intent(inout) :: string  ! units label to be tested
    integer, intent(out) :: imatch
 
    integer i,j,str_length
    character*32 cbuf
 
    imatch=-1
    if(ct%n.le.0) return
    if(ct%n.gt.ctable_max) return
 
    imatch=0
 
    j=0
    cbuf=' '
    do i=1,str_length(string)
       if((string(i:i).ne.' ').and.(string(i:i).ne.char(0)).and. &
            & string(i:i).ne.char(9)) then
          j=j+1
          cbuf(j:j)=string(i:i)
       endif
    enddo
    call uupper(cbuf)
 
    do i=1,ct%n
       if(cbuf.eq.ct%units_string(i)) then
          imatch=i
          exit
       endif
    enddo
 
  end subroutine ctable_match
 
end module geqdsk_aux
