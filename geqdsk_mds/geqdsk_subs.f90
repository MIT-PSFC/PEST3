!------------------------------------------------------
  subroutine geq_ckrzc(r,z,n,ndim)
 
    !  squeeze out duplicate points, make sure end pts connect
 
    implicit NONE
 
    integer, parameter :: r8 = selected_real_kind(12,100)
    integer, intent(in) :: ndim     ! array dims
    integer, intent(inout) :: n     ! actual number in use
    real(r8), dimension(ndim), intent(inout) :: r,z    ! the contour
 
    integer i,j
    real(r8) zmax,rmax,dl,dlmin
 
    ! --------------------------------------
 
    ! squeeze out consecutive duplicate points
 
    rmax=maxval(abs(r(1:n)))
    zmax=maxval(abs(z(1:n)))
    dlmin=5.0e-6_r8*max(rmax,zmax)
 
    i=0
    do
       i=i+1
       if(i.ge.n) exit
       dl=sqrt((r(i+1)-r(i))**2+(z(i+1)-z(i))**2)
       if(dl.lt.dlmin) then
          do j=i+1,n-1
             r(j)=r(j+1)
             z(j)=z(j+1)
          enddo
          n=n-1
          i=i-1
       endif
    enddo
 
    ! make sure end pts are connected
 
    if((r(1).ne.r(n)).or.(z(1).ne.z(n))) then
       n=n+1
       r(n)=r(1)
       z(n)=z(1)
    endif
 
    ! make sure no NaNs left at high end
 
    do i=n+1,ndim
       r(i)=r(n)
       z(i)=z(n)
    enddo
 
  end subroutine geq_ckrzc
!------------------------------------------------------
  subroutine geq_units_echeck_r4(expr,zdatr4,ndata,ierr)
    !
    !  R4 version of geq_units_echeck
    !
    character*(*), intent(in) :: expr    ! MDS+ expression
    integer, intent(in) :: ndata         ! #data pts
    real, intent(inout) :: zdatr4(ndata) ! data pts (may be rescaled)
    integer, intent(out) :: ierr         ! completion code (out)

    integer, parameter :: r8 = selected_real_kind(12,100)
    real(r8) :: zdata(ndata)

    !--------------------------

    zdata = zdatr4
    call geq_units_echeck(expr,zdata,ndata,ierr)
    zdatr4 = zdata

  end subroutine geq_units_echeck_r4
!------------------------------------------------------
  subroutine geq_units_echeck(expr,zdata,ndata,ierr)
    !
    !  geq_units_echeck -- check/convert units for data corresponding
    !    to certain known expressions
    !
    use geqdsk_aux
    implicit NONE

    integer, parameter :: r8 = selected_real_kind(12,100)

    character*(*), intent(in) :: expr    ! MDS+ expression
    integer, intent(in) :: ndata         ! #data pts
    real(r8), intent(inout) :: zdata(ndata) ! data pts (may be rescaled)
    integer, intent(out) :: ierr         ! completion code (out)

    !------------

    character*1 bksl,zcln
    integer :: ii,idelim,ilen,str_length
    character*20 ztest

    !-------------
    !  check units corresponding to certain known MDS+ expressions
    !  callable without "use geqdsk_aux".

    ierr = 0
    zcln = ':'
    bksl = char(92)  ! back slash

    ilen=index(expr,'[')-1    ! avoid subscript expressions
    if(ilen.le.0) ilen=str_length(expr)
    do
       if(expr(ilen:ilen).eq.')') then
          ilen=ilen-1
       else
          exit
       endif
    enddo

    idelim=0
    do ii=ilen,1,-1
       if((expr(ii:ii).eq.bksl).or.(expr(ii:ii).eq.zcln)) then
          idelim = ii
          exit
       endif
    enddo

    ztest = expr(idelim+1:ilen)
    call uupper(ztest)

    if((ztest.eq.'GTIME').or.(ztest.eq.'ATIME')) then
       call geq_units_check(expr,zdata,ndata,time_units)
    else if(ztest.eq.'CPASMA') then
       call geq_units_check(expr,zdata,ndata,current_units)
    else if(ztest.eq.'WMHD') then
       call geq_units_check(expr,zdata,ndata,pressure_units)
    else
       ierr=1
       write(lunmsg,*) ' %geq_units_echeck: not recognized for units check:', &
            ' ',trim(expr)
    endif
  end subroutine geq_units_echeck

!------------------------------------------------------
  subroutine geq_units_check(expr,zdata,ndata,data_units)
    !
    ! fetch units label; check this and convert data if necessary
    ! using the conversion table passed
    !
 
    use geqdsk_aux
    implicit NONE
 
    integer, parameter :: r8 = selected_real_kind(12,100)
    character*(*), intent(in) :: expr
    integer, intent(in) :: ndata
    real(r8), intent(inout) :: zdata(ndata)
    type(ctable), intent(in) :: data_units
 
 ! local:
 
    integer iflag,ilen,str_length
    character*5 ztest
    character*32 actual_units
    integer idescr_cstring
    integer istat,mds_value,iret
    real(auxr8) factor
 
 ! executable code:
 
    call ctable_init

    ilen=index(expr,'[')-1    ! avoid subscript expressions
    if(ilen.le.0) ilen=str_length(expr)
 
    iflag=0
    if(ilen.le.7) then
       iflag=1
    else
       ztest=expr(1:5)
       call uupper(ztest)
       if(ztest.ne.'DATA(') iflag=1
    endif
 
    if(iflag.eq.1) then
       write(lunmsg,*) ' %geqdsk_mod(geq_units_check):  expr = ',expr(1:ilen)
       write(lunmsg,*) '  not a data expression -- likely programming error.'
       return
    endif
 
    istat = mds_value('UNITS('//expr(6:ilen),idescr_cstring(actual_units), &
         & iret)
 
    if(mod(istat,2).ne.1) then
       call mdserr(lunmsg,'%geq_units_check: i/o: UNITS('//expr(6:ilen),istat)
       return
    endif
 
    call ctable_check(data_units, actual_units, factor)
    zdata = factor*zdata
 
  end subroutine geq_units_check
!------------------------------------------------------
  subroutine geq_mds_sub1(zbuf,nbuf,zsav,nsav,indx,nodename, &
       & treename,expr,data_units,istat)
!
!  read and copy a floating scalar fcn of time node, saving valid points only
!
    use geqdsk_aux
    implicit NONE
 
    integer, parameter :: r8 = selected_real_kind(12,100)
    integer, intent(in) :: nbuf,nsav
    real(r8), intent(inout) :: zbuf(nbuf)
    real(r8), intent(out) :: zsav(nsav)
    integer, intent(in) :: indx(nbuf)
 
    character*(*), intent(in) :: treename
    character*(*), intent(in) :: nodename
    character*(*), intent(out) :: expr
 
    type(ctable), intent(in) :: data_units   ! known units labels & conversions
 
    integer, intent(out) :: istat
 
! local
 
    integer iln,ilt,str_length, iok
    integer idims(2),iret,mds_value,idescr_floatarr
    integer i
    character*1 bksl
 
    real, dimension(:), allocatable:: r4buf
 
! executable code
 
    call ctable_init

    allocate(r4buf(nbuf), stat=iok)
    if(iok/=0) then
       write(lunmsg,*) '***failed to allocate r4buf in geq_mds_sub1'
       istat = 9999
       return
    endif
 
    bksl = char(92)   ! backslash
 
    idims(1)=nbuf
    idims(2)=0
 
    iln=str_length(nodename)
    ilt=str_length(treename)
 
    call geq_mk_expr('data(',nodename(1:iln),')',expr)
    istat = mds_value(expr,idescr_floatarr(r4buf,idims,1),iret)
    zbuf=r4buf
    if(mod(istat,2).eq.1) then
       do i=1,nbuf
          if(indx(i).ne.0) zsav(indx(i))=zbuf(i)
       enddo
       call geq_units_check(expr,zsav,nsav,data_units)
    endif
 
    deallocate(r4buf)
    return
  end subroutine geq_mds_sub1
!------------------------------------------------------
  subroutine geq_mds_subi(ibuf,nbuf,isav,nsav,indx,nodename, &
       & treename,expr,istat)
!
!  read and copy an integer scalar fcn of time node, saving valid points only
!
    implicit NONE
 
    integer, parameter :: r8 = selected_real_kind(12,100)
    integer, intent(in) :: nbuf,nsav
    integer, intent(inout) :: ibuf(nbuf)
    integer, intent(out) :: isav(nsav)
    integer, intent(in) :: indx(nbuf)
 
    character*(*), intent(in) :: treename
    character*(*), intent(in) :: nodename
    character*(*), intent(out) :: expr
 
    integer, intent(out) :: istat
 
! local
 
    integer iln,ilt,str_length
    integer idims(2),iret,mds_value,idescr_longarr
    integer i
    character*1 bksl
 
! executable code
 
    bksl = char(92)   ! backslash
 
    idims(1)=nbuf
    idims(2)=0
 
    iln=str_length(nodename)
    ilt=str_length(treename)
 
    call geq_mk_expr('data(',nodename(1:iln),')',expr)
    istat = mds_value(expr,idescr_longarr(ibuf,idims,1),iret)
    if(mod(istat,2).eq.1) then
       do i=1,nbuf
          if(indx(i).ne.0) isav(indx(i))=ibuf(i)
       enddo
    endif
 
    return
  end subroutine geq_mds_subi
!------------------------------------------------------
  subroutine geq_mds_sub2(nd1,nt_orig,zsav,nsav,ktdim,indx,nodename, &
       & treename,expr,data_units,istat)
!
!  read and copy a floating profile fcn of time node, saving valid points only
!  *** the code assumes f(x,t) is stored with x variation contiguous ***
!
    use geqdsk_aux
    implicit NONE
 
    integer, parameter :: r8 = selected_real_kind(12,100)
    integer, intent(in) :: nd1,nt_orig,nsav
    real(r8), intent(out) :: zsav(nd1,nsav)
    integer, intent(in) :: indx(nt_orig)
    integer, intent(in) :: ktdim            ! specifies which dim. is time
 
    character*(*), intent(in) :: treename
    character*(*), intent(in) :: nodename
    character*(*), intent(out) :: expr
 
    type(ctable), intent(in) :: data_units   ! known units labels & conversions
 
    integer, intent(out) :: istat
 
! local
 
    integer iln,ilt,str_length
    integer idims(2),iret,mds_value,idescr_floatarr,idescr_cstring
    integer i
    character*1 bksl
    integer imatch, iok
    character*20 units_string
 
    real, dimension(:,:), allocatable :: zbuf
 
! executable code
 
    call ctable_init

    bksl = char(92)   ! backslash
!
! determine if "time" is 1st or 2nd dimension -- ck 1st dimension...
!
    iln=str_length(nodename)
    ilt=str_length(treename)
 
    if(ktdim.eq.2) then
       idims(1)=nd1
       idims(2)=nt_orig
    else
       idims(1)=nt_orig
       idims(2)=nd1
    endif
 
    allocate(zbuf(idims(1), idims(2)), stat=iok)
    if(iok/=0) then
       write(lunmsg,*) '***failed to allocate zbuf in geq_mds_sub2'
       istat = 9999
       return
    endif
 
    call geq_mk_expr('data(',nodename(1:iln),')',expr)
    istat = mds_value(expr,idescr_floatarr(zbuf,idims,2),iret)
    if(mod(istat,2).eq.1) then
       do i=1,nt_orig
          if(indx(i).ne.0) then
             if(ktdim.eq.2) then
                zsav(1:nd1,indx(i))=zbuf(1:nd1,i)
             else
                zsav(1:nd1,indx(i))=zbuf(i,1:nd1)
             endif
          endif
       enddo
       call geq_units_check(expr,zsav,nd1*nsav,data_units)
    endif
 
    deallocate(zbuf)
 
    return
  end subroutine geq_mds_sub2
 
!------------------------------------------------------
  subroutine geq_mds_sub3(nd1,nd2,nt_orig,zsav,nsav,ktdim,indx,nodename, &
       & treename,expr,data_units,istat)
!
!  read and copy a floating profile f(R,Z,t) node, saving valid points only
!  *** the code assumes R is 1st dimension (contiguous storage), Z is 2nd
!      dimension, t is 3rd dimension
!
    use geqdsk_aux
    implicit NONE
 
    integer, parameter :: r8 = selected_real_kind(12,100)
    integer, intent(in) :: nd1,nd2,nt_orig,nsav
    real(r8), intent(out) :: zsav(nd1,nd2,nsav)
    integer, intent(in) :: indx(nt_orig)
    integer, intent(in) :: ktdim             ! specifies which dim. is "time"
 
    character*(*), intent(in) :: treename
    character*(*), intent(in) :: nodename
    character*(*), intent(out) :: expr
 
    type(ctable), intent(in) :: data_units   ! known units labels & conversions
 
    integer, intent(out) :: istat
 
! local
 
    integer iln,ilt,str_length
    integer idims(3),iret,mds_value,idescr_floatarr,idescr_cstring
    integer i, iok
    character*1 bksl
    integer imatch
    character*20 units_string
 
    real, dimension(:,:,:), allocatable :: zbuf
 
! executable code
 
    call ctable_init

    bksl = char(92)   ! backslash
!
! determine if "time" is 1st or 2nd dimension -- ck 1st dimension...
!
    iln=str_length(nodename)
    ilt=str_length(treename)
 
    if(ktdim.eq.3) then
       idims(1)=nd1
       idims(2)=nd2
       idims(3)=nt_orig
    else if(ktdim.eq.1) then
       idims(1)=nt_orig
       idims(2)=nd1
       idims(3)=nd2
    else
       idims(1)=nd1
       idims(2)=nt_orig
       idims(3)=nd2
    endif
 
    allocate(zbuf(idims(1), idims(2), idims(3)), stat=iok)
    if(iok/=0) then
       write(lunmsg,*) '***failed to allocate zbuf in geq_mds_sub3'
       istat = 9999
       return
    endif
 
    call geq_mk_expr('data(',nodename(1:iln),')',expr)
    istat = mds_value(expr,idescr_floatarr(zbuf,idims,3),iret)
    if(mod(istat,2).eq.1) then
       do i=1,nt_orig
          if(indx(i).ne.0) then
             if(ktdim.eq.3) then
                zsav(1:nd1,1:nd2,indx(i))=zbuf(1:nd1,1:nd2,i)
             else if(ktdim.eq.1) then
                zsav(1:nd1,1:nd2,indx(i))=zbuf(i,1:nd1,1:nd2)
             else
                zsav(1:nd1,1:nd2,indx(i))=zbuf(1:nd1,i,1:nd2)
             endif
          endif
       enddo
       call geq_units_check(expr,zsav,nd1*nd2*nsav,data_units)
    endif
 
    deallocate(zbuf)
 
    return
  end subroutine geq_mds_sub3
!----------------------------------------------------------
  subroutine geq_mk_expr(prefix,nodename,suffix,expr)
!
    use cgeq_paths
!
!  generate a G-EQDISK MDSplus (TDI) expression
!
!    prefix//bksl//'efit_geqdsk:'//nodename//suffix
!
    character*(*), intent(in) :: prefix
    character*(*), intent(in) :: nodename
    character*(*), intent(in) :: suffix
 
    character*(*), intent(out) :: expr
!
    character*1 bksl
    integer str_length,ilpre,ilnod,ilsuf
!
    character*30 ztest
!
!---------------------------
!
    bksl=char(92)   ! backslash
!
    ilpre=str_length(prefix)
    ilnod=str_length(nodename)
    ilsuf=str_length(suffix)
!
    ztest=nodename(1:ilnod)
    call uupper(ztest)
!
    if(ztest.eq.'RBCENT') then
       expr = prefix(1:ilpre)//bksl//trim(cgeq_apath)//nodename(1:ilnod)// &
            & suffix(1:ilsuf)
    else
       expr = prefix(1:ilpre)//bksl//trim(cgeq_gpath)//nodename(1:ilnod)// &
            & suffix(1:ilsuf)
    endif
!
    return
    end
