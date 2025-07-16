!  record of revisions...
! may,12,1982 changed m1 to mf and m2 to ml to avoid conflicts.
!  fixed the writes in saver to treat multiple surfaces correctly.
!
! modified by aplet for up down asymmetric case.  13.12_r8 .93
! need to modify eigenvalue determination later on.
!
!rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
!******************************************************************
      SUBROUTINE pstent23
!******************************************************************
!      entry point for overlay 23
!
!
 USE pstcom

 USE combla

 USE l34com
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER ISET
      INTEGER IOK
      INTEGER IS
      INTEGER N1


!
      data iset /0/
!
      ipolp2 = ipolp
      nfmg   = nfm
      nsfg   = nsf
      nfeg   = nfe
      nfmga  = ipolp2 * nfmg+1
      nfmgb  = ipolp2 * 2 * nfmga
      nad2   = nfmgb *( nfmgb + 1 )/ 2
      nmtg2  = nfmga * nfeg
      nfmgb2 = nfmgb
      nfmgc  = 3*nfmg
!
 ALLOCATE( a(nad2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate a(nad2)'
    stop 'ERROR: cannot allocate in ent23'
 endif
!
 b => a
!
 ALLOCATE( x(nmtg2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate c(nmtg2)'
    stop 'ERROR: cannot allocate in ent23'
 endif
!
  vb => x
 
   c => x

 ALLOCATE( ut(nmtg2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate ut(nmtg2)'
    stop 'ERROR: cannot allocate in ent23'
 endif
!
 ALLOCATE( vt(nmtg2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate vt(nmtg2)'
    stop 'ERROR: cannot allocate in ent23'
 endif
!
 ALLOCATE( xt(nmtg2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate xt(nmtg2)'
    stop 'ERROR: cannot allocate in ent23'
 endif
!
 ALLOCATE( u(nfmgb2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate u(nfmgb2)'
    stop 'ERROR: cannot allocate in ent23'
 endif
!
 ALLOCATE( va(nfmgb2), stat=iok )
 if (iok /= 0 ) then
    print *,'ERROR: cannot allocate va(nfmgb2)'
    stop 'ERROR: cannot allocate in ent23'
 endif

 if ( .not. lmarg  ) then
    ALLOCATE( ashift((nfe+1)*nfm*(2*nfm+1)), stat=iok )
    if (iok /= 0 ) then
       print *,'ERROR: cannot allocate ashift((nfe+1)*nfm*(2*nfm+1))'
       stop 'ERROR: cannot allocate in ent23'
    endif
 end if

!
      iset = iset + 1
!
!      assign disk files ...
!
      is = 1
      n1 = 2  *  jsub(1)
      n1 = n1 * ( n1 + 1 )
      n1 = n1 / 2
!
      if ( iset  >  1 )  go to 10
!
   10 continue
!
!      solve the eigenvalue problem
!
      if ( .not. lmarg ) CALL psteigens
!
!      solve for delta'..
!
      if ( lmarg )  CALL pstdeltap
!
 deallocate(a, stat = iok)
 if ( iok /= 0 ) then
    print *,' ERROR: cannot deallocate a '
    stop 'ERROR: cannot allocate in ent23'
 end if
!
 deallocate(x, stat = iok)
 if ( iok /= 0 ) then
    print *,' ERROR: cannot deallocate x '
    stop 'ERROR: cannot allocate in ent23'
 end if
!
 deallocate(ut, stat = iok)
 if ( iok /= 0 ) then
    print *,' ERROR: cannot deallocate ut '
    stop 'ERROR: cannot allocate in ent23'
 end if
!
 deallocate(vt, stat = iok)
 if ( iok /= 0 ) then
    print *,' ERROR: cannot deallocate vt '
    stop 'ERROR: cannot allocate in ent23'
 end if
!
 deallocate(xt, stat = iok)
 if ( iok /= 0 ) then
    print *,' ERROR: cannot deallocate xt '
    stop 'ERROR: cannot allocate in ent23'
 end if
!
 deallocate(u, stat = iok)
 if ( iok /= 0 ) then
    print *,' ERROR: cannot deallocate u '
    stop 'ERROR: cannot allocate in ent23'
 end if
!
 deallocate(va, stat = iok)
 if ( iok /= 0 ) then
    print *,' ERROR: cannot deallocate va '
    stop 'ERROR: cannot allocate in ent23'
 end if
!
 if ( .not. lmarg  ) then
    deallocate( ashift, stat=iok )
    if (iok /= 0 ) then
       print *,'ERROR: cannot deallocate ashift '
       stop 'ERROR: cannot allocate in ent23'
    endif
 end if

!
   99 return
 7000 CALL psterrmes(outpst,'ent23',is)
      end


