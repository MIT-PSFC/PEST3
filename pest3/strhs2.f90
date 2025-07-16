!
      subroutine pststrhs2(kms,kparit)
!
! load right-hand side into xt...
!
 USE pstcom

 USE combla

 USE comggf
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER KMS
      INTEGER KPARIT
      INTEGER JBEG
      INTEGER J1
      INTEGER J
      INTEGER NTYPE
      real*8 ALAM
      INTEGER NEG
      INTEGER NIT


!
! kms # of rational surface
! kparit=1,2 parity of the solution: kparit=1 odd, kparit=2 even
! case ibas.ne.1 not implemented.
!
      if(ibas /= 1) then
      write(itty,*)'ibas.ne.1 not implemented in pststrhs2'
      stop 'ERROR in strhs: Only ibas=1 is supported'
      end if
!
      xt =  0.0_r8 
!
      mbeg = msing(kms) + 1
      mend = msing(kms+2) - 1
!
      if( kparit == 1 ) then

         do m1 = mbeg,mend
            jbeg = jtot(m1)
            jmax1 = lmax(1) - lmin(1) + 1
            do j1 = 1,jmax1
               j = jbeg + (j1-1)*ibas
               xt(j) = w0l1ou(m1,j1,kms)
            end do
         end do
         !
      else
         do m1 = mbeg,mend
            jbeg = jtot(m1)
            jmax1 = lmax(1) - lmin(1) + 1
            do j1 = 1,jmax1
               j = jbeg + (j1-1)*ibas
               xt(j) = w0l1eu(m1,j1,kms)
            end do
         end do
         !
      end if
if (check1) then
   print *,'in pststrhs2 kparit = ', kparit
   print *,xt
end if

      return
      end subroutine pststrhs2
!
      SUBROUTINE pstivi(ntype,alam,neg,nit)
!--------------------------------------
! up-down compatibility check ap  17.10_r8 .95
!
 USE pstcom
 IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER KMS
      INTEGER KPARIT
      INTEGER JBEG
      INTEGER J1
      INTEGER J
      INTEGER NTYPE
      real*8 ALAM
      INTEGER NEG
      INTEGER NIT


 REAL*8, DIMENSION(10) 	 :: eigvag
 REAL*8, DIMENSION(3) 	 :: al
 REAL*8	 	 	 :: epscon
 REAL*8	 	 	 :: epsmac
 INTEGER, DIMENSION(14) 	 :: npin
 INTEGER, DIMENSION(4) 	 :: npout
 COMMON /c2f90eigvag/ eigvag 
 COMMON /c2f90al/ al 
 COMMON /c2f90epscon/ epscon 
 COMMON /c2f90epsmac/ epsmac 
 COMMON /c2f90npin/ npin 
 COMMON /c2f90npout/ npout 
      al(1)=alam
      npin(1) = ntype
      CALL pstvekit(al,epscon,epsmac,npin,npout)
      neg = npout(1)
      nit = npout(2)
      return
      end


