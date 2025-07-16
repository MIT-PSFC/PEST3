C-----------------------------------------------------------------------
C  CATFUSLIS -- PROVIDE DATA ON FUSION BY FUSION PRODUCT
C
C  modified (dmc Feb 2002 -- catfuslis.for -- now based on arguments
C  instead of COMMON).  Coordinate any changes with comput/catfus*.*
C
C  FUSION REACTION INDEXING IS AS DESCRIBED IN COMMENTS OF SUBROUTINE
C  FUSION (***** see nubeam/fusion.for comments *****).
C  see also comput/catfus*.*
C
C  Caution:  Feb 2002 (dmc):  although catfus has a notion of a "D"
C  RF minority, it is not supported here.
C
C  THIS ROUTINE RECEIVES AS ARGUMENT THE INDEX TO A FAST ION SPECIES
C  USED TO MODEL FUSION PRODUCT FAST IONS.  WHICH SPECIES INDICATES
C  THE FUSION PRODUCT AND THUS ONE OR MORE FUSION REACTIONS INVOLVED.
C
C  THIS ROUTINE GIVES AN ESTIMATE OF THE INITIAL ENERGY OF THE FUSION
C  PRODUCT ** DMC SEPT 1990 ** THIS SHOULD BE CONSISTENT WITH THE MAX
C  ENERGY SET IN S.R. DATCHK FOR THE FUSION PRODUCT SPECIES' COMPUTED
C  DISTRIBUTION FUNCTION.
C
C  ENERGY DATA IS FROM THE NRL PLASMA FORMULARY ... DMC SEPT 1990
C
      SUBROUTINE CATFUSLIS_r8(IDENT,ILREACT,ZEREACT,idim,
     >   nlftrit,nlfhe3,nlrftrit,nlrfhe3,ierr)
C
      use catfus_mod
      implicit NONE
C
C  input:
      character*(*) ident               ! fusion product ID
C             "P" "T" "He3" "He4" "Alpha"
C
      integer idim                      ! dimension of output arrays
C
      logical nlftrit                   ! TRUE if fusion product T present
      logical nlfhe3                    ! TRUE if fusion product He3 present
C
      logical nlrftrit                  ! TRUE if rf minority T present
      logical nlrfhe3                   ! TRUE if rf minority He3 present
C
C  OUTPUT:
C
C  SET .TRUE. FOR REACTIONS THAT PRODUCE A DESIRED FUSION PRODUCT:
      LOGICAL ILREACT(idim)
C
C  INITIAL ENERGY (EV) OF FUSION PRODUCT SO PRODUCED:
      REAL*8 ZEREACT(idim)
C
      integer ierr                      ! completion code: 0=OK
C
C-------------------------------------------------------------
      character*6 test
C-------------------------------------------------------------
C
C  INPUT:
C
C    IDENT = FUSION PRODUCT INDEX (MUST BE NON-ZERO)
C
      test=ident
      call uupper(test)                 ! uppercase conversion
C
      ierr=0
      ilreact=.FALSE.
      zereact=0.0d0
C
      if(idim.lt.nreact) then
         write(6,*) ' ?CatFusLis:  idim = ',idim,' too small.'
         write(6,*) '  idim .ge. nreact = ',nreact,' was expected.'
         ierr=1
         return
      endif
C
      IF(test.eq.'T') THEN
C
C  TRITIUM COMES FROM THE D+D REACTION (PROTON BRANCH)
C
         ILREACT(3)=.TRUE.
         ZEREACT(3)=1.01D6
C
      Else if (test.eq.'P') then
C
C  Proton COMES FROM THE D+D REACTION (PROTON BRANCH)
C
         ilreact(3) = .true.
         zereact(3) = 3.02d6
C
C  D+HE3 AND HE3+D
 
         ilreact(2) = .true.
         zereact(2) = 14.7d6
         ilreact(8) = .true.
         zereact(8) = 14.7d6
 
         IF(NLFHE3) THEN
C  fusion product He3 -- He3+D burnup reactions
            ILREACT(10)=.TRUE.
            ZEREACT(10)=14.7d6          ! proton energy
         ENDIF
C
         IF(NLRFHE3) then
C  RF minority He3
C  He3 + D target
            ILREACT(12)=.TRUE.
            ZEREACT(12)=14.7d6          ! proton energy
         ENDIF
C
      ELSE IF(test.eq.'HE3') THEN
C
C  HELIUM-3 COMES FROM THE D+D REACTION (NEUTRON BRANCH)
C
         ILREACT(4)=.TRUE.
         ZEREACT(4)=0.82D6
C
      ELSE IF((test.eq.'HE4').or.(test.eq.'ALPHA')) THEN
C
C  HELIUM-4 COMES FROM D+T, D+HE3, AND T+T REACTIONS
C
C  D+T AND T+D --
         ILREACT(1)=.TRUE.
         ZEREACT(1)=3.5D6
         ILREACT(7)=.TRUE.
         ZEREACT(7)=3.5D6
C
C  D+HE3 AND HE3+D
         ILREACT(2)=.TRUE.
         ZEREACT(2)=3.6D6
         ILREACT(8)=.TRUE.
         ZEREACT(8)=3.6D6
C
C  T+T (ACTUALLY I DON'T KNOW WHAT ENERGY TO GIVE ...)
C
         ILREACT(5)=.TRUE.
         ZEREACT(5)=3.75D6
C
         IF(NLFTRIT) THEN
C  fusion product T -- T+D burnup reactions
            ILREACT(9)=.TRUE.
            ZEREACT(9)=3.5D6
         ENDIF
C
         IF(NLFHE3) THEN
C  fusion product He3 -- He3+D burnup reactions
            ILREACT(10)=.TRUE.
            ZEREACT(10)=3.6D6
         ENDIF
C
         IF(NLRFTRIT) THEN
C  RF minority T -- T+D burnup reactions
            ILREACT(11)=.TRUE.
            ZEREACT(11)=3.5D6
         ENDIF
C
         IF(NLRFHE3) THEN
C  RF minority He3 -- He3+D burnup reactions
            ILREACT(12)=.TRUE.
            ZEREACT(12)=3.6D6
         ENDIF
C
      else
         ierr=2
         write(6,*) ' ?CatFusLis:  ident="',ident,'" not recognized.'
      ENDIF
C
      RETURN
      END
C--------------------------------------------
C  REAL interface
C
      subroutine catfuslis(IDENT,ILREACT,ZEREACT,idim,
     >   nlftrit,nlfhe3,nlrftrit,nlrfhe3,ierr)
C
      implicit NONE
C
C  input:
      character*(*) ident               ! fusion product ID
C             "P" "T" "He3" "He4" "Alpha"
C
      integer idim                      ! dimension of output arrays
C
      logical nlftrit                   ! TRUE if fusion product T present
      logical nlfhe3                    ! TRUE if fusion product He3 present
C
      logical nlrftrit                  ! TRUE if rf minority T present
      logical nlrfhe3                   ! TRUE if rf minority He3 present
C
C  OUTPUT:
C
C  SET .TRUE. FOR REACTIONS THAT PRODUCE A DESIRED FUSION PRODUCT:
      LOGICAL ILREACT(idim)
C
C  INITIAL ENERGY (EV) OF FUSION PRODUCT SO PRODUCED:
      REAL ZEREACT(idim)
C
      integer ierr                      ! completion code: 0=OK
C
C-------------------------------------------------
      real*8 zereact_r8(idim)
C-------------------------------------------------
C
      zereact_r8=0.0d0
C
      call catfuslis_r8(IDENT,ILREACT,ZEREACT_R8,idim,
     >   nlftrit,nlfhe3,nlrftrit,nlrfhe3,ierr)
 
      zereact = zereact_r8
C
      return
      end
