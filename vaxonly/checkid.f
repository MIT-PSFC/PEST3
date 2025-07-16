C-----------------------------------------------------------------------
C  CHECKID -- PARSE OUT TRANSP RUN ID AND TOKAMAK FROM COMBINED INPUT
C
C History:
C 10/16/96 CAL : Support RunId with 4 - 6 digit shot number
C                                   2 - 3 digit run  number
C 01/07/97 CAL : Fix lrid
C 01/06/98 CAL : Fix extracting shot number
C
C
      SUBROUTINE CHECKID(STRING, lrid, lshot)
C
      implicit NONE
C
      CHARACTER*(*) STRING	! INPUT STRING E.G. "12345A01PBXM"
C				! COULD BE OLD STYLE "1234TFTR"
C                               ! e.g. 123456A123TFTR
      INTEGER  lrid, lshot              ! return values: length of RunId/Shot
C
      CHARACTER*10 ZDIGIT
C
C-------------------------
      integer ilen,idig,ic,i1,i2
C-------------------------
C
C  TWO RUN ID "STYLES" ARE SUPPORTED, THE OLD 4 DIGIT RUN NUMBER STYLE
C  AND THE NEW 8 CHARACTER nnnnnAmm SHOT-TRY NOMENCLATURE
C
C------------------------------------
C
      lshot = 0
      ZDIGIT='0123456789'
C
C
      ILEN=INDEX(STRING,' ')-1
      IF(ILEN.LE.0) ILEN=LEN(STRING)
C	type 1,  string(1:ilen), ilen
C 1	format(1x,a20,i3)
      if (ilen .eq. 4 ) go to 1200
      if (ilen .eq. 5 ) go to 1300
C
C  TEST FOR OLD VS. NEW STYLE IS WHETHER THE 5TH CHARACTER OF THE INPUT
C  STRING IS A DIGIT:  DIGIT INDICATES NEW STYLE
C
C
C  CHECK valid Shot Number
C
          idig = 1
        ic = 0
        DO 20 while (idig .ne. 0 .and. ic .lt. ilen)
          ic = ic+1
          IDIG = INDEX(ZDIGIT,STRING(IC:IC))
 20     CONTINUE
        i1 = ic - 1
C	  if (i1 .eq. 4 .and. i1 .eq. ILEN) go to 1100
        if (ic .eq. 4 ) go to 1200
        if (i1 .eq. 4 ) go to 1100
        if (i1 .lt. 5 .or. i1 .gt. 6) go to 1000
          if (i1 .eq. 6 .and. string(1:1) .eq. '0') go to 1000
 
C
C Check valid Run Number
        idig = 1
        DO 21 while (idig .ne. 0 .and. ic .lt. ilen)
          ic = ic+1
          IDIG = INDEX(ZDIGIT,STRING(IC:IC))
 21     CONTINUE
        if (ic .eq. ilen .and. idig .ne. 0) ic = ic+1
        i2 = ic - i1 - 2
        if (i2 .lt. 2 .or. i2 .gt. 3)  go to 1000
C
 1100 lrid = ic - 1
      lshot = i1
        return
 
 1200 lrid = 4
      lshot = 4
        return

 1300 lrid = 5
      lshot = 5
        return
 
 1000 lrid = 0
      RETURN
      END
 
 
 
 
 
