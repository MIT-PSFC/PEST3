C-----------------------------------------------------------------------
C  PARSID -- PARSE OUT TRANSP RUN ID AND TOKAMAK FROM COMBINED INPUT
C
C History:
C 10/16/96 CAL : Support RunId with 4 - 6 digit shot number
C                                   2 - 3 digit run  number
C
      SUBROUTINE PARSID(STRING,NRUN,TOK)
C
      CHARACTER*(*) STRING	! INPUT STRING E.G. "12345A01PBXM"
C				! COULD BE OLD STYLE "1234TFTR"
C                               ! e.g. 123456A123TFTR
      CHARACTER*(*) NRUN	! OUTPUT RUN ID, "12345A01" OR "1234"
      CHARACTER*(*)  TOK	! OUTPUT TOKAMAKD ID, "PBXM" OR "TFTR"
C
      CHARACTER*10 ZDIGIT
C
C  TWO RUN ID "STYLES" ARE SUPPORTED, THE OLD 4 DIGIT RUN NUMBER STYLE
C  AND THE NEW 8 CHARACTER nnnnnAmm SHOT-TRY NOMENCLATURE
C
C------------------------------------
C
      TOK='????'
C
      if(len(nrun).lt.10) then
         ilen=len(nrun)
         write(6,9901) ilen
 9901    format(
     >       ' ?parsid:  len(rundid)=',i2,
     >       ' incompatible with 6 digit shot number.')
         nrun='len_error'
         return
      endif
C
      ZDIGIT='0123456789'
C
      NRUN='ERROR!!!'
C
      ILEN=INDEX(STRING,' ')-1
      IF(ILEN.LE.0) ILEN=LEN(STRING)
C
C  TEST FOR OLD VS. NEW STYLE IS WHETHER THE 5TH CHARACTER OF THE INPUT
C  STRING IS A DIGIT:  DIGIT INDICATES NEW STYLE
C
      IF(ILEN.LT.5) GO TO 1000
C
      IDIG=INDEX(ZDIGIT,STRING(5:5))
C
      IF(IDIG.EQ.0) THEN
C
C  OLD STYLE
C
C  CHECK LENGTH
C
        IF((ILEN.LT.7).OR.(ILEN.GT.8)) GO TO 1000
C
C  OK, CHECK THAT 1ST FOUR CHARS ARE DIGITS
C
        DO 10 IC=1,4
          IF(INDEX(ZDIGIT,STRING(IC:IC)).EQ.0) GO TO 1000
 10     CONTINUE
C
C  GOOD
C
        NRUN=STRING(1:4)
        TOK=STRING(5:ILEN)
C
      ELSE
C
C  NEW STYLE
C
C  CHECK valid Shot Number
C
        ic = 0
        DO 20 while (idig .ne. 0)
          ic = ic+1
          IDIG = INDEX(ZDIGIT,STRING(IC:IC))
 20     CONTINUE
        i1 = ic - 1
        if (i1 .lt. 5 .or. i1 .gt. 6) go to 1000
          if (i1 .eq. 6 .and. string(1:1) .eq. '0') go to 1000
 
C
C Check valid Run Number
          idig = 1
        DO 21 while (idig .ne. 0)
          ic = ic+1
          IDIG = INDEX(ZDIGIT,STRING(IC:IC))
 21     CONTINUE
        i2 = ic - i1 - 2
        if (i2 .lt. 2 .or. i2 .gt. 3)  go to 1000
 
        NRUN=STRING(1:ic-1)
        TOK=STRING(ic:ILEN)
C
      ENDIF
C
 1000 CONTINUE
      RETURN
      END
