	SUBROUTINE PARSE_LINE(LINE, NAME, VALUE, EQ_POS)
      CHARACTER*100 LINE
      CHARACTER*20 NAME
      CHARACTER*20 VALUE
      INTEGER EQ_POS
      INTEGER I

      NAME = ' '
      VALUE = ' '
      EQ_POS = 0

      DO 30 I = 1, 100
          IF (LINE(I:I) .EQ. '=') THEN
              EQ_POS = I
              GOTO 40
          END IF
30    CONTINUE

40    CONTINUE
      IF (EQ_POS .GT. 0) THEN
          NAME = LINE(1:EQ_POS-1)
          VALUE = LINE(EQ_POS+1:)
          CALL TRIM_STR(NAME)
          CALL TRIM_STR(VALUE)
      END IF
      RETURN
      END

      SUBROUTINE TRIM_STR(STR)
      CHARACTER*(*) STR
      INTEGER LEN, I

      LEN = LEN_TRIM(STR)
      IF (LEN .EQ. 0) RETURN

      DO 50 I = 1, LEN
          IF (STR(I:I) .NE. ' ') THEN
              STR = STR(I:)
              RETURN
          END IF
50    CONTINUE
      RETURN
      END

      INTEGER FUNCTION LEN_TRIM(STR)
      CHARACTER*(*) STR
      INTEGER I

      LEN_TRIM = LEN(STR)
      DO 60 I = LEN(STR), 1, -1
          IF (STR(I:I) .NE. ' ') THEN
              LEN_TRIM = I
              RETURN
          END IF
60    CONTINUE
      LEN_TRIM = 0
      RETURN
      END
