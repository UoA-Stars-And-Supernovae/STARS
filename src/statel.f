**==STATEL.FOR
      SUBROUTINE STATEL(LL, FL, TL, ISTAR)
      IMPLICIT REAL*8(A-H, O-Z)
      COMMON /STAT2 / VF(60)
      DIMENSION IL(24), SF(60,2)
      DATA IL/2, 0, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 
C     :              2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2/
     :     2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2/
C Subtract 15 from L if using 2 stars
      IF (LL.GT.15) THEN
         L = LL - 15
      ELSE
         L = LL
      END IF
      I = IL(L+2)
      IF ( I.NE.1 ) THEN
         CALL STATEF(FL, TL)
         CALL NUCRAT(TL)
         IF ( I.GT.0 ) RETURN
         DO J = 1, 60
            SF(J,ISTAR) = VF(J)
         END DO
         RETURN
      END IF
      DO J = 1, 60
         VF(J) = SF(J,ISTAR)
      END DO
      RETURN
      END
