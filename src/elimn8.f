**==ELIMN8.FOR
      SUBROUTINE ELIMN8(J2, J3, J4, J5, J10, J12, KG)
      IMPLICIT REAL*8(A-H, O-Z)
      PARAMETER (MAXMSH = 2000)
      COMMON /SOLV  / C(MAXMSH+1,51,51),S(50,151),ER(60),D,KH,NE,N2,
     :                N5, N6, N10, N12, N14, MM, IW(101), JW(4)
      IF ( KH.GE.1 ) CALL PRINTC(KG)
      IF ( J2.GT.J4 .OR. 1.GT.J3 .OR. J10.GT.J12 ) RETURN
      NN = J10 - 1
      NP = N12 - N14
      DO K = J2, J4
         M = K - J2 + J5
         DO I = 1, J3
            VX = S(I, K)
            DO J = J10, J12
               S(I, J) = S(I, J) - VX*C(KG, M, J-NN)
            END DO
            DO J = N10, N12
               S(I, J) = S(I, J) - VX*C(KG, M, J-NP)
            END DO
         END DO
      END DO
      RETURN
      END
