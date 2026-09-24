!     Writes debug information to fort.24
!     run_bs calls this "useless debugging information" so should we keep it?
      SUBROUTINE PRINTC(K)

      IMPLICIT NONE

      REAL*8 ER, D, SOLV, S, C
      INTEGER N6, N12, J, JW, I, IW, K, KH
      INTEGER N14, N2, MM, I5, N5, N10, NE, MAXMSH

      PARAMETER (MAXMSH = 2000)
      COMMON /SOLV  / C(MAXMSH+1,51,51),S(50,151),ER(60),D,KH,NE,N2,
     :                N5, N6, N10, N12, N14, MM, IW(101), JW(4)
      IF ( MM.NE.2 ) THEN
            DO I = 1, NE
                  WRITE(24, 99001) K, (S(I,J), J=1, NE)
            END DO

            WRITE(24, 99001) K
      END IF

      DO I = 1, NE
            WRITE(24, 99001) K, (S(I,J), J=N2, N5)
      END DO

      WRITE(24, 99001) K

      DO I = 1, NE
            WRITE(24, 99001) K, (S(I,J), J=N6, N12)
      END DO

      WRITE(24, 99001) K

99001 FORMAT (I5, 1P, 11E11.3)
      RETURN
      END
