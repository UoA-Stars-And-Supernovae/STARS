**==DIVIDE.FOR
      SUBROUTINE DIVIDE(J1, J2, J3, J4, J6, K)
      IMPLICIT REAL*8(A-H, O-Z)
      PARAMETER (MAXMSH = 2000)
      COMMON /SOLV  / C(MAXMSH+1,51,51),S(50,151),ER(60),D,KH,NE,N2,N5, 
     :                N6, N10, N12, N14, MM, IW(101), JW(4)
      DIMENSION KM(51), LM(51), NM(51)
      JC = 1
      IF ( J1.GT.J3 .OR. J2.GT.J4 .OR. J6.GT.N12 ) RETURN
      IF ( KH.GE.1 ) CALL PRINTC(K)
      DO J = J1, J3
         NM(J) = 1
         KM(J) = 5
      END DO
 100  CONTINUE
      JC = JC + 1
      DO J = J1, J3
         IF ( KM(J).LT.3 ) KM(J) = KM(J) + 2
         IF ( KM(J).LT.5 ) NM(LM(J)) = 0
      END DO
      VT = 0.0
C Locate most significant remaining row, significant = ratio of largest
C |element|=VM to sum of remaining |elements| = VS
C RJS 9/6/09 -- set an initial ML so that if a most significant element
C isn't found, at least ML still has a value
      ML = J1
C RJS 9/6/09 -- set JM so that it's never zero
      JM = J1
      DO J = J1, J3
         IF ( KM(J).GE.5 ) THEN
            VM = 0.0
            DO LL = J2, J4
               L = LL - J2 + J1
C Find most significant element
               IF ( NM(L).NE.0 ) THEN
                  VX = DABS(S(J,LL))
                  IF ( VX.GE.VM ) VM = VX
                  IF ( VX.EQ.VM ) ML = L
               END IF
               VS = 0.0D0
            END DO
            LM(J) = ML
C Sum remaining elements
            DO LL = J2, J4
               L = LL - J2 + J1
               IF ( L.NE.ML ) VS = VS + DABS(S(J,LL))
            END DO
            IF ( VS.EQ.0.0D0 ) KM(J) = 1
            IF ( VM.EQ.0.0D0 ) KM(J) = 6
            IF ( VS.EQ.0.0D0 .AND. VM.NE.0.0D0 ) NM(ML) = 0
            IF ( VM.EQ.0.0D0 ) VX = 0.0D0
            IF ( VS.EQ.0.0D0 .AND. VM.GT.0.0D0 ) VX = 2.0D0
            IF ( VS.GT.0.0D0 ) VX = VM/(VM+VS)
            IF ( VX.GE.VT ) VT = VX
            IF ( VX.EQ.VT ) JM = J
         END IF
      END DO
      IF ( KM(JM).EQ.5 ) KM(JM) = 2
C Eliminate elements above and below largest element of most significant row
      DO I = J1, J3
         IM = KM(I)
         IF ( IM.LT.3 ) THEN
            ML = LM(I) + J2 - J1
            LA = J2
            IF ( IM.EQ.1 ) LA = J6
            VX = 1.0D0/S(I, ML)
C What does D do? It is set to zero in each loop, but doesn't seem to be used.
            D = D - DLOG(DABS(VX))
            DO L = LA, N12
               S(I, L) = VX*S(I, L)
            END DO
            S(I, ML) = 1.0D0
            DO J = J1, J3
               IF ( KM(J).GT.3 ) THEN
                  VX = S(J, ML)
                  DO L = LA, N12
                     S(J, L) = S(J, L) - VX*S(I, L)
                  END DO
                  S(J, ML) = 0.0D0
               END IF
            END DO
         END IF
         NM(I) = 1
      END DO
      IF (JC.GT.10000) GOTO 200
      DO J = J1, J3
         IF ( KM(J).EQ.5 .OR. KM(J).LE.2 ) GO TO 100
      END DO
 200  DO L = J6, N12
         DO J = J1, J3
            C(K, LM(J), L-N12+N14) = S(J, L)
         END DO
      END DO
      IF ( KH.EQ.0 ) RETURN
      WRITE(24, 99001) K, KM, LM, NM
99001 FORMAT (32I3)
      M = J6 - N12 + N14
      DO J = 1, NE
         WRITE(24, 99002) (C(K,J,L), L=M, N14)
      END DO
99002 FORMAT (1P, 10E13.5)
      RETURN
      END
