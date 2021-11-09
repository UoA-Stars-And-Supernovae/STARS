**==FDIRAC.F
      SUBROUTINE FDIRAC (F, G)
      IMPLICIT REAL*8(A-H, O-Z)
      COMMON /STATFD/ D(3,3), DTT, DFT, DFF, DTTT, DFTT, DFFT, DFFF
      DIMENSION FF(4), GG(4), C(4,4,3), VW(4,4), VX(4,4)
      DATA NC, C/4, 2.315472, 7.128660, 7.504998, 2.665350, 7.837752, 
     :     23.507934, 23.311317, 7.987465, 9.215560, 26.834068, 
     :     25.082745, 8.020509, 3.693280, 10.333176, 9.168960, 2.668248, 
     :     2.315472, 6.748104, 6.564912, 2.132280, 7.837752, 21.439740, 
     :     19.080088, 5.478100, 9.215560, 23.551504, 19.015888, 
     :     4.679944, 3.693280, 8.859868, 6.500712, 1.334124, 1.157736, 
     :     3.770676, 4.015224, 1.402284, 8.283420, 26.184486, 28.211372, 
     :     10.310306, 14.755480, 45.031658, 46.909420, 16.633242, 
     :     7.386560, 22.159680, 22.438048, 7.664928/
* Evaluate Fermi-Dirac integrals (Eggleton, Faulkner & Flannery 1973).
* 1st row of matrix D contains rho*, P* and Q*; the 2nd row the logarithmic
* derivatives w.r.t. T; the 3rd row logarithmic derivatives w.r.t. f
      VF = 1.0/(1.0+F)
      VG = 1.0/(1.0+G)
      UF = F*VF
      UG = G*VG
      FDF = G*G + G
      FDF = UF*FDF*SQRT(FDF)
      FF(1) = 1.0
      GG(1) = 1.0
      DO I = 2, NC
         FF(I) = F*FF(I-1)
         GG(I) = G*GG(I-1)
         FDF = FDF*VF*VG
      END DO
      ND = 4
      DO I = 1, 3
         IF(I.GT.1) ND = 2
         DO L = 1, ND
            DO M = 1, ND + 1 - L
               VX(L,M) = 0.0
            END DO
         END DO
         DO J = 1, NC
            DO K = 1, NC
               VW(1,1) = C(K,J,I)*GG(J)*FF(K)
               DO L = 1, ND - 1
                  VW(L+1,1) = (J-1)*VW(L,1)
                  DO M = 1, ND - L
                     VW(L,M+1) = (K-1)*VW(L,M)
                  END DO
               END DO
               DO L = 1, ND
                  DO M = 1, ND + 1 - L
                     VX(L,M) = VX(L,M) + VW(L,M)
                  END DO
               END DO
            END DO
         END DO
         WV = 1.0D0/VX(1,1)
         DO L = 1, ND
            DO M = 1, ND + 1 - L
               IF(L.NE.1 .OR. M.NE.1) VX(L,M) = VX(L,M)*WV
            END DO
         END DO
         D(I,1) = FDF*VX(1,1)
         D(I,2) = VX(2,1) + 1.5 + (2.5 - NC)*UG
         D(I,3) = VX(1,2) + 1.0 + (0.5*D(I,2) - NC)*UF
         IF(I.EQ.1) THEN
* second and third order derivatives of density, needed in PRESSI
* !!! Programming should be checked!!!
            DTT = VX(3,1) - VX(2,1)**2 + (2.5 - NC)*UG*VG
            DFT = VX(2,2) - VX(2,1)*VX(1,2) + 0.5*UF*DTT
            DFF = VX(1,3) - VX(1,2)**2 + UF*(DFT + (0.5*D(I,2) - NC)*VF
     &           - 0.25*DTT*UF)
            DTTT = VX(4,1) + VX(2,1)*(2.0*VX(2,1)**2 - 3.0*VX(3,1)) +
     &           (2.5 - NC)*(1.0-G)*UG*VG**2
            DFTT = VX(3,2) - 2.0*VX(2,1)*VX(2,2) - VX(1,2)*(VX(3,1) -
     &           2.0*VX(2,1)**2) + 0.5*UF*DTTT
            DFFT = VX(2,3) - 2.0*VX(2,2)*VX(1,2) - VX(2,1)*(VX(1,3) -
     &           2.0*VX(1,2)**2) + UF*(DFTT + 0.5*VF*DTT - 0.25*UF*DTTT)
            DFFF = VX(1,4) + VX(1,2)*(2.0*VX(1,2)**2 - 3.0*VX(1,3)) +
     &           UF*(1.5*(DFFT + VF*DFT) - UF*(0.75*(DFTT + VF*DTT) - 
     &           0.125*UF*DTTT) + (0.5*D(1,2) - NC)*(1.0-F)*VF**2)
         END IF
      END DO
      RETURN
      END
