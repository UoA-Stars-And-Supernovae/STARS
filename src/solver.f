      SUBROUTINE SOLVER(NOC, ITER, JTER, ERR, IE, NWRT5)

      IMPLICIT NONE

      REAL*8 ER, ELSEIF, SGTHFAC, TOTMP, VX, RLFC, DABS, HAS
      REAL*8 FACSG, BECOME, TRC2, VLHP, HNUC, TOLERANCE, ABS, TRC1
      REAL*8 BOOSTING, DLOG, EXCEEDED, FRR, ERMAX, VZ, VLCP, RESTARTING
      REAL*8 DH, VLEC, C, VLCC, FAC, RLFP, GT, DTCONT
      REAL*8 AT, EPS, ERT, S, EQ, VARIABLE, F7
      REAL*8 DH0, SNAFU, DT2, DMIN1, DEL, VLEP, DHNUC, ERRPR
      REAL*8 DD, ERR, D, VLHC, SOLV, DT1, GE, H
      REAL*8 DMAX1, FACSGMIN, TOTMC, STAR
      INTEGER JIN, N12, J2, NB, N8, KMAX, NE2, IH
      INTEGER KK, N14, K1, MWT, KMESH, N5, IN, N10
      INTEGER N6, NE3, I4, IA, K2, MIXFUD, KTER, MM
      INTEGER JMAX, NOC, NEV, N1, NAN, NE, J1, NUCLEOSYNTHESIS
      INTEGER NV, NMESH, IE, IW, N16, NE4, K, KH
      INTEGER MIXING, ID, N7, LE, NWRT5, NW, N3, MAX
      INTEGER N11, NN, NUCMAT, MAXMSH, I2, ISGFAC, J, N4
      INTEGER N15, ITER, NF, N9, I, L, MWTS, NE1
      INTEGER N13, MIN0, MK, N2, ISTAR, MESH, IVMC, JTER
      INTEGER LT, IVMS, JH, MESHPOINT

      PARAMETER (MAXMSH = 2000)
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,DEL,DH0,NMESH,JIN,IW(200)
      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      COMMON /SOLV  / C(MAXMSH+1,51,51),S(50,151),ER(60),D,KH,NE,N2,N5,
     :                N6, N10, N12, N14, MM, IA(1), NE1, NE2, NE3, NB, 
     :                NEV, NF, J1, J2, IH, JH, ID(90), K1, K2, NV, NW
      COMMON /DTCONT/ VLHP(2), VLEP(2), VLCP(2), RLFP(2), TOTMP(2), VLHC(2),
     :     VLEC(2), VLCC(2), RLFC(2), TOTMC(2)
      COMMON /MIXFUD/ SGTHFAC, FACSGMIN, FACSG, ISGFAC
      COMMON /MESH  / TRC1,TRC2,DD,DT1,DT2,MWT,MWTS,IVMC,IVMS
      DIMENSION MK(60), ERT(60), IE(100)
      LOGICAL REDUCE
      REDUCE = .FALSE.
C EG-style mixing convergence trick
      IF (NOC.EQ.1.AND.ISGFAC.EQ.1) THEN
         FACSG = FACSGMIN
      ELSE
         FACSG = 1d0
      END IF
C 1/9/09 RJS -- the following section of code has been removed. It was originally there
C to stop problems that occurred during matrix inversion due to the use of some elements
C from previous iterations (i.e. the nucleosynthesis). This massively slowed the code in
C non-nucleosynthesis mode, and no longer seems to be a problem. If you find you need it,
C add it back in again :)
C Clear C and S
C        DO J = 1, 51
C           DO I = 1,51
C              DO K = 1,NMESH+1 !MAXMSH+1
C                 C(K,I,J) = 0d0
C              END DO
C           END DO
C        END DO
C        DO J = 1,151
C           DO I = 1,50
C              S(I,J) = 0d0
C           END DO
C        END DO
C Replace loop with written out stuff
      NE1 = IE(1)
      NE2 = IE(2)
      NE3 = IE(3)
      NB = IE(4)
      NEV = IE(5)
      NF = IE(6)
      J1 = IE(7)
      J2 = IE(8)
      IH = IE(9)
      JH = IE(10)
      DO I=1,90 !45
         ID(I) = IE(10+I)
      END DO
      NE = NE1 + NE2
      ERR = 0.0D0
      IF ( NE.EQ.0 ) RETURN
      NV = NE + NEV
      NW = NV + NV
      N1 = NB + 1
      N2 = NE + 1
      N3 = NE + NB
      N4 = N3 + 1
      N5 = NE + NE
      N6 = N5 + 1
      N7 = N5 + NB
      N8 = N7 + 1
      N9 = N5 + NE
      N10 = N9 + 1
      N11 = N9 + NEV
      N12 = N10 + NEV
      N13 = NE - NB
      N15 = N13 + 1
      N16 = N13 + NEV
      N14 = N16 + 1
      NE4 = NE3 + NE2
C SURFACE IS AT K1'TH MESHPOINT, CENTRE AT K2'TH
      K1 = 1 + 0.01D0*J1*(NMESH-1)
      K2 = NMESH - 0.01D0*J2*(NMESH-1)
      KMESH = K2 - K1 + 1
C DETERMINE 'TYPICAL' VALUES FOR EACH  VARIABLE
      IF (NOC.LE.1) THEN
         DO J = 1, NV
            L = ID(J)
            ER(L) = 1.0D0
            IF (L.EQ.14.OR.L.EQ.29) ER(L) = DMAX1(1d-4,H(L,1))
            DO K = K1, K2
               ER(L) = DMAX1(ER(L), DABS(H(L,K)))
            END DO
         END DO
      ELSE
         ISTAR = NOC - 1
         IF (ISTAR.EQ.1) THEN ! Can be done cleaner...
            WRITE(32,*) "Star 1 nucleosynthesis"
         ELSEIF (ISTAR.EQ.2) THEN
            WRITE(32,*) "Star 2 nucleosynthesis"
         END IF
C Scale minor elements relative to most abundant
         DO J = 1,NV
            IF (J.LE.30) THEN
               L = ID(J)
            ELSE
               L = ID(30) + (J-30)
            END IF
            ER(L) = 1d-8
            ER(1) = 1d0
            ER(2) = 1d0
C            ER(28) = 1d0
            DO K = K1, K2
               ER(L) = DMAX1(ER(L), DABS(HNUC(L+50*(ISTAR-1),K)))
            END DO
         END DO
      END IF
      MM = 2 - MIN0(1, NE4)
C BEGIN ITERATIVE LOOP
 600  CONTINUE
      DO KTER = 1, ITER
         D = 0.0D0
         K = K1
         KH = IH
C EVALUATE FUNCTIONS AT SURFACE MESHPOINT
         CALL DIFRNS(K, NOC, ITER, N15, NE, 60+NEV)
         CALL DIVIDE(N15, N6, NE, N7, N8, K)
         K = K1 + 1
         IF ( K.GT.K2 ) GO TO 100
C DITTO NEXT-TO-SURFACE, ELIMINATING SOME UNKNOWNS
         CALL DIFRNS(K, NOC, ITER, 1, NE, 30)
         CALL ELIMN8(N2, NE, N3, N15, N4, N5, K-1)
         CALL DIVIDE(1, N4, NE, N7, N8, K)
 50      CONTINUE
         K = K + 1
         IF ( K.LE.K2 ) THEN
            IF ( K.EQ.K1+3 .OR. K.EQ.JH-1 .OR. K.EQ.JH+2 .OR. 
     :           K.EQ.K2-1 ) KH = IH - KH
C DITTO REMAINING MESHPOINTS
            CALL DIFRNS(K, NOC, ITER, 1, NE, 30)
            CALL ELIMN8(1, NE4, NB, N15, N1, NE, K-2)
            CALL ELIMN8(N1, NE4, NE, 1, N4, N5, K-1)
            CALL ELIMN8(N2, NE, N3, N15, N4, N5, K-1)
            CALL DIVIDE(1, N4, NE, N7, N8, K)
            GO TO 50
         END IF
         CALL DIFRNS(K, NOC, ITER, 1, N16, 60)
         CALL ELIMN8(N2, N16, N3, N15, N4, N5, K-2)
         CALL ELIMN8(N4, N16, N5, 1, N8, N9, K-1)
         CALL ELIMN8(N6, N16, N7, N15, N8, N9, K-1)
C SOLVE FOR CORRECTIONS AT CENTRAL MESHPOINT
         CALL DIVIDE(1, N8, N16, N11, N12, K)
C BY BACK-SUBSTITUTION FIND CORRECTIONS THROUGHOUT
         NN = 1
 60      CONTINUE
         K = K - 1
         IF ( K.NE.K1 .OR. NB.NE.0 ) THEN
            IF ( K.EQ.K1 ) NN = N15
            KK = K + 1
            DO J = 1, N16
               IF ( J.GE.N15 ) KK = K2 + 1
               DO I = NN, NE
                  C(K, I, N14) = C(K, I, N14) - C(K, I, J)
     :                              *C(KK, J, N14)
               END DO
            END DO
            IF ( K.GT.K1 ) GO TO 60
         END IF
         DO K = K1, K2
C Put corrections into first column
            IF ( NEV.NE.0 ) THEN
               DO I = N2, NV
                  C(K, I, 1) = C(K2+1, I-NB, N14)
               END DO
            END IF
            IF ( NB.NE.0 ) THEN
               DO I = 1, NB
                  C(K, I, 1) = C(K, I+N13, N14)
               END DO
            END IF
            DO I = 1, N13
               C(K, I+NB, 1) = C(K+1, I, N14)
            END DO
         END DO
 100     CONTINUE
         IF ( IH.NE.0 ) THEN
            DO K = K1, K2
               WRITE(24, 99001) (C(K,I,1), I=1, NV)
            END DO
99001       FORMAT (10E12.4)
         END IF
         ERR = 0.0D0
         ERMAX = 0.0D0
         DO J = 1, NV
C ESTIMATE ACCURACY OF ITERATION
            VX = 0.0D0
            DO K = K1, K2
               VZ = C(K, J, 1)
               IF ( DABS(VZ).GT.DABS(VX) ) VX = VZ
               IF ( VX.EQ.VZ ) MK(J) = K
               ERR = ERR + DABS(VZ)
            END DO
            ERT(J) = VX
            IF (ABS(ERT(J)).GT.ERMAX) JMAX = J
            IF (ABS(ERT(J)).GT.ERMAX) KMAX = MK(J)
            ERMAX = MAX(ERMAX, ABS(ERT(J)))
         END DO
         ERR = ERR/(NV*KMESH)
*
* Extra code to reduce FAC when in difficulty.
*
         IF (KTER .GT. 3 .AND. ERR .GE. 0.99D0*ERRPR) REDUCE = .TRUE.
         IF (REDUCE) THEN
            IF (ERRPR .GT. ERR) THEN
               FAC = 1.2D0*FAC
               IF (FAC .GE. 1D0) THEN
                  FAC = 1D0
                  REDUCE = .FALSE.
               END IF
            ELSE
               FAC = 0.8D0*FAC
            END IF
         ELSE
            FAC = DEL/DMAX1(ERR, DEL)
         END IF
         FRR = DLOG(ERR)*0.43429D0
         IF (NOC.GE.2.AND.FRR.GT.-5.85) FAC = 0.9d0
         IF (NV.LE.12) THEN
            IF (KTER.GT.NWRT5 .OR. ERR.GE.1D-2 .OR. ERMAX.GE.1D1)
     :           WRITE(32, 99002) KTER, FRR, FAC, (MK(J), ERT(J), J=1, NV)
         ELSE
            IF (KTER.GT.NWRT5 .OR. ERR.GE.1D-2 .OR. ERMAX.GE.1D1)
     :           WRITE(32, 99003) KTER, FRR, FAC, (MK(J), ERT(J), J=1, NV)
         END IF
99002    FORMAT (1X, I2, 2F7.2, 12(I4,F7.4))
99003    FORMAT (1X, I2, 2F7.2, 12(I4,F7.4),3(/,17X,12(I4,F7.4)))
C RJS 10/11/02 - Added lines to print out reason for failure
         IF ((KTER.LT.6).OR.ITER.GT.30) THEN !NOC.GE.2.AND.
            IF (ERR.GE.1d-1) ERR = 1d-2
            IF (ERMAX.GT.1d1) ERMAX = 1d0
         END IF
         IF (ERR.GE.1d-1) THEN 
            WRITE(32,*) "ERR tolerance exceeded:",ERR
         END IF
         IF (ERMAX.GE.1D1) THEN 
            WRITE(32,*) "ERMAX tolerance exceeded in variable ", JMAX
            WRITE(32,*) "at meshpoint ", KMAX, " ERMAX=", ERMAX 
         END IF
C NaN trap
         IF (ERR.NE.ERR) THEN
            WRITE(32,*) "ERR has become NaN (SNAFU), restarting"
            IF (DD.GT.0.1) THEN
               DD = DD / 2.0
            END IF

            ERR = 1d1
         END IF
         CALL FLUSH (32)
         IF ( IH.GT.0 ) IH = IH - 1
         IE(9) = IH
C APPLY CORRECTIONS, SCALED DOWN IF TOO LARGE
         IF (NOC.LE.1) THEN
            DO I = 1, NV
               J = ID(I)
               DO K = K1, K2
                  DH(J, K) = DH(J, K) - C(K, I, 1)*FAC*ER(J)
               END DO
            END DO
         ELSE
            DO I = 1, NV
               IF (I.LE.30) THEN
                  J = ID(I)
               ELSE
                  J = ID(30) + (I-30)
               END IF
               DO K = K1, K2
                  DHNUC(J+50*(ISTAR-1),K) = DHNUC(J+50*(ISTAR-1),K)
     :                 - C(K, I, 1)*FAC*ER(J)
               END DO
            END DO
         END IF
*
* Don't let the error get too large.
*
         JTER = KTER
         IF (ERR.LT.EPS.AND.FACSG.NE.1d0) THEN
            FACSG = 10.0*FACSG
            FACSG = DMIN1(1d0,FACSG)
            WRITE(32,*) "Boosting mixing", FACSG
            ERR = 1d-3
C Yes, I know I shouldn't do this, but haven't worked out another way yet
            GOTO 600
         END IF
C Sort out neutrons once convergence ok
         IF (ERR.LT.EPS.AND.NOC.GE.2) CALL NEUTRON(ISTAR)
         IF (ERR .LT. EPS .OR. ERR.GE.1d-1 .OR. ERMAX.GE.1D1) RETURN
 200     ERRPR = ERR
C CONTINUE ITERATING IF NOT YET ACCURATE ENOUGH
      END DO
      RETURN
      END