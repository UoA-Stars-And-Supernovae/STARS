      SUBROUTINE DIFRNS(K, NOC, ITER, NX, NY, NZ)

      IMPLICIT NONE

      REAL*8 ER, DINE, TEMPSTORE, OUTE, FNSTORE, VX, DABS, DOUTE
      REAL*8 DFN1, HNUC, EQU, DVAR, DOUTF, VARN, BC4, VAR
      REAL*8 DH, C, GT, EVMODE, YUK, OUTF, FN1, EPS
      REAL*8 DINF, DX, BC5, S, EQ, BC3, DH0, DVARN
      REAL*8 BC2, BC1, DEL, DHNUC, DINFN, D, SOLV, DEQU
      REAL*8 FN2, DVX, GE, H, XD, DFN2
      INTEGER NZ, L4, JIN, II, N12, J2, NB, NE2
      INTEGER IH, N14, K1, NY, IMODE, N5, N10, JI
      INTEGER N6, M, INFN, NE3, I4, IA, K2, MM
      INTEGER NOC, NEV, NE, J1, NV, NMESH, IW, K
      INTEGER KH, INE, ID, LE, NW, NUCMAT, MAXMSH, J
      INTEGER ITER, NF, I, L, NE1, L3, N2, ISTAR
      INTEGER NX, JJ, INF, JH

      PARAMETER (MAXMSH = 2000)
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,DEL,DH0,NMESH,JIN,IW(200)
      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      COMMON /SOLV  / C(MAXMSH+1,51,51),S(50,151),ER(60),D,KH,NE,N2,N5, 
     :                N6,N10, N12, N14, MM, IA(1), NE1, NE2, NE3, NB,  
     :                NEV, NF, J1, J2, IH, JH, ID(90), K1, K2, NV, NW
      COMMON /INF   / VAR(60)
      COMMON /DINF  / DVAR(60)
      COMMON /INFN  / VARN(50)
      COMMON /DINFN / DVARN(50)
      COMMON /OUTF  / FN1(154)
      COMMON /DOUTF / DFN1(50, 154)
      COMMON /INE   / FN2(3, 154)
      COMMON /DINE  / DFN2(3, 50, 154)
      COMMON /OUTE  / EQU(50)
      COMMON /DOUTE / DEQU(50, 3, 50)
      COMMON /YUK   / XD(MAXMSH)
      COMMON /EVMODE/ IMODE
      DIMENSION TEMPSTORE(154), FNSTORE(3,154)

      IF ( K.LE.K2 .OR. NOC.NE.1 ) THEN
            IF ( K.GT.K2 .AND. NOC.GE.2 ) THEN
                  GO TO 100
            END IF
C REDEFINE PREVIOUS CURRENT MESH POINT AS CURRENT PREVIOUS MESH POINT
            DO J = MM, 2
                  JJ = J + 1
                  DO I = 1, NF
                        FN2(J, I) = FN2(JJ, I)
                  END DO
                  DO L = 1, NV
                        DO I = 1, NF
                              DFN2(J, L, I) = DFN2(JJ, L, I)
                        END DO
                  END DO
            END DO
            DO L = 1, 30
                  DVAR(L) = DH(L, K)
                  VAR(L) = H(L, K) + DVAR(L)
            END DO
            IF ( NOC.GE.2 ) THEN
C Call funcs1 to get new values of all physical variables
                  IF (IMODE.EQ.2) THEN
C place star 2 data into INF
                        DO L = 1,15
                              DVAR(L) = DH(L+15, K)
                              VAR(L) = H(L+15, K) + DVAR(L+15)
                              DVAR(L+15) = DH(L,K)
                              VAR(L+15) = H(L,K) + DH(L,K)
                        END DO
                        CALL FUNCS1(0, K, K1, K2, 2)
C Restore VAR stuff 
                        DO L = 1, 30
                              DVAR(L) = DH(L, K)
                              VAR(L) = H(L, K) + DVAR(L)
                        END DO
                  END IF
C call funcs1 for star 1
                  CALL FUNCS1(0, K, K1, K2, 1)
C Call massloss routine to get BC right
                  CALL MASSLOSS(0,BC1,BC2,BC3,BC4,BC5)
C evaluate nucleosynthesis for this star
                  ISTAR = NOC - 1

                  DO L = 1, NV
                        DVARN(L) = DHNUC(L+50*(ISTAR-1), K)
                        VARN(L) = HNUC(L+50*(ISTAR-1), K) + DVARN(L)
                  END DO

                  CALL FUNCS2(K, K1, K2, ITER, ISTAR)

                  DO J = 1, NV
                        DO I = 1, NF
                              DFN2(3, J, I) = DFN1(J, I)
                        END DO
                  END DO

                  DO I = 1, NF
                        FN2(3, I) = FN1(I)
                  END DO

                  GO TO 100
            ELSE
C EVALUATE ARGUMENT LIST OF VARIABLES AND THEIR INCREMENTS AT CURRENT ME
                  DO L = 1, NW
                        XD(L) = XD(L+NV)
                  END DO
C EVALUATE AND STORE THE REQUIRED FUNCTIONS OF THE INDEPENDENT VARIABLES
                  IF (IMODE.EQ.2) THEN
C place star 2 data into INF
                        DO L = 1,15
                              DVAR(L) = DH(L+15, K)
                              VAR(L) = H(L+15, K) + DVAR(L+15)
                              DVAR(L+15) = DH(L,K)
                              VAR(L+15) = H(L,K) + DH(L,K)
                        END DO

                        CALL FUNCS1(0, K, K1, K2, 2)
C Put output into the latter half of INE
                        DO I = 1,51
                              FN1(I+51) = FN1(I)
                        END DO
C Restore VAR stuff 
                        DO L = 1, 60
                              DVAR(L) = DH(L, K)
                              VAR(L) = H(L, K) + DVAR(L)
                        END DO
                  END IF
C call funcs1 for star 1
                  CALL FUNCS1(0, K, K1, K2, 1)
C call boundary condition if this is surface MP
C will need to ammend if new things are passed to equns
                  IF (K.EQ.K1) THEN
                        CALL MASSLOSS(0,BC1,BC2,BC3,BC4,BC5)
                        FN1(1) = BC1
                        FN1(1+51) = BC2
                        FN1(4) = BC3
                        FN1(5) = BC4
                        FN1(5+51) = BC5
                  END IF

                  IF ( KH.GE.1 ) THEN
                        WRITE(24, 99001) VAR, FN1
                  END IF

99001             FORMAT (/, 20(/,1P8E16.8))

                  DO I = 1, NF
                        FN2(3, I) = FN1(I)
                  END DO
C BY VARYING EACH INDEPENDENT VARIABLE IN TURN, EVALUATE NUMERIC DERIVAT
                  DO L = 1, NV
                        M = ID(L)
                        VX = VAR(M)
                        DVX = DVAR(M)
                        DX = DH0*(DABS(VX)+ER(M))
                        XD(L+NW) = ER(M)/DX
                        VAR(M) = VX + DX
                        DVAR(M) = DVX + DX

                        IF (IMODE.EQ.2) THEN
C Call for star 2
                              DO I = 1,30
C Store variations to be restored later
                                    TEMPSTORE(I) = VAR(I)
                                    TEMPSTORE(I+30) = DVAR(I)
                              END DO

                              DO I = 1,15
                                    VAR(I) = VAR(I+15)
                                    DVAR(I) = DVAR(I+15)
                                    VAR(I+15) = TEMPSTORE(I)
                                    DVAR(I+15) = TEMPSTORE(I+30)
                              END DO
C               IF (K.EQ.199) write (*,*) L, M, ID(L)
                              CALL FUNCS1(M, K, K1, K2, 2)
C               IF (K.EQ.199) write (*,*) L, M, ID(L)
                              IF ( KH.GE.1 ) THEN
                                    WRITE(24, 99001) VAR, FN1
                              END IF

                              DO I = 1, 51  !NF
                                    DFN2(3, L, I+51) = FN1(I)
                              END DO
C Restore order of variations
                              DO I = 1,30
                                    VAR(I) = TEMPSTORE(I)
                                    DVAR(I) = TEMPSTORE(I+30)
                              END DO
                        END IF
C Call for star 1
                        CALL FUNCS1(M, K, K1, K2, 1)
C call boundary condition if this is surface MP
C will need to ammend if new things are passed to equns
                        IF (K.EQ.K1) THEN
                              CALL MASSLOSS(M,BC1,BC2,BC3,BC4,BC5)

                              FN1(1) = BC1
                              DFN2(3,L,1+51) = BC2
                              FN1(4) = BC3
                              FN1(5) = BC4
                              DFN2(3, L, 5+51) = BC5
                        END IF

                        IF ( KH.GE.1 ) THEN
                              WRITE(24, 99001) VAR, FN1
                        END IF

                        DO I = 1, 51     !NF
                              DFN2(3, L, I) = FN1(I)
                        END DO

                        DVAR(M) = DVX
                        VAR(M) = VX
                  END DO
            END IF
      END IF
C EVALUATE AND STORE THE DIFFERENCE RELATIONS WHICH ARE TO BE SATISFIED
      IF (IMODE.EQ.2) THEN
C Store FN2 prior to shunting numbers around
            DO J = 1,3
                  DO I = 1,102
                        FNSTORE(J,I) = FN2(J,I)
                  END DO
            END DO
C Move star 2 stuff to front of FN2
            DO J = 1,3
                  DO I = 1, 51
                        FN2(J,I) = FN2(J,I+51)
                  END DO
            END DO
C Call equns1 for star 2
            CALL EQUNS1(K, K1, K2, 2, 0)
C Store output in latter half of OUTE
            DO I = 1,15
                  EQU(I+15) = EQU(I)
            END DO
C Restore original order of FN2
            DO J = 1,3
                  DO I = 1,102
                        FN2(J,I) = FNSTORE(J,I)
                  END DO
            END DO
      END IF
C Call equns1 for star 1
      CALL EQUNS1(K, K1, K2, 1, 0)

      DO J = NX, NY
            M = ID(J+NZ)
            S(J, N12) = EQU(M)
      END DO
C BY VARYING EACH INDEPENDENT VARIABLE IN TURN, EVALUATE NUMERIC DERIVAT
      L3 = MM*NV - NV

      DO II = MM, 3
            DO JI = 1, NV
                  L3 = L3 + 1

                  IF ( JI.LE.NE ) THEN
                        L4 = NE*(II-1) + JI
                  ELSE
                        L4 = N5 + JI
                  END IF

                  DO I = 1, NF
                        FN1(I) = FN2(II, I)
                        FN2(II, I) = DFN2(II, JI, I)
                  END DO

                  IF (IMODE.EQ.2) THEN
C Store FN2 prior to shunting numbers around
                        DO J = 1,3
                              DO I = 1,102
                                    FNSTORE(J,I) = FN2(J,I)
                              END DO
                        END DO
C Move star 2 stuff to front of FN2
                        DO J = 1,3
                              DO I = 1, 51
                                    FN2(J,I) = FN2(J,I+51)
                              END DO
                        END DO
C Stuff to get composition accretion right
                        L = ID(JI)
C Call equns1 for star 2
                        CALL EQUNS1(K, K1, K2, 2, L)
C Store output in latter half of OUTE
                        DO I = 1,15
                              EQU(I+15) = EQU(I)
                        END DO
C Restore original order of FN2
                        DO J = 1,3
                              DO I = 1,102
                                    FN2(J,I) = FNSTORE(J,I)
                              END DO
                        END DO
                  END IF
C Call equns1 for star 1
                  CALL EQUNS1(K, K1, K2, 1, L)

                  DO I = 1, NF
                        FN2(II, I) = FN1(I)
                  END DO

                  DO J = NX, NY
                        M = ID(J+NZ)
                        IF ( JI.LE.NE .OR. II.EQ.MM ) THEN
                              S(J, L4) = (EQU(M)-S(J,N12))*XD(L3)
                        END IF

                        IF ( JI.GT.NE .AND. II.GT.MM ) THEN
                              S(J, L4) = (EQU(M)-S(J,N12))*XD(L3) + S(J, L4)
                        END IF
                  END DO
            END DO
      END DO

      RETURN

 100  CONTINUE

 501  FORMAT (50e15.6)
 502  FORMAT (I4,50e15.6)

      CALL EQUNS2(K, K1, K2, NE)

      DO J = NX, NY
C 29/9/03 RJS - no longer using ID
            IF (J.LE.30) THEN
                  M = ID(J+NZ)        ! finds equation no.
            ELSE
                  M = ID(30+NZ) + (J - 30)
            END IF

            DO I = 1,NE
                  IF (I.GT.30) THEN
                        L = ID(30)+(I - 30)
                  ELSE
                        L = ID(I)        ! finds location in H(J,K)
                  END IF

                  S(J, I) = DEQU(I, 1, M)*ER(L)
                  S(J, I+NE) = DEQU(I, 2, M)*ER(L)
                  S(J, I+N5) = DEQU(I, 3, M)*ER(L)
            END DO
            IF ( NEV.NE.0 ) THEN
                  DO I = N2, NV
                        S(J, I+N5) = (DEQU(I,1,M)+DEQU(I,2,M)+DEQU(I,3,M))
     :                             * ER(I+12)
                  END DO
            END IF

            S(J, N12) = EQU(M)
      END DO

      RETURN
      END