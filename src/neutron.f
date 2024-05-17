**==NEUTRON.F
C Sorts out where n's get used up - I hope!
C NEEDS TO BE OVERHAULED!
      SUBROUTINE NEUTRON(ISTAR)
      IMPLICIT NONE
      REAL*8 RANE22, RANE, RPC14, RANE22N, W6, RCO, ETH, ABUND
      REAL*8 WW2, RANA23NT, W5, RPBE, RPLI, RAN15G, TRB, RPNE21
      REAL*8 RPNA23A, XA2, R33, RPMG26, RPB11, W1, RANE21, RNNA22A
      REAL*8 RBE, HNUC, RPNE22, RN26TP, RAC13N, RN26GP, RPMG24, RP26G
      REAL*8 RPF19A, STAT2, RN26TA, PROTONS, RPNA23, RPMG26MN, RPNA23N, RPO18A
      REAL*8 RPNG, RAMG25P, RAO18N, RN26GA, RANA23NM, RPMG25T, RPNA22, RPSI30
      REAL*8 DT, GRAD, RANA23NG, EGR, RPMG25G, RPF19, RAC14, RAMG25
      REAL*8 RG26MP, DH, RAO, R34, RCC, RPNE20, RAMG26P, CSI
      REAL*8 AF, BARYN, AT, RPSI28, V, RP26T, VLEDD, RPAL27
      REAL*8 ROO, RPAL27A, W7, RPO, RPSI29, XA, EPS, VM
      REAL*8 ABUND2, RANE21N, RAF19, RPMG26GN, OP, RPP, RP26M, RAO18
      REAL*8 W8, RNNA22P, RGMG, R, RAN, RAMG25N, RPN15, RPO17A
      REAL*8 RCCG, RAMG26, R3A, WR, DHNUC, CFE, RLIA, REBE
      REAL*8 RN26MP, CMG, RPMG25M, RAAL27N, RAMG26N, RG26TP, RPN, RLF
      REAL*8 RPC, GR, RPD, H, RPO18, RAO17N, RATETOT, RPO17
      REAL*8 RG26GP, RAC, RGNE, RPC13, RATES, RAMG24, ZS, RN26MA
      REAL*8 FRAC, RBP, SODDS, RPN15A, RPMG26TN, QQ, WW
      INTEGER JIN, NE, NMESH, IE, IW, K, ID, NUCMAT
      INTEGER MAXMSH, I, NEUTRO, ISTAR, LT

      PARAMETER (MAXMSH = 2000)
      COMMON /NEUTRO/ FRAC(50,MAXMSH)
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,V(2),NMESH,JIN,ID(100),IE(100)
      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      COMMON /SODDS / W5(2), CMG, CSI, CFE, W8(11), DT, W6(23), IW(50), TRB
      COMMON /STAT2 / W1(19), RPP, R33, R34, RBE, RBP, 
     :                RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO,
     :                ROO, RGNE, RGMG, RCCG, RPNG, WR(21)
      COMMON /RATES / W7(1),RPD,RPBE,RPLI,REBE,RPC13,RPN15A,RPO17a,
     :     RAC13n, RAN15g, RAO17n, RLIA, RPB11, RPC14, RAC14,
     :     RPO18,RPO18a,RAO18,RAO18n,RPF19,RPF19A,RAF19,RPNE21,
     :     RANE21,RANE21n,RPNE22,RANE22,RANE22n,RnNA22p,RnNA22a,RPNA22,
     :     RPNA23,RPNA23n,RPNA23a,RPMG24,RAMG24,RPMG25T,RPMG25M,RPMG25G,
     :     RAMG25,RAMG25n,RAMG25p,RPMG26,RPMG26Tn,RPMG26Mn,RPMG26Gn,
     :     RAMG26,RAMG26n,RAMG26p,Rg26Tp,Rn26Tp,Rn26Ta,Rp26T,Rg26Mp,
     :     Rn26Mp,Rn26Ma,Rp26M,Rg26Gp,Rn26Gp,Rn26Ga,Rp26G,RPAL27,
     :     RPAL27a,RAAL27n,RANA23nT,RANA23nM,RANA23nG,RPSI28,RPSI29,
     :     RPSI30,RPN15,RPO17, RPNE20
      COMMON /ABUND / XA(10), WW(14)
      COMMON /ABUND2/ XA2(50), WW2(50)
      COMMON /OP    / ZS, VLEDD, VM, GR, GRAD, ETH, RLF, EGR, R, QQ
      DIMENSION BARYN(50)
      DIMENSION PROTONS(MAXMSH)
      DATA BARYN/62.0, 1.0, 2.0, 3.0, 7.0, 7.0, 11.0, 13.0, 14.0,  15.0,
     :     17.0, 18.0, 19.0, 21.0, 22.0, 22.0, 23.0, 24.0, 25.0, 26.0,
     :     26.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,56.0,57.0,
     :     58.0,59.0,60.0,59.0,58.0,59.0,60.0,61.0,1.0,4.0,12.0,14.0,16.0,
     :     20.0,1.0,1.0,1.0,1.0/ 
      DO K=1,NMESH
         DO I=1,50
            IF (I.NE.2.AND.I.NE.41) THEN
C Don't change proton and neutron fractions - they are dealt with below
               DHNUC(I+50*(ISTAR-1),K) = DHNUC(I+50*(ISTAR-1),K) - 
     :              FRAC(I,K)*BARYN(I)*DHNUC(2+50*(ISTAR-1),K)
            END IF
         END DO
         PROTONS(K) = - FRAC(41,K)*BARYN(41)*DHNUC(2+50*(ISTAR-1),K)
C         DH(53,K) = DH(53,K) + FRAC(41,K)*BARYN(41)*DH(14,K)
         DHNUC(2+50*(ISTAR-1),K) = 0d0
      END DO
C Sort out remaining protons - RJS 27/8/04
      DO K=1,NMESH
C Call to statel to get electron density & no. abundances.
         AF = H(1+15*(ISTAR-1),K) + DH(1+15*(ISTAR-1),K)
         AT = H(2+15*(ISTAR-1),K) + DH(2+15*(ISTAR-1),K)
         CALL STATEL (0, AF, AT, ISTAR)
C Blank any -ve abundances. There shouldn't be any anyway...
C RJS 19/1/04
         XA(1) = H(5+15*(ISTAR-1),K) + DH(5+15*(ISTAR-1),K)
         XA(2) = H(9+15*(ISTAR-1),K) + DH(9+15*(ISTAR-1),K)
         XA(3) = H(10+15*(ISTAR-1),K) + DH(10+15*(ISTAR-1),K)
         XA(4) = H(12+15*(ISTAR-1),K) + DH(12+15*(ISTAR-1),K)
         XA(5) = H(11+15*(ISTAR-1),K) + DH(11+15*(ISTAR-1),K)
         XA(6) = H(3+15*(ISTAR-1),K) + DH(3+15*(ISTAR-1),K)
         XA(7) = CMG*ZS
         XA(8) = CSI*ZS
         XA(9) = CFE*ZS
         DO I = 1, 10
            IF (XA(I).LT.0d0) THEN
               XA(I) = 0d0
               WW(I) = 0d0
            END IF
         END DO
         DO I = 1, 50
            XA2(I) = HNUC(I+50*(ISTAR-1),K) + DHNUC(I+50*(ISTAR-1),K)
            IF (XA2(I).LT.0d0) XA2(I) = 0d0
            XA2(2) = 0d0
         END DO
C Pipe excess protons to nucrat as only source of them
         XA2(41) = PROTONS(K)
C Work out abundances for nucrat2 - RJS
         DO I = 1, 50
            WW2(I) = XA2(I)/BARYN(I)
         END DO
         CALL NUCRAT2(AT)
         RATETOT = 2d0*RPP + RPC + RPN + RPO + RPNG + RPD + RPBE + RPLI
     :        + RPC13 + RPN15A + RPO17a + RPB11 + RPC14 + RPO18 + RPO18a 
     :        + RPF19 + RPF19A + RPNE21 + RPNE22 + RPNA22 + RPNA23
     :        + RPNA23n + RPNA23a + RPMG24 + RPMG25M + RPMG25G + RPMG26
     :        + RPMG26Mn + RPMG26Gn + Rp26M + Rp26G + RPAL27 + RPAL27a
     :        + RPSI28 + RPSI29 + RPSI30 + RPN15 + RPO17 + RPNE20
         IF (RATETOT.NE.0d0) THEN
            DO I = 1,50
               FRAC(I,K) = 0d0
            END DO
            FRAC(3 ,K) = - RPP + RPD
            FRAC(4 ,K) = - RPD
            FRAC(5 ,K) = RPLI
            FRAC(6 ,K) = RPBE
            FRAC(7 ,K) = RPB11
            FRAC(8 ,K) = RPC13 - RPC
            FRAC(9 ,K) = RPC14
            FRAC(10,K) = RPN15A - RPN - RPO18a + RPN15
            FRAC(11,K) = - RPO + RPO17a + RPO17
            FRAC(12,K) = RPO18 + RPO18a - RPO17
            FRAC(13,K) = RPF19 + RPF19A - RPO18
            FRAC(14,K) = RPNE21 - RPNE20
            FRAC(15,K) = RPNE22
            FRAC(16,K) = RPNA22 - RPNE21
            FRAC(17,K) = RPNA23 + RPNA23n + RPNA23a - RPNE22
            FRAC(18,K) = RPMG24 - RPNA23 - RPAL27a
            FRAC(19,K) = RPMG25M + RPMG25G
            FRAC(20,K) = RPMG26 + RPMG26Mn + RPMG26Gn
            FRAC(21,K) = Rp26M - RPMG25M - RPMG26Mn
            FRAC(22,K) = Rp26G - RPMG25G - RPMG26Gn
            FRAC(23,K) = RPAL27 + RPAL27a - RPMG26
            FRAC(24,K) = RPSI28 - RPAL27
            FRAC(25,K) = RPSI29 - RPSI28
            FRAC(26,K) = RPSI30 - RPSI29
            FRAC(27,K) = - RPSI30
C Why didn't I do the major isotopes as well?
            DO I=1,50
               FRAC(I,K) = FRAC(I,K)/RATETOT
               DHNUC(I+50*(ISTAR-1),K) = DHNUC(I+50*(ISTAR-1),K)
     :              - FRAC(I,K)*BARYN(I)*PROTONS(K)
            END DO
C            DHNUC(41+50*(ISTAR-1),K) = 0d0
         END IF
      END DO
      RETURN   
      END
