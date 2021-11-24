**==PRINTB.FOR
      SUBROUTINE PRINTB(DTY,PER,NWRT1,NWRT2,NWRT3,NWRT4,NMONT,NMOD,IEND)
      IMPLICIT REAL*8(A-H, L, M, O-Z)
      SAVE
      INTEGER MAXMSH
      PARAMETER (MAXMSH = 2000)
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,DEL,DH0,NMESH,JIN,IW(200)
      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      COMMON /TRANS / HT(26,MAXMSH,2)
      COMMON /SODDS / ALPHA, RML, CMG, CSI, CFE, CT(10), AGE, DT, M1,
     :  EC, BM, ANG, CM, MTA, MTB, TM(2), T0, M0, TC(2), OS, AC, RCD,
     :  RMG, RHL, XF, DR, AK1 ,RMT, AK2, ITH, IX, IY, IZ, IB, ISX(45),
     :  TRB
* extra common for mesh-spacing
      COMMON /PMESH / PMH(2), PME(2), IAGB
      COMMON /INF   / Q(60)
      COMMON /DINF  / QD(60)
      COMMON /OP    / ZS, LEDD, VM, GR, GRAD, ETH, RLF, EGR, R, QQ
      COMMON /ABUND / WA1(10), YA(10), WA(3), AVM
c GR is ov value, XX is Schwarzschild value, of gradr-grada
      COMMON /STAT2 / PL, RL, U, P, RHO, FK, T, SF, ST, ZT, GRADA, CP,
c     :                CH, S, PR, PG, PF, PT, EN, RR(20), EX, WX, EXX(3),
c     :                RW(16)
     :                CH, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP,
     :                RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO,
     :                ROO, RGNE, RGMG, RCCG, RPNG, EX, ENX, EXX(3), AAWT,
     :                WT, RW(14)
      COMMON /YUK1  / PX(34), WMH, WMHE, VMH, VME, VMC, VMG, BE, VLH,
     :                VLE, VLC, VLN, VLT, MCB(12),WWW(100)
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CL, CD, CG, CR, CC, CEVB,
     &                CEN, CPL, CMEVMU, CSECYR, LSUN, MSUN, RSUN, TSUNYR
      COMMON /ATDATA/ DH2(4), CHI(26,9), OMG(27), AM(10), BN(10), JZ(10)
      COMMON /NCDATA/ QPP, Q33, Q34, QBE, QBP, QPC, QPNA, QPO, Q3A, QAC,
     :                QAN, QAO, QANE, QCCA, QCO, QOO, QGNE, QGMG, QCCG,
     :                QPNG, CNUC(390)
      COMMON /WT    / PWT,SG,MK,MT2,LQP,LKP
*
* Extra COMMON for main-sequence evolution.
*
      COMMON /ZAMS  / TKH(2)
      COMMON /MESH  / TRC1,TRC2,DD,DT1,DT2,MWT,MWTS,IVMC,IVMS
      COMMON /MIX   / KBICZ,KICZ
      COMMON /DHBLOC/ IDREDGE
      COMMON /WIND  / WIND,SWIND,PERIOD
      COMMON /EVMODE/ IMODE
      COMMON /INERTI/ VI(2)
      COMMON /MONTAG/ WCV, WL
C Extra common for extra timestep control stuff - P for previous, C for current
      COMMON /DTCONT/ VLHP(2), VLEP(2), VLCP(2), RLFP(2), TOTMP(2), VLHC(2),
     :     VLEC(2), VLCC(2), RLFC(2), TOTMC(2)
      COMMON /WINDS / WINDML(2), FAKEWIND(2), BE2
      COMMON /TIDES / MENV(2), RENV(2)
      PS(VX) = 0.5D0*(VX+DABS(VX))
      DIMENSION XA(10), TCB(12), RCB(12)
      DIMENSION SX(45,MAXMSH+1)
      DIMENSION XINIT(44),YIELD(44)
      CHARACTER*5 AX(34), BX(16), CX(16), DX(17)
      INTEGER KCB(12), KMX(3), KEX(12), KENV(2)
      REAL*8 EMAX(3), VMX(3), MEX(12), THB(12) ! EXX(3),
      INTEGER IPX(25)
      DATA IPX/9, 17, 2, 3, 4, 5, 6, 8, 10, 11, 12, 13, 14, 15, 16,
     &     18, 19, 20, 21, 28, 27, 7, 24, 25, 26/
C Data required to compute yields - Anders & Grevasse '89 abunds
C Used A&G for H,He,C,N,O,Ne rather than input for convenience
      DATA XINIT/4.801d-5,2.929d-5,9.353d-9,0d0,4.725d-9,3.650d-5,0d0,
     :     4.363d-6,3.887d-6,2.167d-5,4.051d-7,4.127d-6,1.302d-4,0d0,
     :     3.339d-5,5.148d-4,6.766d-5,7.760d-5,0d0,0d0,5.798d-5,6.530d-4,
     :     3.426d-5,2.352d-5,8.155d-5,3.958d-4,3.222d-6,1.866d-5,1.169d-3,
     :     2.855d-5,3.697d-6,0d0,0d0,3.358d-6,4.944d-5,0d0,1.958d-5,8.594d-7,
     :     7.057d-1,2.752d-1,3.032d-3,1.105d-3,9.592d-3,1.619d-3/
      DATA YIELD/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,
     :     0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,
     :     0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/
      DATA IDREDGE/0/
      DATA AX/' psi ','  P  ',' rho ','  T  ','kappa','grada',' grad',
     :'gr-ga','  m  ','  H1 ',' He4 ',' C12 ',' N14 ',' O16 ',' Ne20',
     :' He3 ','  r  ','  L  ',' Eth ',' Enuc',' Eneu','  dm ',' k**2',
     :'n/n+1','lr/lp','lm/lp',' U   ',' S   ','L/Edd',' mu  ','  mu ',
     :' LK ','Ebind','beta '/
      DATA BX/'  g  ','  n  ','  D  ',' He3 ',' Li7 ',' Be7 ',' B11 ',
     :     ' C13 ',' C14 ',' N15 ',' O17 ',' O18 ',' F19 ',' Ne21',
     :     ' Ne22','  m  '/
      DATA CX/' Na22',' Na23',' Mg24',' Mg25',' Mg26','Al26M',
     :     'Al26G',' Al27',' Si28',' Si29',' Si30',' P31 ',' S32 ',
     :     ' S33 ',' S34 ','  m  '/
      DATA DX/' Fe56',' Fe57',' Fe58',' Fe59',' Fe60',' Co59',
     :     ' Ni58',' Ni59',' Ni60',' Ni61', '  p  ',' He4 ',' C12 ',
     :     ' N14 ',' O16 ',' Ne20','  m  '/
99001 FORMAT (/, '  K', 17(4X, A5, 2X),/)
99002 FORMAT (I4, 1P, 15(1X,E10.3))
99003 FORMAT (/, 1X, I7, '  dt/age/MH/MHe  tn/tKH/Mb  P/rlf/dM LH/LHe/',
     :'LC Lth/Lnu/m  H1      He4     C12     N14     O16     Ne20    ',
     :'He3       psi    dens    temp', 2(/, F8.4, 1P, E16.9, 4E10.3,
     :0P, 7F8.5, F8.3, 2F8.4, 2X, A4), /, 3F8.4, F10.4, 1P, 2E10.3, 0P,
     :F10.4, 7F8.5, F8.3, 2F8.4, 2X, A4, /, F8.4, 12F8.3, 5F8.4, I6,
     :/,2I4,12I8,/,15F8.4,1P,3E10.3,0P,/,15I8)
99004 FORMAT (I13, 1P, 9E13.5, 3(/, 10E13.5))
99005 FORMAT (/, '  K', 15(6X, A5, 3X),/)
99006 FORMAT (I4, 1P, 11E14.7)
C Store previous values of extra timestep control variables
      VMHP = VMH
      DO I=1,2
         VLHP(I) = VLHC(I)
         VLEP(I) = VLEC(I)
         VLCP(I) = VLCC(I)
         RLFP(I) = RLFC(I)
         TOTMP(I) = TOTMC(I)
      END DO
            DO ISTAR = 1,IMODE
               JC = 1
               JTHB = 1
               JE = 1
               KTM = 0
               KCE = 1
               KMH = 0
               KME = 0
               TMAX = 0.0D0
               VI(ISTAR) = 0d0
               DO K=1,NMESH - 1
                  VI(ISTAR) = VI(ISTAR) + 2.0/3.0*(DEXP(H(4+15*(ISTAR-1),K)+DH(4+15*(ISTAR-1),K))
     : - DEXP(H(4+15*(ISTAR-1),K+1)+DH(4+15*(ISTAR-1),K)))*DEXP(H(7+15*(ISTAR-1),K)
     :                 +DH(7+15*(ISTAR-1),K))**2.0
               END DO
               VI(ISTAR) = VI(ISTAR) + 2.0/5.0*DEXP(H(4+15*(ISTAR-1),NMESH)+DH(4+15*(ISTAR-1),NMESH))*
     :              DEXP(H(7+15*(ISTAR-1),NMESH)+DH(7+15*(ISTAR-1),NMESH))**2.0
               DO J = 1, 3
                  KMX(J) = 1
                  VMX(J) = 0.0D0
C     EXX(J) = 0.0D0
                  EMAX(J) = 0.0D0
               END DO
               EXLIM = 10.0*H(8,1)/EXP(H(4,1))
               DO J = 1, 34
                  PX(J) = 0.0D0
                  SX(J,1) = 1.0D-8
                  IF (J .LE. 12) KCB(J) = 0
                  IF (J .LE. 12) MCB(J) = 0d0
                  IF (J .LE. 12) KEX(J) = 0
                  IF (J .LE. 12) MEX(J) = 0.0D0
                  IF (J .LE. 12) THB(J) = 0.0D0
               END DO
               WMH = 0d0
               WMHE = 0d0
               VMH = 0d0
               VME = 0d0
               VMC = 0d0
               VMG = 0d0
               BE = 0d0
               VLH = 0d0
               VLE = 0d0
               VLC = 0d0
               VLN = 0d0
               VLT = 0d0
               NPRINT = 0
C     Decide whether to print the interior or not
               IF ( MOD(NMOD,NWRT1).EQ.0 ) NPRINT = 1
               IF ( NPRINT.EQ.1 ) WRITE(32+20*(ISTAR-1), 99001) (AX(ISX(J)), J=1,15)
               IF (NPRINT.EQ.1.AND.IW(102).NE.0) WRITE (36+20*(ISTAR-1),*) NMOD
               IF (NPRINT.EQ.1.AND.IW(102).NE.0) WRITE (37+20*(ISTAR-1),*) NMOD
               IF (NPRINT.EQ.1.AND.IW(102).NE.0) WRITE (38+20*(ISTAR-1),*) NMOD
      IF ( NPRINT.EQ.1.AND.IW(102).NE.0 ) WRITE(36+20*(ISTAR-1), 99001) (BX(J), J=1,16)
      IF ( NPRINT.EQ.1.AND.IW(102).NE.0 ) WRITE(37+20*(ISTAR-1), 99001) (CX(J), J=1,16)
      IF ( NPRINT.EQ.1.AND.IW(102).NE.0 ) WRITE(38+20*(ISTAR-1), 99001) (DX(J), J=1,17)
         DO KK = 1, NMESH
            K = NMESH + 1 - KK
            DO J = 1, 20
               IF ( IEND.LT.0 ) QD(J) = 0.0D0
               IF ( IEND.LT.0 ) Q(J) = H(J+15*(ISTAR-1), K)
               IF ( IEND.GE.0 ) QD(J) = DH(J+15*(ISTAR-1), K)
               IF ( IEND.GE.0 ) Q(J) = H(J+15*(ISTAR-1), K) + QD(J)
            END DO
            IF (ISTAR.EQ.2) Q(28) = H(13, K)
            CALL FUNCS1 (-1, K, 1, NMESH, ISTAR)
C call mass-loss to get surface conditions
            IF (ISTAR.EQ.1.AND.IMODE.EQ.2.AND.K.EQ.1) THEN
               DO J=1,30
                  IF ( IEND.LT.0 ) QD(J) = 0.0D0
                  IF ( IEND.LT.0 ) Q(J) = H(J, K)
                  IF ( IEND.GE.0 ) QD(J) = DH(J, K)
                  IF ( IEND.GE.0 ) Q(J) = H(J, K) + QD(J)
               END DO
               CALL MASSLOSS(-1,A1,A2,A3,A4,A5)
            END IF
*
* Correct abundances for actual atomic masses
*
            DO I = 1, 10
               XA(I) = AM(I)*YA(I)/AVM
            END DO
            XH = XA(1)
            XHE = XA(2)
            XC = XA(3)
            XN = XA(4)
            XO = XA(5)
            XNE = XA(6)
            XMG = XA(7)
            XSI = XA(8)
            XFE = XA(9)
            XHE3 = XA(10)
C Evaluate the functions to be printed
            WF=DSQRT(1.0D0+DEXP(Q(1)))
            PX(1) = 2.0D0*(WF-DLOG(WF+1.0D0))+Q(1)
            PX(2) = P
            PX(3) = RHO
            PX(4) = T
            TMAX = MAX(TMAX, SX(4,KK))
            IF ( PX(4).GE.TMAX ) KTM = KK + 1
            PX(5) = FK
            PX(6) = GRADA
            PX(7) = GRAD
            PX(8) = DMIN1(9.99D0, GR)
            PX(9) = VM/MSUN
C locate convective/radiative boundaries (GR=0)
            IF (KK.LT.2.OR.JC.GT.12.OR.SX(8,KK)*PX(8).GT.0.0D0
     &           .OR. AAWT.LT.1.0D-1) GO TO 1
C            IF (K.GT.KBICZ+1) GOTO 1
               MCB(JC) = (SX(9,KK)*PX(8)-PX(9)*SX(8,KK))/(PX(8)-SX(8,KK))
               KCB(JC) = K
               TCB(JC) = (SX(4,KK)*PX(8)-PX(4)*SX(8,KK))/(PX(8)-SX(8,KK)) !T
               RCB(JC) = (SX(3,KK)*PX(8)-PX(3)*SX(8,KK))/(PX(8)-SX(8,KK)) !RHO
C     Find T at base of conv. envelope
            IF (K.GT.50.AND.K.LT.300) THEN
               TBCE = (SX(4,KK)*PX(8)-PX(4)*SX(8,KK))/(PX(8)-SX(8,KK))
            END IF
C Find T and rho of convective pulse extrema
C         IF (K.LT.750) THEN
C            TCB(JC) = T
C            RCB(JC) = RHO
C         END IF
            IF ( PX(4).GT.1.0D5 ) KCE = KK
            JC = JC + 1
C locate convective/semiconvective boundaries (loosely, GR=0.01)
 1          IF ((PX(8)-DR)*(SX(8,KK)-DR).GT.0.0D0.OR.KK.LT.2.OR.JC.GT.12
     &           .OR. AAWT.LT.1.0D-1) GO TO 2
C            IF (K.GT.KBICZ+1) GOTO 2
               MCB(JC)=-(SX(9,KK)*(PX(8)-DR)-PX(9)*(SX(8,KK)-DR))
     :              /(PX(8)-SX(8,KK))
               KCB(JC) = -K
               JC=JC + 1
C locate burning shell boundaries
 2          IF (XH.GT.XF.AND.SX(10,KK).LT.XF) THEN
               VMH = (PX(9)*(XF-SX(10,KK))+SX(9,KK)*(XH-XF))/(XH-SX(10,KK))
               KMH = K
               PMH(ISTAR) = (PX(2)*(XF-SX(10,KK))+SX(2,KK)*(XH-XF))/(XH-SX(10,KK))
            END IF
            IF (XHE.GT.XF.AND.SX(11,KK).LT.XF.AND.K.GT.NMESH/2) THEN
               VME = (PX(9)*(XF-SX(11,KK))+SX(9,KK)*(XHE-XF))/
     &              (XHE-SX(11,KK))
               KME = K
               PME(ISTAR) = (PX(2)*(XF-SX(11,KK))+SX(2,KK)*(XHE-XF))/
     &              (XHE-SX(11,KK))
            END IF
            IF (XC.GT.XF.AND.SX(12,KK).LT.XF) VMC = (PX(9)*
     :           (XF-SX(12,KK))+SX(9,KK)*(XC-XF))/(XC-SX(12,KK))
C N abundance in H-burning ash
            IF (XH.GT.1d-10.AND.SX(10,KK).LT.1d-10) XASH = XN
C     Find T_hmid (XH = 0.345)
            IF (XH.GT.0.345.AND.SX(10,KK).LT.0.345) THEN
               THMID = (PX(4)*(0.345-SX(10,KK))+SX(4,KK)*(XH-0.345))/
     :              (XH-SX(10,KK))
            END IF
C Find T_hbase (XH = 0.05)
            IF (XH.GT.0.05.AND.SX(10,KK).LT.0.05) THEN
               THBASE = (PX(4)*(0.05-SX(10,KK))+SX(4,KK)*(XH-0.05))/
     :              (XH-SX(10,KK))
            END IF
C Find T_intmid (Xhe = 0.45)
            IF (XHE.GT.0.45.AND.SX(11,KK).LT.0.45) THEN
               TINTMID = (PX(4)*(0.45-SX(11,KK))+SX(4,KK)*(XHE-0.45))/
     :              (XHE-SX(11,KK))
            END IF
C Find T_intbase (Xhe = 0.05)
            IF (XHE.GT.0.05.AND.SX(11,KK).LT.0.05) THEN
               TINTBASE = (PX(4)*(0.05-SX(11,KK))+SX(4,KK)*(XHE-0.05))/
     :              (XHE-SX(11,KK))
            END IF
C find points of max energy generation
            EXX(1) = CMEVMU*((QPP + 0.5*Q33)*RPP + QPC*RPC + QPO*RPO
     :           + QPNA*RPN + QPNG*RPNG)
            EXX(2) = CMEVMU*(Q3A*R3A + QAC*RAC + QAN*RAN + QAO*RAO
     :           + QANE*RANE)
            EXX(3) = CMEVMU*(QCCA*RCC + QCCG*RCCG + QCO*RCO +QOO*ROO
     :           + QGNE*RGNE + QGMG*RGMG)
            DO J = 1, 3
               IF (EXX(J).GE.EMAX(J)) THEN
                  IF (EXX(J).GE.EXLIM) KMX(J) = K
                  IF (EXX(J).GE.EXLIM) VMX(J) = PX(9)
                  EMAX(J) = EXX(J)
               END IF
            END DO
            PX(10) = XH
            PX(11) = XHE
            PX(12) = XC
            PX(13) = XN
            PX(14) = XO
            PX(15) = XNE
            PX(16) = XHE3 !XMG
            PX(17) = R/RSUN
            PX(18) = Q(8)/LSUN
            PX(19) = ETH
            PX(20) = EX
C locate boundaries of burning regions (EX > EXLIM)
            IF((PX(20)-EXLIM)*(SX(20,KK)-EXLIM).LE.0.0D0.AND.JE.LT.12.AND.
     &           KK.GE.2) THEN
               MEX(JE)= (SX(9,KK)*(PX(20)-EXLIM)-PX(9)*(SX(20,KK)-EXLIM))
     :              /(PX(20)-SX(20,KK))
               KEX(JE) = K
               JE = JE  + 1
            END IF
            PX(21) = -EN
            PX(22) = PX(9)-SX(9,KK)
            WF = SX(17,KK)/PX(17)
            PX(23) = (SX(9,KK)*SX(23,KK)*WF**2+PX(22)*0.4D0*(1.0D0-WF**5)
     :           /(1.0D0-WF**3))/PX(9)
            WF = 1.0D0/DLOG10(PX(2)/SX(2,KK))
            IF ( KK.EQ.1 ) WF = 0.0D0
            PX(24) = WF*DLOG10(PX(3)/SX(3,KK))
            PX(25) = -WF*DLOG10(PX(17)/SX(17,KK))
            PX(26) = -WF*DLOG10(PX(9)/SX(9,KK))
            PX(27) = U
            PX(28) = S
            PX(29) = Q(8)/LEDD
            TADJ = FK*CP*T*(1.0D11*PX(22)/(4*PI4*R*R))**2/(3*CL*PR)/CSECYR
            PX(30) = LKP !TADJ
            PX(31) = CR*RHO*T/PG
C         PX(32) = HT(2,K)*HT(4,K)*1.0D66
            PX(32) = MT2*PS(PX(30) - SX(30,KK))! Thermohaline mixing coeff ! LKP
C RJS 4/4/07 Locate TH boundaries
            IF (KK.GT.2.AND.JTHB.LT.12.AND.((PX(32).LE.-1d-20.AND
     :   .SX(32,KK).GT.-1d-20).OR.(PX(32).GT.-1d-20.AND.SX(32,KK).LE.-1d-20))) THEN
               THB(JTHB) = (SX(9,KK)*PX(32)-PX(9)*SX(32,KK))/(PX(32)-SX(32,KK))
               JTHB = JTHB + 1
            END IF
C Print the interior details on first `page', if required
            IF ( NPRINT.EQ.1.AND.MOD(K-1,NWRT2).EQ.0 ) WRITE(32+20*(ISTAR-1),99002) K,
     :           (PX(ISX(J)), J=1,15)
99007       FORMAT (I4, 1P, 17(1X,E10.3))
            IF ( NPRINT.EQ.1.AND.MOD(K,10*NWRT2).EQ.0 ) WRITE(32+20*(ISTAR-1),99002)
            IF ( NPRINT.EQ.1.AND.MOD(K-1,NWRT2).EQ.0.AND.IW(102).NE.0)
     :           WRITE(36+20*(ISTAR-1),99007) K,
     :           HNUC(1+50*(ISTAR-1),K), HNUC(2+50*(ISTAR-1),K), HNUC(3+50*(ISTAR-1),K),
     :           HNUC(4+50*(ISTAR-1),K), HNUC(5+50*(ISTAR-1),K), HNUC(6+50*(ISTAR-1),K),
     :           HNUC(7+50*(ISTAR-1),K), HNUC(8+50*(ISTAR-1),K), HNUC(9+50*(ISTAR-1),K),
     :           HNUC(10+50*(ISTAR-1),K),HNUC(11+50*(ISTAR-1),K),HNUC(12+50*(ISTAR-1),K),
     :           HNUC(13+50*(ISTAR-1),K),HNUC(14+50*(ISTAR-1),K),HNUC(15+50*(ISTAR-1),K),
     :           PX(9)
            IF ( NPRINT.EQ.1.AND.MOD(K,10*NWRT2).EQ.0.AND.IW(102).NE.0)
     :           WRITE(36+20*(ISTAR-1),99002)
            IF ( NPRINT.EQ.1.AND.MOD(K-1,NWRT2).EQ.0.AND.IW(102).NE.0)
     :           WRITE(37+20*(ISTAR-1),99007) K,
     :           HNUC(16+50*(ISTAR-1),K),HNUC(17+50*(ISTAR-1),K),HNUC(18+50*(ISTAR-1),K),
     :           HNUC(19+50*(ISTAR-1),K),HNUC(20+50*(ISTAR-1),K),HNUC(21+50*(ISTAR-1),K),
     :           HNUC(22+50*(ISTAR-1),K),HNUC(23+50*(ISTAR-1),K),HNUC(24+50*(ISTAR-1),K),
     :           HNUC(25+50*(ISTAR-1),K),HNUC(26+50*(ISTAR-1),K),HNUC(27+50*(ISTAR-1),K),
     :           HNUC(28+50*(ISTAR-1),K),HNUC(29+50*(ISTAR-1),K),HNUC(30+50*(ISTAR-1),K),
     :           PX(9)
            IF ( NPRINT.EQ.1.AND.MOD(K,10*NWRT2).EQ.0.AND.IW(102).NE.0)
     :           WRITE(37+20*(ISTAR-1),99002)
            IF ( NPRINT.EQ.1.AND.MOD(K-1,NWRT2).EQ.0.AND.IW(102).NE.0 )
     :           WRITE(38+20*(ISTAR-1),99007) K,
     :           HNUC(31+50*(ISTAR-1),K),HNUC(32+50*(ISTAR-1),K),HNUC(33+50*(ISTAR-1),K),
     :           HNUC(34+50*(ISTAR-1),K),HNUC(35+50*(ISTAR-1),K),HNUC(36+50*(ISTAR-1),K),
     :           HNUC(37+50*(ISTAR-1),K),HNUC(38+50*(ISTAR-1),K),HNUC(39+50*(ISTAR-1),K),
     :           HNUC(40+50*(ISTAR-1),K),HNUC(41+50*(ISTAR-1),K),HNUC(42+50*(ISTAR-1),K),
     :           HNUC(43+50*(ISTAR-1),K),HNUC(44+50*(ISTAR-1),K),HNUC(45+50*(ISTAR-1),K),
     :           HNUC(46+50*(ISTAR-1),K), PX(9)
            IF ( NPRINT.EQ.1.AND.MOD(K,10*NWRT2).EQ.0.AND.IW(102).NE.0)
     :           WRITE(38+20*(ISTAR-1),99002)
            DO J=1,34
               SX(J,KK+1)=PX(J)
            END DO
C Some integrated quantities, to be printed in the short summary
            WMH = WMH + PX(22)*XH
            WMHE = WMHE + PX(22)*XHE
            IF (PX(9).GT.VMHP) THEN
               BE = BE - PX(22)*EGR
            ELSE
               BE = 0d0
            END IF
            BE2 = BE+0.01
         PX(33) = BE
         PX(34) = PR/(PG+PR)
            WF = PX(22)*MSUN/LSUN
            VLH = VLH + WF*EXX(1)
            VLE = VLE + WF*EXX(2)
            VLC = VLC + WF*EXX(3)
            VLN = VLN - WF*EN
            VLT = VLT + WF*ETH
            IF (NMONT.NE.0) THEN
               IF (MOD(NMOD,NMONT).EQ.0) THEN
C Output for use in MONTAGE -- may need smart output control RJS 29/5/08
99008             FORMAT (1P,12(1X,E14.7))
99018             FORMAT (I4, 1X, I6, 1X, 1P, E20.13, 2(1X,E12.5))
                  IF (K.EQ.NMESH) THEN
                     WRITE (42+20*(ISTAR-1),99018) NMESH, NMOD, AGE*CSECYR, VLH, VLE
                  END IF
                  WRITE (42+20*(ISTAR-1),99008) DLOG10(PX(4)), DLOG10(PX(3)), DLOG10(WL),
     :                   DLOG10(DMAX1(WCV,1d-30)), DLOG10(PX(17)*RSUN*1d11), PX(9),
     :                   PX(10), PX(16), PX(11), PX(12), PX(13), PX(14)
               END IF
            END IF
         END DO
C sneplot file here?????? - SMR + JJE
99019    FORMAT (I6,E1.9,0P)
         LOGG = DLOG10(6.67E-11*PX(9)/PX(17)**2)
         WRITE(49+(ISTAR-1), 99019) NMOD, LOGG
C End SNEPLOT
         IF (NWRT3.EQ.1) GO TO 3
C Write further `pages' for each detailed model, if required
         DO I=2, NWRT3
            IJ = 15*I-15
            IF ( NPRINT.EQ.1 ) WRITE(32+20*(ISTAR-1), 99001) (AX(ISX(J+IJ)), J=1,15)
            DO KK = 1, NMESH
               K = NMESH + 1 - KK
C Have to choose the SX's that are wanted
               IF ( NPRINT.EQ.1.AND.MOD(K-1,NWRT2).EQ.0 ) WRITE(32+20*(ISTAR-1),99002)
     :              K, (SX(ISX(J+IJ),KK+1), J=1,15)
               IF ( NPRINT.EQ.1.AND.MOD(K,10*NWRT2).EQ.0 ) WRITE(32+20*(ISTAR-1),99002)
            END DO
         END DO
 3       TM(ISTAR) = VM
C         PER = BM/3.55223D0*(ANG/(TM(1)*(BM-TM(1))))**3
         IF (ISTAR.EQ.1) THEN
            PER = BM/3.55223D0*(Q(13)/(TM(1)*(BM-TM(1))))**3
         ELSE
C            PER = BM/3.55223D0*(H(13,1)/(TM(1)*(BM-TM(1))))**3
            PER = (TM(1)+TM(2))/3.55223D0*(H(13,1)/(TM(1)*TM(2)))**3
         END IF
         DMT = PX(9)*DH(4+15*(ISTAR-1), 1)/DTY
         TKH(ISTAR) = DABS(BE)/(Q(8)*CSECYR)
         TN = 4.0D10*PX(9)/Q(8)
C convective mixing timescale needed in FUNCS1
         TC(ISTAR) = RCD/(CSECYR*TN*DSQRT(VM))
*
* Write to the numerical data storage unit for plotting purposes.
*
C Separation
         IF (IMODE.EQ.2) THEN
            SEP = (TM(1)+TM(2))*(H(13,1)/(TM(1)*TM(2)))**2.0
         ELSE
            SEP = BM*(H(13,1)/(TM(1)*(BM-TM(1))))**2.0
         END IF
C Total angular momentum
         HTOT = H(14,1) + H(29,1) + H(13,1)
C Total spin angular momentum
         HSPINTOT = H(14,1) + H(29,1)
C Orbital angular velocity
         OORB = H(13,1)*(TM(1)+TM(2))/(TM(1)*TM(2)*SEP**2.0)
         VIORB = (TM(1)*TM(2)/(TM(1)+TM(2)))*SEP**2.0
         OSPIN = Q(14)/VI(ISTAR)
         OCRIT = DSQRT(6.67d-11*PX(9)*2d30/(6.96d8*PX(17))**3.0) !/DSQRT(CG)
         SEP = SEP/RSUN
         WRITE (33+20*(ISTAR-1),115) NMOD,AGE,LOG10(PX(17)),LOG10(PX(4)),
     &        LOG10(PX(18)),PX(9),VMH,VME,LOG10(MAX(VLH,1.01D-10)),
     &        LOG10(MAX(VLE,1.01D-10)), LOG10(MAX(VLC,1.01D-10)),
     &        MCB,VMX(1),VMX(2),LOG10(FK),DTY,(PX(JJ),JJ=10,14), PX(16),
     &        RLF,Q(14),PER,SEP,BM/MSUN,H(13,1),HSPINTOT,HTOT,
     &        OORB*DSQRT(CG), OSPIN*DSQRT(CG), VI(ISTAR),OCRIT, DMT,MEX,THB,MENV(ISTAR)/MSUN,
     &        RENV(ISTAR)/RSUN, DLOG10(SX(3,2)), DLOG10(SX(4,2))
C There are 99 things in the format statement and 74 have been used.
 115     FORMAT (I6,1P,E16.9,0P,24F10.5,1P,3E13.6,18(1X,E12.5),0P,52F9.5)
         CALL FLUSH(33+20*(ISTAR-1))
C Write output for nucleosynthesis stuff
         write (41+20*(ISTAR - 1),116) NMOD,AGE,VMH,VME,VMH - VME,PX(9),XASH,(TCB(I)/1d8,
     :        I=1,12),RCB,
     :        TBCE,THMID,THBASE,TINTMID,TINTBASE,PERIOD,DABS(DMT)
 116     format (I6,1P,E16.9,0P,17(1X,E13.6),1P,19(1X,E13.6),15I5)
         CALL FLUSH(41+20*(ISTAR - 1))
         IF (IW(102).NE.0) THEN
C Write surface/centre abundances, yields -- but only if using NS code!
            write (39+20*(ISTAR - 1),117) NMOD,AGE,(HNUC(I+50*(ISTAR-1),1),I=1,50),
     :           WINDML(ISTAR)/MSUN*CSECYR, FAKEWIND(ISTAR)/MSUN*CSECYR, DTY
 117        format (I6,1P,E16.9,54(1X,E13.6))
            CALL FLUSH(39+20*(ISTAR - 1))
            write (40+20*(ISTAR - 1),117) NMOD,AGE,(HNUC(I+50*(ISTAR-1),NMESH),I=1,50)
     :           , H(5+15*(ISTAR - 1),NMESH)
            call flush(40+20*(ISTAR - 1))
         END IF
*
* Write detailed model data for plotting purposes. -- Note new unit number needed! RJS
*
C      NWRT5 = 21
C      IF (IEND.LT.0) WRITE (11,'(2I6)') NMESH, NWRT5
C      IF (NPRINT.EQ.1) THEN
C         WRITE (11,10001) NMOD, AGE
C         DO KK = 1, NMESH
C            WRITE (11,10002) (SX(IPX(J),KK+1), J=1, NWRT5)
C         END DO
C      END IF
10001    FORMAT (I6,1P,E17.9,0P)
10002    FORMAT (1P,E13.6,4E11.4,25E11.3,0P)
*
C Print the short summary, for every NWRT4'th model
         IF (MOD(NMOD,NWRT4).NE.0) RETURN
         BMS = BM/MSUN
         IF (IMODE.EQ.1) THEN
            BM = BM - (WINDML(1)+FAKEWIND(1))*DT
            BMS = BM/MSUN
         END IF
         FR = DLOG10(PX(17))
         FL = DLOG10(PX(18))
* print at max He energy generation, rather than Tmax
* Restore to Tmax
         KHE = NMESH + 2 - KMX(2)
         SDC = DLOG10(SX(3,2))
         SDM = DLOG10(SX(3,KTM))
         SDS = DLOG10(PX(3))
         STC = DLOG10(SX(4,2))
         STM = DLOG10(SX(4,KTM))
         STS = DLOG10(PX(4))
         PER0 =  1.5D0*FR - 0.5D0*DLOG10(8.157D0*PX(9))
         WMI = DLOG10(SX(9,KCE)*SX(17,KCE)**2*SX(23,KCE))
C RJS - Viscous mesh
         IF (IVMC.EQ.1) THEN
C            MWT = DMAX1(0d0, 0.75*(DTY/1d-4 - 1))
            MWT = DMAX1(0d0, 0.25*(DTY/1d-4 - 1))
            MWT = DEXP(-DSQRT(MWT))
            MWT = DMAX1(0.1d0,DMIN1(1d0,MWT))
         ELSE
            MWT = 0d0
         END IF
C Timestep control for TP-AGB
C Switch on mechanism for dealing with dredge up when necessary
         VLEC(ISTAR) = VLE
         LHEF = (VLEC(ISTAR) - VLEP(ISTAR))/VLEP(ISTAR)
         IF (IAGB.EQ.1) THEN
            IF (IDREDGE.EQ.0.AND.VLEC(ISTAR).GT.5d5.AND.LHEF.LT.0d0) IDREDGE = 1
C Find out where next to outer conv boundary is
            IF (IDREDGE.EQ.1) THEN
               MAXB = 0d0
               DO II = 1,12
                  IF (MCB(II).GT.MAXB.AND.MCB(II).LT.1d0) THEN
                     MAXB = MCB(II)
                  END IF
               END DO
               IDREDGE = 2
            END IF
C If boundary moves in, reduce DD
            IF (IDREDGE.EQ.2.AND.VLEC(ISTAR).LT.1d5) THEN
               MAXB2 = 0d0
               DO II = 1,12
                  IF (MCB(II).GT.MAXB2.AND.MCB(II).LT.1d0) MAXB2 = MCB(II)
               END DO
               IF (MAXB2.LT.MAXB) IDREDGE = 3
            END IF
            IF (IDREDGE.EQ.3) THEN
               IF (VLEC(ISTAR).GT.2d4) THEN
                  DD = 5.0
               ELSE
                  DD = 1.0
               END IF
            END IF
            IF (IDREDGE.EQ.3.AND.VLEC(ISTAR).LT.3d3) THEN
               DD = 5.0
               IDREDGE = 0
            END IF
            KICZ = KCB(3)
            KBICZ = KMX(2)
         END IF
         WRITE(32+20*(ISTAR-1),99003) NMOD,PX(9),DTY,TN,PER,VLH,VLT,SX(10,2),SX(11,2),
     :        SX(12,2),SX(13,2),SX(14,2),SX(15,2),SX(16,2),SX(1,2),SDC,STC,
     :        'cntr',VMH,AGE,TKH(ISTAR),RLF,VLE,VLN,XH,XHE,XC,XN,XO,XNE,XHE3,PX(1),
     :        SDS,STS,'srfc',VME,WMH,WMHE,BMS,DMT,VLC,SX(9,KTM),SX(10,KTM),
     :        SX(11,KTM),SX(12,KTM),SX(13,KTM),SX(14,KTM),SX(15,KTM),SX(16,KTM),
     :        SX(1,KTM),SDM,STM,'Tmax',VMC,MCB,WMI,PX(23),PER0,FR,FL,NMOD,
     :        KMH,KME,KCB,VMX,MEX,EMAX,KMX,KEX
         CALL FLUSH (32+20*(ISTAR-1))
         CALL FLUSH (36+20*(ISTAR-1))
         CALL FLUSH (37+20*(ISTAR-1))
         CALL FLUSH (38+20*(ISTAR-1))
C Find convective envelope base
         MENV(ISTAR) = 0d0
         RENV(ISTAR) = 0d0
        DO K=NMESH-1,1,-1 !surf to centre
            IF (SX(8,K+1).GT.0d0.AND.SX(8,K).LT.0d0.AND.SX(9,K).GT.MENV(ISTAR)/MSUN
     :          .AND.K.LT.NMESH-20) THEN
               MENV(ISTAR) = MSUN*(SX(9,K)*SX(8,K+1)-SX(9,K+1)*SX(8,K))/(SX(8,K+1)-SX(8,K))
               RENV(ISTAR) = RSUN*(SX(17,K)*SX(8,K+1)-SX(17,K+1)*SX(8,K))/(SX(8,K+1)-SX(8,K))
               KENV(ISTAR) = NMESH - K
            END IF
         END DO
C Update current values for extra timestep control variables
         VLHC(ISTAR) = VLH
         VLEC(ISTAR) = VLE
         VLCC(ISTAR) = VLC
         RLFC(ISTAR) = RLF
         TOTMC(ISTAR) = TM(ISTAR)
      END DO
      RETURN
      END
