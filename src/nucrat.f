**==NUCRAT.FOR
      SUBROUTINE NUCRAT(TL)
* Compute rates of (at present) 20 nuclear reactions, and the corresponding
* energy and neutrino release
      IMPLICIT REAL*8(A-H, N-Z)
      COMMON /AUXIN / ICL, ION, JW, IOP, INUC, IBC, ICN, IML(2), ISGTH,
     :     IMO, IDIFF
      COMMON /STAT1 / CAT(119230), CRT(200, 92), JCSX
      COMMON /STAT2 / W1(4), RHO, W2(4), ZT, W3(8), 
     :                RRT(21), EX, ENX, EXH, EXHE, EXC, WR(16)
      COMMON /ABUND / XA(10), N(10),
     &                NE, NI, NZZ, AVM
      COMMON /CNSTS / CPI, PI4, CLN10, CDUM(9), CPL, CMEVMU, CDUM2(5)
      COMMON /NCDATA/ QRT(20), QNT(20), CZA(91), CZB(91), CZC(91),
     &                CZD(91), VZ(10)
      dimension czw(20)
      CBRT(VX) = DEXP(DLOG(VX)/3.0D0)
      DATA CSA, CSB, CSC, CSD, CXD /0.624, 0.316, 0.460, 0.38, 0.86/
      data czw /2,8,8,0,8,12,14,16,16,24,28,32,40,72,96,128,0,0,0,0/
* RHB is 'baryon density': 1 amu * number of baryons per cm3
      RHB = RHO/AVM
* Electron screening theory from Graboske, DeWitt, Grossman & Cooper (1973),
* for strong (ZA, ZB, ZC) are intermediate screening (ZD). The reaction
* dependent charge parameters are stored in CZA ... CZD.
      WC = 0.0
      DO I = 1, 10
         WC = WC + N(I)*VZ(I)
      END DO
C Ionization details for He3 not worked out - RJS
      N1 = N(1)
      N4 = N(2)
      N12 = N(3)
      N14 = N(4)
      N16 = N(5)
      N20 = N(6)
      N24 = N(7)
      N28 = N(8)
      N56 = N(9)
      N3 = N(10)
      WC = WC/NI
      WB = NE/NI
      WA = ZT*ZT/(WB*WB)
      XB = CBRT(ABS(WB))
      VL = CPL*SQRT(ABS(NI)*RHB*EXP(-3.0D0*TL))
      ZA = CSA*XB*EXP(LOG(VL)/1.5D0)
      ZB = CSB*XB*ZA
      ZC = CSC/(XB*XB)
      ZD = CSD*WC*WA*EXP(LOG(VL/(WA*ZT))*CXD)
* weak screening
      zw = 0.5*zt*vl
* Reaction rates interpolated in T, mostly from Caughlan & Fowler (1988)
      TF = TL/CLN10
      DO J = 1, 19
         RN = 0.0D0
         TT = 50.0*(TF - 6.0) + 1.0
         IF ( TT .GE. 1.0D0 ) THEN
            IT = MAX(1, MIN(199, INT(TT)))
            TT = TT - IT
            TU = 1.0D0 - TT
            RR = TU*CRT(IT, J) + TT*CRT(IT+1, J)
            IF ( RR .GE. -50.0D0 ) THEN
               SCRN = ZD*CZD(J)
               STRN = ZA*CZA(J) + ZB*CZB(J)
               DSTR = ZC*CZC(J)
               IF (DSTR .LT. 0.29*STRN) SCRN = MIN(SCRN, STRN - DSTR)
* weak screening
               if (inuc.ge.10) scrn = zw*czw(j)
               RN = EXP(CLN10*(RR + 20.0D0) + SCRN)*1.0D-20
            END IF
         END IF
         RRT(J+1) = RN
      END DO
C Sort out rates
      RPP = RRT(2)
      R33 = RRT(3)
      R34 = RRT(4)
      RBE = RRT(5)
      RBP = RRT(6)
      RPC = RRT(7)
      RPN = RRT(8)
      RPO = RRT(9)
      R3A = RRT(10)
      RAC = RRT(11)
      RAN = RRT(12)
      RAO = RRT(13)
      RANE = RRT(14)
      RCC = RRT(15)
      RCO = RRT(16)
      ROO = RRT(17)
      RGNE = RRT(18)
      RGMG = RRT(19)
      RCCG = RRT(20)
      RPNG = RRT(21)
      if(mod(inuc,10).eq.1)then
* correct rates to simulate Bahcall (1992) cross sections
         rpp = rpp*0.9828
         r33 = r33*0.971
         r34 = r34*0.987
         t6r = 10.0**(0.5*tf - 3.0)
         rbe = 5.54d-9/t6r*(0.936 + 0.004*t6r*t6r)
         rbp = rbp*0.933
      elseif(mod(inuc,10).eq.2)then
* idem, for Bahcall (1995) cross sections
         rpp = rpp*0.9557
         r33 = r33*0.969
         r34 = r34*0.970
         t6r = 10.0**(0.5*tf - 3.0)
         rbe = 5.54d-9/t6r*(0.936 + 0.004*t6r*t6r)
         rbp = rbp*0.933
         rpn = rpn*0.991
      end if
* Multiply with density and abundances to get rates per baryon per second,
* note that abundances of He3 and Be7 are not needed in equilibrium
      RPP = RHB*N1*N1*RPP/2.0
      R33 = RHB*N3*N3*R33/2.0 !RHB*R33/2.0
      R34 = RHB*N3*N4*R34 !RHB*N4*R34
      RBE = RHB*NE*RBE
      RBP = RHB*N1*RBP
      RPC = RHB*N1*N12*RPC
      RPN = RHB*N1*N14*RPN
      RPO = RHB*N1*N16*RPO
      R3A = RHB*RHB*N4*N4*N4*R3A/6.0
      RAC = RHB*N4*N12*RAC
      RAN = RHB*N4*N14*RAN
      RAO = RHB*N4*N16*RAO
      RANE = RHB*N4*N20*RANE
      RCC = RHB*N12*N12*RCC/2.0
      RCO = RHB*N12*N16*RCO
      ROO = RHB*N16*N16*ROO/2.0
      RGNE = N20*RGNE
      RGMG = N24*RGMG
* Branching of pN and CC reactions
      FPNG = 8.0D-4
      RPNA = (1.0 - FPNG)*RPN
      RPNG = FPNG*RPN
      RPN = RPNA
      FCCG = RCCG
      RCCA = (1.0 - FCCG)*RCC
      RCCG = FCCG*RCC
      RCC = RCCA
* PP chain in equilibrium, RPP becomes effective rate of 2 H1 -> 0.5 He4
C      F34 = 0.0
C      IF (R34 .GT. 1.0D-20)
C     :     F34 = 2.0/(1.0 + SQRT(1.0 + 8.0*RPP*R33/(R34*R34)))
C      RPP = RPP*(1.0 + F34)
C      PP2 = 1.0
C      IF (RBE+RBP .GT. 1.0D-20) PP2 = RBE/(RBE + RBP)
C      PP3 = 1.0 - PP2
C      QPP = QRT(1) + 0.5*QRT(2)
C      QNPP = (QNT(1) + F34*(QNT(4)*PP2 + QNT(5)*PP3))/(1.0 + F34)
C Put rates back to RRT
      RRT(2) = RPP
      RRT(3) = R33
      RRT(4) = R34
      RRT(5) = RBE
      RRT(6) = RBP
      RRT(7) = RPC
      RRT(8) = RPN
      RRT(9) = RPO
      RRT(10) = R3A
      RRT(11) = RAC
      RRT(12) = RAN
      RRT(13) = RAO
      RRT(14) = RANE
      RRT(15) = RCC
      RRT(16) = RCO
      RRT(17) = ROO
      RRT(18) = RGNE
      RRT(19) = RGMG
      RRT(20) = RCCG
      RRT(21) = RPNG
* calculate energy release and neutrino loss, in erg/gram/sec
      EX = 0d0
      ENX = 0d0
      DO J = 1, 3
         EX = EX + QRT(J)*RRT(J+1)
         ENX = ENX + QNT(J)*RRT(J+1)
      END DO
C      EX = QPP*RPP
c      EXH = CMEVMU*((QPP + 0.5*Q33)*RPP + QPC*RPC + QPO*RPO + QPNA*RPN 
c     :     + QPNG*RPNG)/AVM
c      EXHE = CMEVMU*(Q3A*R3A + QAC*RAC + QAN*RAN + QAO*RAO + QANE*RANE)
c     :     /AVM
c      EXC = CMEVMU*(QCCA*RCC + QCCG*RCCG + QCO*RCO + QOO*ROO + QGNE*RGNE
c     :     + QGMG*RGMG)/AVM
c      EX = EXH + EXHE + EXC
C      ENX = QNPP*RPP
      DO J = 6, 20
         EX = EX + QRT(J)*RRT(J+1)
         ENX = ENX + QNT(J)*RRT(J+1)
      END DO
      EX = CMEVMU*EX/AVM
      ENX = -CMEVMU*ENX/AVM
      RETURN
      END
