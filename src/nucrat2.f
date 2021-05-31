**==NUCRAT2.FOR
      SUBROUTINE NUCRAT2(TL)
C RJS 14/8/03 - funcs2 compatible version of nucrat.f
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
      COMMON /NCDATA/ QRT(20), QNT(20), CZA(92), CZB(92), CZC(92),
     &                CZD(92), VZ(10)
      COMMON /ABUND2/ XA2(50),NNG, Nn, NN2,NN3,NNL7,NNB7,NN11,NN13,NN14,
     : NN15,NN17,NN18,NN19,NNE21,NNE22,NNA22,NNA23,NMG24,NMG25,NMG26,
     : N26M,N26G,NAL27,NSI28,NSI29,NSI30,NP31,NS32,NS33,NS34,
     : NFE56,NFE57,NFE58,NFE59,NFE60,NCO59,NNI58,NNI59,NNI60,NNI61,NN1,
     : NN4,NN12,NNN14,NN16,NN20,W(4)
      COMMON /RATES / RRT2(73)
      COMMON /NDATA / NRT(200,45)
      COMMON /NRATES/ NRATE(45)
      COMMON /DECAY / RDAL26G,RDNA22,RD26M,RDFe59,RDFe60,RDNi59,RDn,RDC14
      CBRT(VX) = DEXP(DLOG(VX)/3.0D0)
      dimension czw(92)
      DATA CSA, CSB, CSC, CSD, CXD /0.624, 0.316, 0.460, 0.38, 0.86/
      data czw /2,8,8,0,8,12,14,16,16,24,28,32,40,72,96,128,0,0,0,0,
     :     2,6,16,0,12,14,16,24,28,32,12,10,12,48,16,16,64,64,18,18,
     :     72,20,80,80,20,80,80,0,0,22,22,22,22,24,96,24,24,24,96,96,
     :     96,24,24,24,24,96,96,96,0,0,0,26,0,0,0,26,0,0,0,26,26,26,104,
     :     88,88,88,28,28,28,14,16,20/
* RHB is 'baryon density': 1 amu * number of baryons per cm3
      RHB = RHO/AVM
* Electron screening theory from Graboske, DeWitt, Grossman & Cooper (1973),
* for strong (ZA, ZB, ZC) are intermediate screening (ZD). The reaction
* dependent charge parameters are stored in CZA ... CZD.
      WC = 0.0
      DO I = 1, 10
         WC = WC + N(I)*VZ(I)
      END DO
      N1 = N(1)
      N4 = N(2)
      N12 = N(3)
      N14 = N(4)
      N16 = N(5)
      N20 = N(6)
      N24 = N(7)
      N28 = N(8)
      N56 = N(9)
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
      DO J = 1, 92
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
         IF (J+1.LE.20) THEN
            RRT(J+1) = RN
         ELSE
            RRT2(J-19) = RN
         END IF
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
      RPP = RHB*NN1*NN1*RPP/2.0
C      RPP = RHB*RPP/2.0 ! Why divide by 2?
      R33 = RHB*NN3*NN3*R33/2.0
      R34 = RHB*NN3*NN4*R34
      RBE = RHB*NE*RBE
      RBP = RHB*NN1*RBP
      RPC = RHB*NN1*NN12*RPC
      RPN = RHB*NN1*NNN14*RPN
      RPO = RHB*NN1*NN16*RPO
      R3A = RHB*RHB*NN4*NN4*NN4*R3A/6.0
      RAC = RHB*NN4*NN12*RAC
      RAN = RHB*NN4*NNN14*RAN
      RAO = RHB*NN4*NN16*RAO
      RANE = RHB*NN4*NN20*RANE
      RCC = RHB*NN12*NN12*RCC/2.0
      RCO = RHB*NN12*NN16*RCO
      ROO = RHB*NN16*NN16*ROO/2.0
      RGNE = NN20*RGNE
C      RGMG = N24*RGMG
* Branching of pN and CC reactions
      FPNG = 8.0D-4
      RPNA = (1.0 - FPNG)*RPN
      RPNG = FPNG*RPN
      RPN = RPNA
      FCCG = RCCG
      RCCA = (1.0 - FCCG)*RCC
      RCCG = FCCG*RCC
      RCC = RCCA
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
C Minor variable reaction rates - 3/9/03 RJS
      RRT2(1) = 0d0
      RRT2(2) = RHB*NN1*NN2*RRT2(2)
      RRT2(3) = RHB*NN1*NNB7*RRT2(3)
      RRT2(4) = RHB*NN1*NNL7*RRT2(4)
      RRT2(5) = RHB*NE*NNB7*RRT2(5)
      RRT2(6) = RHB*NN1*NN13*RRT2(6)
      RRT2(7) = RHB*NN1*NN15*RRT2(7)
      RRT2(8) = RHB*NN1*NN17*RRT2(8) 
      RRT2(9) = RHB*NN4*NN13*RRT2(9)
      RRT2(10) = RHB*NN4*NN15*RRT2(10)
      RRT2(11) = RHB*NN4*NN17*RRT2(11)
      RRT2(12) = RHB*NN4*NNL7*RRT2(12)
      RRT2(13) = RHB*NN1*NN11*RRT2(13)
C C14 reactions
      RRT2(14) = RHB*NN1*NN14*RRT2(14)
      RRT2(15) = RHB*NN4*NN14*RRT2(15)
C O18 reactions
      RRT2(16) = RHB*NN1*NN18*RRT2(16)
      RRT2(17) = RHB*NN1*NN18*RRT2(17)
      RRT2(18) = RHB*NN4*NN18*RRT2(18)
      RRT2(19) = RHB*NN4*NN18*RRT2(19)
C F19 reactions
      RRT2(20) = RHB*NN1*NN19*RRT2(20)
      RRT2(21) = RHB*NN1*NN19*RRT2(21)
      RRT2(22) = RHB*NN4*NN19*RRT2(22)
C Ne21
      RRT2(23) = RHB*NN1*NNE21*RRT2(23)
      RRT2(24) = RHB*NN4*NNE21*RRT2(24)
      RRT2(25) = RHB*NN4*NNE21*RRT2(25)
C Ne22
      RRT2(26) = RHB*NN1*NNE22*RRT2(26)
      RRT2(27) = RHB*NN4*NNE22*RRT2(27)
      RRT2(28) = RHB*NN4*NNE22*RRT2(28)
C Na22
      RRT2(29) = RHB*Nn*NNA22*RRT2(29)
      RRT2(30) = RHB*Nn*NNA22*RRT2(30)
      RRT2(31) = RHB*NN1*NNA22*RRT2(31)
C Na23
      RRT2(32) = RHB*NN1*NNA23*RRT2(32)
      RRT2(33) = RHB*NN1*NNA23*RRT2(33)
      RRT2(34) = RHB*NN1*NNA23*RRT2(34)
C Mg24
      RRT2(35) = RHB*NN1*NMG24*RRT2(35)
      RRT2(36) = RHB*NN4*NMG24*RRT2(36)
C Mg25
      RRT2(37) = RHB*NN1*NMG25*RRT2(37)
      RRT2(38) = RHB*NN1*NMG25*RRT2(38)
      RRT2(39) = RHB*NN1*NMG25*RRT2(39)
      RRT2(40) = RHB*NN4*NMG25*RRT2(40)
      RRT2(41) = RHB*NN4*NMG25*RRT2(41)
      RRT2(42) = RHB*NN4*NMG25*RRT2(42)
C Mg26
      RRT2(43) = RHB*NN1*NMG26*RRT2(43)
      RRT2(44) = RHB*NN1*NMG26*RRT2(44)
      RRT2(45) = RHB*NN1*NMG26*RRT2(45)
      RRT2(46) = RHB*NN1*NMG26*RRT2(46)
      RRT2(47) = RHB*NN4*NMG26*RRT2(47)
      RRT2(48) = RHB*NN4*NMG26*RRT2(48)
      RRT2(49) = RHB*NN4*NMG26*RRT2(49)
C Al26T
C      RRT2(50) = RHB*N26T*RRT2(50)
C      RRT2(51) = RHB*Nn*N26T*RRT2(51)
C      RRT2(52) = RHB*Nn*N26T*RRT2(52)
C      RRT2(53) = RHB*NN1*N26T*RRT2(53)
      RRT2(50) = 0d0
      RRT2(51) = 0d0
      RRT2(52) = 0d0
      RRT2(53) = 0d0
C Al26M
      RRT2(54) = RHB*N26M*RRT2(54)
      RRT2(55) = RHB*Nn*N26M*RRT2(55)
      RRT2(56) = RHB*Nn*N26M*RRT2(56)
      RRT2(57) = RHB*NN1*N26M*RRT2(57)
C Al26G
      RRT2(58) = RHB*N26G*RRT2(58)
      RRT2(59) = RHB*Nn*N26G*RRT2(59)
      RRT2(60) = RHB*Nn*N26G*RRT2(60)
      RRT2(61) = RHB*NN1*N26G*RRT2(61)
C Al27
      RRT2(62) = RHB*NN1*NAL27*RRT2(62)
      RRT2(63) = RHB*NN1*NAL27*RRT2(63)
      RRT2(64) = RHB*NN4*NAL27*RRT2(64)
C Na23(a,n)Al26TGM
      RRT2(65) = RHB*NN4*NNA23*RRT2(65)
      RRT2(66) = RHB*NN4*NNA23*RRT2(66)
      RRT2(67) = RHB*NN4*NNA23*RRT2(67)
C Si reactions
      RRT2(68) = RHB*NN1*NSI28*RRT2(68)
      RRT2(69) = RHB*NN1*NSI29*RRT2(69)
      RRT2(70) = RHB*NN1*NSI30*RRT2(70)
C N15(p,gamma)
      RRT2(71) = RHB*NN1*NN15*RRT2(71)
C O17(p,gamma)F18(beta+)O18
      RRT2(72) = RHB*NN1*NN17*RRT2(72)
C Ne20(p,gamma)Na21(beta+)Ne21
      RRT2(73) = RHB*NN1*NN20*RRT2(73)
C Unstable particle decays
C Al26G t_0.5 = 0.72 Myr
      CLN2 = 0.69314718
      RDAL26G = (RHB/26.0)*N26G*CLN2/2.27d13
C C14 t_0.5 = 5730 yr
      RDC14 = (RHB/14.0)*NN14*CLN2/1.81d11
C Na22 t_0.5 = 2.6 yr
      RDNA22 = (RHB/22.0)*NNA22*CLN2/8.199d7
C Al26M t_0.5 = 6 s Should just shortcut the network...
      RD26M = (RHB/26.0)*N26M*CLN2/6.0
C Fe59 t_0.5 = 44.6 d
      RDFe59 = (RHB/59.0)*NFE59*CLN2/3.85d6
C Fe60 t_0.5 = 1.5 Myr
      RDFe60 = (RHB/60.0)*NFE60*CLN2/4.73d13
C Ni59 t_0.5 = 0.075 Myr
      RDNi59 = (RHB/59.0)*NNi59*CLN2/2.365d12
C Free n t_0.5 = 10.3 min
      RDn = (RHB/1.0)*Nn*CLN2/6.18d2
C (n,g) reactions
      TF = TL/CLN10
      DO J = 1, 45
         RN = 0.0D0
         TT = 50.0*(TF - 6.0) + 1.0
         IF ( TT .GE. 1.0D0 ) THEN
            IT = MAX(1, MIN(199, INT(TT)))
            TT = TT - IT
            TU = 1.0D0 - TT
            RR = TU*NRT(IT, J) + TT*NRT(IT+1, J)
            IF (RR.GE.-50.0d0) THEN
               RN = EXP(CLN10*(RR + 20.0D0))*1.0D-20
            END IF
         END IF
         NRATE(J) = RN
      END DO
      NRATE(1) = RHB*Nn*NFE56*NRATE(1)
      NRATE(2) = RHB*Nn*NFE57*NRATE(2)
      NRATE(3) = RHB*Nn*NFE58*NRATE(3)
      NRATE(4) = RHB*Nn*NCO59*NRATE(4)
      NRATE(5) = RHB*Nn*NNI58*NRATE(5)
      NRATE(6) = RHB*Nn*NNI59*NRATE(6)
      NRATE(7) = RHB*Nn*NNI60*NRATE(7)
      NRATE(8) = RHB*Nn*NN1*NRATE(8)
      NRATE(9) = RHB*Nn*NN3*NRATE(9)
      NRATE(10) = RHB*Nn*NNL7*NRATE(10)
      NRATE(11) = RHB*Nn*NN12*NRATE(11)
      NRATE(12) = RHB*Nn*NN13*NRATE(12)
      NRATE(13) = RHB*Nn*NN14*NRATE(13)
      NRATE(14) = RHB*Nn*NNN14*NRATE(14)
      NRATE(15) = RHB*Nn*NN15*NRATE(15)
      NRATE(16) = RHB*Nn*NN16*NRATE(16)
      NRATE(17) = RHB*Nn*NN18*NRATE(17)
      NRATE(18) = RHB*Nn*NN19*NRATE(18)
      NRATE(19) = RHB*Nn*NN20*NRATE(19)
      NRATE(20) = RHB*Nn*NNE21*NRATE(20)
      NRATE(21) = RHB*Nn*NNE22*NRATE(21)
      NRATE(22) = RHB*Nn*NNA23*NRATE(22)
      NRATE(23) = RHB*Nn*NMG24*NRATE(23)
      NRATE(24) = RHB*Nn*NMG25*NRATE(24)
      NRATE(25) = RHB*Nn*NMG26*NRATE(25)
      NRATE(26) = RHB*Nn*NAL27*NRATE(26)
      NRATE(27) = RHB*Nn*NSI28*NRATE(27)
      NRATE(28) = RHB*Nn*NSI29*NRATE(28)
      NRATE(29) = RHB*Nn*NSI30*NRATE(29)
      NRATE(30) = RHB*Nn*NP31*NRATE(30)
      NRATE(31) = RHB*Nn*NS32*NRATE(31)
      NRATE(32) = RHB*Nn*NS33*NRATE(32)
      NRATE(33) = RHB*Nn*NS34*NRATE(33)
      NRATE(34) = RHB*Nn*N26G*NRATE(34)
      NRATE(35) = RHB*Nn*NS33*NRATE(35)
      NRATE(36) = RHB*Nn*NNN14*NRATE(36)
      NRATE(37) = RHB*Nn*NNI59*NRATE(37)
      NRATE(38) = RHB*Nn*NNI59*NRATE(38)
      NRATE(39) = RHB*Nn*NN17*NRATE(39)
      NRATE(40) = RHB*Nn*N26G*NRATE(40)
      NRATE(41) = RHB*Nn*NS32*NRATE(41)
      NRATE(42) = RHB*Nn*NFE59*NRATE(42)
      NRATE(43) = RHB*Nn*NFE60*NRATE(43)
      NRATE(44) = RHB*Nn*NS34*NRATE(44)
      NRATE(45) = RHB*Nn*NNI61*NRATE(45)
C Check for negative rates
      DO I=1,45
         IF (NRATE(I).LT.0d0) THEN
            write (32,*) "-ve rate in ",I
            NRATE(I) = 0d0
         END IF
      END DO
      RETURN
      END
