**==FUNCS2.FOR
      SUBROUTINE FUNCS2(K, K1, K2, ITER, ISTAR)
      IMPLICIT REAL*8(A-H, O-Z)
      PARAMETER (MAXMSH = 2000)
      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      COMMON /AUXIN / ICL, ION, JW, IOP, INUC, IBC, ICN, IML(2), ISGTH,
     :     IMO, IDIFF
      COMMON /SODDS / W5(2), CMG, CSI, CFE, W8(11), DT, W6(23), IW(50), TRB
      COMMON /STAT2 / W1(4), RHO, W2(4), ZT, W3(9), RPP, R33, R34, RBE, RBP, 
     :                RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO,
     :                ROO, RGNE, RGMG, RCCG, RPNG, WR(21)
      COMMON /INF   / W(60)
      COMMON /DINF  / DW(60)
      COMMON /INFN  / WN(50)
      COMMON /DINFN / DWN(50)
C      COMMON /CEE   / MHC(2), MENVC(2), DSEP, ICE, ICEP, ALPHACE
      COMMON /OUTF  / YDOT(50), Y(50), SIG, DMT, SGTH, DMU, DSG(50)
      COMMON /DOUTF / YLIN(50, 154)
      COMMON /ABUND / XA(10), WW(14)
C       COMMON /YUK1  / PX(34), WMH, WMHE, VMH, VME, VMC, VMG, BE, VLH,
C      :                VLE, VLC, VLN, VLT, MCB(12),WWW(100)
      COMMON /ABUND2/ XA2(50), WW2(50)
      COMMON /TRANS / HT(28,MAXMSH,2)
      COMMON /SURFCO/ SURFCOMPOS(50), ISTAROTHER
      COMMON /RATES / W7(1),RPD,RPBE,RPLI,REBE,RPC13,RPN15A,RPO17a,
     :     RAC13n, RAN15g, RAO17n, RLIA, RPB11, RPC14, RAC14,
     :     RPO18,RPO18a,RAO18,RAO18n,RPF19,RPF19A,RAF19,RPNE21,
     :     RANE21,RANE21n,RPNE22,RANE22,RANE22n,RnNA22p,RnNA22a,RPNA22,
     :     RPNA23,RPNA23n,RPNA23a,RPMG24,RAMG24,RPMG25T,RPMG25M,RPMG25G,
     :     RAMG25,RAMG25n,RAMG25p,RPMG26,RPMG26Tn,RPMG26Mn,RPMG26Gn,
     :     RAMG26,RAMG26n,RAMG26p,Rg26Tp,Rn26Tp,Rn26Ta,Rp26T,Rg26Mp,
     :     Rn26Mp,Rn26Ma,Rp26M,Rg26Gp,Rn26Gp,Rn26Ga,Rp26G,RPAL27,
     :     RPAL27a,RAAL27n,RANA23nT,RANA23nM,RANA23nG,RPSI28,RPSI29,
     :     RPSI30,RPN15,RPO17,RPNE20
      COMMON /NRATES/ RnFe56,RnFe57,RnFe58,RnCo59,RnNi58,RnNi59,RnNi60,
     :     Rnp,RnHe3,RnLi7,RnC12,RnC13,RnC14,RnN14,RnN15,RnO16,RnO18,
     :     RnF19,RnNe20,RnNe21,RnNe22,RnNa23,RnMg24,RnMg25,RnMg26,
     :     RnAl27,RnSi28,RnSi29,RnSi30,RnP31,RnS32,RnS33,RnS34,Rn26G,
     :     RnS33a,RnN14p,RnNi59p,RnNi59a,RnO17a,Rn26Gad,RnS32a,RnFe59,
     :     RnFe60,RnS34s,RnNi61s
      COMMON /DIFFU2/ D(50), A12(50)
      COMMON /DECAY / RDAL26G,RDNA22,RD26M,RDFe59,RDFe60,RDNi59,RDn,RDC14
      COMMON /NEUTRO/ FRAC(50,MAXMSH)
      COMMON /OP    / ZS, VLEDD, VM, GR, GRAD, ETH, RLF, EGR, R, QQ
      DIMENSION BARYN(50), IZZ(50), M(2), R(2)
      DATA BARYN/1.0, 1.0, 2.0, 3.0, 7.0, 7.0, 11.0, 13.0, 14.0,  15.0,
     : 17.0, 18.0, 19.0, 21.0, 22.0, 22.0, 23.0, 24.0, 25.0, 26.0,
     : 26.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,56.0,57.0,58.0,
     : 59.0,60.0,59.0,58.0,59.0,60.0,61.0,1.0,4.0,12.0,14.0,16.0,20.0,1.0,1.0,1.0,1.0/ 
      DATA IZZ/30,0,1,2,3,4,5,6,6,7,8,8,9,10,10,11,11,12,12,12,13,13,13,14,14,14,15,
     :16,16,16,26,26,26,26,26,27,28,28,28,28,1,2,6,7,8,10,0,0,0,0/
C Set up boundary condition data
      IF (K.LE.K1) THEN
         IF (ISTAR.EQ.1) THEN
C Check for RLOF/Wind accretion and set BC appropriately
            ISTAROTHER = 2
            IF (HT(24,1,2).GT.0d0.OR.HT(23,1,1).LT.0d0) THEN
               DO I = 1,50
                  SURFCOMPOS(I) = HNUC(50+I,1)
               END DO
            ELSE
               DO I = 1,50
                  SURFCOMPOS(I) = WN(I)
               END DO
            END IF
         ELSE
            ISTAROTHER = 1
            IF (HT(24,1,1).GT.0d0.OR.HT(23,1,2).LT.0d0) THEN
               DO I = 1,50
                  SURFCOMPOS(I) = HNUC(I,1)
               END DO
            ELSE
               DO I = 1,50
                  SURFCOMPOS(I) = WN(I)
               END DO
            END IF
C NB - NO CONDITION FOR COMMON ENVELOPE
         END IF

         ! OLD new CEE prescription was here
      END IF
C Set up composition for reaction calculation
      XA(1) = W(5+15*(ISTAR - 1))
      XA(2) = W(9+15*(ISTAR - 1))
      XA(3) = W(10+15*(ISTAR - 1))
      XA(4) = W(12+15*(ISTAR - 1)) 
      XA(5) = W(3+15*(ISTAR - 1))
      XA(6) = W(11+15*(ISTAR - 1)) 
      XA(7) = CMG*ZS
      XA(8) = CSI*ZS
      XA(9) = CFE*ZS
C Set up thermodynamic variables
      AF = HT(21, K, ISTAR)
      AT = HT(22, K, ISTAR)
      DO I = 1,14
         WW(I) = HT(4+I, K, ISTAR)
      END DO
C Take SGTH and MU from evolution model
      IF (ISGTH.EQ.1) THEN
         SGTH = HT(19, K, ISTAR)
      ELSE
         SGTH = 0d0
      END IF
      DMU = HT(20, K, ISTAR)
C Blank any -ve abundances. There shouldn't be any anyway...
C RJS 19/1/04
      DO I = 1, 9
         IF (XA(I).LT.0d0) THEN
            XA(I) = 0d0
            WW(I) = 0d0
         END IF
      END DO
      DO I = 1, 50
         XA2(I) = WN(I)
         IF (XA2(I).LT.0d0) XA2(I) = 0d0
C Neutron kill
         XA2(2) = 0d0
      END DO
C Proton fudging
      IF (XA(1).LT.1d-10) XA2(41) = 0d0
C Work out abundances for nucrat2 - RJS
      DO I = 1, 50
         WW2(I) = XA2(I)/BARYN(I)
      END DO
      DO I = 1, 50
         DO J = 1, 50
            YLIN(I, J) = 0d0
         END DO
      END DO
      DAM = DW(4+15*(ISTAR - 1)) 
      RHO = HT(1, K, ISTAR)
      SIG = HT(2, K, ISTAR)
      ZT = HT(3, K, ISTAR)
      AM = W(4+15*(ISTAR - 1))
      DM = DEXP(AM)
C Pierre's mass mod.
      DMT = (DM - DEXP(AM-DAM))/DT
      DMK = HT(4, K, ISTAR)
      W1(4) = 0.0
      CALL NUCRAT2(AT)
      IF (IDIFF.EQ.1) THEN
C Compute gravitational settling coefficients 19/9/08 RJS
         T = DEXP(AT)
C Factors need to compute settling coefficients
         GKT = HT(25, K, ISTAR)
         PIRHOR2 = HT(26, K, ISTAR)
         DMUREAL = HT(27, K, ISTAR)
         ATM = HT(28, K, ISTAR)
         CALL DIFFUSION2(RHO,T)
C     n's aren't charged and don't diffuse
         D(2) = 0d0
         A12(1) = 0d0
         A12(2) = 0d0
C can't diffuse through H if we don't have any...
         DO II=1,46
            IF (XA2(II).EQ.0d0) D(II) = 0d0
            IF (XA2(41).LT.0.1) D(II) = 0d0
            IF (XA2(41).LT.0.1) A12(II) = 0d0
            IF (D(II).NE.D(II)) THEN
               WRITE(*,*) "NaN in element", II
               STOP
            END IF
            IF (A12(II).NE.A12(II)) THEN
               WRITE(*,*) "NaN in element, A12 ", II, K
               STOP
            END IF
C Assuming all species full ionized -- not this is not consistent with the
C structure code
            A12(II) = A12(II)*D(II)*ATM*PIRHOR2
            D(II) = D(II)*GKT*(BARYN(II)/(1d0+IZZ(II)) - DMU)
            D(II) = D(II) - A12(II)
            D(II) = D(II)*PIRHOR2
            DSG(II) = D(II)
         END DO
      ELSE
         DO II = 1,50
            DSG(II) = 0d0
         END DO
      END IF
C  EVALUATE COMPOSITION CHANGE
C 1/9/03 Minor element reactions from Clayton 1963
      YDOT(1) = 0d0
      YDOT(2) =  RDn - RAC13n - RAO17n - RANE21n - RAO18n  - RANE22n  + 
     :     RnNA22p + RnNA22a - RPNA23n - RAMG25n - RPMG26Mn - RPMG26Gn
     :     - RAMG26n + Rn26Mp + Rn26Ma + Rn26G
     :     + Rn26Gp + Rn26Ga - RAAL27n - RANA23nM - RANA23nG
     :     + RnFe56 + RnFe57 + RnFe58 + RnCo59 + RnNi58 + RnNi59
     :     + RnNi60 + Rnp + RnHe3 + RnLi7 + RnC12 + RnC13 + RnC14
     :     + RnN14 + RnN15 + RnO16 + RnO18 + RnF19 + RnNe20 + RnNe21
     :     + RnNe22 + RnNa23 + RnMg24 + RnMg25 + RnMg26 + RnAl27
     :     + RnSi28 + RnSi29 + RnSi30 + RnP31 + RnS32 + RnS33 + RnS34
     :     + RnS33a + RnN14p + RnNi59p + RnNi59a + RnO17a + RnS32a
      YDOT(3) = RPD - RPP - Rnp
      YDOT(4) = - RPD + 2.0*R33 + R34 + RnHe3
      YDOT(5) =  RPLI - REBE + RLIA + RnLi7
      YDOT(6) = -R34 + RPBE + REBE
      YDOT(7) = RPB11 - RLIA
      YDOT(8) = RPC13 - RPC + RAC13n + RnC13 - RnC12
      YDOT(9) = RPC14 + RAC14 + RnC14 - RnN14p - RnO17a - RnC13 + RDC14
      YDOT(10)= RPN15A - RPN + RAN15g - RPO18a + RnN15 -  RnC14 - RnN14
     :     + RPN15
      YDOT(11)= - RPO + RPO17a + RAO17n - RnO16 + RnO17a + RPO17
      YDOT(12)= - RAN - RAC14 + RPO18 + RPO18a + RAO18 + RAO18n + RnO18
     :     - RPO17
      YDOT(13)= - RAN15g + RPF19 + RPF19A + RAF19 - RnNA22a - RnO18 
     :     + RnF19 - RPO18
      YDOT(14)= - RAO18n + RPNE21 + RANE21 + RANE21n + RnNe21 - RnNe20
     :     - RPNE20
      YDOT(15)= - RAO18 - RAF19 + RPNE22 + RANE22 + RANE22n - RnNA22p
     :     - RDNA22 - RnNe21 + RnNe22
      YDOT(16)= RnNA22p + RnNA22a + RPNA22 - RPNE21 + RDNA22
      YDOT(17)= RPNA23 + RPNA23n + RPNA23a - RPNE22 - Rn26Ga 
     :     - Rn26Ma + RANA23nM + RANA23nG - RnNe22 + RnNa23
      YDOT(18)= RPMG24 + RAMG24 - RCC - RANE - RANE21n - RPNA23 
     :     - RPAL27a - RnNa23 + RnMg24
      YDOT(19)= RPMG25M + RPMG25G + RAMG25 + RAMG25n + RAMG25p
     :     - RANE21 - RANE22n - Rg26Mp - Rg26Gp + RnMg25 - RnMg24
      YDOT(20)= RPMG26 + RPMG26Mn + RPMG26Gn  
     :     + RAMG26 + RAMG26n + RAMG26p - RANE22 - Rn26Mp 
     :     - Rn26Gp - RDAL26G - RD26M + RnMg26 - RnMg25
      YDOT(21)= Rg26Mp + Rn26Mp + Rn26Ma + Rp26M - RPMG25M - RPMG26Mn
     :     - RANA23nM + RD26M
      YDOT(22)= Rg26Gp + Rn26Gp + Rn26Ga + Rp26G - RPMG25G - RPMG26Gn
     :     - RANA23nG + RDAL26G + Rn26G
      YDOT(23)= RPAL27 + RPAL27a + RAAL27n - RPMG26 - RnMg26 + RnAl27
     :     - Rn26G
C assume decay of Si28(p,g)P29 -> Si29
      YDOT(24)= RPSI28 - RCO - RAMG24 - RAMG25n - RPAL27 - RnAl27 + 
     :     RnSi28
C assume decay of Si29(p,g)P30 -> Si30
      YDOT(25)= RPSI29 - RAMG25 - RAMG26n - RPSI28 - RnSi28 + RnSi29
     :          - RnS32a
      YDOT(26)= RPSI30 - RAMG26 - RPSI29 - RnSi29 + RnSi30 - RnS33a
      YDOT(27)= - RPSI30 - RnSi30 + RnP31
      YDOT(28)= - RnP31 + RnS32 + RnS32a
      YDOT(29)= - RnS32 + RnS33 + RnS33a
      YDOT(30)= - RnS33 + RnS34 + RnS34s
      YDOT(31)= RnFe56 - RnNi59a
      YDOT(32)= RnFe57 - RnFe56
      YDOT(33)= RnFe58 - RnFe57
      YDOT(34)= - RnFe58 + RnFe59 + RDFe59
      YDOT(35)= - RnFe59 + RnFe60 + RDFe60
      YDOT(36)= RnCo59 - RDFe59 - RnNi59p - RDNi59
      YDOT(37)= RnNi58 
      YDOT(38)= RnNi59 - RnNi58 + RnNi59a + RnNi59p + RDNi59
      YDOT(39)= RnNi60 - RnNi59
      YDOT(40)= RnNi61s - RnNi60
C Protons duplicate
      YDOT(41)= 2.0*RPP + RPD + Rnp - 2.0*R33 + RPLI + RPBE + RPB11
     :     + RPC + RPC13 + RPC14 + RPN + RPN15A + RPO18a + RPO
     :     + RPO17a - RDn + RPO18 + RPF19 + RPF19A + RPNE21 
     :     + RPNE22 - RnNa22p + RPNA22 + RPNA23 + RPNA23n + RPNA23a
     :     + RPMG24 + RPAL27a + RPMG25M + RPMG25G - RAMG25p + RPMG26 
     :     + RPMG26Mn + RPMG26Gn + RPN15
     :     + Rp26M + Rp26G + RPAL27 + RPSI28 + RPSI29 + RPSI30
     :     + RPO17 + RPNE20
C Helium duplicate - RPO17a is O17(p,a)N14
      YDOT(42)= - R33 + R34 - 2.0*RPLI - 2.0*RPBE + RLIA + RAC13n + RAN
     :     - RPN15A - RPO18a + RAC14 + RAN15g + RAO17n + RAO + 3.0*R3A 
     :     + RAC + RAO18 + RAO18n - RPO17a - RPF19A + RAF19 - RnNa22a !- RCC
     :     + RANE21 + RANE21n + RANE22 + RANE22n - RnNa22a - RPNA23a 
     :     - Rn26Ga - Rn26Ma + RANA23nM + RANA23nG + RANE + RAMG24 
     :     - RPAL27a + RAMG25 + RAMG25n + RAMG25p + RAMG26 + RAMG26n 
     :     + RAMG26p + RAAL27n - RnS32a - RnS33a - RnNi59a
C Carbon duplicate - RCC to alpha and Ne20 - should be ~50-50 with 
C this and Na23 + p ! CHEATED! sending all RCC to Mg24
      YDOT(43)= - RPB11 + RPC - RPN15A - R3A + RAC + 2.0*RCC + RnC12 + RCO     
C Nitrogen duplicate
      YDOT(44)= - RPC13 + RnN14p + RPN + RAN - RPO17a
C Oxygen duplicate
      YDOT(45)= - RAC13n + RPO + RAO - RAC + RnO16 - RPF19A + RCO-RPN15
C Neon duplicate
      YDOT(46)= - RAO - RPF19 - RnF19 + RnNe20 + RANE - RPNA23a + RPNE20! - RCC
C  EVALUATE DERIVATIVES
C YLIN(derived by,equation)
C RJS 1/9/03 - Added species for PP-I,II chain from Clayton 1963
C 1st eq gallinoes

C 2nd eq neutrons
      YLIN(2, 2) = (RDn + RnNA22p + RnNA22a + Rn26Mp + Rn26Ma
     :              + Rn26Gp + Rn26Ga + RnFe56 + RnFe57 + RnFe58
     :              + RnCo59 + RnNi58 + RnNi59 + RnNi60 + Rnp
     :              + RnHe3 + RnLi7 + RnC12 + RnC13 + RnC14 + RnN14
     :              + RnN15 + RnO16 + RnO18 + RnF19 + RnNe20 + RnNe21
     :              + RnNe22 + RnNa23 + RnMg24 + RnMg25 + RnMg26
     :              + RnAl27 + RnSi28 + RnSi29 + RnSi30 + RnP31 + RnS32
     :              + RnS33 + RnS34 + Rn26G+ RnS33a + RnN14p + RnNi59p
     :              + RnNi59a + RnO17a + RnS32a)/XA2(2)
      YLIN(4, 2) = RnHe3/XA2(4)
      YLIN(5, 2) = RnLi7/XA2(5)
      YLIN(8, 2) = (RnC13- RAC13n)/XA2(8) 
      YLIN(9, 2) = RnC14/XA2(9)
      YLIN(10,2) = RnN15/XA2(10)
      YLIN(11,2) = (RnO17a - RAO17n)/XA2(11)
      YLIN(12,2) = (RnO18 - RAO18n)/XA2(12)
      YLIN(13,2) = RnF19/XA2(13)
      YLIN(14,2) = (RnNe21 - RANE21n)/XA2(14) 
      YLIN(15,2) = (RnNe22 - RANE22n)/XA2(15)
      YLIN(16,2) = (RnNA22p + RnNA22a)/XA2(16)
      YLIN(17,2) = (RnNa23 - RPNA23n - RANA23nM - RANA23nG)/XA2(17)
      YLIN(18,2) = RnMg24/XA2(18)
      YLIN(19,2) = (RnMg25 - RAMG25n)/XA2(19)
      YLIN(20,2) = (RnMg26 - RAMG26n - RPMG26Mn - RPMG26Gn)/XA2(20) 
      YLIN(21,2) = (Rn26Mp + Rn26Ma)/XA2(21)
      YLIN(22,2) = (Rn26Gp + Rn26Ga + Rn26G)/XA2(22)
      YLIN(23,2) = (- RAAL27n + RnAl27)/XA2(23)
      YLIN(24,2) = RnSi28/XA2(24)
      YLIN(25,2) = RnSi29/XA2(25)
      YLIN(26,2) = RnSi30/XA2(26)
      YLIN(27,2) = RnP31/XA2(27)
      YLIN(28,2) = (RnS32 + RnS32a)/XA2(28)
      YLIN(29,2) = (RnS33 + RnS33a)/XA2(29)
      YLIN(30,2) = RnS34/XA2(30)
      YLIN(31,2) = RnFe56/XA2(31)
      YLIN(32,2) = RnFe57/XA2(32)
      YLIN(33,2) = RnFe58/XA2(33)
      YLIN(36,2) = RnCo59/XA2(36)
      YLIN(37,2) = RnNi58/XA2(37)
      YLIN(38,2) = (RnNi59 + RnNi59p + RnNi59a)/XA2(38)
      YLIN(39,2) = RnNi60/XA2(39)
      YLIN(41,2) = (Rnp - RPNA23n - RPMG26Mn - RPMG26Gn)/XA2(41)
      YLIN(42,2) = - (RAC13n + RAO17n + RAO18n + RANE21n + RANE22n 
     :     + RANA23nM + RANA23nG + RAMG25n + RAMG26n + RAAL27n)/XA2(42)
      YLIN(43,2) = RnC12/XA2(43)
      YLIN(45,2) = RnO16/XA2(45)
      YLIN(46,2) = RnNe20/XA2(46)
C 3rd eq D
      YLIN(2, 3) = - Rnp/XA2(2)
      YLIN(3, 3) = RPD/XA2(3)
      YLIN(41,3) = (RPD - Rnp - 2.0*RPP)/XA2(41)
C 4th eq He-3
      YLIN(2, 4) = RnHe3/XA2(2)
      YLIN(3, 4) = -RPD/XA2(3)
      YLIN(4, 4) = (4.0*R33+R34+RnHe3)/XA2(4)
      YLIN(41,4) = - RPD/XA2(41)
      YLIN(42,4) = R34/XA2(42)
C 5th eq Li-7
      YLIN(2, 5) = RnLi7/XA2(2)
      YLIN(5, 5) = (RPLI+RLIA+RnLi7)/XA2(5)
      YLIN(6, 5) = -REBE/XA2(6)
      YLIN(41,5) = RPLI/XA2(41)
      YLIN(42,5) = RLIA/XA2(42)
C 6th eq Be-7
      YLIN(4, 6) = -R34/XA2(4)
      YLIN(6, 6) = (RPBE + REBE)/XA2(6)
      YLIN(41,6) = RPBE/XA2(41)
      YLIN(42,6) = - R34/XA2(42)
C 7th eq B-11
      YLIN(5, 7) = - RLIA/XA2(5)
      YLIN(7, 7) = RPB11/XA2(7)
      YLIN(41,7) = RPB11/XA2(41)
      YLIN(42,7) = - RLIA/XA2(42)
C CNO elements. Will assume beta-decay is instantaneous (which is ok on all but the very shortest timescales)
C 8th eq C-13
      YLIN(2, 8) = (RnC13 - RnC12)/XA2(2)
      YLIN(8, 8) = (RPC13+RAC13n+RnC13)/XA2(8)
      YLIN(41,8) = (RPC13 - RPC)/XA2(41)
      YLIN(42,8) = RAC13n/XA2(42)
      YLIN(43,8) = (- RPC - RnC12)/XA2(43)
C 9th eq C-14
      YLIN(2, 9) = (RnC14 - RnC13 - RnN14p - RnO17a)/XA2(2)
      YLIN(8, 9) = - RnC13/XA2(8)
      YLIN(9, 9) = (RPC14 + RAC14 + RnC14 + RDC14)/XA2(9)
      YLIN(11,9) = - RnO17a/XA2(11)
      YLIN(41,9) = RPC14/XA2(41)
      YLIN(42,9) = RAC14/XA2(42)
      YLIN(44,9) = - RnN14p/XA2(44)
C 10th eq N-15
C Assume instantaneous C15 -> N15 beta decay
      YLIN(2 ,10) = (RnN15 - RnC14 - RnN14)/XA2(2)
      YLIN(9 ,10) = -RnC14/XA2(9)
      YLIN(10,10) = (RPN15A+RAN15g+RnN15+RPN15)/XA2(10)
      YLIN(12,10) = - RPO18a/XA2(12)
      YLIN(41,10) = (RPN15A - RPN - RPO18a + RPN15)/XA2(41)
      YLIN(42,10) = RAN15g/XA2(42)
      YLIN(44,10) = (- RnN14 - RPN)/XA2(44)
C 11th eq O-17
C Assume O17(p,g)->O18 via beta decay of F18
      YLIN(2 ,11) = (RnO17a - RnO16)/XA2(2)
      YLIN(11,11) = (RPO17a+RAO17n+RnO17a+RPO17)/XA2(11)
      YLIN(41,11) = (RPO17a+RPO17)/XA2(41)
      YLIN(42,11) = RAO17n/XA2(42)
      YLIN(45,11) = - RnO16/XA2(45)
C 12th eq O18
      YLIN(2, 12) = RnO18/XA2(2)
      YLIN(9 ,12) = - RAC14/XA2(9)
      YLIN(11,12) = - RPO17/XA2(11)
      YLIN(12,12) = (RPO18 + RPO18a + RAO18 + RAO18n + RnO18)/XA2(12)
      YLIN(41,12) = (RPO18 + RPO18a - RPO17)/XA2(41)
      YLIN(42,12) = (RAO18 + RAO18n - RAN)/XA2(42)
      YLIN(44,12) = - RAN/XA2(44)
C 13th eq F19
      YLIN(2 ,13) = (RnF19 - RnNA22a - RnO18)/XA2(2)
      YLIN(10,13) = - RAN15g/XA2(10)
      YLIN(12,13) = - (RnO18 + RPO18)/XA2(12)
      YLIN(13,13) = (RPF19 + RPF19A + RAF19 + RnF19)/XA2(13)
      YLIN(16,13) = - RnNA22a/XA2(16)
      YLIN(41,13) = (RPF19 + RPF19A - RPO18)/XA2(41)
      YLIN(42,13) = (RAF19 - RAN15g)/XA2(42)
C 14th eq Ne21
      YLIN(2, 14) = (RnNe21-RnNe20)/XA2(2)
      YLIN(12,14) = - RAO18n/XA2(12)
      YLIN(14,14) = (RPNE21 + RANE21 + RANE21n + RnNe21)/XA2(14)
      YLIN(41,14) = (RPNE21 - RPNE20)/XA2(41)
      YLIN(42,14) = (RANE21 + RANE21n - RAO18n)/XA2(42)
      YLIN(46,14) = -RPNE20/XA2(46)
C 15th eq Ne22
      YLIN(2 ,15) = (RnNe22 - RnNe21 - RnNA22p)/XA2(2)
      YLIN(12,15) = - RAO18/XA2(12)
      YLIN(13,15) = - RAF19/XA2(13)
      YLIN(14,15) = -RnNe21/XA2(14)
      YLIN(15,15) = (RPNE22 + RANE22 + RANE22n + RnNe22)/XA2(15)
      YLIN(16,15) = - (RnNA22p + RDNA22)/XA2(16)
      YLIN(41,15) = RPNE22/XA2(41)
      YLIN(42,15) = (RANE22 + RANE22n - RAO18 - RAF19)/XA2(42)
C 16th eq Na22
      YLIN(2 ,16) = (RnNA22p + RnNA22a)/XA2(2)
      YLIN(14,16) = - RPNE21/XA2(14)
      YLIN(16,16) = (RnNA22p + RnNA22a + RPNA22 + RDNA22)/XA2(16)
      YLIN(41,16) = (RPNA22 - RPNE21)/XA2(41)
C 17th eq Na23
      YLIN(2 ,17) = (RnNa23 - Rn26Ma - Rn26Ga - RnNe22)/XA2(2)
      YLIN(15,17) = - (RPNE22+RnNe22)/XA2(15)
      YLIN(17,17) = (RPNA23 + RPNA23n + RPNA23a + RANA23nM
     :     + RANA23nG + RnNa23)/XA2(17)
      YLIN(21,17) = - Rn26Ma/XA2(21)
      YLIN(22,17) = - Rn26Ga/XA2(22)
      YLIN(41,17) = (RPNA23 + RPNA23n + RPNA23a - RPNE22)/XA2(41)
      YLIN(42,17) = (RANA23nM + RANA23nG)/XA2(42)
C 18th eq Mg24
      YLIN(2 ,18) = (RnMg24 - RnNa23)/XA2(2)
      YLIN(14,18) = - RANE21n/XA2(14)
      YLIN(17,18) = (- RnNa23 - RPNA23)/XA2(17)
      YLIN(18,18) = (RPMG24 + RAMG24 + RnMg24)/XA2(18)
      YLIN(23,18) = - RPAL27a/XA2(23)
      YLIN(41,18) = (RPMG24 - RPAL27a - RPNA23)/XA2(41)
      YLIN(42,18) = (RAMG24 - RANE - RANE21n)/XA2(42)
      YLIN(46,18) = - RANE/XA2(46)
C 19th eq Mg25
      YLIN(2 ,19) = (RnMg25 - RnMg24)/XA2(2) 
      YLIN(14,19) = - RANE21/XA2(14)
      YLIN(15,19) = - RANE22n/XA2(15)
      YLIN(18,19) = - RnMg24/XA2(18)
      YLIN(19,19) = (RPMG25M + RPMG25G + RAMG25 + RAMG25n 
     :      + RAMG25p + RnMg25)/XA2(19)
      YLIN(21,19) = - Rg26Mp/XA2(21)
      YLIN(22,19) = - Rg26Gp/XA2(22)
      YLIN(41,19) = (RPMG25M + RPMG25G)/XA2(41)
      YLIN(42,19) = (RAMG25 + RAMG25n + RAMG25p)/XA2(42)
C 20th eq Mg26
      YLIN(2 ,20) = (RnMg26 - RnMg25 - Rn26Mp - Rn26Gp)/XA2(2)
      YLIN(15,20) = - RANE22/XA2(15)
      YLIN(19,20) = - RnMg25/XA2(19)
      YLIN(20,20) = (RPMG26 + RPMG26Mn + RPMG26Gn  
     :     + RAMG26 + RAMG26n + RAMG26p + RnMg26)/XA2(20)
      YLIN(21,20) = - (Rn26Mp + RD26M)/XA2(21)
      YLIN(22,20) = - (Rn26Gp + RDAL26G)/XA2(22)
      YLIN(41,20) = (RPMG26 + RPMG26Mn + RPMG26Gn)/XA2(41)
      YLIN(42,20) = (RAMG26 + RAMG26n + RAMG26p)/XA2(42)
C 21st eq Al26M
      YLIN(2 ,21) = (Rn26Mp + Rn26Ma)/XA2(2)
      YLIN(17,21) = - RANA23nM/XA2(17)
      YLIN(19,21) = - RPMG25M/XA2(19)
      YLIN(20,21) = - RPMG26Mn/XA2(20)
      YLIN(21,21) = (Rg26Mp + Rn26Mp + Rn26Ma + Rp26M + RD26M)/XA2(21)
      YLIN(41,21) = (Rp26M - RPMG25M - RPMG26Mn)/XA2(41)
      YLIN(42,21) = - RANA23nM/XA2(42)
C 22nd eq Al26G
      YLIN(2 ,22) = (Rn26Gp + Rn26Ga + Rn26G)/XA2(2)
      YLIN(17,22) = - RANA23nG/XA2(17)
      YLIN(19,22) = - RPMG25G/XA2(19) 
      YLIN(20,22) = - RPMG26Gn/XA2(20)
      YLIN(22,22) = (Rg26Gp + Rn26Gp + Rn26Ga + Rp26G + RDAL26G
     :     + Rn26G)/XA2(22)
      YLIN(41,22) = (Rp26G - RPMG25G - RPMG26Gn)/XA2(41)
      YLIN(42,22) = - RANA23nG/XA2(42)
C 23rd eq Al27
      YLIN(2 ,23) = (RnAl27 - RnMg26 - Rn26G)/XA2(2)
      YLIN(20,23) = (-RnMg26 - RPMG26)/XA2(20)
      YLIN(22,23) = -Rn26G/XA2(22)
      YLIN(23,23) = (RPAL27 + RPAL27a + RAAL27n + RnAl27)/XA2(23)
      YLIN(41,23) = (RPAL27 + RPAL27a - RPMG26 - Rp26G - Rp26M)/XA2(41)
      YLIN(42,23) = RAAL27n/XA2(42)
C 24th eq Si28
      YLIN(2, 24) = (RnSi28 - RnAl27)/XA2(2)
      YLIN(18,24) = - RAMG24/XA2(18) 
      YLIN(19,24) = - RAMG25n/XA2(19)
      YLIN(23,24) = (-RnAl27 - RPAL27)/XA2(23)
      YLIN(24,24) = (RnSi28 + RPSI28)/XA2(24) 
      YLIN(41,24) = (RPSI28 - RPAL27)/XA2(41)
      YLIN(42,24) = (- RAMG24 - RAMG25n)/XA2(42)
      YLIN(43,24) = - RCO/XA2(43)
C 25th eq Si29
      YLIN(2 ,25) = (RnSi29 - RnSi28 - RnS32a)/XA2(2)
      YLIN(19,25) = - RAMG25/XA2(19)
      YLIN(20,25) = - RAMG26n/XA2(20)
      YLIN(24,25) = - (RnSi28 + RPSI28)/XA2(24)
      YLIN(25,25) = (RPSI29 + RnSi29)/XA2(25)
      YLIN(28,25) = - RnS32a/XA2(28)
      YLIN(41,25) = (RPSI29 - RPSI28)/XA2(41)
      YLIN(42,25) = (- RAMG25 - RAMG26n)/XA2(42)
C 26th eq Si30
      YLIN(2 ,26) = (RnSi30 - RnSi29 - RnS33a)/XA2(2)
      YLIN(20,26) = - RAMG26/XA2(20)
      YLIN(25,26) = - (RPSI29 + RnSi29)/XA2(25)
      YLIN(26,26) = (RPSI30+RnSi30)/XA2(26)
      YLIN(29,26) = - RnS33a/XA2(29)
      YLIN(41,26) = (RPSI30 - RPSI29)/XA2(41)
      YLIN(42,26) = - RAMG26/XA2(42)
C 27th eq P31
      YLIN(2 ,27) = (RnP31 - RnSi30)/XA2(2)
      YLIN(26,27) = - (RPSI30 + RnSi30)/XA2(26)
      YLIN(27,27) = RnP31/XA2(27)
      YLIN(41,27) = - RPSI30/XA2(41)
C 28th eq S32
      YLIN(2 ,28) = (RnS32 - RnP31 + RnS32a)/XA2(2)
      YLIN(27,28) = - RnP31/XA2(27)
      YLIN(28,28) = (RnS32 + RnS32a)/XA2(28)
C 29th eq S33
      YLIN(2 ,29) = (RnS33 - RnS32 + RnS33a)/XA2(2)
      YLIN(28,29) = - RnS32/XA2(28)
      YLIN(29,29) = (RnS33 + RnS33a)/XA2(29) 
C 30th eq S34
      YLIN(2 ,30) = (RnS34 - RnS33 + RnS34s)/XA2(2)
      YLIN(29,30) = - RnS33/XA2(29)
      YLIN(30,30) = (RnS34 + RnS34s)/XA2(30) 
C 31st eq Fe56
      YLIN(2 ,31) = (RnFe56 - RnNi59a)/XA2(2)
      YLIN(31,31) = RnFe56/XA2(31)
      YLIN(38,31) = -RnNi59a/XA2(38)
C 32nd eq Fe57
      YLIN(2 ,32) = (RnFe57 - RnFe56)/XA2(2)
      YLIN(31,32) = - RnFe56/XA2(31)
      YLIN(32,32) = RnFe57/XA2(32)
C 33rd eq Fe58
      YLIN(2 ,33) = (RnFe58-RnFe57)/XA2(2)
      YLIN(32,33) = - RnFe57/XA2(32)
      YLIN(33,33) = RnFe58/XA2(33)
C 34th eq Fe59
      YLIN(2 ,34) = (RnFe59 - RnFe58)/XA2(2)
      YLIN(33,34) = - RnFe58/XA2(33)
      YLIN(34,34) = (RnFe59 + RDFe59)/XA2(34)
C 35th eq Fe60
      YLIN(2 ,35) = (RnFe60 - RnFe59)/XA2(2)
      YLIN(34,35) = - RnFe59/XA2(34)
      YLIN(35,35) = (RnFe60 + RDFe60)/XA2(35)
C 36th eq Co59
      YLIN(2 ,36) = (RnCo59 - RnNi59p)/XA2(2)
      YLIN(34,36) = - RDFe59/XA2(34)
      YLIN(36,36) = RnCo59/XA2(36)
      YLIN(38,36) = (-RnNi59p - RDNi59)/XA2(38)
C 37th eq Ni58
      YLIN(2 ,37) = RnNi58/XA2(2)
      YLIN(37,37) = RnNi58/XA2(37)
C 38th eq Ni59
      YLIN(2 ,38) = (RnNi59 - RnNi58 + RnNi59p + RnNi59a)/XA2(2)
      YLIN(37,38) = - RnNi58/XA2(37)
      YLIN(38,38) = (RnNi59+ RnNi59p + RnNi59a + RDNi59)/XA2(38)
C 39th eq Ni60
      YLIN(2 ,39) = (RnNi60 - RnNi59)/XA2(2)
      YLIN(38,39) = - RnNi59/XA2(38)
      YLIN(39,39) = RnNi60/XA2(39)
C 40th eq Ni61
      YLIN(2 ,40) = - RnNi60/XA2(2)
      YLIN(39,40) = - RnNi60/XA2(39)
C 41st eq H - duplicate
      YLIN(2 ,41) = (Rnp - RnN14p - RDn - Rn26Mp - Rn26Gp - RnNi59p)/XA2(2)
      YLIN(3 ,41) = RPD/XA2(3)
      YLIN(4 ,41) = - 4.0*R33/XA2(4)
      YLIN(5 ,41) = RPLI/XA2(5)
      YLIN(6 ,41) = RPBE/XA2(6)
      YLIN(7 ,41) = RPB11/XA2(7)
      YLIN(8 ,41) = RPC13/XA2(8)
      YLIN(9 ,41) = RPC14/XA2(9)
      YLIN(10,41) = (RPN15A+RPN15)/XA2(10)
      YLIN(11,41) = (RPO17a+RPO17)/XA2(11)
      YLIN(12,41) = RPO18a/XA2(12)
      YLIN(13,41) = (RPF19 + RPF19A - RAF19)/XA2(13)
      YLIN(14,41) = RPNE21/XA2(14)
      YLIN(15,41) = RPNE22/XA2(15)
      YLIN(16,41) = RPNA22/XA2(16)
      YLIN(17,41) = (RPNA23 + RPNA23n + RPNA23a)/XA2(17)
      YLIN(18,41) = RPMG24/XA2(18)
      YLIN(19,41) = (RPMG25M + RPMG25G - RAMG25p)/XA2(19)
      YLIN(20,41) = (RPMG26 + RPMG26Mn + RPMG26Gn - RAMG26p)/XA2(20)
      YLIN(21,41) = (Rp26M - Rn26Mp - Rg26Mp)/XA2(21)
      YLIN(22,41) = (Rp26G - Rn26Gp - Rg26Gp)/XA2(22)
      YLIN(23,41) = (RPAL27 + RPAL27a)/XA2(23)
      YLIN(24,41) = RPSI28/XA2(24)
      YLIN(25,41) = RPSI29/XA2(25)
      YLIN(26,41) = RPSI30/XA2(26)
      YLIN(38,41) = - RnNi59p/XA2(38)
      YLIN(41,41) = (4.0*RPP + RPD + Rnp + RPLI + RPBE + RPB11
     :     + RPC13 + RPC14 + RPN + RPN15A + RPO18a + RPO + RPO17a + RPF19 
     :     + RPF19A + RPNE21 + RPNE22 + RPNA22 + RPNA23 + RPNA23n 
     :     + RPNA23a + RPMG24 + RPAL27a + RPMG25M + RPMG25G + RPMG26  
     :     + RPMG26Mn + RPMG26Gn + Rp26M + Rp26G + RPAL27 + RPSI28
     :     + RPSI29 + RPSI30 + RPN15 + RPO17 + RPNE20)/XA2(41)
      YLIN(42,41) = - (RAF19 + RAMG25p + RAMG26p)/XA2(42)
      YLIN(43,41) = RPC/XA2(43)
      YLIN(44,41) = (RPN - RnN14p)/XA2(44)
      YLIN(45,41) = RPO/XA2(45)
      YLIN(46,41) = RPNE20/XA2(46)
C 42nd eq He4 - duplicate
      YLIN(2 ,42) = (- Rn26Ga - Rn26Ma - RnNa22a - RnS32a - RnS33a
     :     - RnNi59a)/XA2(2)
      YLIN(4 ,42) = (- 4.0*R33 + R34)/XA2(4)
      YLIN(5 ,42) = (RLIA - 2.0*RPLI)/XA2(5)
      YLIN(6 ,42) = -2.0*RPBE/XA2(6)
      YLIN(8 ,42) = RAC13n/XA2(8)
      YLIN(9 ,42) = RAC14/XA2(9)
      YLIN(10,42) = (RAN15g - RPN15A)/XA2(10)
      YLIN(11,42) = (RAO17n - RPO17a)/XA2(11)
      YLIN(12,42) = - RPO18a/XA2(12)
      YLIN(13,42) = (RAF19 - RPF19A)/XA2(13)
      YLIN(14,42) = (RANE21 + RANE21n)/XA2(14)
      YLIN(15,42) = (RANE22 + RANE22n)/XA2(15)
      YLIN(16,42) = - RnNa22a/XA2(16)
      YLIN(17,42) = (RANA23nM + RANA23nG)/XA2(17)
      YLIN(18,42) = RAMG24/XA2(18)
      YLIN(19,42) = (RAMG25 + RAMG25n + RAMG25p)/XA2(19)
      YLIN(20,42) = (RAMG26 + RAMG26n + RAMG26p)/XA2(20)
      YLIN(21,42) = - Rn26Ma/XA2(21)
      YLIN(22,42) = - Rn26Ga/XA2(22)
      YLIN(23,42) = (RAAL27n - RPAL27a)/XA2(23)
      YLIN(28,42) = - RnS32a/XA2(28)
      YLIN(29,42) = - RnS33a/XA2(29)
      YLIN(38,42) = - RnNi59a/XA2(38)
      YLIN(41,42) = (-2.0*RPLI - 2.0*RPBE - RPN15A - RPO18a - RPO17a
     :     - RPF19A - RPNA23a - RPAL27a)/XA2(41)
      YLIN(42,42) = (R34+RLIA+RAC13n+RAN + RAC14 + RAN15g + RAO17n
     :     + RAO + 9.0*R3A + RAC + RAF19 + RANE21 + RANE21n + RANE22
     :     + RANE22n + RANA23nM + RANA23nG + RANE + RAMG24 + RAMG25 
     :     + RAMG25n + RAMG25p + RAMG26 + RAMG26n + RAMG26p + RAAL27n)
     :     /XA2(42)
      YLIN(43,42) = (RAC)/XA2(43) ! - 2.0*RCC
      YLIN(44,42) = RAN/XA2(44)
      YLIN(45,42) = RAO/XA2(45)
      YLIN(46,42) = RANE/XA2(46)
C 43rd eq C12 - duplicate
      YLIN(7 ,43) = - RPB11/XA2(7)
      YLIN(10,43) = - RPN15A/XA2(10)
      YLIN(41,43) = (RPC - RPB11 - RPN15A)/XA2(41)
      YLIN(42,43) = (RAC - 9.0*R3A)/XA2(42)
      YLIN(43,43) = (RPC + RAC + 4.0*RCC + RCO)/XA2(43)
      YLIN(45,43) = RCO/XA2(45)
C 44th eq N14 - duplicate
      YLIN(2 ,44) = RnN14p/XA2(2)
      YLIN(8 ,44) = - RPC13/XA2(8)
      YLIN(11,44) = - RPO17a/XA2(11)
      YLIN(41,44) = (RPN - RPC13 - RPO17a)/XA2(41)
      YLIN(44,44) = (RPN + RAN + RnN14p)/XA2(44)
C 45th eq O16 - duplicate
      YLIN(2, 45) = RnO16/XA2(2)
      YLIN(8 ,45) = - RAC13n/XA2(8)
      YLIN(10,45) = - RPN15/XA2(10)
      YLIN(13,45) = - RPF19A/XA2(13)
      YLIN(41,45) = (RPO - RPF19A - RPN15)/XA2(41)
      YLIN(42,45) = (RAO - RAC13n - RAC)/XA2(42)
      YLIN(43,45) = (RCO - RAC)/XA2(43)
      YLIN(45,45) = (RPO + RAO + RnO16 + RCO)/XA2(45)
C 46th eq Ne20 - duplicate
      YLIN(2 ,46) = (RnNe20 - RnF19)/XA2(2)
      YLIN(13,46) = (-RPF19 - RnF19)/XA2(13)
      YLIN(17,46) = - RPNA23a/XA2(17)
      YLIN(41,46) = (RPNE20 - RPF19- RPNA23a)/XA2(41)
      YLIN(42,46) = (- RAO + RANE)/XA2(42)
      YLIN(43,46) = 0d0 !- 2.0*RCC/XA2(43)
      YLIN(45,46) = - RAO/XA2(45)
      YLIN(46,46) = (RnNe20 + RANE + RPNE20)/XA2(46)
C Remove any NaN issues, by blanking rates if abundance = 0
      DO I=1, 50
         IF (XA2(I).EQ.0d0) THEN
            DO J=1, 50
               YLIN(I,J) = 0d0
            END DO
         END IF
      END DO
C Kill all reaction rates if doing common envelope evolution
C      IF (ICE.EQ.1) THEN
C         DO I = 1, 50
C            YDOT(I) = 0d0
C            DO J=1,50
C               YLIN(J,I) = 0d0
C            END DO
C         END DO
C      END IF
      DO I = 1, 50
         DO J = 1, 50
            YLIN(J, I) = YLIN(J, I)*BARYN(I)*DMK
         END DO
         Y(I) = WN(I)
         YLIN(I, I) = YLIN(I, I) + DMK/DT
         YDOT(I) = (YDOT(I)*BARYN(I)+DWN(I)/DT)*DMK
      END DO
C Calculate neutron out rates for use later
      DO I = 1, 50
         XA2(I) = WN(I)
         IF (XA2(I).LT.0d0) XA2(I) = 0d0
C Neutron kill
C         IF (XA2(2).LT.1d-10) XA2(2) = 0d0
      END DO
C Work out abundances for nucrat2 - RJS
      DO I = 1, 50
         WW2(I) = XA2(I)/BARYN(I)
      END DO
      CALL NUCRAT2(AT)
      RATETOT = RDn + RnNA22p + RnNA22a + Rn26Mp + Rn26Ma + Rn26G
     :     + Rn26Gp + Rn26Ga
     :     + RnFe56 + RnFe57 + RnFe58 + RnCo59 + RnNi58 + RnNi59
     :     + RnNi60 + Rnp + RnHe3 + RnLi7 + RnC12 + RnC13 + RnC14
     :     + RnN14 + RnN15 + RnO16 + RnO18 + RnF19 + RnNe20 + RnNe21
     :     + RnNe22 + RnNa23 + RnMg24 + RnMg25 + RnMg26 + RnAl27
     :     + RnSi28 + RnSi29 + RnSi30 + RnP31 + RnS32 + RnS33 + RnS34
     :     + RnS33a + RnN14p + RnNi59p + RnNi59a + RnO17a + RnS32a
     :     + RnS34s + RnFe59 + RnFe60 + RnNi61s
      IF (RATETOT.EQ.0d0) RATETOT = 1d0
      FRAC(1 ,K) = - RnNi61s
      FRAC(2 ,K) = 0d0
      FRAC(3 ,K) = - Rnp
      FRAC(4 ,K) = RnHe3
      FRAC(5 ,K) = RnLi7
      FRAC(6 ,K) = 0d0
      FRAC(7 ,K) = 0d0
      FRAC(8 ,K) = RnC13 - RnC12
      FRAC(9 ,K) = RnC14 - RnN14p - RnO17a - RnC13
      FRAC(10,K) = RnN15 -  RnC14 - RnN14
      FRAC(11,K) = - RnO16 + RnO17a
      FRAC(12,K) = RnO18
      FRAC(13,K) = - RnNA22a - RnO18 + RnF19 
      FRAC(14,K) = RnNe21 - RnNe20
      FRAC(15,K) = - RnNe21 + RnNe22 - RnNA22p
      FRAC(16,K) = RnNA22p + RnNA22a
      FRAC(17,K) = - Rn26Ga - Rn26Ma - RnNe22 + RnNa23
      FRAC(18,K) = - RnNa23 + RnMg24
      FRAC(19,K) = RnMg25 - RnMg24
      FRAC(20,K) = - Rn26Mp - Rn26Gp + RnMg26 - RnMg25
      FRAC(21,K) = Rn26Mp + Rn26Ma
      FRAC(22,K) = Rn26Gp + Rn26Ga + Rn26G
      FRAC(23,K) = - RnMg26 + RnAl27 - Rn26G
      FRAC(24,K) = - RnAl27 + RnSi28
      FRAC(25,K) = - RnSi28 + RnSi29 - RnS32a
      FRAC(26,K) = - RnSi29 + RnSi30 - RnS33a
      FRAC(27,K) = - RnSi30 + RnP31
      FRAC(28,K) = - RnP31 + RnS32 + RnS32a
      FRAC(29,K) = - RnS32 + RnS33 + RnS33a
      FRAC(30,K) = - RnS33 + RnS34 + RnS34s
      FRAC(31,K) = RnFe56 - RnNi59a
      FRAC(32,K) = RnFe57 - RnFe56
      FRAC(33,K) = RnFe58 - RnFe57
      FRAC(34,K) = RnFe59 - RnFe58
      FRAC(35,K) = - RnFe59 + RnFe60
      FRAC(36,K) = RnCo59 - RnNi59p
      FRAC(37,K) = RnNi58 
      FRAC(38,K) = RnNi59 - RnNi58 + RnNi59a + RnNi59p
      FRAC(39,K) = RnNi60 - RnNi59
      FRAC(40,K) = RnNi61s - RnNi60
      FRAC(41,K) = Rnp - RDn - Rn26Mp - Rn26Gp - RnNi59p - RnN14p
      FRAC(42,K) = - Rn26Ga - Rn26Ma - RnS32a - RnNa22a - RnS33a
     :                - RnNi59a! NEEDS RnHe4
      FRAC(43,K) = RnC12
      FRAC(44,K) = RnN14p
      FRAC(45,K) = RnO16
      FRAC(46,K) = - RnF19
      DO I = 1,50
         FRAC(I,K) = FRAC(I,K)/RATETOT
      END DO
      RETURN
      END
