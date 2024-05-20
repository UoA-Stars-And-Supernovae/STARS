      SUBROUTINE STATEF(FL, TL)

      IMPLICIT NONE

      REAL*8 CG, D2TI, H2T, XH, DZH2TT, OPDAT, H2BT, D3
      REAL*8 PNLAST, COCOMPOS, DPBF, PI, THETA, EU, DL, PE
      REAL*8 CNSTS, ABUND, SQRT, U, T4, AM, STATFD, UF
      REAL*8 NCO, NEF, DSBF, DELTAN, VX, DABS, PSI, DUB
      REAL*8 ENB, SF, TT, WF, P, DE, D2, FKCN
      REAL*8 DUA, HG, F, GAMMA, DC, TH4, CP, AVM
      REAL*8 DST, NOH, PET, VA, DZH2T, BF, DH2, FR
      REAL*8 RET, TL, DV, PC, PF, HI, NXF, FDD
      REAL*8 DSF, P0, FT, NA, SI, REF, ZH2T, TST
      REAL*8 CNIU, TH2, NIO, DKCN, CEN, RT, ABS, POLAST
      REAL*8 STAT2, DSAF, CEVB, Q, DKN2, ZET, HA, H2F
      REAL*8 RU, ZH2S, NH2O, UI, TFAKE, CLN10, PH2, DKH2O
      REAL*8 OMG, DPAT, OPR, HGF, DPB, FKH2O, ZH2TT, CH
      REAL*8 NEO, CB, CW, DVF, FY, STAT1, DPBT, D1
      REAL*8 PCLAST, ZT, D1TI, SEF, H2, PC2, FC, PT
      REAL*8 NXT, DSAT, CT, DSBT, PH2O, HIT, RL, NE
      REAL*8 DSA, XT, HIF, XU, CR, FKOH, PI4, SCHA
      REAL*8 EXP, SET, GT, AA, AF, FKL, FK, TI
      REAL*8 SHA, CS, AT, RF, PL, W, CAMU, NET
      REAL*8 DSB, DKOH, RE, XHE, PN2, RHO, ZH2, PHI
      REAL*8 PO, PHG, FXP, QET, ST, DPA, BN, XA
      REAL*8 EPS, VM, SJHA, PR, FL, PEF, T, PG
      REAL*8 OPT, DVT, NN2, QE, FKH, S, XCOMPOS, HGT
      REAL*8 XF, PN, EQ, FX, DPAF, CD, OBASE, PSUBO
      REAL*8 DKCO, DKC2, UE, CSX, AD, AE, FKTOT, H2AT
      REAL*8 TKAPPA, AB, H2A, SE, SIT, CC, WR, DH2TI
      REAL*8 TH3, FKCO, NZZ, PSUBN, G, SIF, DELTAC, POH
      REAL*8 DB, ET, CA, NSTORE2, DELTAO, FKC2, EN, BT
      REAL*8 QEF, NCN, CPL, GE, PCO, ATDATA, D3TI
      REAL*8 DMAX1, CPI, NSTORE, NX, CHI, TF, PSUBC
      REAL*8 ENP, TCR, AUXIN, NC2, ZS, FO, CBASE, PCN
      REAL*8 GRADA, DEXP, B, TU, NI, AC
      INTEGER INT, JR, JW, IZ, IO, IML, JT, IONISE
      INTEGER IN, JCSX, IMO, IM, K2, ICN, IX, IOP
      INTEGER IDIFF, ISGTH, IT, LOG, INUC, K, MIN, LE
      INTEGER IBC, IC, MAX, JK, JX, J, ION, I
      INTEGER ICL, ILOOP, LOG10, LT, IR

      PARAMETER (EPS=1.0D-30)
      COMMON /AUXIN / ICL, ION, JW, IOP, INUC, IBC, ICN, IML(2), ISGTH,
     :     IMO, IDIFF
      COMMON /STAT1 / CSX(10), CS(90,127,10), CNIU(60,41,2), W(18400), 
     :                JCSX
      COMMON /STAT2 / PL, RL, U, P, RHO, FK, T, SF, ST, ZT, GRADA, CP, 
     :                CH, S, PR, PG, PF, PT, EN, WR(41)
      COMMON /STATFD/ RE, PE, QE, RET, PET, QET, REF, PEF, QEF, FDD(7)
      COMMON /ABUND / XA(10), NA(10), NEO, NIO, NZZ, AVM
      COMMON /ATDATA/ DH2, D1, D2, D3, CHI(26,9), OMG(27), AM(10), BN(10),
     &                IZ(10)
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CC, CD, CG, CR, CT, CEVB,
     &                CEN, CPL, CW(6)
      COMMON /IONISE/ NSTORE(8,5), NSTORE2(5)
      COMMON /OPDAT / cbase,obase,opT(141),opR(31),ZS
      FXP(VX) = DEXP(MAX(-50.0D0,MIN(50.0D0,VX)))
      DIMENSION HA(26), VA(26), tkappa(2,2),Xcompos(5),COcompos(8)
      data Xcompos /0.0d0,0.03d0,0.1d0,0.35d0,0.7d0/
      data COcompos /0.0d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,0.6d0,1.0d0/
      DATA JT, JX /2, 2/
* Evaluate Fermi-Dirac integrals according to Eggleton, Faulkner &
* Flannery (1973): RE, PE, SE and UE correspond to rho*, P*, S* and U*. 
* PSI is the usual degeneracy parameter
      F = EXP(FL)
      T = EXP(TL)
      UF = F/(1.0+F)
      WF = SQRT(1.0+F)
      PSI = FL + 2.0*(WF-LOG(1.0+WF))
      G = CT*T*WF
      CALL FDIRAC(F, G)
      PE = G*PE
      PET = PET + 1.0
      PEF = PEF + 0.5*UF
      QE = QE/(RE*WF)
      SE = QE + 2.0*WF - PSI
      SEF = QE*(QEF-REF-0.5*UF) - 1.0/WF
      SET = QE*(QET-RET)
      UE = SE + PSI - PE/(RE*CT*T)
* Some quantities that do not depend on the state of ionization:
* the NA are the element number densities (per baryon); AVM is the average
* mass per baryon (in amu); NEO and NIO are the numbers of electrons and
* ions assuming complete ionization.
      NEO = 0.0
      NIO = 0.0
      NZZ = 0.0
      AVM = 0.0
      DO I = 1, 10
         NA(I) = XA(I)/BN(I)
         AVM = AVM + AM(I)*NA(I)
         NIO = NIO + NA(I)
         NEO = NEO + IZ(I)*NA(I)
         NZZ = NZZ + IZ(I)*IZ(I)*NA(I)
      END DO
* PRESSI gives a crude model for pressure ionization and (if ICL=1) 
* a model for Coulomb interactions, returning corrections to the electron
* chemical potential, pressure, entropy and internal energy.
* TI is 1eV/kT, DE is 1 amu * number of electrons/cm3
      TI = CEVB/T
      DE = RE*CD
      CALL PRESSI(ICL, TI, DE, REF, RET, F, DC, DVT, DVF, DPA, DPAT,
     &     DPAF, DSA, DSAT, DSAF, DUA)
      DV = DC - PSI
      DVF = DVF - WF
* Contributions of the completely ionized species
      NE = 0.0
      NEF = 0.0
      NET = 0.0
      SI = 0.0
      SIF = 0.0
      SIT = 0.0
      UI = 0.0
      DO I = ION+1, 9
         NE = NE + IZ(I)*NA(I)
         VM = AM(I)*SQRT(AM(I))
         SI = SI - NA(I)*LOG(NA(I)/VM + EPS)
      END DO
* Calculate ionization of the first ION elements.
      DO I = ION, 1, -1
         SHA = 1.0
         SJHA = 0.0
         SCHA = 0.0
* compute potentials VA and number ratios HA of ionization state J
* relative to the ground state
         DO J = 1, IZ(I)
            VA(J) = -CHI(J,I)*TI + J*DV
            IF(J.EQ.1) THEN
               HA(J) = FXP(VA(J))*OMG(IZ(I))/OMG(IZ(I)+1)
            ELSE
               HA(J) = HA(J-1)*FXP(VA(J) - VA(J-1))*OMG(IZ(I)+1-J)
     &              /OMG(IZ(I)+2-J)
            END IF
            SHA = SHA + HA(J)
            SJHA = SJHA + J*HA(J)
            SCHA = SCHA + CHI(J,I)*TI*HA(J)
            NSTORE(J,I) = HA(J)
         END DO
         NSTORE2(I) = SHA
         VM = AM(I)*SQRT(AM(I))
         SI = SI + NA(I)*LOG(VM)
         IF(I.GT.1) THEN
* contributions to electron number density NE, entropy SI and
* internal energy UI for Helium and heavier
            VX = NA(I)/SHA
            SI = SI - VX*LOG(VX/OMG(IZ(I)+1) + EPS)
            DO J = 1, IZ(I)
               NX = HA(J)*VX
               NXF = NX*DVF*(J - SJHA/SHA)
               NXT = NXF*DVT/DVF + NX*(CHI(J,I)*TI - SCHA/SHA)
               NE = NE + J*NX
               NEF = NEF + J*NXF
               NET = NET + J*NXT
               SIF = SIF - VA(J)*NXF
               SIT = SIT - VA(J)*NXT
               SI = SI - NX*LOG(NX/OMG(IZ(I)+1-J) + EPS)
               UI = UI + CHI(J,I)*NX
            END DO
         END IF
      END DO
* Ionization and molecular dissciation of Hydrogen.
* partition function for H2 from Vardya (1960), Webbink (1975)
      DH2TI = DH2*TI
      D1TI = D1*TI
      D2TI = (D2*TI)**2
      D3TI = (D3*TI)**3
      ZET = 1.0 - (1.0 + DH2TI)*EXP(-DH2TI)
      DZH2T = -DH2TI**2*EXP(-DH2TI)/ZET
      DZH2TT = (DH2TI - 2.0 - DZH2T)*DZH2T
      ZH2 = 6608.8*ZET*DH2TI**(-2.5)*EXP(-D1TI - D2TI - D3TI)
      ZH2T = 2.5 + D1TI + 2.0*D2TI + 3.0*D3TI + DZH2T
      ZH2TT = -D1TI - 4.0*D2TI - 9.0*D3TI + DZH2TT
      ZH2S = ZH2*SQRT(8D0)/VM
      H2A = CEN*(ZH2S/4.0)*DE/(T*SQRT(T))*EXP(DH2TI)
      H2BT = DH2TI + 1.5 - ZH2T
      H2AT = RET - H2BT
* solve for densities of H+, H, and H2
      AA = 2*H2A + HA(1)*(1.0 + HA(1))
      AB = NE + HA(1)*(NE - NA(1))
      AC = NA(1)*NE
      HG = 2*AC/(SQRT(AB*AB + 4*AA*AC) + AB)
      HI = HA(1)*HG
      NE = NE + HI
      EN = 1.0/NE
      H2 = H2A*HG*HG*EN
      NI = NIO - H2
* derivatives w.r.t. F and T
      AA = NE + 4*HG*H2A
      AB = HA(1)*(NE - 2*H2)
      AD = 1/(AA + AB)
      AC = 2*H2*AD
      AF = (NEF - NE*REF)*AC
      AT = (NET - NE*H2AT)*AC
      AE = HG*AD
      BF = DVF*AE
      BT = (CHI(1,1)*TI + DVT)*AE
      HGF = AF - AB*BF
      HGT = AT - AB*BT
      HIF = HA(1)*(AF + AA*BF)
      HIT = HA(1)*(AT + AA*BT)
      NEF = NEF + HIF
      NET = NET + HIT
      H2F = H2*REF + EN*(2*H2A*HG*HGF - H2*NEF)
      H2T = H2*H2AT + EN*(2*H2A*HG*HGT - H2*NET)
* hydrogen contribution to entropy, internal energy
      SIF = SIF - VA(1)*HIF - H2BT*H2F
      SIT = SIT - VA(1)*HIT - H2BT*H2T + H2*(ZH2T + ZH2TT)
      SI = SI - HI*LOG(HI/OMG(1) + EPS) - HG*LOG(HG/OMG(2) + EPS)
     &        - H2*(LOG(H2/ZH2S + EPS) + ZH2T)
      UI = UI + CHI(1,1)*HI + 0.5*DH2*(HI + HG)
* DB is 1 amu * number of baryons/cm3; RHO is the mass density in g/cm3
      DB = EN*DE
      DL = LOG(DB)
      RHO = DB*AVM
      RL = LOG(RHO)
      RT = RET - EN*NET
      RF = REF - EN*NEF
* second call to PRESSI compensates for spurious pressure and entropy terms
      DE = DB*NEO
      CALL PRESSI(0, TI, DE, RF, RT, F, DC, DVT, DVF, DPB, DPBT, DPBF, 
     :            DSB, DSBT, DSBF, DUB)
* pressure terms
      PE = CB*PE
      TCR = T*CR
      P0 = TCR*DB
      PI = NI*P0
      T4 = T*T*T*T
      PR = CA*T4/3.0
      B = 4.0*PR/P0
      PG = PE + PI + TCR*(DPA-DPB)
      P = PG + PR
      PF = (PE*PEF      + PI*RF - H2F*P0 + TCR*(DPAF-DPBF))/P
      PT = (PE*PET + PI + PI*RT - H2T*P0 + TCR*(DPAT-DPBT) + PR*4.0)/P
      PL = LOG(P)
* entropy, in erg/g/K
      DSF = NEF*DSA + NE*DSAF - NEO*DSBF - RF*B
      DST = NET*DSA + NE*DSAT - NEO*DSBT - (RT - 3.0)*B
      SF = CR*(-NI*RF       + NEF*SE + NE*SEF + SIF + DSF)/AVM
      ST = CR*( NI*(1.5-RT) + NET*SE + NE*SET + SIT + DST)/AVM
      S = CR*(SE*NE + DSA*NE - DSB*NEO + B + (1.5*TL - DL + 2.5
     :        - LOG(CEN))*NI + SI)/AVM
* internal energy, in erg/g
      U = TCR*(UE*NE + DUA*NE - DUB*NEO + 1.5*NI + ZH2T*H2 + 0.75*B
     :         + TI*UI)/AVM
* other thermodynamic quantities
      Q = PT*SF - PF*ST
      CP = -Q/PF
      GRADA = SF/Q
      GAMMA = Q/(RT*SF-RF*ST)
      ZT = SQRT(ABS((NE*REF/WF + NZZ)/NI))
* TST ought to be zero, if all the above programming is correct
      TST = SF/CR - P*(RT*PF-RF*PT)/(TCR*RHO)
*** END OF THERMODYNAMIC CALCULATION. BEGINNING OF TABLE-BASED CALCULATIONS
      FR = RL/CLN10
      TF = TL/CLN10
* Opacity tables from Alexander & Ferguson (1994; molecular), Itoh (1983;
* electron conduction) and Iglesias & Rogers (1992; the rest)
C      XF = 2.0*X + Y + 1.0
C RJS -- Attempt at estimating molecular opacities
C RJS 14/7/07 - Compute pressures for HG, HI and H2
      IF (T.LT.1d4.AND.IMO.EQ.1) THEN
C Low T's cause trouble...
         TFAKE = DMAX1(1d3,T)
         PHG = RHO*HG*(CR*TFAKE)
         PHI = RHO*HI*(CR*TFAKE)       ! Note that P(ressure)HI spells PHI!!!!
         PH2 = RHO*H2*(CR*TFAKE)
         PC = RHO*NA(3)*(CR*TFAKE)
         PN = RHO*NA(4)*(CR*TFAKE)
         PO = RHO*NA(5)*(CR*TFAKE)
C Compute dissociation coefficients - from Rossi & Maciel '83
         THETA = 5040d0/TFAKE
         TH2 = THETA*THETA
         TH3 = TH2*THETA
         TH4 = TH3*THETA
         DKCN = 1d1**(1.2699d1-8.0228*THETA+1.8920d-2*TH2-6.1584e-4*TH3
     :        -1.1931e-6*TH4)
         DKCO = 1d1**(1.3590d1-1.1523d1*THETA+6.6519e-2*TH2-6.3343d-3*TH3
     :        +2.3781e-4*TH4)
         DKC2 = 1d1**(1.2654d1-6.4287*THETA+3.4061d-2*TH2-3.0069d-3*TH3
     :        +2.0995d-4*TH4)
         DKN2 = 1d1**(1.3236d1-1.0177d1*THETA+6.3664d-2*TH2-5.9114d-3*TH3
     :        +2.1845d-4*TH4)
         DKOH = 1d1**(1.2198d1-4.8641*THETA+6.6674d-2*TH2-5.6697d-3*TH3
     :        +1.9546d-4*TH4)
         DKH2O = 1d1**(2.5315d1-1.0343d1*THETA+1.0807d-1*TH2-8.7408d-3*TH3
     :        +2.9465d-4*TH4)
C Iterate computation of these things
         PsubN = PN
         PsubC = PC
         PsubO = PO
         IF (NA(3).GT.NA(5)) THEN
            PsubO = 0d0
            PsubC = PC
         ELSE
            PsubC = 0d0
            PsubO = PO
         END IF
         PClast = 0d0
         PNlast = 0d0
         POlast = 0d0
         ILOOP = 0
 10      continue
         ILOOP = ILOOP + 1
         PsubO = PO/(1+PsubC/DKCO+PHG/DKOH+PHG*PHG/(DKOH*DKH2O))
         PsubC = PC/(1+PsubO/DKCO+PsubC/DKC2+PsubN/DKCN)
         PsubN = -(1+PsubC/DKCN)
         PsubN = PsubN + SQRT(PsubN*PsubN+4*PN/DKN2)
         PsubN = DKN2/2d0*PsubN
         deltaC = dabs((PsubC-PClast)/PsubC)
         deltaN = dabs((PsubN-PNlast)/PsubN)
         deltaO = dabs((PsubO-POlast)/PsubO)
         PClast = PsubC
         PNlast = PsubN
         POlast = PsubO
         IF (deltaC.lt.1d-3.and.deltaN.lt.1d-3.and.deltaO.lt.1d-3) THEN
            CONTINUE
         ELSE
            IF (ILOOP.GT.5000) THEN
C               write (*,*) "Convergence failed", T, TFAKE, deltaC, deltaN, deltaO
            ELSE
               GOTO 10
            END IF
         END IF
C Compute molecular partial pressures
         PCO = PsubC*PsubO/DKCO
         PC2 = PsubC*PsubC/DKC2
         PN2 = PsubN*PsubN/DKN2
         PCN = PsubC*PsubN/DKCN
         POH = PsubO*PHG/DKOH
         PH2O = PHG*PHG*PsubO/(DKOH*DKH2O)
C Compute number abundances from pressures - note these are number fractions
         NCO = PCO/(CR*TFAKE*RHO)
         NC2 = PC2/(CR*TFAKE*RHO)
         NN2 = PN2/(CR*TFAKE*RHO)
         NCN = PCN/(CR*TFAKE*RHO)
         NOH = POH/(CR*TFAKE*RHO)
         NH2O = PH2O/(CR*TFAKE*RHO)
C Fits to opacity as a function of temp:
C CN per molecule - Scalo & Ulrich 1975
         CAMU = 1.6605402D-24
         FKCN = 1d1**(-19.12+2.22479*THETA-2.8069*TH2+0.76000*TH3-0.078384*TH4)
C Assume C2 as opaque as CN
         FKC2 = FKCN*NC2/(CAMU)
         FKCN = FKCN*NCN/(CAMU)
C CO, OH and H2O from Keeley 1970
         FKCO = 2.75d-26*NCO/CAMU
         FKOH = 17*(NOH/CAMU)*(1.4d-21*(T/1d4)**6d0)/(0.1+(T/1d4)**6d0)
C H2O, second term modified as in Marigo 2002
         FKH2O = (2.6d-27/(4.23d-4 + (T/1d4)**4d0)) + 
     :        (9.72d-18*DEXP(-3.2552/(T/1d4))/(1+3.78d3*(T/1d4)**1d1))
         FKH2O = 18*NH2O*FKH2O/CAMU
C Note that I have been unable to reproduce Marigo's plots with this, CN seems
C about right though and that's the important one
         FKTOT = FKCO+FKC2+FKCN+FKOH+FKH2O
      END IF
      if(IOP.LE.1) THEN
         XH = NA(1)*AM(1)/AVM
         XHE = NA(2)*AM(2)/AVM
         XF = 2.0*XH + XHE + 1.0
         IF ( JX.LT.2 ) JX = 2
         IF ( JX.GT.JCSX ) JX = JCSX
         JX = JX - 2
 1       CONTINUE
         JX = JX + 1
         IF ( XF.LT.CSX(JX+1) .AND. JX.LT.JCSX-1 ) GO TO 1
         IF ( XF.GE.CSX(JX) .AND. JX.GT.1 ) THEN
            JX = JX - 2
            GO TO 1
         ELSE
            XT = (XF-CSX(JX))/(CSX(JX+1)-CSX(JX))
            XU = 1.0 - XT
            IF (IOP.EQ.0) THEN
*     linear interpolation in opacity table
               TT = 20.0*(TF - 2.95)
               JT = MAX(1, MIN(126, INT(TT)))
               TT = TT - JT
               TU = 1.0 - TT
               RT = 4.0*(FR + 12.25)
               JR = MAX(1, MIN(89, INT(RT)))
               RT = MAX(0.0, RT - JR)
               RU = 1.0 - RT
               DO J = 1, 2
                  K = JX + J-1
                  IF (J .EQ. 1) FKL = TT*(RT*CS(JR+1,JT+1,K)
     &                 + RU*CS(JR,JT+1,K)) + TU*(RT*CS(JR+1,JT,K) 
     &                 + RU*CS(JR,JT,K))
                  IF (J .EQ. 2) FKH = TT*(RT*CS(JR+1,JT+1,K)
     &                 + RU*CS(JR,JT+1,K)) + TU*(RT*CS(JR+1,JT,K) 
     &                 + RU*CS(JR,JT,K))
               END DO
            ELSE
*     bicubic spline interpolation
               CALL OPACTY(JX,TF,FR,FKL,FKH)
            END IF
            FK = XT*10D0**FKH + XU*10D0**FKL
            CH = 4D0*CC*PR/(FK*RHO*RHO*CP*T)
*     Neutrino loss rates from Itoh et al (1983-1992)
            EN = 0.0
            IF ( TF .GE. 7.0 .AND. FR .GE. 0.0 ) THEN
               TT = 20.0*(TF - 6.95)
               IT = MAX(1, MIN(59, INT(TT)))
               TT = TT - IT
               TU = 1.0 - TT
               RT = 4.0*(FR + 0.25)
               IR = MAX(1, MIN(40, INT(RT)))
               RT = RT - IR
               RU = 1.0 - RT
               ENP = TT*(RT*CNIU(IT+1,IR+1,1) + RU*CNIU(IT+1,IR,1))
     :              + TU*(RT*CNIU(IT,IR+1,1) + RU*CNIU(IT,IR,1))
               ENB = TT*(RT*CNIU(IT+1,IR+1,2) + RU*CNIU(IT+1,IR,2))
     :              + TU*(RT*CNIU(IT,IR+1,2) + RU*CNIU(IT,IR,2))
               EN = -10.0D0**ENP - NZZ*10.0D0**ENB
            END IF
         END IF
      ELSE
!! DO MY NEW OPACITY STUFF - JJ 5/11/02
         fT=LOG10(EXP(TL))
C RJS - note fR is now log rho/(T6^3) not log rho as it was above
         fR=LOG10(RHO)-3d0*(fT-6d0)
         fX = NA(1)*AM(1)/AVM
         fY = NA(2)*AM(2)/AVM
         if(IOP.eq.2) then
            fC=0d0
            fO=0d0
         else
            fC = (NA(3)*AM(3)/AVM-cbase)/(1d0-ZS-fX)
            fO = (NA(5)*AM(5)/AVM-obase)/(1d0-ZS-fX)
         endif
         
         if(fX.lt.1d-10.and.IOP.ge.5) then
            fO=1d0-fX-fY-ZS-fC
         endif
!     add in N,Ne,Mg,Si
         if(fC.le.1d-6.and.fO.le.1d-6) then
                                !Pure X only Spline Opacity
            FKL=0d0
            FKH=0d0
            if(fX.ge.0.1) then
               if(fX.ge.0.35) then
                  call Xopacty(1+3*61,fT,fR,FKL,FKH)
                  FK=((0.7-fX)*10d0**FKL+(fX-0.35)*10d0**FKH)/0.35
               else
                  call Xopacty(1+2*61,fT,fR,FKL,FKH)
                  FK=((0.35-fX)*10d0**FKL+(fX-0.1)*10d0**FKH)/0.25
               endif
            else
               if(fX.ge.0.03) then
                  call Xopacty(1+61,fT,fR,FKL,FKH)
                  FK=((0.1-fX)*10d0**FKL+(fX-0.03)*10d0**FKH)/0.07
               else
                  call Xopacty(1,fT,fR,FKL,FKH)
                  FK=((0.03-fX)*10d0**FKL+(fX)*10d0**FKH)/0.03       
               endif
            endif
         else
            if(fC.lt.1d-20) fC=0d0
            if(fO.lt.1d-20) fO=0d0
                                !Full 5D interpolate
            if(fX.ge.0.1) then
               if(fX.ge.0.35) then
                  iX=3
               else
                  iX=2
               endif
            else
               if(fX.ge.0.03) then
                  iX=1
               else
                  iX=0
               endif
            endif
            do K2=1,8
               if(COcompos(K2).le.fC) iC=K2
               if(COcompos(K2).le.fO) iO=K2
            enddo
            if(iC.le.7) then
               JX=iO+(iC-1)*8+iX*61
            else
               JX=6*8+7+iX*61+iO
            endif
            do im=0,1
               do in=0,1
                  if(iO.le.6) then
                     JK=JX+im*8+in
                  else
                     JK=JX+im*7+in
                  endif
                  call Xopacty(JK,fT,fR,FKL,FKH)
                  FK=((Xcompos(iX+2)-fX)*10d0**FKL+(fX-Xcompos(iX+1))
     :                 *10d0**FKH)/(Xcompos(iX+2)-Xcompos(iX+1))
                  tkappa(im+1,in+1)=FK               
               enddo
            enddo
C RJS 19/09/05
C Fudge to keep things inside the array - I have no idea what effect this has...
            IF (iC.GT.7) iC = 7
            IF (iO.GT.7) iO = 7
            et=(fC-COcompos(iC))/(COcompos(iC+1)-COcompos(iC))
            eu=(fO-COcompos(iO))/(COcompos(iO+1)-COcompos(iO))
            FK=(1d0-et)*(1d0-eu)*tkappa(1,1)+et*(1d0-eu)*tkappa(2,1)
     :           +et*eu*tkappa(2,2)+(1d0-et)*eu*tkappa(1,2)
            
         endif
         CH = 4D0*CC*PR/(FK*RHO*RHO*CP*T)
C     Restore FR to log rho
         FR = RL/CLN10
         TF = TL/CLN10
* Neutrino loss rates from Itoh et al (1983-1992)
         EN = 0.0
         IF ( TF .GE. 7.0 .AND. FR .GE. 0.0 ) THEN
            TT = 20.0*(TF - 6.95)
            IT = MAX(1, MIN(59, INT(TT)))
            TT = TT - IT
            TU = 1.0 - TT
            RT = 4.0*(FR + 0.25)
            IR = MAX(1, MIN(40, INT(RT)))
            RT = RT - IR
            RU = 1.0 - RT
            ENP = TT*(RT*CNIU(IT+1,IR+1,1) + RU*CNIU(IT+1,IR,1))
     :           + TU*(RT*CNIU(IT,IR+1,1) + RU*CNIU(IT,IR,1))
            ENB = TT*(RT*CNIU(IT+1,IR+1,2) + RU*CNIU(IT+1,IR,2))
     :           + TU*(RT*CNIU(IT,IR+1,2) + RU*CNIU(IT,IR,2))
            EN = -10.0D0**ENP - NZZ*10.0D0**ENB
         END IF
      END IF
C Add molecular opacities to tabulated ones
      IF (T.LT.1d4.AND.IMO.EQ.1) THEN
         FK = FKTOT + FK
      END IF
      RETURN
      END
