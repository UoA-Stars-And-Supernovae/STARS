**==FUNCS1.FOR
      SUBROUTINE FUNCS1(I, K, K1, K2, ISTAR)
      IMPLICIT REAL*8(A-H, L, M, O-Z)
      INTEGER MAXMSH
      PARAMETER (MAXMSH = 2000)
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,DEL,DH0,NMESH,JIN,IW(200)
      COMMON /TRANS / HT(28,MAXMSH,2)
      COMMON /AUXIN / ICL, ION, JW, IOP, INUC, IBC, ICN, IML(2), ISGTH,
     :     IMO, IDIFF
      COMMON /SODDS / ALPHA, RML, CMG, CSI, CFE, CT(10), AGE, DT, M1,
     :  EC, BM, ANG, CM, MTA, MTB, TM(2), T0, M0, TC(2), OS, AC, RCD,
     :  RMG, RHL, XF, DR, AK1 ,RMT, AK2, ITH, IX, IY, IZ, IB, ISX(45),
     :  TRB
      COMMON /ATDATA/ DH2(4), CCHI(26,9), OMG(27), AMASS(10), BN(10), IZZ(10)
* extra common for mesh-spacing
      COMMON /PMESH / PMH(2), PME(2), IAGB
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, CP,
     :                CHI, QP, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP,
     :                RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO,
     :                ROO, RGNE, RGMG, RCCG, RPNG, EX, ENX, EXX(3), AAWT,
     :                AWT, WR(14)
C NB - if new variables are added, will need to change the WI(n) numbers in the
C rest of this subroutine. RJS 5/7/06
      COMMON /INF   / AF, AT, VX16, AM, VX1, VQK, AR, L, VX4, VX12,
     :                VX20, VX14, HORB, HSPIN, VX3, WI(45)
      COMMON /DINF  / DAF, DAT, DX16, DAM, DX1, V1(3), DX4, DX12, DX20,
     :                DX14, DHORB, DHSPIN, DX3, DWI(45)
      COMMON /OUTF  / BC1, BC2, BC3, BCHORB, BCHSPIN, VP, VPK, VR, VRK,
     :                VT, VTK, VL, LK,
     :                LQ, GTMA, MT, VM, VMK, QK, SG, WT, X1, X1T, X16,
     :                X16T, X4, X4T, X12, X12T, X20, X20T, X14, X14T,
     :                X24, SGTH, MU, X3, X3T, DLEV, DMIX(6),
     :                D4, D12, D14, D16, D20, D3, WX(103)
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XHE3, XW(14)
      COMMON /OP    / ZS, LEDD, M, DG, GRADT, ETH, RLF, EGR, R, Q
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CC, CD, CG, CR(2), CEVB,
     &                CEN, CPL, CMEVMU, CSECYR, LSUN, MSUN, RSUN, TSUNYR
      COMMON /WT    / PWT,SIG,MK,MT2,LQP,LKP
      COMMON /MIX   / KBICZ,KICZ
      COMMON /MASLOS/ AIJ(6,5),baseN
      COMMON /ACCRET/ ACCOMPOS(7,31, 2)
      COMMON /EVMODE/ IMODE
      COMMON /INERTI/ VI(2)
      COMMON /MONTAG/ MIXV, MIXL
      COMMON /MIXFUD/ SGTHFAC, FACSGMIN, FACSG, ISGFAC
      COMMON /IONISE/ MSTORE(8,5), MSTORE2(5)
      COMMON /DIFFUS/ D(10), A12(10)
      DIMENSION XSPEC(6), DDMIX(10), MUXX(6)
      CBRT(VX) = DEXP(DLOG(VX)/3.0D0)
      PS(VX) = 0.5D0*(VX+DABS(VX))
      RLOBE(VX) = 0.49D0*VX*VX/(0.6D0*VX*VX+DLOG(1.0D0+VX))
      M = DEXP(AM)
c Check for age of star and also timestep and stop if either too large - JJE - 2/5/2021
      IF(AGE.GE.1d15) THEN
         write(*,*) "Age of star beyond 1000 TYrs - Stopping - ",AGE
         STOP
      ENDIF
      IF(DT/CSECYR.GE.1d12) THEN
         write(*,*) "DT in years is over 1 Tyrs - stopping",DT/CSECYR
         STOP
      ENDIF
C Spin period to equns - this isn't the best way of doing this...
      BCHSPIN = HSPIN
C Compute binary mass if in binary mode
      IF (IMODE.EQ.2) THEN
         BM = M + DEXP(WI(4))
      END IF
      VM3 = CBRT(M)
C Pierre's mass mod.
      MT = (M - DEXP(AM-DAM))/DT
C      MT = M*DAM/DT
      MT2 = MT
C set up the composition variables
      XH = VX1
      XHE = VX4
      XC = VX12
      XN = VX14
      XO = VX16
      XNE = VX20
      XMG = CMG*ZS
      XSI = CSI*ZS
      XFE = CFE*ZS
      XHE3 = VX3
* Put everything else in Mg24.
      XMG = 1.0D0 - XH - XHE - XC - XO - XNE - XN - XSI - XFE - XHE3
      IF (XMG .LT. 0.0D0) XMG = 0.0D0
* following if N14 is a variable rather than Ne20:
C      XNE = 1.0D0 - XH - XHE - XC - XN - XO - XMG - XSI - XFE
C      IF (XNE .LT. 0.0D0) XNE = 0.0D0
      CALL STATEL (I, AF, AT, ISTAR)
      R = DEXP(AR)
      R2 = R*R
      R2M = 2.0D0/(PI4*RHO*R)
C pressure gradient equation; gravity corrected for corotation
C ANGMOM is the orbital angular momentum
      GRAV = 1.0D11*CG*M/R2
      SILLYM = DMAX1(M,0.95D0*TM(ISTAR))
      IF (IB.EQ.1) GO TO 1
      FDT=AGE*3.1557D7+DT-T0
      M1=M0-PS(-FDT)*MTA+PS(FDT)*MTB
    1 SEP = BM*(ANG/(SILLYM*(BM-SILLYM)))**2
C Angular momentum in G=1 units, according to RPC
C Replacing old separation with new one based on eigenvalue HORB
      IF (ISTAR.EQ.2) HORB = WI(13)
      IF (IMODE.EQ.2) THEN
         SEP = (TM(1)+TM(2))*(HORB/(TM(1)*TM(2)))**2.0
      ELSE
         SEP = BM*(HORB/(TM(1)*(BM - TM(1))))**2.0
      END IF
*     GRAV = GRAV - GRAV*BM*R*R2/(1.5D0*M*SEP**3)
C gravity corrected for rotation
      OSPIN = (H(14+15*(ISTAR - 1),1) + DH(14+15*(ISTAR - 1),1))/VI(ISTAR)
C Need to mess about with units to get this right
      OSPIN = OSPIN*SQRT(CG) ! now in s^-1
      GRAV = GRAV - 1d11*2.0/3.0*OSPIN**2.0*R
      IF (1.0D11*CG*M/R2.LT.1d11*2.0/3.0*OSPIN**2.0*R.AND.I.EQ.-1) THEN
         WRITE (32+(ISTAR-1),*) "break-up velocity reached!"
         write (32+(ISTAR-1),*) ISTAR, 1.0D11*CG*M/R2, 1d11*2.0/3.0*OSPIN**2.0*R, VI(ISTAR)
         write (32+(ISTAR-1),*) OSPIN, R, H(14+15*(ISTAR - 1),1) + DH(14+15*(ISTAR - 1),1)
         STOP
      END IF
      APM = -1.0D11*GRAV/(PI4*R2*P)
C temperature gradient equation
      GRADR = 1.0D11*FK*P*L/(4.0D0*PI4*CC*GRAV*R2*PR)
      DG = GRADR - GRADA
      S2 = GRADA*CP*T
      HP = P/(GRAV*RHO)
C mixing-length theory
      W0 = 0.5D0*S2*(ALPHA*ALPHA*HP/(9.0D0*CHI))**2
      W2 = 546.75D0*W0*DMAX1(0.0D0,DG) + 73.0D0
      W3 = CBRT(W2+DSQRT(W2*W2+12167.0D0))
      W5 = DMAX1(W3-23.0D0/W3-2.0D0,0.0D0)
      WCV = W5*CHI/(3.0D0*ALPHA*HP)
      WL = ALPHA*HP*WCV
      MIXL = ALPHA*HP
      MIXV = WCV
      IF (DG.LT.0d0) MIXV = 0d0
C if GRADR < GRADA, we get WCV = 0 and GRADT = GRADR
      IF (KBICZ.GT.600.AND.K.GT.KBICZ.AND.IAGB.EQ.1) WCV = 0d0
      GRADT = GRADR-4.0D0*HP*WCV**3/(ALPHA*S2*CHI)
      GTMA = GRADT - GRADA
      ATM = GRADT*APM
C Old MSF 2000
      CT1 = 1.0D1**(1.0D1*CT(1)) !!! should change input format !!!
      VP = CT(4)*AP + CT(5)*LOG((P+CT(9))/(P+CT1)) +
     &     CT(2)*LOG((P+CT(10))/(P+CT1))
      VPP = CT(4) + P/(P+CT1)*(CT(5)*(CT1-CT(9))/(P+CT(9))
     &     + CT(2)*(CT1-CT(10))/(P+CT(10)))
      CT10 = 2.0D4              !!! fixed, should change input format !!!
      VT = CT(7)*DLOG(T/(T+CT10))
      VTT = CT(7)*CT10/(T+CT10)
C mesh spacing equation, with modified pressure, temperature gradient equations
      IF (IAGB.EQ.1) THEN
         CP1 = CT(1)*PME(ISTAR)
         CP2 = CT(10)*PME(ISTAR)
         CP3 = CT(9)*PMH(ISTAR)
         VP = CT(4)*AP + CT(5)*LOG((P+CP3)/(P+CP1)) +
     &        CT(2)*LOG((P+CP2)/(P+CP1))
         VPP = CT(4) + P/(P+CP1)*(CT(5)*(CP1-CP3)/(P+CP3)
     &        + CT(2)*(CP1-CP2)/(P+CP2))
         CT10 = 2.0D4           !!! fixed, should change input format !!!
         VT = CT(7)*DLOG(T/(T+CT10))
         VTT = CT(7)*CT10/(T+CT10)
      END IF
C VR and VM must go to zero at centre like r**2, m**(2/3)
      VM = CT(6)*CBRT(TM(ISTAR)*TM(ISTAR))
      VMM = VM + VM3*VM3
      VM = DLOG(VM/VMM)
      VMM = -0.66666667D0/(VMM*VM3)
      VR = -CT(3)*DLOG(R2/CT(8)+1.0D0)
      R21 = (3D0*M/(PI4*RHO))**(2D0/3D0)
      VR1 = -CT(3)*DLOG(R21/CT(8)+1.0D0)
      VRR = -CT(3)/(R2+CT(8))
C Q is the quantity that the meshpoints should be at equal intervals of
      Q = VP + VT + VM + VR
      QK = VQK
      QM = VPP*APM + VTT*ATM + VMM + VRR*R2M
      MK = VQK/QM
      VMK = VMM*MK
      VPK = VPP*APM*MK
      VTK = VTT*ATM*MK
      VRK = VRR*R2M*MK
      WT = 1D-7*R2*CHI*DT/(VRK*VRK)
      PWT = WT
C energy equation
      VL = L
      EN = EN + ENX
*** assume thermal eq. in central region if timestep small
*      ITHC = ITH
*      IF(K.GT.19*NMESH/20 .AND. DT .LT. 1D-1*CSECYR) ITHC = 0
*      LK = (EX + EN + EC - ITHC*T*(SF*DAF+ST*DAT)/DT)*MK
*      LQ = ITHC*CP*T*APM*MK
*** weight by a function of WT that drops from 1 to 0 for WT between 10^10
*** and 10^9
      IF (K.GT.39*NMESH/40) THEN
         AWT = 2D-16*R2*CHI*DT/(R2M*MK)**2
         AWT = 2D-15*R2*CHI*DT/(R2M*MK)**2
         AWT4 = AWT**4
         AAWT = AWT4/(1D0 + AWT4)
      ELSE
         AAWT = 1.0D0
      END IF
      LK = (EX + EN + EC - ITH*T*(SF*DAF+ST*DAT)/DT)*MK
      LKP = EX + EN + EC - ITH*T*(SF*DAF+ST*DAT)/DT
      SIG = SF*DAF
      LQ = ITH*CP*T*APM*MK
      LQP = CP*T*APM*GTMA
      L1 = M*(EX + EN + EC - ITH*T*(SF*DAF+ST*DAT)/DT)
C radius equation
      R2K = MK/(0.5D0*PI4*RHO*R)
C composition equations
* simplified to exclude reactions involving Mg24
C hydrogen equation with MS baryon correction when ICN = 1
      X1 = VX1
*      X1T = (2.0*IX*((1-ICN)*RPP + RPC + RPNG + (1-2*ICN)*(RPN + RPO))
      X1T = (IX*((1-ICN)*RPP*3 + RPC*2 + RPNG*2 -R33*2 +R34  +2* (1-2*ICN)*(RPN + RPO))
     :     + DX1/DT)*MK
*     x1t = (-5d-7*(0.76d0-x1)/csecyr + dx1/dt)*mk
C helium equation
      X4 = VX4
*      X4T = (4.0*(-IX*(0.5*RPP + RPN + RPO)*(1-ICN)
      X4T = (4.0*(-IX*(RPN + RPO + R33 +R34)*(1-ICN)
*     :     + IY*(3.0*R3A + RAC + 1.5*RAN + RAO + RANE)
*     :     - IZ*(RCC + RCO + 2.0*ROO + RGNE + RGMG)) + DX4/DT)*MK
     :     + IY*(3.0*R3A + RAC + 1.5*RAN + RAO)
     :     - IZ*(RCC + RGNE)) + DX4/DT)*MK
*     x4t = (5d-7*(0.76-x1)/csecyr + dx4/dt)*mk
C carbon equation
      X12 = VX12
      X12T = (12.0*(IX*(RPC - RPN) - IY*(R3A - RAC)
*     :            + IZ*(2.0*(RCC + RCCG) + RCO)) + DX12/DT)*MK
     :            + IZ*2.0*RCC) + DX12/DT)*MK
C nitrogen equation
      X14 = VX14
      X14T = (14.0*(IX*(RPN + RPNG - RPC - RPO) + IY*RAN) + DX14/DT)*MK
C oxygen equation
      X16 = VX16
      X16T = (16.0*(IX*(RPO - RPNG) - IY*(RAC - RAO)
*     :            + IZ*(RCO + 2.0*ROO - RGNE)) + DX16/DT)*MK
     :            - IZ*RGNE) + DX16/DT)*MK
C neon equation
      X20 = VX20
*      X20T = (20.0*(IY*(RANE - RAN - RAO) + IZ*(RGNE - RGMG - RCC))
      X20T = (20.0*(-IY*(RAN + RAO) + IZ*(RGNE - RCC))
     :            + DX20/DT)*MK
C He3 equation
      X3 = VX3
       X3T = (3.0*(-IX*RPP+IX*(2d0*R33+R34))+DX3/DT)*MK
C fudged convective diffusion coefficient: OS > 0 gives overshooting
      B = PR/PG
      AMAP=0.1D0/DABS(APM*M)
      EG = DG+OS/(2.5D0+B*(2.0D1+1.6D1*B))/(AMAP+1.0D0)
      UG = PS(EG)
      SG = UG*UG*TC(ISTAR)/(QM*QK)

      ! Always mix the outermost layers of the star...
      IF ( K.LT.15 ) SG = TC(ISTAR)/(QM*QK)
      IF(K.GT.NMESH-10) SG = TC(ISTAR) / (QM*QK)

      SIG = SG
C modified diffusion coefficient according to MLT
      SGMLT = (WL*1.0D-22*(PI4*RHO*R2)**2)/(3.0*MK)
      VG = UG/GRADR
      IF (IAGB.EQ.1) THEN
         SG = SGMLT
         IF (XH.GT.1.0D-6) THEN ! H-rich envelope...
            SG = SGMLT * AK1*VG/(1.0 - (1.0-AK1)*VG)
         ELSE                   ! interior...
            SG = SGMLT * AK2*VG/(1.0 - (1.0-AK2)*VG)
         END IF
      END IF
      IF (KBICZ.GT.600.AND.K.GT.KBICZ.AND.IAGB.EQ.1) SG = 0d0
C Mixing fudge
      SG = FACSG*SG
C Compute mean molecular weight, inc. for partial ionisation
C Should have IOP equal to 5
      XH1 = H(5+15*(ISTAR - 1),K)
      XHE1 =  H(9+15*(ISTAR - 1),K)
      XC1 =  H(10+15*(ISTAR - 1),K)
      XN1 =  H(12+15*(ISTAR - 1),K)
      XO1 = H(3+15*(ISTAR - 1),K)
      XNE1 =  H(11+15*(ISTAR - 1),K)
      XHE31 = H(15+15*(ISTAR - 1),K)
      LKP = 0d0
      XSPEC(1) = XH
      XSPEC(2) = XHE
      XSPEC(3) = XC
      XSPEC(4) = XN
      XSPEC(5) = XO
      XSPEC(6) = XHE3
      DO II = 1,ION
         MUX = 1d0
         DO JJ = 1, IZZ(II)
            MUX = MUX + JJ*MSTORE(JJ,II)/MSTORE2(II)
         END DO
C Store mu of each species
         MUXX(II) = 1/(MUX/AMASS(II))
         IF (II.EQ.2) MUXX(6) = 1/(MUX/AMASS(10))
         LKP = LKP + MUX*XSPEC(II)/AMASS(II)
C pretend He3 is like He4...
         IF (II.EQ.2) LKP = LKP + MUX*XSPEC(6)/AMASS(10)
      END DO
      LKP = 1/LKP
      MU = LKP
      MUREAL = MU
C diffusion requires the use of the actual mu, not one based on full
C ionisation. Note that it is the actual mu that is passed to printb
C Need to test what this does to the TH boundary output stuff...
      IF (IDIFF.EQ.1) THEN
C Paquette diffusion stuff
         CALL DIFFUSION(RHO,T)
         DO II = 1, 10
C can't diffuse through H if we don't have any...
            IF (XH.LT.0.1) D(II) = 0d0
            DDMIX(II) = D(II)*(1.0D-11*(PI4*RHO*R2))**2d0/MK
            A12(II) = A12(II)*D(II)*ATM*(1.0D-11*(PI4*RHO*R2))
            IF (II.LE.5) THEN
               D(II) = D(II)*GRAV/(CR(1)*T)*(MUXX(II) - MU)
            ELSE
               D(II) = D(II)*GRAV/(CR(1)*T)*(AMASS(II)/(1d0+IZZ(II)) - MU)
               IF (II.EQ.10) D(II) = D(II)*GRAV/(CR(1)*T)*(MUXX(6) - MU)
            END IF
            D(II) = D(II) - A12(II)
            D(II) = D(II)*(1.0D-11*(PI4*RHO*R2))
         END DO
         D4 = D(2)
         D12 = D(3)
         D14 = D(4)
         D16 = D(5)
         D20 = D(6)
         D3 = D(10)
         DMIX(1) = DDMIX(2)
         DMIX(2) = DDMIX(3)
         DMIX(3) = DDMIX(4)
         DMIX(4) = DDMIX(5)
         DMIX(5) = DDMIX(6)
         DMIX(6) = DDMIX(10)
      ELSE
         DO II = 1,6
            DMIX(II) = 0d0
         END DO
         D4 = 0d0
         D12 = 0d0
         D14 = 0d0
         D16 = 0d0
         D20 = 0d0
         D3 = 0d0
      END IF
      IF (ISGTH.EQ.1) THEN
         MU = 1.0/(2.0*XH1 + 0.75*XHE1 + 7.0/12.0*XC1 + 8.0/14.0*XN1 + 9.0/16.0*XO1
     :        + 11.0/20.0*XNE1 + XHE31) ! + 13.0/24.0*XMG + 15.0/28.0*XSI + 27.0/56.0*XFE)
C Full implicit mu -- is much less stable!
C         MU = 1.0/(2.0*XH + 0.75*XHE + 7.0/12.0*XC + 8.0/14.0*XN + 9.0/16.0*XO
C     :        + 11.0/20.0*XNE + XHE3) ! + 13.0/24.0*XMG + 15.0/28.0*XSI + 27.0/56.0*XFE)
         CHI1 = 4.0D0*CC*PR/(FK*RHO*RHO*CP*T)
         DTH = 3.0D0*CHI1/(MU*(-GTMA)*APM*MK)
         DTH = DMAX1(0d0,DTH)
         SGTH = SGTHFAC*DTH*1.0D-22*(PI4*R2*RHO)**2/MK
C     IF (XHE31.LT.1d-12.AND.K.LT.400) SGTH = 0d0
         IF (K.LE.4) THEN
            SGTH = 0d0
         END IF
      ELSE
         SGTH = 0d0
      END IF
      MT2 = SGTH
C the outermost 25 meshpoints are always mixed, even if not convective;
C      IF ( K.LT.5 ) SG = TC(ISTAR)/(QM*QK)
      IF ( K.LT.0.03*NMESH ) SG = 1d2*TC(ISTAR)/(QM*QK)
      IF ( K.GT.NMESH-0.01*NMESH) SG = TC(ISTAR)/(QM*QK)
      IF ( I.LE.0 ) THEN
C sundry numbers for PRINTB, FUNCS2
         ETH = -T*(SF*DAF+ST*DAT)/DT + CP*T*APM*GTMA*MT
         EGR = 1.0D22*CG*M/R - U
         LEDD = PI4*CC*CG*M/FK
         HT(1, K, ISTAR) = RHO
         HT(2, K, ISTAR) = SG
         HT(3, K, ISTAR) = ZT
         HT(4, K, ISTAR) = MK
         DO J = 1,14
            HT(4+J, K, ISTAR) = XW(J)
         END DO
C Store thermohaline stuff
         HT(19, K, ISTAR) = SGTH
         HT(20, K, ISTAR) = MU
C Store temp and ln f
         HT(21, K, ISTAR) = AF
         HT(22, K, ISTAR) = AT
C HT(23-24,K) are for mass loss, computed in massloss.f
C Store diffusion coeff related stuff
         HT(25, K, ISTAR) = GRAV/(CR(1)*T)
         HT(26, K, ISTAR) = 1.0D-11*(PI4*RHO*R2)
         HT(27, K, ISTAR) = MUREAL
         HT(28, K, ISTAR) = ATM
      END IF
C Calc RLF for printb - only if I=-1
      IF (K.EQ.1.AND.I.EQ.-1) THEN
C surface boundary conditions
         IF (IMODE.EQ.2) THEN
            BM = TM(1) + TM(2)
         END IF
         RAT = M/(BM-M)
         RLF = AR - DLOG(SEP*RLOBE(CBRT(RAT)))
         IF (AGE.LT.1d4) THEN
            RLF = -1d-1
         END IF
      END IF
      IF ( K.NE.K1 ) RETURN
      PRB = CA*TRB**4/3D0
      BC2 = DLOG(FK/GRAV*(1.5D0*PG + 0.75D0*PR)*(PR-2D0*PRB)/PR)
      BC3 = 1.0D11*L - (0.75D0*PI4*CC*R2*(PR-2D0*PRB))
      IF (I.GE.0) THEN
         IF (ISTAR.EQ.1) ISTAROTHER = 2
         IF (ISTAR.EQ.2) ISTAROTHER = 1
         ACCOMPOS(1,I+1,ISTAR) = 0.5*(H(5+15*(ISTAR-1),1) + H(5+15*(ISTAROTHER-1),1)) !X1
         ACCOMPOS(2,I+1,ISTAR) = 0.5*(H(9+15*(ISTAR-1),1) + H(9+15*(ISTAROTHER-1),1)) !X4
         ACCOMPOS(3,I+1,ISTAR) = 0.5*(H(10+15*(ISTAR-1),1) + H(10+15*(ISTAROTHER-1),1)) !X12
         ACCOMPOS(4,I+1,ISTAR) = 0.5*(H(12+15*(ISTAR-1),1) + H(12+15*(ISTAROTHER-1),1)) !X14
         ACCOMPOS(5,I+1,ISTAR) = 0.5*(H(3+15*(ISTAR-1),1) + H(3+15*(ISTAROTHER-1),1)) !X16
         ACCOMPOS(6,I+1,ISTAR) = 0.5*(H(11+15*(ISTAR-1),1) + H(11+15*(ISTAROTHER-1),1)) !X20
         ACCOMPOS(7,I+1,ISTAR) = 0.5*(H(15+15*(ISTAR-1),1) + H(15+15*(ISTAROTHER-1),1))
      END IF
      RETURN
      END