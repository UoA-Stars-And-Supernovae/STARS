C     Computes effects on the EOS of a free energy contribution
C     delta(F) = -k*T*N_e*G(XI,TI), where XI = n_e*m_u and TI = 1eV/(k*T).
      SUBROUTINE PRESSI(ICL, TI, XI, XF, XT, F, DC, DVT, DVF, DP, DPT,
     :                  DPF, DS, DST, DSF, DU)

      IMPLICIT NONE

      REAL*8 CDUM, EEG, CNSTS, XFFT, THX, XFTT, ABUND, SQRT
      REAL*8 THFT, STATFD, VX, DGDG, DPF, WF, BEGG, CA1
      REAL*8 F, DC, DD1, GAMY, ZI, DST, AVM, GAMXY
      REAL*8 EE, CF2, DD2, THF, DSF, XFT, NA, BB
      REAL*8 CEN, ABS, BEG, DGDXY, CEVB, DGDGG, CP1, CP4
      REAL*8 THXX, DP, DLOG, AX, THT, THYY, BXX, GAMXX
      REAL*8 EEX, BYY, DVF, GAMX, THC, DU, DS, XFFF
      REAL*8 BX, THFF, GI, NE, XT, WX, CA2, EXP
      REAL*8 AA, TI, LOG, BXY, TH, CX, FX, DGDYY
      REAL*8 FF, CP3, CA3, CF1, CBRT, GEG, BY, GAMYY
      REAL*8 XA, XFF, CXX, DVT, DGDXX, XF, DGDY, EQ
      REAL*8 CDUM2, CP2, BE, DGDX, CC, NZZ, YI, XTTT
      REAL*8 XTT, DD, DPT, GAM, AXX, FF1, THY, CPL
      REAL*8 GE, FF2, XI, THXY, THTT, FD, DEXP, NI
      INTEGER ICL

      COMMON /STATFD/ FD(9), XTT, XFT, XFF, XTTT, XFTT, XFFT, XFFF
      COMMON /ABUND / XA(10), NA(10), NE, NI, NZZ, AVM
      COMMON /CNSTS / CDUM(10), CEVB, CEN, CPL, CDUM2(6)

      CBRT(VX) = DEXP(DLOG(VX)/3.0D0)
      DATA CA1, CA2, CA3 /0.89752, 0.768, 0.208/
      DATA CF1, CF2 /1.07654, 0.56956/
      DATA CP1, CP2, CP3, CP4 /3.0, 0.25, 2.0, 3.0D-2/
* Pressure ionization. 
* FF and CC are approximations of the degeneracy parameters F and PSI
      YI = 13.6*TI
      ZI = CF1*XI*YI*SQRT(YI)
      DD = CF2*ZI
      DD1 = 1.0 + DD
      DD2 = DD1 + DD/3.0
      FF = ZI*CBRT(DD1)
      FX = DD2/DD1
      CC = 2*SQRT(FF)
      CX = CC*0.5*FX
      CXX = CX*(0.5*FX + DD/(3.0*DD2*DD1))
      EE = 1.0/(1.0 + XI/CP4)
      EEX = CP3*(XI/CP4)*EE
      BB = YI + CC - CP3*LOG(EE)
      BX = CX + EEX
      BY = YI + 1.5*CX
      BXX = CXX + EE*EEX
      BXY = 1.5*CXX
      BYY = YI + 1.5*BXY
      WX = (CP1/XI)**CP2
      AA = EXP(-WX)
      AX = CP2*WX*AA
      AXX = CP2*AX*(WX - 1.0)
      GI = AA*BB
      DGDX = AA*BX + AX*BB
      DGDY = AA*BY
      DGDXX = AA*BXX + 2.0*AX*BX + AXX*BB
      DGDXY = AA*BXY + AX*BY
      DGDYY = AA*BYY

      IF (ICL.EQ.1) THEN
* Coulomb interaction.
* the input variable F is used, rather than the above approximation, 
* resulting in quite a bit of programming to get the derivatives of 
* TH = d log(n_e)/d(psi). XFF, XFFF etc are 2nd and 3rd derivatives of n_e.
            FF1 = F + 1.0
            WF = SQRT(FF1)
            FF2 = -0.5*F/FF1
            TH = XF/WF
            THT = XFT/WF
            THF = XFF/WF + FF2*TH
            THX = THF/XF
            THY = THX*XT - THT
            THTT = XFTT/WF
            THFT = XFFT/WF + FF2*THT
            THFF = XFFF/WF + FF2*(2.0*THF - (2.0 + F)*FF2*TH/F)
            THXX = (THFF - THX*XFF)/XF**2
            THXY = THXX*XT - (THFT - THX*XFT)/XF
            THYY = THXX*XT**2 + 2*XT/XF*(THX*XFT - THFT) + THTT - THX*XTT
* GAM is the plasma interaction parameter. Note that THC = ZT**2*NI/NE
            THC = TH + ABS(NZZ/NE)
            GAM = CBRT(XI*(CPL*NE/NI)**2/3.0)*TI/CEVB*THC
            BB = (CA1*SQRT(3.0/GAM))**(1.0/CA2)
            EE = (GAM/(GAM + CA3))**(1.0/CA2)
            BE = EE + BB
            EEG = CA3/(GAM + CA3)
            BEG = (EEG*EE - 0.5*BB)/BE
            BEGG = (EEG*EEG*(1.0 - GAM*CA2/CA3)*EE + 0.25*BB)/BE
            GE = (NI/NE)*CA1*GAM/(BE**CA2)
            GEG = 1.0 - BEG
            DGDG = GE*GEG
            DGDGG = GE*(GEG*GEG + (BEG*BEG - BEGG)/CA2)
            GAMX = 1.0/3.0 + THX/THC
            GAMY = 1.0 + THY/THC
            GAMXX = THXX/THC - (THX/THC)**2
            GAMXY = THXY/THC - THX*THY/THC**2
            GAMYY = THYY/THC - (THY/THC)**2
            GI = GI + GE
            DGDX = DGDX + DGDG*GAMX
            DGDY = DGDY + DGDG*GAMY
            DGDXX = DGDXX + DGDGG*GAMX**2 + DGDG*GAMXX
            DGDXY = DGDXY + DGDGG*GAMX*GAMY + DGDG*GAMXY
            DGDYY = DGDYY + DGDGG*GAMY**2 + DGDG*GAMYY
      END IF
* evaluate changes to ionization potential (DC), pressure (DP), and
* entropy (DS), and their derivatives w.r.t. log(f) and log(T)
      DC = DGDX + GI
      DVT = DGDXX + DGDX
      DVF = DVT*XF
      DVT = DVT*XT
      DP = -XI*DGDX
      DPF = -XI*DVF
      DPT = XI*(DGDXY - DVT) + DP
      DVT = DVT - DGDXY - DGDY
      DS = GI - DGDY
      DST = DGDX - DGDXY
      DSF = DST*XF
      DST = DST*XT - DGDY + DGDYY
      DU = -DGDY

      RETURN
      END
