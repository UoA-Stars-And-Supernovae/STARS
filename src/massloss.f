C Program to sort out the slightly complicated issue of binary star mass-loss
C ML1 and ML2 are BC1 for stars 1 and 2 - Assumes input is set up with star 1 data
C first
      SUBROUTINE MASSLOSS(I,ML1,ML2, BCHORB, BCHSPIN1, BCHSPIN2)
      IMPLICIT REAL*8(A-H, L, M, O-Z)
      INTEGER MAXMSH
      PARAMETER (MAXMSH = 2000)
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,DEL,DH0,NMESH,JIN,IW(200)
      COMMON /TRANS / HT(26,MAXMSH,2)
      COMMON /AUXIN / ICL, ION, JW, IOP, INUC, IBC, ICN, IML(2),ISGTH,
     :     IMO, IDIFF
      COMMON /SODDS / ALPHA, RML, CMG, CSI, CFE, CT(10), AGE, DT, M1,
     :  EC, BM, ANG, CM, MTA, MTB, TM(2), T0, M0, TC(2), OS, AC, RCD,
     :  RMG, RHL, XF, DR, AK1 ,RMT, AK2, ITH, IX, IY, IZ, IB, ISX(45),
     :  TRB
      COMMON /INF   / VIN(60)
      COMMON /DINF  / DVIN(60)
      COMMON /OP    / ZS, LEDD, MM, DG, GRADT, ETH, RRLF, EGR, RR, Q
      COMMON /MASLOS/ AIJ(6,5),baseN
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CC, CD, CG, CR(2), CEVB,
     &                CEN, CPL, CMEVMU, CSECYR, LSUN, MSUN, RSUN, TSUNYR,
     &                STEFBOLTZ
      COMMON /EVMODE/ IMODE
      COMMON /INERTI/ VI(2)
      COMMON /ANGMOM/ VROT1, VROT2, FMAC, FAM, IRAM, IRS1, IRS2
      COMMON /WINDS / WINDML(2), FAKEWIND(2), BE
      COMMON /TIDES / MENV(2), RENV(2)
      COMMON /ZAMS  / TKH(2)
      DIMENSION T(2), AM(2), M(2), MT(2), AR(2), R(2), DAR(2), L(2), XH(2),
     :     XHE(2), XC(2), XO(2), HSPIN(2), DHSPIN(2), ML(2),
     :     RAT(2), RLF(2), DAM(2), WINDACC(2), MOMINER(2),
     :     BCHSPIN(2), DMOMINER(2), HSPINDT(2), R2O(2), TF(2), HTF(2), HSTF(2),
     :     OSPIN(2), ACCLIMIT(2), OSC(2)
      CBRT(VX) = DEXP(DLOG(VX)/3.0D0)
      PS(VX) = 0.5D0*(VX+DABS(VX))
      RLOBE(VX) = 0.49D0*VX*VX/(0.6D0*VX*VX+DLOG(1.0D0+VX))
C Set important variables into convenient pairs and calculate other
C necessary quantities
      DO ISTAR = 1,IMODE!2
         T(ISTAR) = DEXP(VIN(2+15*(ISTAR - 1)))
         AM(ISTAR) = VIN(4+15*(ISTAR - 1))
         DAM(ISTAR) = DVIN(4+15*(ISTAR - 1))
         M(ISTAR) = DEXP(AM(ISTAR))
         MT(ISTAR) = (M(ISTAR) - DEXP(AM(ISTAR) - DAM(ISTAR)))/DT
         AR(ISTAR) = VIN(7+15*(ISTAR - 1))
         R(ISTAR) = DEXP(AR(ISTAR))
         DAR(ISTAR) = DVIN(7+15*(ISTAR - 1))
         L(ISTAR) = VIN(8+15*(ISTAR - 1))
         XH(ISTAR) = VIN(5+15*(ISTAR - 1))
         XHE(ISTAR) = VIN(9+15*(ISTAR - 1))
         XC(ISTAR) = VIN(10+15*(ISTAR - 1))
         XO(ISTAR) = VIN(3+15*(ISTAR - 1))
         HSPIN(ISTAR) = VIN(14+15*(ISTAR - 1))
         DHSPIN(ISTAR) = DVIN(14+15*(ISTAR - 1))
      END DO
      HORB = VIN(13)
      DHORB = DVIN(13)
      IF (IMODE.EQ.2) THEN
         BM = M(1) + M(2)
      ELSE
         M(2) = BM - M(1)
      END IF
      SEP = (M(1)+M(2))*(HORB/(M(1)*M(2)))**2.0
C Orbital angular velocity
      OORB = HORB*BM/(M(1)*M(2)*SEP**2.0)
C surface boundary conditions
      DO ISTAR = 1,IMODE!2
         IF (ISTAR.EQ.1) ISTAROTHER = 2
         IF (ISTAR.EQ.2) ISTAROTHER = 1
         RAT(ISTAR) = M(ISTAR)/M(ISTAROTHER)
         RLF(ISTAR) = AR(ISTAR) - DLOG(SEP*RLOBE(CBRT(RAT(ISTAR))))
C FUDGE TO AVOID PRE-MS RLOF
         IF (AGE.LT.1d3) THEN
            RLF(ISTAR) = -1d-1
         END IF
C Save RLF data for use in funcs2
         HT(24, 1, ISTAR) = RLF(ISTAR)
C Different surface mass bc for *1 or *2 of binary
         IF (IB.EQ.1) THEN
C Select mass loss - RJS 24/6/03
            IF (IML(ISTAR).EQ.0) BC1 = 0d0
            IF (IML(ISTAR).EQ.1) THEN
               BC1 = RML*L(ISTAR)*R(ISTAR)/M(ISTAR)
            END IF
            IF (IML(ISTAR).EQ.2) THEN
               BC1 = RML*L(ISTAR)*R(ISTAR)/M(ISTAR)
C Convert Reimers to Blocker mass loss - need to sort out how to do M_ZAMS = 1.5
               BC1 = 4.83d-9*(1.5d0)**(-2.1)*(L(ISTAR)/LSUN)**2.7*BC1
            END IF
            IF (IML(ISTAR).EQ.3) THEN
               PERIOD = -2.07 + 1.94*DLOG10(1.4577*R(ISTAR)) - 0.9*DLOG10(M(ISTAR)/2.0)
               PERIOD = 10**PERIOD
               BC1 = 10**(-11.4 + 0.0123*PERIOD)
               BC1 = MSUN/CSECYR * BC1
               WIND = BC1
C Superwind, limited to vexp = 15 kms^-1 for P=>500 days
C            VEXP = DMAX1(15d0,-13.5d0+5.6d-2*PERIOD)
               VEXP = DMAX1(3.0d0,DMIN1(15d0,-13.5d0+5.6d-2*PERIOD))
               SWIND = L(ISTAR)*1d26/(CC*1d-2*VEXP*1d3)
C            SWIND = 2d-12*CSECYR/(MSUN*1d30)*SWIND
               SWIND = SWIND/(MSUN*1d30)
               BC1 = DMIN1(BC1,SWIND)
            END IF
            IF (IML(ISTAR).EQ.4) THEN
C Different surface mass bc for *1 or *2 of binary
c Stuff taken from thesis of L.Dray 2003
               COHe=(XC(ISTAR)/3.0+XO(ISTAR)/4.0)/XHE(ISTAR)
               SURFXH=XH(ISTAR)
               SURFXHe=XHE(ISTAR)/4.0
c     if(XHE.lt.0.1d0) THEN
c         if(PME.lt.3d18) THEN
c            CP3 = CT(9)*3d18
c         endif
c      endif


c*****NL***************************************
cc     c     pre-WR: de Jager 1998
               zml1=(log10(T(ISTAR))-4.05d0)/0.75d0
               zml2=(log10(L(ISTAR)/LSUN)-4.6d0)/2.1d0
               BC1=0d0
               do n2=0,5
                  do i2=0,n2
                     j2=n2-i2
C RJS 16/12/05 - prevent array out of bounds
                     IF (J2+1.LE.5) THEN
                        BC1=BC1-AIJ(i2+1,j2+1)*dcos(i2*dacos(zml1))
     :                       *dcos(j2*dacos(zml2))
                     END IF
                  enddo
               enddo
               BC1=(SQRT(ZS*50d0))*(10d0**BC1)/CSECYR
               if(zml1.ge.1d0.or.zml2.ge.1d0) BC1=0d0 !-5d90
               if(zml1.le.-1d0.or.zml2.le.-1d0) BC1=0d0 !-5d90
               BC1 = 2.0*BC1*MSUN
C RJS 9/6/06 - Added some bits to make smooth transitions between WR mass loss
               BCPREWR = BC1
c     !!!!WOLF RAYET MASS-LOSS!!!!!
CWhen XH(surface)<0.4 and log T > 4.0 the star is in the WNL phase:
CUse a constant rate of: 8e-5 M(sun) yr^-1
               BCWNL = 8d-5*MSUN/CSECYR
               IF (SURFXH.LT.0.4.AND.log10(T(ISTAR)).GT.3.9) THEN
                  BC1 = 1d1*(BCWNL - BCPREWR)*(log10(T(ISTAR)) - 3.9) + BCPREWR
               END IF
               IF (SURFXH.LT.0.4.AND.log10(T(ISTAR)).GT.4.0) BC1 = 8d-5*MSUN/CSECYR
               BCWC = 1d-7*(M(ISTAR)/MSUN)**2.5*MSUN/CSECYR
               IF (SURFXH.LT.3d-3.AND.log10(T(ISTAR)).GT.4.0) THEN
                  BC1 = 5d2*(BCWNL - BCWC)*(SURFXH - 1d-3) + BCWC
               END IF
CWhen XH(surface)<1e-3 and log T > 4.0 the star is in the WNE,WC or WO
Cphase:
CFor WNE use: 1.0e-7 (M(WR)/M(sun))**2.5 M(sun) yr^-1
               IF (SURFXH.LT.1d-3.AND.log10(T(ISTAR)).GT.4.0) THEN
                  BC1 = 1d-7*(M(ISTAR)/MSUN)**2.5*MSUN/CSECYR
                  IF (COHe.GT.3d-2) BC1 = 0.6d-7*(M(ISTAR)/MSUN)**2.5*MSUN/CSECYR
               END IF
CFor WC/WO use: 0.6e-7 (M(WR)/M(sun))**2.5 M(sun) yr^-1
CStars become WC when (C+O)/He > 3e-2
            END IF
            IF (IML(ISTAR).EQ.5) THEN
C - so this is the same as IML 5 but changing so that it fixes a few bugs
               COHe=(XC(ISTAR)/3.0+XO(ISTAR)/4.0)/XHE(ISTAR)
               SURFXH=XH(ISTAR)
               SURFXHe=XHE(ISTAR)/4.0
cc     c     pre-WR: de Jager 1988
               zml1=(log10(T(ISTAR))-4.05d0)/0.75d0
               zml2=(log10(L(ISTAR)/LSUN)-4.6d0)/2.1d0
               if(zml2.ge.1d0) zml2=1d0 !if more luminous that bounds, set to max
               if(zml1.ge.1d0) zml1=1d0 !if more hot that bounds, set to max
               if(zml1.le.-1d0) zml1=-1d0 !if more cool that bounds, set to max
               BC1=0d0
               do n2=0,5
                  do i2=0,n2
                     j2=n2-i2
C prevent array out of bounds
                     IF (J2+1.LE.5) THEN
                        BC1=BC1-AIJ(i2+1,j2+1)*dcos(i2*dacos(zml1))
     :                       *dcos(j2*dacos(zml2))
                     END IF
                  enddo
               enddo
               BC1=(SQRT(ZS*50d0))*(10d0**BC1)*MSUN/CSECYR
               if(zml2.le.-1d0) BC1=0d0 !so when the star becomes too faint no mass loss
               BCPREWR1 = BC1

C Vink et al. 2001 rates for OB stars - taken from JJE's code
               RVIN=(ZS*50d0)**0.13d0
C RMVA is for hot side of stability jump  27500 - 50000 K
               RMVA=-6.697+2.194*LOG10(L(ISTAR)/(LSUN*1d5))-1.313*LOG10(M(ISTAR)/(MSUN*30d0))
     :  -1.226*LOG10(RVIN*2.6/2.0)+0.933*LOG10(T(ISTAR)/4d4)-10.92*(LOG10(T(ISTAR)/4d4))**2
     :  +0.85*LOG10(50d0*ZS)
C RMVB is for cool side of jump  12500 - 22500 K
               RMVB=-6.688+2.210*LOG10(L(ISTAR)/(LSUN*1d5))-1.339*LOG10(M(ISTAR)/(MSUN*30d0))
     :              -1.601*LOG10(RVIN*1.3/2.0)+1.07*LOG10(T(ISTAR)/2d4)+0.85*LOG10(50d0*ZS)
C Am I missing an MSUN out of all this? I think so...
               IF(T(ISTAR).LE.5d4.AND.T(ISTAR).GT.2.75d4) BC1 =(10d0**RMVA)*MSUN/CSECYR
               IF(T(ISTAR).LE.2.25d4.AND.T(ISTAR).GE.1.25d4) BC1 =(10d0**RMVB)*MSUN/CSECYR
               IF(T(ISTAR).GT.2.25d4.AND.T(ISTAR).LE.2.75d4) THEN
                  BC1=((T(ISTAR)-2.25d4)*(10d0**RMVA)+(2.75d4-T(ISTAR))*(10d0**RMVB))*MSUN/(CSECYR*5d3)
               ENDIF
               IF(T(ISTAR).LT.1.25d4.AND.T(ISTAR).GE.1d4) THEN
                  BC1=((T(ISTAR)-1d4)*(10d0**RMVB)*MSUN/CSECYR+(1.25d4-T(ISTAR))*BCPREWR1)/2.5d3
               ENDIF
               IF(T(ISTAR).GT.5d4.AND.T(ISTAR).LT.6d4) THEN
                  BC1=((6d4-T(ISTAR))*(10d0**RMVA)*MSUN/CSECYR+(T(ISTAR)-5d4)*(BCPREWR1))/1d4
               ENDIF

               BCPREWR = BC1
C Note this smooths out the bistability jump... RJS
CWhen XH(surface)<0.4 and log T > 4.0 the star is in the WNL phase:
               BCWNL = -13.6 + 1.63*log10(L(ISTAR)/LSUN)+2.22*log10(XHE(ISTAR))
C Now has metallicity scaling
               BCWNL = ((ZS*50d0)**0.5d0)*(10.0**BCWNL)*MSUN/CSECYR
               IF (SURFXH.LT.0.4.AND.log10(T(ISTAR)).GT.3.9) THEN
                  BC1 = 1d1*(BCWNL - BCPREWR)*(log10(T(ISTAR)) - 3.9) + BCPREWR
               END IF
               IF (SURFXH.LT.0.4.AND.log10(T(ISTAR)).GE.4.0) BC1 = BCWNL
CWhen XH(surface)<1e-3 and log T > 4.0 the star is in the WNE,WC or WO
Cphase:
               BCWC = -8.3 + 0.84*log10(L(ISTAR)/LSUN) + 2.04*log10(XHE(ISTAR)) +
     :              1.04*log10(1d0-XHE(ISTAR))
C Also now metallicity scaled
               BCWC = ((ZS/0.02)**0.5d0)*(10.0**BCWC)*MSUN/CSECYR
               IF (SURFXH.LT.1d-3.AND.log10(T(ISTAR)).GE.4.0) THEN
C     WC rate
                  IF (COHe.GT.2d-2) THEN
                     BC1 = 1d2*(BCWC-BCWNL)*(COHe - 2d-2) + BCWNL
                  END IF
                  IF (COHe.GT.3d-2) BC1 = BCWC
C               IF (COHe.GT.1d0) BC1 = 1.9d-5*MSUN/CSECYR
               END IF
            END IF
C     New Mass loss from Bestenlehner 2020 - WIP
            IF (IML(ISTAR).EQ.7) THEN
C     Bremsstrahlung
               RHO = M(ISTAR)/MSUN
               RADI = R(ISTAR)/RSUN
               VOL = (4d0 / 3d0 * CPI * RADI ** 3d0)
               RHO = RHO / VOL
               KAP = 3.68D22 * GFF * (1d0 - ZS) * (1d0 + XH(ISTAR))
               KAP = KAP * RHO * T(ISTAR) ** (-3.5d0)

C     KAP is the average electron opacity
               VTH = SQRT(4d0 * BOLTZM * T(ISTAR) / AME)

C     VTH is the thermal velocity
               CFAC = -2d0 / CG**1.5d0
               CFAC = CFAC * (3 * CL / (CPI * STEFBOLTZ)) ** 0.5d0
               CFAC = CFAC * CR(1)**2D0 * (-2.01824D0)

C     ok now compute the mass loss
               BC1 = CFAC * 4*CPI*CG / (KAP * VTH)
               BC1 = BC1 * (2d0*XH(ISTAR)+0.75d0*XHE(ISTAR)+0.5*ZS) ** 2d0
               BC1 = BC1 * 0.0032421531d0 * (L(ISTAR) / LEDD) ** (27d0/14d0)
               BC1 = BC1 / (1-L(ISTAR)/LEDD)**(1d0/9d0)
            END IF
c     New mass loss to get to target mass by JJE - 2/5/2021
            IF (IML(ISTAR).EQ.9) THEN
               BC1 = 0.1d0*(M(ISTAR)/MSUN-RML*2.5d12*CSECYR*RSUN*LSUN/(MSUN**2d0))
     :              *MSUN/CSECYR/TKH(ISTAR)
C               write(*,*) BC1,M(ISTAR)/MSUN,RML*2.5d12*CSECYR*RSUN*LSUN/(MSUN**2d0) !- debug line
!mass loss rate is from target mass - current mass limited to 10% of thermal timesvale
!note - the funny constants around RML are because the value is modified to code units in printa.f, the simplest way to adjust the code is therefore to undo it here in this line without modifying the rest of the code to not to this when IML=9
            ENDIF
C Save wind mass loss for orbital angular momentum calculation
C Enhance wind mass loss by 1/(1-Omega/Omega_crit)
            OSPIN(ISTAR) = HSPIN(ISTAR)/VI(ISTAR)*SQRT(CG)
            OCRIT = DSQRT(6.67d-11*(M(ISTAR)*1d30)/((1d9*R(ISTAR))**3.0))
            OSC(ISTAR) = OSPIN(ISTAR)/OCRIT
            BC1 = BC1/DMAX1(1d-2,(1-OSC(ISTAR)/0.8))
            WINDML(ISTAR) = BC1
         END IF
         ML(ISTAR) = BC1
      END DO
C Now decide on what the BC really is - needs to be outside the loop to know about
C both Roche Lobes
C Messed up Common Envelope stuff
      FAKEWIND(1) = 0d0
      FAKEWIND(2) = 0d0
      MASSLIMIT = 1d-2 !2.5d-4
C Set limit for mass accretion at M/kelvin-helmholtz timescale
      DO ISTAR = 1,IMODE
         ACCLIMIT(ISTAR) = M(ISTAR)/TKH(ISTAR)*DMAX1(0d0,(1-OSC(ISTAR)/0.8))
      END DO
C Fakewind deals with mass-loss in CE systems - will interfer with normal evolution
      DO ISTAR=1,IMODE
         IF (ISTAR.EQ.1) ISTAROTHER = 2
         IF (ISTAR.EQ.2) ISTAROTHER = 1
         ML(ISTAR) = MT(ISTAR) - RMG*M(ISTAR) + ML(ISTAR)
     :        + DMIN1(RMT*(PS(RLF(ISTAR)))**3,MASSLIMIT*MSUN/CSECYR)
         IF (IMODE.EQ.2) THEN
C Add (1-omega/omega_crit) to reduce accretion rate
            ML(ISTAR) = ML(ISTAR) - DMIN1(FMAC*DMIN1(RMT*(PS(RLF(ISTAROTHER)))**3.0
     :           ,MASSLIMIT*MSUN/CSECYR),ACCLIMIT(ISTAR)*MSUN/CSECYR)
         END IF
      END DO
      FAKEWIND(1) = DMAX1(0d0,DMIN1(RMT*(PS(RLF(1)))**3,MASSLIMIT*MSUN/CSECYR) -
C     :     DMIN1(FMAC*RMT*(PS(RLF(1)))**3.0,ACCLIMIT(1)*MSUN/CSECYR))
     :     DMIN1(FMAC*RMT*(PS(RLF(1)))**3.0,ACCLIMIT(2)*MSUN/CSECYR))
C Should ACCLIMIT be the other star? I think so...
C Fakewind is supposed to be total mass lost by star one minus mass accreted by 2
      FAKEWIND(2) = DMAX1(0d0,DMIN1(RMT*(PS(RLF(2)))**3,MASSLIMIT*MSUN/CSECYR) -
C     :     DMIN1(FMAC*RMT*(PS(RLF(2)))**3.0,ACCLIMIT(2)*MSUN/CSECYR))
     :     DMIN1(FMAC*RMT*(PS(RLF(2)))**3.0,ACCLIMIT(1)*MSUN/CSECYR))
C Ignore wind accretion if RLOF occurs
      WINDACC(1) = 0d0
      WINDACC(2) = 0d0
      IF ((PS(RLF(1)))**3.0.EQ.0d0.AND.(PS(RLF(2))).EQ.0d0) THEN
C This assumes you're only changing the mass by winds with no RLF
C Wind accretion from Hurley et al. 2002
C Accrete material from the greater mass loser
         IF (WINDML(1).GT.WINDML(2)) THEN
            IDONOR = 1
            IACC = 2
         ELSE
            IDONOR = 2
            IACC = 1
         END IF
         VORB2 = CG*BM/SEP
C BetaW depends on spectral type = 7 for O-type stars
         BETAW = 7.0
         VESC2 = 2.0*BETAW*CG*M(IDONOR)/R(IDONOR)
         V2 = VORB2/VESC2
C Have assumed a circular orbit - with alpha_w = 3/2
         MTACC = (CG*M(IACC)/VESC2)**2.0*(3.0/(4.0*SEP**2.0))/(1+V2)**(3.0/2.0)
         MTACC = MTACC*WINDML(IDONOR)
         MTACC = DMIN1(MTACC,0.8*WINDML(IDONOR))
         WINDACC(IACC) = MTACC
         ML(IACC) = ML(IACC) - MTACC
      END IF
C Store WINDML + WINDACC for use in determining when variable composition accretion
C is needed - only if sum of one of them is -ve. Only do this if there's no RLOF
      IF ((PS(RLF(1)))**3.0.EQ.0d0.AND.(PS(RLF(2))).EQ.0d0) THEN
         HT(23,1,1) = 0d0 !WINDML(1) - WINDACC(1)
         HT(23,1,2) = 0d0 !WINDML(2) - WINDACC(2)
      ELSE
         HT(23,1,1) = 0d0
         HT(23,1,2) = 0d0
      END IF
C Fudge CE stuff
C      IF(RLF(1).GT.0d0) THEN
C         ML(1) = 1d-2*MSUN/CSECYR + ML(1)
C         WINDML(1) = 1d-2*MSUN/CSECYR
C         IF (M(2)/MSUN.LT.17) THEN
C            ML(2) = -5d-3*MSUN/CSECYR + ML(2)
C         END IF
C      END IF
      ML1 = ML(1)
      ML2 = ML(2)
C Tidal friction stuff. From Hut (1981), using Zahn (1977) for frictional timescale
      DO ISTAR = 1,IMODE
C TF in years
C         TF(ISTAR) = 3.5*(M(ISTAR)/MSUN*(R(ISTAR)/RSUN)**2.0
C     :        /(L(ISTAR)/LSUN))**1.0/3.0*CSECYR
C         DkT = 0.1/(R(ISTAR)**3.0/(CG*M(ISTAR)*TF(ISTAR)))
         QQ = 1/RAT(ISTAR)
         FQ = QQ*(1+QQ)
         RADSEP = R(ISTAR)/SEP
C RJS 14/1/08
C Alternative tidal prescription from Hurley, Trout & Pols (2002)
C Convective envelope stars
C Convective turnover time
         ME = DMAX1(0d0,M(ISTAR) - MENV(ISTAR))
         RE = DMAX1(0d0,R(ISTAR) - RENV(ISTAR))
         TCONV = ME*RE*(R(ISTAR) - 0.5*RE)
         TCONV = 0.4311*(TCONV/(3d0*L(ISTAR)))**(1d0/3d0)*CSECYR
         TCONV = DMAX1(TCONV,1d0)
C Tidal pumping timescale
         PID = 1d0/DABS(OORB - OSPIN(ISTAR))/SQRT(CG)
C Numerical factor supressing tide
         FCONV = DMIN1(1d0,(0.5d0*PID/TCONV)**2d0) ! check units match, should be ok now...
C Form of k/T
         DkT = (2d0/21d0)*(FCONV/TCONV)*(ME/M(ISTAR))
C Dynamical tide with radiative damping, also from HTP 2002
C Don't understand their formula, but from eq 28 & 41
         DkTR = SQRT(CG*M(ISTAR)/R(ISTAR)**3d0)*(1d0+QQ)**(5d0/6d0)
C This is a fit to tabulated data for MS stars from Zahn '75 -- how valid is it
C for different metallicities or non-MS stars?
         E2 = 1.592d-9*(M(ISTAR)/MSUN)**2.84d0
         DkTR = DkTR*E2*RADSEP**(5d0/2d0)
C         write (*,*) DkT, DkTR
         DkT = DkT + DkTR
C HTF is rate of change of orbital angular momentum
C Neglecting eccentricity - for now...
         OSPIN(ISTAR) = HSPIN(ISTAR)/VI(ISTAR)
         HTF(ISTAR) = -3.0*DKT*FQ*RADSEP**8.0*HORB*(1 - (OSPIN(ISTAR)/OORB))
C Spin angular momentum - HSTF.
         HSTF(ISTAR) = 3.0*DKT*QQ**2.0*RADSEP**6.0*M(ISTAR)*R(ISTAR)**2.0
     :        *OORB*(1 - (OSPIN(ISTAR)/OORB))
      END DO
      IF (AGE.LT.1d3) THEN
         HTF(1) = 0d0
         HTF(2) = 0d0
         HSTF(1) = 0d0
         HSTF(2) = 0d0
      END IF
C      IF (RLF(1).GT.0d0) THEN
C      IF (R(1)/SEP.GT.0.45) THEN
C         FACAM = DMAX1(0d0, (0.55 - (R(1)/SEP))/0.1)
C         HTF(1) = 0d0
C         HTF(2) = 0d0
C         HSTF(1) = 0d0
C         HSTF(2) = 0d0
C      END IF
C Boundary condition for orbit
C Evolution of orbital Ang. Mom. from Hurley et al (2002)
      IF (IMODE.EQ.2) THEN
         BCHORB = DHORB/DT + (HORB/BM)*((WINDML(1)+FAKEWIND(1))*M(2)/M(1) - WINDACC(1)
C     :        + (WINDML(2)+FAKEWIND(2))*M(2)/M(1) - WINDACC(2)) - HTF(1) - HTF(2)
     :        + (WINDML(2)+FAKEWIND(2))*M(1)/M(2) - WINDACC(2)) - HTF(1) - HTF(2)
      ELSE
         BCHORB = DHORB/DT + HORB/BM*((WINDML(1)+FAKEWIND(1))*M(2)/M(1)) - HTF(1)
      END IF
C      IF (RLF(1).GT.0d0) THEN
C         BCHORB = DHORB/DT + 1.5*HORB*WINDML(1)/(M(1)+M(2))
C      END IF
C Boundary condition for spin period
C Assuming solid body rotation
      R2O(1) = HSPIN(1)/VI(1)*R(1)**2.0
      R2O(2) = HSPIN(2)/VI(2)*R(2)**2.0
      DO ISTAR = 1,IMODE
         IF (ISTAR.EQ.1) ISTAROTHER = 2
         IF (ISTAR.EQ.2) ISTAROTHER = 1
C Assume mass lost/gained as a shell - hence factor of 2/3. Assume mu_w = 1
C Loss from star
C This has to be made consistent with mass loss & accretion limits above
C         HSPINDT(ISTAR) = 2.0/3.0*(WINDML(ISTAR) + RMT*(PS(RLF(ISTAR)))**3)*R2O(ISTAR)
         HSPINDT(ISTAR) = 2.0/3.0*(WINDML(ISTAR)/DMAX1(1d-2,(1-OSC(ISTAR)/0.8))
     :        + DMIN1(RMT*(PS(RLF(ISTAR)))**3
     :        ,MASSLIMIT*MSUN/CSECYR))*R2O(ISTAR)
         IF (IMODE.EQ.2) THEN
            HSPINDT(ISTAR) =  HSPINDT(ISTAR) - 2.0/3.0*(WINDACC(ISTAROTHER) +
     :           FAM*DMIN1(FMAC*DMIN1(RMT*(PS(RLF(ISTAROTHER)))**3.0
     :           ,MASSLIMIT*MSUN/CSECYR),ACCLIMIT(ISTAR)*MSUN/CSECYR))*R2O(ISTAROTHER)
         END IF
         BCHSPIN(ISTAR) = DHSPIN(ISTAR)/DT + HSPINDT(ISTAR) - HSTF(ISTAR)
      END DO
      BCHSPIN1 = BCHSPIN(1)
      BCHSPIN2 = BCHSPIN(2)
      RETURN
      END