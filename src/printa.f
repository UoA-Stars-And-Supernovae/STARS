**==PRINTA.FOR
      SUBROUTINE PRINTA(IEND, NP, IT1, IT2, NT5)
      IMPLICIT REAL*8(A-H, L, M, O-Z)
      real*8 MAT(4,141),Xcompos(3,305),COcompos(8)
      SAVE
      INTEGER MAXMSH
      PARAMETER (MAXMSH = 2000)
      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      COMMON /PREVM / HPR(60,MAXMSH),DHPR(60,MAXMSH),PR(14),MS(9999),
     &                ST(9999),NPR,KPR, HPPR(60,MAXMSH),
     &                HNUCPR(100,MAXMSH), DHNUCPR(100, MAXMSH)
      COMMON H(60,MAXMSH), DH(60,MAXMSH), EP(3), NH, JIN, ID(200)
      COMMON /AUXIN / ICL, ION, IAM, IOP, INUC, IBC, ICN, IML(2), ISGTH,
     :     IMO, IDIFF
      COMMON /SODDS / ALPHA, RML, CMG, CSI, CFE, CT(10), AGE, DT, M1, 
     :  EC, BM, ANG, CM, MTA, MTB, TM(2), T0, M0, TC(2), OS, AC, RCD,
     :  RMG, RHL, XF, DR, AK1 ,RMT, AK2, IZ(4), IB, ISX(45),
     :  TRB
      COMMON /STAT1 / CSX(10), CS(90, 127, 10), HAT(23320), NCSX
      COMMON /OP    / ZS, LEDD, VM, GR, GRAD, ETH, RLF, EGR, R, QQ
      COMMON /ATDATA/ DH2(4), CHI(26,9), OMG(27), AM(10), BN(10), JZ(10)
      COMMON /NDATA / RATEN(9000)
      COMMON /CNSTS / CPI, PI4, CLN10, CDUM(11), CSECYR, LSUN, MSUN,
     &                RSUN, TSUNYR
      COMMON /YUK1  / PX(34), WMH, WMHE, VMH, VME, VMC, VMG, BE, VLH, 
     :                VLE, VLC, VLN, VLT, MCB(12),WWW(100)
      COMMON /OPDAT / cbase,obase,opT(141),opR(31),fZ
      COMMON /XOPDAT/ opac(4,4,141,31,5)
      COMMON /COPDAT/ opacCO(4,4,141,31,305)
      COMMON /EVMODE/ IMODE
      COMMON /MIXFUD/ SGTHFAC, FACSGMIN, FACSG, ISGFAC
C extra common for mesh-spacing
      COMMON /PMESH / PMH(2), PME(2), IAGB
C first guess of pressure at H, He-burning shell (should be in input file!)
C      data pmh, pme /1.0e17, 7.5e19/
C
C Extra COMMON for main-sequence evolution.
C
      COMMON /ZAMS  / TKH(2)
      COMMON /MESH  / TRC1,TRC2,DD,DT1,DT2,MWT,MWTS,IVMC,IVMS
      COMMON /DHBLOC/ IDREDGE
      COMMON /DTCONT/ VLHP(2), VLEP(2), VLCP(2), RLFP(2), TOTMP(2), VLHC(2),
     :     VLEC(2), VLCC(2), RLFC(2), TOTMC(2)
      COMMON /ANGMOM/ VROT1, VROT2, FMAC, FAM, IRAM, IRS1, IRS2
      COMMON /DIFCOE/ DC(50,4,3), DCD(50,4)
      CBRT(VX) = DEXP(DLOG(VX)/3.0D0)
      RLOBE(VX) = 0.49D0*VX*VX/(0.6D0*VX*VX+DLOG(1.0D0+VX)) 
      DIMENSION WW(16),WX(52),DELDAT(22)
      data COcompos /0.0d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,0.6d0,1.0d0/
C99002 FORMAT (1X, 1P, 12E13.5, 0P)
99002 FORMAT (1P, 50E15.8, 0P)
99003 FORMAT (12I4,/,12I4,/,7I4,/,1P,5E8.1,0P,/,2(10I3,/,3(30I3,/)),3(15I
     :3,/), 9F5.2, 1P, 3E8.1,
     :/, E9.2, 0P, 9F6.3, /, 1P, 2(7E9.2, /), 0P, I2, 2(I2,1X,E8.2),2(1X,F4.2)
     : ,/, I2,F6.1,I2,F6.1, 1X, F4.2, I2, I2, 2(1X, E8.2),
     :/,I2,E8.1)
99004 FORMAT (1X, 10F7.3)
99005 FORMAT (1X, 1P, 2E14.6, E17.9, 3E14.6, 0P, 4I6, 1P, 2E11.3)
      IF ( IEND.NE.-1 ) GO TO 30
C Initialize physical constants
      CALL CONSTS
C Read opacity, nuclear reaction and neutrino loss rate data
      READ (11,'(I4)') NCSX
      READ (11,99004) CSX
      DO 1 N = 1, NCSX
    1    READ (11,99004) ((CS(JR,JT,N),JR=1,90),JT=1,127)
      READ (11,99004) HAT
      READ (13,99004) RATEN
C RJS 18/4/08 - read in spline coefficients for diffusion
99006 FORMAT (E12.5,3(1X,E12.5))
      DO K = 1,3
         DO I = 1,50
            READ (14,99006) (DC(I,J,K),J=1,4)
         END DO
      END DO
C d-coefficients
      DO I = 1,50
         READ (14,99006) (DCD(I,J),J=1,4)
      END DO
      DO 2 J = 1, 60
C      CT(J) = 0.0D0
         DO 2 K = 1,MAXMSH
            H(J, K) = 0.0D0
    2       DH(J, K) = 0.0D0
      DO 21 J=1,9999
         MS(J) = 0.0D0
   21    ST(J) = J

C Read in data
      READ  (1,99003) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2),ISGTH,IMO,IDIFF,
     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
     :TRB,
     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
     :ISTART, HKH

C Idiot proofing -- otherwise the logic in solver will fail
      FACSGMIN = DMIN1(1d0, FACSGMIN)

C Read first line of modin
      READ  (30, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)

C Adjust parameters if we are doing an evolution run
      IF (ISTART.EQ.1) THEN
         write(*,*) "Age, DT, NMOD overriden"
         DTY = 3e7/(SM**2d0) * HKH
         AGE = 0d0
         NMOD = 0
      END IF

      WRITE (32,99003) NH2,IT1,IT2,JIN,JOUT,NCH,JP,IZ,IMODE,
     :ICL,ION,IAM,IOP,IBC,INUC,ICN,IML(1),IML(2), ISGTH, IMO, IDIFF,
     :NT1,NT2,NT3,NT4,NT5,NSV,NMONT,
     :EP,DT3,DD,ID,ISX,DT1,DT2,CT,ZS,ALPHA,CH,CC,CN,CO, 
     :CNE,CMG,CSI,CFE,RCD,OS,RML,RMG,ECA,XF,DR,RMT,RHL,AC,AK1,AK2,ECT,
     :TRB,
     :IRAM, IRS1, VROT1, IRS2, VROT2, FMAC, FAM,
     :IVMC, TRC1, IVMS, TRC2, MWTS, IAGB, ISGFAC, FACSGMIN, SGTHFAC,
     :ISTART, HKH
      WRITE (32, 99005)
      WRITE (32, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
C Convert RML from eta to coefficient required
      RML = 4d-13*RML
C
C Create the spline interpolation data.
C
      IF (IOP .EQ. 1) CALL OPSPLN
!extra lines for COopac bit
      write(*,*) 'Selection for opacity is:',IOP
      cbase=ZS*0.173  !CC       - Changed by SMR (4/6/21) to allow for different CO abundances
      obase=ZS*0.482  !CO
      fZ=ZS
C READ IN NEW OPACITY DATA and SETUP STUFF - JJ 4/11/02
C     Read in Opal Data
C     Setup format statements
99042 FORMAT (F5.2, 31F7.3)
99043 FORMAT (5F7.3)
99045 FORMAT (3F7.3)
      if(IOP.ge.2) then
         write(*,*) 'Reading in base tables and setting up splines'
         do I=1,141
            opT(I)=3d0+0.05d0*(I-1)
         enddo
         do J=1,31
            opR(J)=-8d0+0.5d0*(J-1)
         enddo
C     Load in CO tables and setup splines
         write(*,*) 'Reading in Variable tables and setting up splines'
C         OPEN(10,FILE='COtables',STATUS='unknown',ACCESS='SEQUENTIAL')
         do K=1,305
            read(10,99045) b3,b1,b2
C            write (*,*) K,b3
            do I=1,141
               read (10,99042) temp,(opacCO(1,1,I,J,K),J=1,31)
            enddo
         enddo
C         CLOSE(10)
C     Bit to add in variable molecular bits from old paper in Marigo
C     Setup composition matrix
         if(IOP.eq.4.or.IOP.eq.6) then
            write(*,*) "Not Available"
         endif
C     Setup CO spline tables
         do K=1,305
C            write(*,*) K
C     Construct splines in T direciton
            do J=1,31
               do I=1,141
                  MAT(1,I)=opacCO(1,1,I,J,K)
               enddo
               call spline(141,opT,MAT)
               do I=1,140
                  opacCO(2,1,I,J,K)=MAT(2,I)
                  opacCO(3,1,I,J,K)=MAT(3,I)
                  opacCO(4,1,I,J,K)=MAT(4,I)
               enddo
            enddo
C     Construct splines in R direction
            do I=1,140
               do IC=1,4
                  do J=1,31
                     MAT(1,J)=opacCO(IC,1,I,J,K)
                  enddo
                  MAT(2,1)=0d0
                  MAT(3,1)=0d0
                  MAT(2,31)=0d0
                  MAT(3,31)=0d0
                  call spline(31,opR,MAT)
                  do J=1,30
                     opacCO(IC,2,I,J,K)=MAT(2,J)
                     opacCO(IC,3,I,J,K)=MAT(3,J)
                     opacCO(IC,4,I,J,K)=MAT(4,J)
                  enddo
               enddo
            enddo
         enddo
C         write (*,*) 'DONE!!'
      endif
                                !! END of new opac tables bit - JJ - 4/11/02
C
C If IAM=0, use integer atomic weights
C
      IF(IAM.EQ.0)THEN
         DO J = 1, 9
            AM(J) = BN(J)
         END DO
      END IF
C Read the initial model         
      DO  K = 1, NH
         READ (30, 99002) (H(J,K), J=1, JIN)
      END DO
C If available, read initial (last converged) changes
      DO K = 1, NH
         READ (30, 99002, END = 61, ERR = 61) (DH(J,K), J=1, JIN)
         DO 15 J = 1,JIN
            DHPR(J,K) = DH(J,K)
   15    CONTINUE
      END DO
 61   CONTINUE
C Read in first line of star 2 - most of this gets ignored
      IF (IMODE.EQ.2) THEN
         READ (50, 99005) SM2, DTY2, AGE2, PER2, BMS2, EC2,NNH2,NP2,NMOD2,IB2,PMH(2),PME(2)
         DO K = 1, NH
            READ (50, 99002) (H(J,K), J=16, JIN+15)
         END DO
      END IF
C Attempt to read in nucleosynthesis input - but don't worry if it doesn't exist.
      DO J = 1, 100
         DO K = 1, NH
            HNUC(J,K) = 0d0
            DHNUC(J,K) = 0d0
         END DO
      END DO
C Star 1 nucleosynthesis data
      READ (31, 99005, ERR = 12, END = 12)
      DO K = 1, NH
         READ(31, 99002, ERR = 12, END = 12) (HNUC(J,K), J=1, 50)
      END DO
      IF (IMODE.EQ.2) THEN
C Star 2 nucleosynthesis data
         READ (51, 99005, ERR = 12, END = 12)
         DO K = 1, NH
            READ(51, 99002, ERR = 12, END = 12) (HNUC(J,K), J=51, 100)
         END DO
      END IF
C Convert some things to `cgs' units: 10**11 cm, 10**33 gm, 10**33 erg/s 
   12 DT = CSECYR*DTY
      TM(1) = MSUN*SM
      TM(2) = MSUN*SM2
      RMG = RMG/CSECYR
      RMT = MSUN*RMT/CSECYR
      RML = RML*MSUN**2/LSUN/RSUN/CSECYR
C Optionally, re-initialise mass
      IF ( NCH.GE.1 ) THEN
         H(4,1) = DLOG(TM(1)) 
         HPR(4,1) = DLOG(TM(1)) 
         IF (IMODE.EQ.2) THEN
            H(19,1) = DLOG(TM(2))
            HPR(19,1) = DLOG(TM(2))
         END IF
      END IF
      IF (IMODE.EQ.2) THEN
         BM = MSUN*(SM + SM2)   !MSUN*BMS
      ELSE
         BM = MSUN*BMS
      END IF
      M0=BM-TM(1)
      ANG = TM(1)*(BM-TM(1))*(3.55223D0*PER/BM)**(1D0/3D0)
 13   CONTINUE
  102 FORMAT(2(F8.4,1PE16.9,4E10.3,0P7F8.5,F8.3,2F8.4,/),3F8.4,F10.4, 
     :1P2E10.3,0PF10.4,7F8.5,F8.3,2F8.4,/,F8.4,12F8.3,5F8.4,/,3F8.4,I6)
C REMESH optionally rezones the model, e.g. for different no. of meshpoints
   14 CALL REMESH ( NH2, NCH, CH, CO, CC, CNE )
      DO K=1,NH
         DO J = 1,JOUT
            HPR(J,K) = H(J,K)
         END DO
C Store nucleosynthesis
         DO J=1,100
            HNUCPR(J,K) = HNUC(J,K)
         END DO
      END DO
C COMPOS puts composition variables to zero if they are very small
      CALL COMPOS
      CALL PRINTB ( DTY, PER, NT1, NT2, NT3, NT4, NMONT, NMOD, IEND )
C
C If initial timestep is negative calculate DT as a fraction of the
C Kelvin-Helmholtz timescale and scale mass loss to evolve up the main
C sequence.
C
C Does this still work??? I never use it...
      IF (ICN .EQ. 1) DT = CSECYR*5D0*TKH(ISTAR)/(SM*SM)
      IF (ICN .EQ. 1) RMG = CLN10/(DT*2D2)
      IHOLD = 4
      CLOSE (17)
      CLOSE (20)
      CLOSE (25)
C Store certain previous values, for possible emergency restart
      PR(1) = AGE
      PR(2) = DT
      PR(3) = M1
      PR(4) = EC
      PR(5) = BM
      PR(6) = ANG      
      PR(7) = CM
      PR(8) = MTA
      PR(9) = MTB
      PR(10) = TM(1)
      PR(11) = T0
      PR(12) = M0
      PR(13) =  TC(1)
      PR(14) = TC(2)
C      DO 16 J = 1,13
C   16    PR(J) = CT(J+10)
      NPR = NMOD
      KPR = KS
      GO TO 40
C Almost end of initial input section. Start of regular update section
   30 IF ( IEND.NE.0 ) GO TO 31
      CALL COMPOS
C Store certain previous values, for possible emergency restart
      PR(1) = AGE
      PR(2) = DT
      PR(3) = M1
      PR(4) = EC
      PR(5) = BM
      PR(6) = ANG      
      PR(7) = CM
      PR(8) = MTA
      PR(9) = MTB
      PR(10) = TM(1)
      PR(11) = T0
      PR(12) = M0
      PR(13) =  TC(1)
      PR(14) = TC(2)
C      DO 10 J = 1,13
C   10    PR(J) = CT(J+10)
      NPR = NMOD
      KPR = KS
C PRINTB prints out every NT2'th meshpoint of every NT1'th model; NT3
C `pages' per printed model; also 4-line summary for every NT4'th model
      DTY = DT/CSECYR
      AGE = AGE + DTY
      NMOD = NMOD + 1
      PPER = PER
      CALL PRINTB ( DTY, PER, NT1, NT2, NT3, NT4, NMONT, NMOD, IEND )
      ANG = ANG/(1.0D0 + RHL*DTY)
      EC = EC*(1.0D0 + DTY*ECT)/(1.0D0 - DTY*ECA*EC) 
C     TRB = TRB*(1.0D0 + DTY*ECT)
C FUDGE TO DEAL WITH KS outside range. THIS ***WILL*** SCREW UP BINARIES!
C This is no longer used, so I don't care whether it works or not. RJS
      KS = 1
      IF ( AGE.GT.ST(KS+1) ) KS = KS+1
      IF ( IB.EQ.2 .AND. (ST(KS+2).EQ.0.0D0.OR.RLF.GT.0.0D0) ) STOP
      DELTA = 0.0D0
      DO 5 K = 1, NH
         DO 4 J = 1, 30 !60
C Don't use L, HORB in delta
            IF (J.NE.8.AND.J.NE.23.AND.J.NE.13.AND.J.NE.28.AND.J.NE.14.AND.J.NE.29) then
               DELTA = DELTA + DABS(DH(J,K))
            end if
            HPPR(J,K) = HPR(J,K)
            HPR(J, K) = H(J, K)
            DHPR(J,K) = DH(J,K)
    4       H(J, K) = H(J, K) + DH(J, K)
    5       IF (DTY.GT.3D-4) DELTA = DELTA + AC*DABS(DH(8,K)/HPR(8,1))
C Update nucleosynthesis matrix
      DO K = 1, NH
         DO J=1,100
            HNUCPR(J, K) = HNUC(J, K)
            DHNUCPR(J, K) = DHNUC(J, K)
            HNUC(J, K) = HNUC(J, K) + DHNUC(J, K)
C Blank DHNUC each time
C            DHNUC(J,K) = 0d0
         END DO
      END DO
      write (32,*) "DELTA =",DELTA, " DD = ", DD
      DTF = DMIN1 (DT2, DD/DELTA)      
C     IF ( IHOLD .LE. 3 ) DTF = 1.0D0
      IF ( IHOLD .LE. 2 ) DTF = 1.0D0
      DTY = DMAX1(DT1,DTF)*DTY
      DO ISTAR = 1,IMODE
         IF (IAGB.EQ.1) THEN
C Reduce timestep if He luminosity is increasing too fast -- useful on AGB
            IF ((VLEC(ISTAR) - VLEP(ISTAR))/VLEP(ISTAR).GT.0.05.AND.
     :           VLEC(ISTAR).GT.1d3) DTY = 0.8*DTY
         END IF
C Extra control mechanisms that can be uncommented as necessary - RJS
C         IF ((VLHC(ISTAR) - VLHP(ISTAR))/VLHP(ISTAR).GT.0.10) DTY = 0.8*DTY
C         IF ((VLCC(ISTAR) - VLCP(ISTAR))/VLCP(ISTAR).GT.0.10) DTY = 0.8*DTY
C         IF (((RLFC(ISTAR)-RLFP(ISTAR))/RLFP(ISTAR)).GT.0.1
CC     :        .AND.RLFP(ISTAR).GT.-1d-1.AND.RLFP(ISTAR).LT.-1d-3) THEN
C     :        .AND.RLFP(ISTAR).GT.-1d-1) THEN
C            DTY = 0.1*DTY
C            write (32,*) "RLF issues - reducing timestep"
C         END IF
      END DO
C      IF (DABS((PER - PPER)/PPER).GT.0.01) DTY = 0.5*DTY
      IF ( IB.EQ.2 ) DTY = DMIN1(DTY,ST(KS+2)-AGE)
      DT=CSECYR*DTY
C      IF (DT1.EQ.1d0) GO TO 6
C      IF (IDREDGE.EQ.3) GO TO 6      
      IF ( (JP.EQ.1 .AND. DTF.GE.DT1).OR.DTY.LT.6d-5 ) GO TO 6
C clear DH in some circumstances
      WRITE (32,*) "Clearing DH..."
      DO 7 K = 1, NH
         DO 7 J = 1, 60
    7       DH(J, K) = 0.0D0
    6 IHOLD = IHOLD + 1 
C CNO equilibrium on the main sequence.
C Does this still work???
      IF (ICN .EQ. 1) DT = MSUN*MSUN*CSECYR*5D0*TKH(ISTAR)/(VM*VM)
      IF (ICN .EQ. 1) RMG = CLN10/(DT*2D2)
      CALL COMPOS
C For *2, some factors relating to accretion from *1. Ignored if this is *1
 40   CONTINUE
C   40 T0 = CSECYR*ST(KS+1)
C      M0 = MSUN*MS(KS+1)
C RJS added to allow -C compile and run
C      IF (KS.NE.0) THEN
C         MTA = MSUN/CSECYR*(MS(KS+1)-MS(KS))/(ST(KS+1)-ST(KS))
C         MTB = MSUN/CSECYR*(MS(KS+2)-MS(KS+1))/(ST(KS+2)-ST(KS+1))
C      END IF
      IF ( MOD(NMOD,NSV).NE.0 .OR. IEND.EQ.-1 ) RETURN
C End of regular update section. Intermediate or final output section
   31 IF (IEND.EQ.2) GO TO 32
      SM = DEXP(H(4,1))/MSUN 
      SM2 = DEXP(H(19,1))/MSUN
      BMS = (SM + SM2) !BM/MSUN
      IF (IMODE.EQ.1) BMS = BM/MSUN
      DTY = DT/CSECYR
      WRITE (34, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      DO K = 1, NH
C Need to do this better - at present I'm writing out blanks
        WRITE (34, 99002) (H(J,K), J=1, 15)
      END DO
      DO K = 1, NH
         WRITE (34, 99002) (DH(J,K), J=1, 15)
      END DO
      CALL FLUSH(34)
      IF (IMODE.EQ.2) THEN
C Write out star 2
         WRITE (54, 99005) SM2, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(2),PME(2)
         DO K = 1, NH
C Copy HORB from star 1
            H(28,K) = H(13,K)
            DH(28,K) = H(28,K)
C Need to do this better - at present I'm writing out blanks
            WRITE (54, 99002) (H(J,K), J=16, 30)
         END DO
         DO K = 1, NH
            WRITE (54, 99002) (DH(J,K), J=16, 30)
         END DO
         CALL FLUSH(54)
      END IF
C write out nucleosynthesis files - star 1 first
      WRITE (35, 99005) SM, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(1)
      DO K = 1, NH
         WRITE (35, 99002) (HNUC(J,K),J=1,50)
      END DO
      CALL FLUSH(35)
      IF (IMODE.EQ.2) THEN
         WRITE (55, 99005) SM2, DTY, AGE, PER, BMS, EC,NH,NP,NMOD,IB,PMH(1),PME(2)
         DO K = 1, NH
            WRITE (55, 99002) (HNUC(J,K),J=51,100)
         END DO
         CALL FLUSH(55)
      END IF
      RETURN
C End of final output section. Start of emergency restart section 
   32 IF ( IHOLD .LE. 0 ) GO TO 34
      AGE = PR(1)
      DT = PR(2)
      M1 = PR(3)
      EC = PR(4)
      BM = PR(5)
      ANG = PR(6)      
      CM = PR(7)
      MTA = PR(8)
      MTB = PR(9)
      TM(1) = PR(10)
      T0 = PR(11)
      M0 = PR(12)
      TC(1) = PR(13)
      TC(2) = PR(14)
C      DO 11 J = 1,13
C   11    CT(J+10) = PR(J)
      NMOD = NPR
   34 DT=0.8*DT
      IF ( DT .LT. 0.01*PR(2) ) STOP
   33 DO 9 K = 1, NH
         DO 9 J = 1, JOUT
            DH(J,K) = DHPR(J,K)
    9       H(J,K) = HPR(J,K) 
            DTOLD = DTOLDP
C Sort out nucleosynthesis for restart
            DO K = 1, NH
               DO J = 1,100
                  DHNUC(J, K) = DHNUCPR(J, K)
                  HNUC(J, K) = HNUCPR(J, K)
               END DO
            END DO
      IHOLD = 0
      RETURN
      END
