      SUBROUTINE OVERFLOW(RMT_MODE, RLOF_MASS_LOSS, M1, M2, RLF, RMT)
      
      ! M1, M2 in Eggleton units (10^33 g)
      
      IMPLICIT NONE
      
      REAL*8 RLOF_MASS_LOSS, OPTICALLY_THIN_MASS_LOSS, OPTICALLY_THICK_MASS_LOSS, M1, M2, RLF, AR
      REAL*8 KB, AMU, V_SOUND, DR_RLOBE, RMT, ACCRETION_MASS_LIMIT, CLAEYS_FACTOR
      REAL*8 PNEXT, G_1, F3, GAMMA_1
      REAL*8 Q, Q_RITTER, F1, GAMMA, HP_RITTER
      REAL*8 H, DH, EPS, DEL, DH0, IW, PS, VX
      REAL*8 AR_LOBE, R2
      REAL*8 HP_0, SX, HPFUNCS
      REAL*8 OVERFLOWRLF
      REAL*8 CPI, PI4, CLN10, CA, CB, CC, CD, CG
      REAL*8 CR, CEVB, CEN, CPL, CMEVMU, CSECYR, L_SUN, M_SUN
      REAL*8 R_SUN, TSUNYR
      REAL*8 R, P, RHO, T, MU, BETA
      REAL*8 R_LOBE, R_LOBE_P, R_LOBE_RHO, R_LOBE_T, R_LOBE_MU, R_LOBE_BETA
      REAL*8 INTERPOLATED_R, INTERPOLATED_P, INTERPOLATED_RHO, INTERPOLATED_T, INTERPOLATED_MU, INTERPOLATED_BETA
      REAL*8 D_R, D_P, D_RHO, D_T, D_MU, D_BETA
      REAL*8 SURFACE_R, SURFACE_P, SURFACE_RHO, SURFACE_T, SURFACE_MU, SURFACE_BETA
      REAL*8 CORE_P, CORE_RHO, CORE_T, CORE_MU, CORE_BETA
      REAL*8 INNER_R, INNER_P, INNER_RHO, INNER_T, INNER_MU, INNER_BETA
      REAL*8 OUTER_R, OUTER_P, OUTER_RHO, OUTER_T, OUTER_MU, OUTER_BETA
      REAL*8 INTERPOLATION_SPACING, INTERPOLATION_DISTANCE
      
      INTEGER MAXMSH, K, KK, R_LOBE_MESHPOINT, N_MESH, JIN, NMOD
      INTEGER INDEX, N_POINTS_INTERPOLATION, OUTER_MESHPOINT, INNER_MESHPOINT, RMT_MODE
      
      PARAMETER (MAXMSH=2000)
      PARAMETER (N_POINTS_INTERPOLATION=100)
      
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,DEL,DH0,N_MESH,JIN,IW(200)
      COMMON /OVRFLW/ HPFUNCS, SX(45,MAXMSH+1), DR_RLOBE, HP_RITTER,
     &                OPTICALLY_THIN_MASS_LOSS, OPTICALLY_THICK_MASS_LOSS,
     &                SURFACE_P, CORE_P, SURFACE_RHO, CORE_RHO, SURFACE_T,
     &                CORE_T, SURFACE_MU, CORE_MU, SURFACE_BETA, CORE_BETA
      COMMON /OVRFW2/ R_LOBE_MESHPOINT
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CC, CD, CG, CR(2), CEVB,
     &                CEN, CPL, CMEVMU, CSECYR, L_SUN, M_SUN, R_SUN, TSUNYR
      COMMON /MISC  / NMOD
      COMMON /CHECK / OVERFLOWRLF
      
      DIMENSION R(N_MESH), P(N_MESH), RHO(N_MESH), T(N_MESH), MU(N_MESH), BETA(N_MESH)
      DIMENSION INTERPOLATED_R(N_POINTS_INTERPOLATION), INTERPOLATED_P(N_POINTS_INTERPOLATION),
     &          INTERPOLATED_RHO(N_POINTS_INTERPOLATION), INTERPOLATED_T(N_POINTS_INTERPOLATION), 
     &          INTERPOLATED_MU(N_POINTS_INTERPOLATION), INTERPOLATED_BETA(N_POINTS_INTERPOLATION)
      
      PS(VX) = 0.5D0*(VX+DABS(VX))
      F3(VX) = (VX**0.5D0) * (2D0/(VX+1D0))**((VX+1)/(2D0*(VX-1)))
      GAMMA_1(VX) = (32D0-24D0*VX-3D0*(VX*VX))/(24D0-21D0*VX)
      
      ! Drop relevant physial quantities into specifically named arrays to make code more human readable
      DO K = 1, N_MESH
            R(K) = DEXP(H(7,K))*1D11  ! cm
            P(K) = SX(2,N_MESH+2-K)   ! g/(cm s2)
            RHO(K) = SX(3,N_MESH+2-K) ! g/cm3
            T(K) = SX(4,N_MESH+2-K)   ! K
            MU(K) = SX(30,N_MESH+2-K)
            BETA(K) = SX(34,N_MESH+2-K)
      END DO
      
      ! Start with optically thin mass loss from Ritter 1988
      OPTICALLY_THIN_MASS_LOSS = 0D0
            
      ! Want Roche lobe radius in solar units
      
      AR = H(7,1) ! Get the radius of the primary in Eggleton units
      SURFACE_R = R(1)            ! cm
      R2 = SURFACE_R*SURFACE_R    ! cm2
      
      AR_LOBE = AR-RLF            ! log(Eggleton)
      R_LOBE = DEXP(AR_LOBE)*1D11 ! cm
      
      DR_RLOBE = SURFACE_R-R_LOBE         ! cm
            
      SURFACE_P = P(1)         ! g/(cm s2)
      SURFACE_RHO = RHO(1)     ! g/cm3
      SURFACE_T = T(1)         ! K
      SURFACE_MU = MU(1)
      SURFACE_BETA = BETA(1)
      
      CORE_P = P(N_MESH)
      CORE_RHO = RHO(N_MESH)
      CORE_T = T(N_MESH)
      CORE_MU = MU(N_MESH)
      CORE_BETA = BETA(N_MESH)
      
      G_1 = GAMMA_1(SURFACE_BETA)
      
      V_SOUND = (G_1*SURFACE_P/SURFACE_RHO)**(0.5D0) ! cm/s
      
      Q = M2/M1 ! Using the version of Q given by Ritter
      
      ! Restrict Q to valid range for Ritter (1988) A9
      Q_RITTER = DMIN1(DMAX1(Q, 0.5D0), 10D0)
      
      F1 = 1.23D0 + 0.5D0 * DLOG10(Q_RITTER)
      
      ! Restrict Q to valid range for Ritter (1988) Eq 7
      Q_RITTER = DMIN1(DMAX1(Q, 4D-2), 20D0)
      
      IF (Q_RITTER.LE.1d0) THEN
            GAMMA = 0.954D0 + 0.025D0 * DLOG10(Q_RITTER) - 0.038D0 * (DLOG10(Q_RITTER))**2D0
      ELSE
            GAMMA = 0.954D0 + 0.039D0 * DLOG10(Q_RITTER) + 0.114D0 * (DLOG10(Q_RITTER))**2D0
      END IF
      
      HP_0 = 0D0
      HP_0 = (8.3145D7 * SURFACE_T * (SURFACE_R)**2D0)/(CG * M1 * 1D33 * SURFACE_MU) ! cm
      
      HP_RITTER = HP_0/GAMMA ! cm
      
C     This is a bit of a hack for now, I want to change this to allow for Ritter mass loss to be toggleable with a separate variable.
      IF (RMT_MODE.GE.2) THEN
            OPTICALLY_THIN_MASS_LOSS = (2D0*CPI/DEXP(0.5D0))*(V_SOUND**3D0)*((R_LOBE**3D0)/(CG*M1*1D33))*SURFACE_RHO*F1 ! g/s
      
            ! Only allow exponential term in optically thin mass transfer when Roche lobe is underfilled (optically thin is saturated 
            ! when Roche lobe filled).
            IF (PS(RLF).EQ.0d0) THEN 
                  OPTICALLY_THIN_MASS_LOSS = OPTICALLY_THIN_MASS_LOSS * DEXP(DR_RLOBE/HP_RITTER) ! g/s
            END IF
      
            ! Rescale back to code units (Eggleton -- 10^33 g/s)
            OPTICALLY_THIN_MASS_LOSS = OPTICALLY_THIN_MASS_LOSS/1D33
      
      END IF
      
      R_LOBE_MESHPOINT = 0
      OPTICALLY_THICK_MASS_LOSS = 0D0

      IF (PS(RLF).NE.0d0) THEN
            ! RMT_MODE=0 -- Hurley; RMT_MODE=2 -- Hurley + Ritter
            IF (RMT_MODE.EQ.0.OR.RMT_MODE.EQ.2) THEN ! HPT (2002) -- old version of RLOF in the code
C Set limit for mass accretion at M/kelvin-helmholtz timescale
                  ACCRETION_MASS_LIMIT = 1d-2
                  
                  OPTICALLY_THICK_MASS_LOSS = DMIN1((RMT*((M1/M_SUN)**2d0)*((PS(RLF))**3d0)), ACCRETION_MASS_LIMIT*M_SUN/CSECYR)                  
            
            ! RMT_MODE=1 -- Claeys; RMT_MODE=3 -- Claeys + Ritter
            ELSE IF (RMT_MODE.EQ.1.OR.RMT_MODE.EQ.3) THEN ! Claeys et al (2014)
C                  RMT = M_SUN*3d-6/CSECYR ! This would overrule the data file, do I want this?
                  CLAEYS_FACTOR = 0d0
                  ACCRETION_MASS_LIMIT = 1d-2

                  IF (Q.LE.1d0) THEN
                        CLAEYS_FACTOR = 1d3 ! Set f to 1000 if Q<1
                  ELSE
                        CLAEYS_FACTOR = DMAX1(1d0, 1d3/Q * DEXP(-0.5 * (DLOG10(Q)/0.15)**2d0)) ! Set f variable in Q if Q > 1
                  END IF

                  OPTICALLY_THICK_MASS_LOSS = DMIN1((RMT*CLAEYS_FACTOR*((M1/M_SUN)**2d0)*((PS(RLF))**3d0)), ACCRETION_MASS_LIMIT*M_SUN/CSECYR)
            
            ELSE IF (RMT_MODE.EQ.4) THEN
                  DO K = 1, N_MESH ! start at surface meshpoint and go inwards
                        IF (R(K).GT.R_LOBE) THEN
                              R_LOBE_MESHPOINT = R_LOBE_MESHPOINT + 1 ! If meshpoint exceeds R_LOBE then increase the meshpoint number 
                                                                      ! of the last meshpoint to exceed the R_LOBE.
                                                                      ! Therefore, the Roche lobe sits somewhere in between 
                                                                      ! R_LOBE_MESHPOINT and R_LOBE_MESHPOINT+1
                        END IF
                  END DO
            
                  ! Find physical quantities at Roche lobe
            
                  ! Find values between meshpoints above and below Roche lobe
                  
                  OUTER_R = R(R_LOBE_MESHPOINT)
                  INNER_R = R(R_LOBE_MESHPOINT+1)
                  
                  OUTER_P = P(R_LOBE_MESHPOINT)
                  INNER_P = P(R_LOBE_MESHPOINT+1)
                  
                  OUTER_RHO = RHO(R_LOBE_MESHPOINT)
                  INNER_RHO = RHO(R_LOBE_MESHPOINT+1)
                  
                  OUTER_T = T(R_LOBE_MESHPOINT)
                  INNER_T = T(R_LOBE_MESHPOINT+1)
                  
                  OUTER_MU = MU(R_LOBE_MESHPOINT)
                  INNER_MU = MU(R_LOBE_MESHPOINT+1)
                  
                  OUTER_BETA = BETA(R_LOBE_MESHPOINT)
                  INNER_BETA = BETA(R_LOBE_MESHPOINT+1)
                  
                  D_R = INNER_R-OUTER_R
                  D_P = INNER_P-OUTER_P
                  D_RHO = INNER_RHO-OUTER_R
                  D_T = INNER_T-OUTER_T
                  D_MU = INNER_MU-OUTER_MU
                  D_BETA = INNER_BETA-OUTER_BETA
            
                  ! Recast w.r.t radius
            
                  D_P = D_P/D_R
                  D_RHO = D_RHO/D_R
                  D_T = D_T/D_R
                  D_MU = D_MU/D_R
                  D_BETA = D_BETA/D_R
                  
                  R_LOBE_P = D_P * (R_LOBE - INNER_R) + INNER_P
                  R_LOBE_RHO = D_RHO * (R_LOBE - INNER_R) + INNER_RHO
                  R_LOBE_T = D_T * (R_LOBE - INNER_R) + INNER_T
                  R_LOBE_MU = D_MU * (R_LOBE - INNER_R) + INNER_MU
                  R_LOBE_BETA = D_BETA * (R_LOBE - INNER_R) + INNER_BETA
                  
                  ! Place Roche lobe quantites into interpolated arrays as entry 1 and surface quantities as entry 
                  ! N_POINTS_INTERPOLATION
                  
                  INTERPOLATED_R(1) = R_LOBE
                  
                  INTERPOLATED_P(1) = R_LOBE_P
                  INTERPOLATED_P(N_POINTS_INTERPOLATION) = SURFACE_P
                  
                  INTERPOLATED_RHO(1) = R_LOBE_RHO
                  INTERPOLATED_RHO(N_POINTS_INTERPOLATION) = SURFACE_RHO
                  
                  INTERPOLATED_T(1) = R_LOBE_T
                  INTERPOLATED_T(N_POINTS_INTERPOLATION) = SURFACE_T
                  
                  INTERPOLATED_MU(1) = R_LOBE_MU
                  INTERPOLATED_MU(N_POINTS_INTERPOLATION) = SURFACE_MU
                  
                  INTERPOLATED_BETA(1) = R_LOBE_BETA
                  INTERPOLATED_BETA(N_POINTS_INTERPOLATION) = SURFACE_BETA
                  
                  ! Interpolate part of mesh above Roche lobe to place N_POINTS_INTERPOLATION points between Roche lobe and surface, 
                  ! evenly spaced in pressure.
                  INTERPOLATION_SPACING = (SURFACE_P-R_LOBE_P)/(N_POINTS_INTERPOLATION-1) ! cm
                  
                  ! Distribute interpolation points
                  DO INDEX = 1, N_POINTS_INTERPOLATION
                        INTERPOLATED_P(INDEX) = R_LOBE_P + (INDEX-1)*INTERPOLATION_SPACING
                  END DO
                  
                  ! Interpolate remaining quantities
                  DO INDEX = 2, N_POINTS_INTERPOLATION
                        OUTER_MESHPOINT = 0
                        DO K = 1, R_LOBE_MESHPOINT
                              IF (P(K).LT.INTERPOLATED_P(INDEX)) THEN    ! Pressure decreases monotonically, so if meshpoint pressure 
                                    OUTER_MESHPOINT = OUTER_MESHPOINT+1  ! is less than interpolated point pressure then the meshpoint
                                                                         ! is still above the interpolated point. OUTER_MESHPOINT is
                                                                         ! therefore the last meshpoint for which the pressure is less
                                                                         ! than the interpolated pressure
                              END IF
                        END DO
                        INNER_MESHPOINT = OUTER_MESHPOINT+1
                        
                        INNER_R = R(INNER_MESHPOINT)
                        OUTER_R = R(OUTER_MESHPOINT)
                        
                        INNER_P = P(INNER_MESHPOINT)
                        OUTER_P = P(OUTER_MESHPOINT)
                        
                        INNER_RHO = RHO(INNER_MESHPOINT)
                        OUTER_RHO = RHO(OUTER_MESHPOINT)
                        
                        INNER_T = T(INNER_MESHPOINT)
                        OUTER_T = T(OUTER_MESHPOINT)
                        
                        INNER_MU = MU(INNER_MESHPOINT)
                        OUTER_MU = MU(OUTER_MESHPOINT)
                        
                        INNER_BETA = BETA(INNER_MESHPOINT)
                        OUTER_BETA = BETA(OUTER_MESHPOINT)
                  
                        ! Obtain gradients between inner and outer meshpoint
                        D_R = (OUTER_R-INNER_R)/(OUTER_P-INNER_P)
                        D_RHO = (OUTER_RHO-INNER_RHO)/(OUTER_P-INNER_P)
                        D_T = (OUTER_T-INNER_T)/(OUTER_P-INNER_P)
                        D_MU = (OUTER_MU-INNER_MU)/(OUTER_P-INNER_P)
                        D_BETA = (OUTER_BETA-INNER_BETA)/(OUTER_P-INNER_P)
                  
                        ! Distance (in pressure) between interpolated point and inner meshpoint pressure
                        INTERPOLATION_DISTANCE = INTERPOLATED_P(INDEX) - INNER_P
                        
                        ! Complete interpolation
                        INTERPOLATED_R(INDEX) = D_R * INTERPOLATION_DISTANCE + INNER_R
                        INTERPOLATED_RHO(INDEX) = D_RHO * INTERPOLATION_DISTANCE + INNER_RHO
                        INTERPOLATED_T(INDEX) = D_T * INTERPOLATION_DISTANCE + INNER_T
                        INTERPOLATED_MU(INDEX) = D_MU * INTERPOLATION_DISTANCE + INNER_MU
                        INTERPOLATED_BETA(INDEX) = D_BETA * INTERPOLATION_DISTANCE + INNER_BETA
                  END DO
                  
            ! Now integrate over interpolated region of the star for optically thick mass loss from Kolb & Ritter 1990
            
            ! Compute integral in Kolb & Ritter A17 for each interpolated meshpoint above the Roche lobe
                  DO INDEX = 1, N_POINTS_INTERPOLATION-1 ! Go from surface to R_LOBE
                        G_1 = GAMMA_1(INTERPOLATED_BETA(INDEX))
                        D_P = INTERPOLATED_P(INDEX)-INTERPOLATED_P(INDEX+1)     ! g/(cm s2)
                        V_SOUND = (G_1*INTERPOLATED_P(INDEX)/INTERPOLATED_RHO(INDEX))**(0.5D0) ! cm/s
                        OPTICALLY_THICK_MASS_LOSS = OPTICALLY_THICK_MASS_LOSS + F3(G_1) * V_SOUND * D_P ! should be in cgs (g/s)?
                  END DO
                  ! Now add factors outside integrand
                  OPTICALLY_THICK_MASS_LOSS = (2D0*CPI*F1)*(R_LOBE**3D0)/(CG*M1*1D33)*OPTICALLY_THICK_MASS_LOSS ! Are these units right?
                  
                  ! Rescale back to code units (Eggleton -- 10^33 g/s)
                  OPTICALLY_THICK_MASS_LOSS = OPTICALLY_THICK_MASS_LOSS/1D33
            END IF
            
      END IF
      
      ! Merge with OPTICALLY_THIN_MASS_LOSS
      RLOF_MASS_LOSS = OPTICALLY_THIN_MASS_LOSS + OPTICALLY_THICK_MASS_LOSS
      
      RETURN
      
      END
