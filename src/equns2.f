**==EQUNS2.FOR
      SUBROUTINE EQUNS2(K, K1, K2, NE)
      IMPLICIT REAL*8(A-H, O-Z)
      PARAMETER (MAXMSH = 2000)
      COMMON /INE   / XT(3,50), X(3,50), SG(3), DMT(3), SGTH(3), DMU(3), DSG(3,50)
      COMMON /DINE  / DXT(3, 50, 154)
      COMMON /OUTE  / EQU(50)
      COMMON /DOUTE / DEQU(50, 3, 50)
      COMMON /TRANS / HT(28,MAXMSH,2)
      COMMON /SURFCO/ SURFCOMPOS(50), ISTAROTHER
      COMMON /EVMODE/ IMODE
      PS(VX) = 0.5*(VX+DABS(VX))
C Get ISTAR from ISTAROTHER
      IF (ISTAROTHER.EQ.1) ISTAR = 2
      IF (ISTAROTHER.EQ.2) ISTAR = 1
C next to surface BC
      IF (K.EQ.K1+1) THEN
            DO J = 1, NE
               DO I = 1, NE
                  DEQU(I, 1, J) = 0.0
                  DEQU(I, 2, J) = 0.0
                  DEQU(I, 3, J) = 0.0
               END DO
               IF (DMT(2).LE.3d-16) THEN
                  DEQU(J, 2, J) = DEQU(J, 2, J) - 1.0
                  DEQU(J, 3, J) = 1.0
                  EQU(J) = X(3,J) - X(2,J)
C               END IF
               ELSE
C Boundary condition for fudged surface accretion
                  DO I = 1, NE
                     DEQU(I, 2, J) = -DXT(2,I,J)
                  END DO
C                  SG2 = 0.5d0*(SG(3)+SG(2)) - PS(-DMT(3))
C                  SG2 = SG2 + 0.5*(SGTH(2)+SGTH(3))*PS((DMU(2)-DMU(3)))
                  SG2 = -1d-5
C                  DEQU(J, 2, J) = DEQU(J, 2, J) - 1.0
C                  DEQU(J, 3, J) = 0.0
                  DEQU(J,2,J) = DEQU(J,2,J) - SG2 + PS(DMT(2))
                  DEQU(J,3,J) = SG2
                  EQU(J) = SG2*(X(3,J)-X(2,J)) + PS(DMT(2))*(X(2,J)-SURFCOMPOS(J)) - XT(2,J)
               END IF
            END DO
      END IF
      IF ( K.GT.K1+1 ) THEN
         IF ( K.LE.K2 ) THEN
C Interior points
            SG1 = 0.5d0*(SG(2)+SG(1)) - PS(DMT(2))
            SG2 = 0.5d0*(SG(3)+SG(2)) - PS(-DMT(3))
C Plus thermohaline mixing
            SG1 = SG1 + 0.5*(SGTH(1)+SGTH(2))*PS((DMU(1) - DMU(2)))
            SG2 = SG2 + 0.5*(SGTH(2)+SGTH(3))*PS((DMU(2)-DMU(3)))
            DO J = 1, NE
               DO I = 1, NE
                  DEQU(I, 2, J) = -DXT(2, I, J)
               END DO
C For all species but H
               IF (J.NE.41) THEN
                  DEQU(J, 1, J) = SG1
                  DEQU(J, 2, J) = DEQU(J, 2, J) - SG1 - SG2 - DSG(2,J) 
                  DEQU(J, 3, J) = SG2 + DSG(3,J)
                  EQU(J) = ((X(3,J)-X(2,J))*SG2 - (X(2,J)-X(1,J))
     :                 *SG1 - XT(2, J)) - DSG(2,J)*X(2,J) + DSG(3,J)*X(3,J)
               ELSE
                  DEQU(J, 1, J) = SG1
                  DEQU(J, 2, J) = (DEQU(J, 2, J) - SG1 - SG2)
                  DEQU(J, 3, J) = SG2
                  EQU(J) = ((X(3,J)-X(2,J))*SG2 - (X(2,J)-X(1,J))
     :                 *SG1 - XT(2, J))
C Diffusion of all other species to H 
                 DO I=1,NE
                     IF (I.NE.41) THEN
                        DEQU(I, 2, J) = DEQU(I,2,J) + DSG(2,I)
                        DEQU(I, 3, J) = -DSG(3,I)*X(3,I)
                        EQU(J) = EQU(J) + DSG(2,I)*X(2,I) - DSG(3,I)*X(3,I)
                     END IF
                  END DO
               END IF
C               DEQU(J, 1, J) = SG1
C               DEQU(J, 2, J) = (DEQU(J, 2, J) - SG1 - SG2)
C               DEQU(J, 3, J) = SG2
C               EQU(J) = ((X(3,J)-X(2,J))*SG2 - (X(2,J)-X(1,J))
C     :              *SG1 - XT(2, J))
            END DO
         ELSE
C Central conditions
            SG2 = 0.5d0*(SG(3)+SG(2))
C Plus thermohaline mixing
            SG2 = SG2 + 0.5*(SGTH(2)+SGTH(3))*PS((DMU(2)-DMU(3)))
            DO J = 1, NE
               DO I = 1, NE
                  DEQU(I, 1, J) = 0.0
                  DEQU(I, 2, J) = 0.0
                  DEQU(I, 3, J) = DXT(3, I, J)
               END DO
               DEQU(J, 2, J) = - SG2
               DEQU(J, 3, J) = (DEQU(J, 3, J) + SG2)     
               EQU(J) = (SG2*(X(3,J)-X(2,J)) + XT(3, J))
            END DO
         END IF
      END IF
 500  FORMAT (I4,17e15.6)
      RETURN
      END
