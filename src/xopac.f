**==XOPAC.FOR
      SUBROUTINE Xopacty (JX,fT,fR,FKL,FKH)
      IMPLICIT REAL*8 (A-H,N-Z)
      PARAMETER (MT = 141,MR = 31)
      COMMON /OPDAT / cbase,obase,TFM(141),FRM(31),fZ
      COMMON /COPDAT/ F(4,4,141,31,305)
      SAVE
*
* Calculate a bicubic spline interpolation fit for the temperature
* and density opacity fit.  Do not stop if the input lies outside the
* array but rather use the previous result.
*
* Check that we are in the table.
*
      IF((fT.LT.TFM(1)).OR.(fT.GE.TFM(MT)).OR.
     &   (fR.LT.FRM(1)).OR.(fR.GE.FRM(MR))) THEN
         FKL = FKLO
         FKH = FKHO
         WRITE (6,100) fT,fR,FKL,FKH
      ELSE
         FKLO = FKL
         FKHO = FKH
*
* Find interval in which target point lies.
*
         I=int((fT-TFM(1))*20.0)+1
         J=int((fR-FRM(1))*2.0)+1
!         I = 1 + (MT-1)*(TF-TFM(1))/(TFM(MT)-TFM(1))
!         J = 1 + (MR-1)*(FR-FRM(1))/(FRM(MR)-FRM(1))
C RJS 29/9/05 - Fudge to stop going out of array
         IF (I.LT.1) I = 1
         IF (J.LT.1) J = 1
         DT = fT-TFM(I)
         DR = (fR-FRM(J))
*
* Evaluate the splines.
*
         FKL = F(1,1,I,J,JX) + DR*(F(1,2,I,J,JX)
     &   + DR*(F(1,3,I,J,JX) + DR*F(1,4,I,J,JX)))
     &   + DT*(F(2,1,I,J,JX) + DR*(F(2,2,I,J,JX)
     &   + DR*(F(2,3,I,J,JX) + DR*F(2,4,I,J,JX)))
     &   + DT*(F(3,1,I,J,JX) + DR*(F(3,2,I,J,JX)
     &   + DR*(F(3,3,I,J,JX) + DR*F(3,4,I,J,JX)))
     &   + DT*(F(4,1,I,J,JX) + DR*(F(4,2,I,J,JX)
     &   + DR*(F(4,3,I,J,JX) + DR*F(4,4,I,J,JX))))))
*
C Stop going outside array!
         IF (JX+61.GT.305) JX = 305 - 61
         FKH = F(1,1,I,J,JX+61) + DR*(F(1,2,I,J,JX+61)
     &   + DR*(F(1,3,I,J,JX+61) + DR*F(1,4,I,J,JX+61)))
     &   + DT*(F(2,1,I,J,JX+61) + DR*(F(2,2,I,J,JX+61)
     &   + DR*(F(2,3,I,J,JX+61) + DR*F(2,4,I,J,JX+61)))
     &   + DT*(F(3,1,I,J,JX+61) + DR*(F(3,2,I,J,JX+61)
     &   + DR*(F(3,3,I,J,JX+61) + DR*F(3,4,I,J,JX+61)))
     &   + DT*(F(4,1,I,J,JX+61) + DR*(F(4,2,I,J,JX+61)
     &   + DR*(F(4,3,I,J,JX+61) + DR*F(4,4,I,J,JX+61))))))

      ENDIF
  100 FORMAT ('OPACITY OUT OF RANGE',' TF FR FKL FKH ',4F9.4)
      RETURN
      END


