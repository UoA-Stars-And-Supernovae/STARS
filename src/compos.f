      SUBROUTINE COMPOS

      IMPLICIT NONE

      REAL*8 DH0, GT, DEL, DHNUC, EPS, H, HNUC, EQ
      REAL*8 DH
      INTEGER JIN, J, NMESH, I, IW, K, LT, NUCMAT
      INTEGER MAXMSH, JX

      PARAMETER (MAXMSH = 2000)

      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,DEL,DH0,NMESH,JIN,IW(200)
      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      DIMENSION JX(15)
      DATA JX/3, 5, 9, 10, 11, 12, 15, 14, 15, 16, 17, 18, 19, 20 ,21/

      DO I = 1, 7
            J = JX(I)
            DO K = 1, NMESH
                  IF (H(J,K) + DH(J,K).GT.1.0D-12) THEN
                        EXIT
                  END IF

                  H(J, K) = 0.0D0
                  DH(J, K) = 0.0D0
         END DO
      END DO
! Composition check for star 2 data
      DO I = 1, 7
            J = JX(I)+ 15
            DO K = 1, NMESH
                  IF (H(J,K) + DH(J,K).GT.1.0D-12) THEN
                        EXIT
                  END IF
                  H(J, K) = 0.0D0
                  DH(J, K) = 0.0D0
            END DO
      END DO
      DO K = 1,NMESH
            IF (H(15,K)+DH(15,K).LT.1d-12) THEN
                  H(15,K) = 0d0
                  DH(15,K) = 0d0
            END IF
      END DO
! Keep angular momentum positive
      DO K = 1,NMESH
            IF (H(13,K)+DH(13,K).LT.0d0) THEN
                  H(13,K) = 0d0
                  DH(13,K) = 0d0
            END IF
            IF (H(14,K)+DH(14,K).LT.0d0) THEN
                  H(14,K) = 0d0
                  DH(14,K) = 0d0
            END IF
            IF (H(29,K)+DH(29,K).LT.0d0) THEN
                  H(29,K) = 0d0
                  DH(29,K) = 0d0
            END IF
      END DO
! Remove He tails
!      DO K = 1,NMESH
!         IF (H(9,K)+DH(9,K).LT.1d-5) THEN
!            H(9,K) = 0d0
!            DH(9,K) = 0d0
!         END IF
!      END DO
! Nucleosynthesis composition check
      DO K = 1, NMESH
            DO J=1,100
                  IF (HNUC(J,K)+DHNUC(J,K).LT.1d-20) THEN
                        HNUC(J,K) = 0d0
                        DHNUC(J,K) = 0d0
                  END IF
            END DO
      END DO
!      DO K = 1,NMESH
!         IF (H(14,K)+DH(14,K).LT.1d-100) THEN
!            H(14,K) = 0d0
!            DH(14,K) = 0d0
!         END IF
!      END DO
! Deal with problem protons in core
      DO K = 1,NMESH
            IF (H(5,K).EQ.0d0) HNUC(41,K) = 0d0
      END DO
! Fudge D, He3 in core...
!      DO K = 350,NMESH
!         IF (H(5,K).LT.1d-3) HNUC(3,K) = 0d0
!         IF (HNUC(3,K).GT.0d0) HNUC(3,K) = 0d0
!         IF (HNUC(4,K).GT.0d0) HNUC(4,K) = 0d0
!      END DO
      RETURN
      END