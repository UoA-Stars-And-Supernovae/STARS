!
!
!  ░▒▓███████▓▒░▒▓████████▓▒░▒▓██████▓▒░░▒▓███████▓▒░ ░▒▓███████▓▒░
! ░▒▓█▓▒░         ░▒▓█▓▒░  ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
! ░▒▓█▓▒░         ░▒▓█▓▒░  ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
!  ░▒▓██████▓▒░   ░▒▓█▓▒░  ░▒▓████████▓▒░▒▓███████▓▒░ ░▒▓██████▓▒░
!        ░▒▓█▓▒░  ░▒▓█▓▒░  ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░
!        ░▒▓█▓▒░  ░▒▓█▓▒░  ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░
! ░▒▓███████▓▒░   ░▒▓█▓▒░  ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░
!
!

      IMPLICIT NONE

      REAL*8 GT, V, STATUS, F10, EVMODE, ERR, DT
      REAL*4 CPU
      REAL*8 EPS, H, EQ, TCPU, DH
      INTEGER JIN, I6, ITER, NTER, NPR, I3, NMESH
      INTEGER IE, I7, NSTEP, IMERGE, NWRT4, ID, NM, KTER
      INTEGER IMODE, NWRT5, ITER2, IDET, LT, ITER1, NE, MAXMSH
      INTEGER SNAFUS, SNAFUNMOD

      PARAMETER (MAXMSH = 2000)

      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,V(2),NMESH,JIN,ID(100),IE(100)
      COMMON /EVMODE/ IMODE
      COMMON /STATUS/ IDET, IMERGE
      COMMON /ERRORS/ SNAFUS, SNAFUNMOD


!       REAL dtime,cpu(2),dt,tcpu
      DIMENSION CPU(2)

      SNAFUS = 0
      SNAFUNMOD = 0

! Read physical data and an initial model
      CALL PRINTA(-1, NSTEP, ITER1, ITER2, NWRT5)

      IF (NSTEP.NE.0) THEN ! If the step size is zero - just stop the program now.
            ITER = ITER1
            nter = 0
            nm = 0 ! Current model number
            npr = nm ! Previous model number
            tcpu = dtime(cpu) ! This is a sneaky line. dtime returns the time elapsed since this was *last* called.
                              ! In future, we want it to be "since execution began". So we call it here to update
                              ! the *next* time it is called - and manually set the execution time to zero below.
            tcpu = 0.0 ! This is the execution time of the script.

! Begin evolutionary loop of NSTEP time steps.
            DO WHILE (NM.NE.NSTEP)
! Solve for structure, mesh, and major composition variables
                  CALL SOLVER(1, ITER, KTER, ERR, ID, NWRT5)

! Store the previous model number, and increase the current one,
! if the error isn't too large.
                  IF (ERR.LT.EPS) THEN
                        npr = nm
                        nm = nm + 1
                  END IF

                  nter = nter + kter
                  dt = dtime(cpu)
                  tcpu = tcpu + dt ! Update the current execution time

                  IF (ERR.GT.EPS) THEN
                        kter = -kter
                  END IF
!         write(61,99000) nm, kter, dt, dt/kter, nter, tcpu, tcpu/nter,
!     &        tcpu/nm
99000             FORMAT(I6,I3,2F8.4,I7,F10.2,2F9.5)
                  CALL flush(31) ! 31 is nucleosynthesis file.
! If model didn't converge, restart from 2 steps back with DT halved
                  IF (ERR.GT.EPS) THEN
                        nm = npr ! nm - 2
                        CALL PRINTA(2, NSTEP, ITER1, ITER2, NWRT5)
                        GO TO 1
                  END IF
! If required, solve for minor composition variables, with mesh, structure fixed.
! First pass for star 1
                  CALL SOLVER ( 2, ITER, KTER, ERR, IE, NWRT5)
! Restart if failed on composition
                  IF (ERR.GT.EPS) THEN
                        nm = npr ! nm - 2
                        CALL PRINTA(2, NSTEP, ITER1, ITER2, NWRT5)
                        GO TO 1
                  END IF

                  IF (IMODE.EQ.2) THEN
! Second nucleosynthesis pass - this time for star 2
                        CALL SOLVER(3, ITER, KTER, ERR, IE, NWRT5)
! Restart if failed on composition
                        IF (ERR.GT.EPS) THEN
                              nm = npr ! nm - 2
                              CALL PRINTA(2, NSTEP, ITER1, ITER2, NWRT5)
                              GO TO 1
                        END IF
                  END IF
! If model didn't converge, give up
                  IF (ERR.GT.EPS) THEN
                        EXIT
                  END IF
                  IF (ERR.LT.EPS) THEN
                        CALL PRINTA(0, NSTEP, ITER1, ITER2, NWRT5)
                  END IF

    1             ITER = ITER2
            END DO
      END IF

! Output the last converged model
      CALL PRINTA(1, NSTEP, ITER1, ITER2, NWRT4)

      WRITE (*,*) "Code termination reached -- stopping."

      STOP
      END
