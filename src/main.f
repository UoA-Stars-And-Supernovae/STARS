C          _______.___________.    ___      .______          _______.
C         /       |           |   /   \     |   _  \        /       |
C        |   (----`---|  |----`  /  ^  \    |  |_)  |      |   (----`
C         \   \       |  |      /  /_\  \   |      /        \   \
C     .----)   |      |  |     /  _____  \  |  |\  \----.----)   |
C     |_______/       |__|    /__/     \__\ | _| `._____|_______/
C
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE
      PARAMETER (MAXMSH = 2000)
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,V(2),NMESH,JIN,ID(100),IE(100)
      COMMON /EVMODE/ IMODE
      COMMON /STATUS/ IDET, IMERGE, IDOKR

      real dtime,cpu(2),dt,tcpu

C Read physical data and an initial model
      CALL PRINTA ( -1, NSTEP, ITER1, ITER2, NWRT5 )
      IF (NSTEP.EQ.0) GO TO 3 ! If the step size is zero - just stop the program now.
      ITER = ITER1
      nter = 0
      nm = 0 ! Current model number
      IDOKR = 0
      npr = nm ! Previous model number
      tcpu = dtime(cpu) ! This is a sneaky line. dtime returns the time elapsed since this was *last* called.
                        ! In future, we want it to be "since execution began". So we call it here to update
                        ! the *next* time it is called - and manually set the execution time to zero below.
      tcpu = 0.0 ! This is the execution time of the script.

C Begin evolutionary loop of NSTEP time steps.
C This is effectively the world's most convoluted while loop.
    2 IF (NM .EQ. NSTEP) GOTO 3 ! 3 is past the end of the "loop".
C Solve for structure, mesh, and major composition variables
         CALL SOLVER ( 1, ITER, KTER, ERR, ID, NWRT5 )

C Store the previous model number, and increase the current one,
C if the error isn't too large.
         IF (ERR.LT.EPS) THEN
            npr = nm
            nm = nm + 1
         END IF

         nter = nter + kter
         dt = dtime(cpu)
         tcpu = tcpu + dt ! Update the current execution time

         IF (ERR.GT.EPS) kter = -kter
C         write(61,99000) nm, kter, dt, dt/kter, nter, tcpu, tcpu/nter,
C     &        tcpu/nm
99000    format(i6,i3,2f8.4,i7,f10.2,2f9.5)
         call flush(31) ! 31 is nucleosynthesis file.
C If model didnt converge, restart from 2 steps back with DT halved
         IF (ERR.GT.EPS) THEN
            nm = npr ! nm - 2
            CALL PRINTA ( 2, NSTEP, ITER1, ITER2, NWRT5 )
            GO TO 1
         END IF
C If required, solve for minor composition vbles, with mesh, structure fixed.
C First pass for star 1
         CALL SOLVER ( 2, ITER, KTER, ERR, IE, NWRT5 )
C Restart if failed on composition
         IF (ERR.GT.EPS) THEN
            nm = npr ! nm - 2
            CALL PRINTA ( 2, NSTEP, ITER1, ITER2, NWRT5 )
            GO TO 1
         END IF

         IF (IMODE.EQ.2) THEN
C Second nucleosynthesis pass - this time for star 2
            CALL SOLVER ( 3, ITER, KTER, ERR, IE, NWRT5 )
C Restart if failed on composition
            IF (ERR.GT.EPS) THEN
               nm = npr ! nm - 2
               CALL PRINTA ( 2, NSTEP, ITER1, ITER2, NWRT5 )
               GO TO 1
            END IF
         END IF
C If model didn't converge, give up
         IF (ERR.GT.EPS) GO TO 3
         IF (ERR.LT.EPS) CALL PRINTA ( 0, NSTEP, ITER1, ITER2, NWRT5 )
    1 ITER = ITER2
      GOTO 2 ! go back to start of "loop"
C Output the last converged model.
    3 CALL PRINTA ( 1, NSTEP, ITER1, ITER2, NWRT4 )
      STOP
      END
