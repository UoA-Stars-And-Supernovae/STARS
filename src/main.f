C Between the comments added to the code and the file 'bswriteup.tex'
C you can hopefully understand how this all works... RJS 9/6/06
C PPEs original comment:
C The working of this programme is described in a file `writeup.tex';
C Consequently I do not include many comments in the body of the programme.
C Following is the main routine
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE
      PARAMETER (MAXMSH = 2000)
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EPS,V(2),NMESH,JIN,ID(100),IE(100)
      COMMON /EVMODE/ IMODE
      real dtime,cpu(2),dt,tcpu
C Read physical data and an initial model
      CALL PRINTA ( -1, NSTEP, ITER1, ITER2, NWRT5 )
      IF (NSTEP.EQ.0) GO TO 3
      ITER = ITER1
      nter = 0
      nm = 0
      npr = nm
      tcpu = dtime(cpu)
      tcpu = 0.0
C Begin evolutionary loop of NSTEP time steps
    2 IF (NM .EQ. NSTEP) GOTO 3
C Solve for structure, mesh, and major composition variables
         CALL SOLVER ( 1, ITER, KTER, ERR, ID, NWRT5 )
         IF (ERR.LT.EPS) npr = nm
         IF (ERR.LT.EPS) nm = nm + 1
         nter = nter + kter
         dt = dtime(cpu)
         tcpu = tcpu + dt
         IF (ERR.GT.EPS) kter = -kter
C         write(61,99000) nm, kter, dt, dt/kter, nter, tcpu, tcpu/nter,
C     &        tcpu/nm
99000    format(i6,i3,2f8.4,i7,f10.2,2f9.5)
         call flush(31)
C If model didnt converge, restart from 2 steps back with DT halved 
         IF (ERR.GT.EPS) nm = npr ! nm - 2
         IF (ERR.GT.EPS) CALL PRINTA ( 2, NSTEP, ITER1, ITER2, NWRT5 ) 
         IF (ERR.GT.EPS) GO TO 1
C If required, solve for minor composition vbles, with mesh, structure fixed.
C First pass for star 1
         CALL SOLVER ( 2, ITER, KTER, ERR, IE, NWRT5 )
C Restart if failed on composition
         IF (ERR.GT.EPS) nm = npr ! nm - 2
         IF (ERR.GT.EPS) CALL PRINTA ( 2, NSTEP, ITER1, ITER2, NWRT5 ) 
         IF (ERR.GT.EPS) GO TO 1
         IF (IMODE.EQ.2) THEN
C Second nucleosynthesis pass - this time for star 2
            CALL SOLVER ( 3, ITER, KTER, ERR, IE, NWRT5 )
C Restart if failed on composition
            IF (ERR.GT.EPS) nm = npr ! nm - 2
            IF (ERR.GT.EPS) CALL PRINTA ( 2, NSTEP, ITER1, ITER2, NWRT5 ) 
            IF (ERR.GT.EPS) GO TO 1
         END IF
C If model didn't converge, give up
         IF (ERR.GT.EPS) GO TO 3
         IF (ERR.LT.EPS) CALL PRINTA ( 0, NSTEP, ITER1, ITER2, NWRT5 )
    1 ITER = ITER2
      GOTO 2
C Output the last converged model. 
    3 CALL PRINTA ( 1, NSTEP, ITER1, ITER2, NWRT4 )
      STOP
      END
