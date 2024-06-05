      SUBROUTINE MDROT
c
c     ROTATE POLYATOMIC MOLECULES RANDOMLY BETWEEN [0,2PI]
c     ATTENTION!!! THIS ROUTINE DOES NOT ROTATE DIATOMIC MOLECULES
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c
      PARAMETER(NCOMPMAX=5)
      PARAMETER(NMOLMAX=256)
      PARAMETER(NSPECMAX=10)
      PARAMETER(NATSPMAX=100) 
      PARAMETER(NPARTCL=500)
c      
      PARAMETER(PI=3.141592654D0)      
c
      COMMON/INTGRS/IN(14)      
c
      COMMON/POSIT /RX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)
      COMMON/PARMOL/NMOL(NMOLMAX),NSPEC(NSPECMAX),
     &NATMSP(NCOMPMAX,NATSPMAX)
      COMMON/ATCOMP/NATMOL(NPARTCL) 
c
      DIMENSION RXOLD(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &          RYOLD(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &          RZOLD(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)
c     
      EQUIVALENCE (IN(1),NCOMP)
c
c     STORE OLD POSITIONS AND PUT LIGAND ATOMS AT THE  SYSTEM OF AXIS WITH ORIGIN AT (0,0,0)
c      
      DO 5 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 15 J=1,JMAX
            DO 25 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 35 L = 1,LMAX
                  RXOLD(I,J,K,L)=RX(I,J,K,L)-RX(I,J,1,1)
                  RYOLD(I,J,K,L)=RY(I,J,K,L)-RY(I,J,1,1)
                  RZOLD(I,J,K,L)=RZ(I,J,K,L)-RZ(I,J,1,1)
   35          CONTINUE
   25       CONTINUE      
   15    CONTINUE
    5 CONTINUE           
c
c     ATTENTION!!! KK = 2 for polyatomic molecules NATMOL > 2 hence only ligands are rotated
c     THE FUNCTION RANF RETURNS A NUMBER IN THE INTERVAL [0-1] -> (0°-360°)
c
      DO 45 I = 1,NCOMP
      IF((NMOL(I).GT.1).AND.(NATMOL(I).GT.2))THEN
c...Polyatomic molecules_rotate only ligands      
      IF(NATMOL(I).GE.3)KK=2
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 55 J=1,JMAX
c           CHOOSE A SPACE FIXED AXIS AT RANDOM
            IAXIS = IDINT(3.0D0*RANF(DUMMY))+1
c           CHOOSE A RANDOM ROTATION [0,2PI]........................... STEP 1/3
            PHI = 2.0D0*PI*RANF(DUMMY)
            DO 65 K = KK,KMAX
               LMAX = NATMSP(I,K)
               DO 75 L = 1,LMAX
c              PERFORM ROTATIONS
                  IF (IAXIS.EQ.1) THEN
c...Rotation around the x-axis
c                 CALCULATE OLD AND NEW ANGLE PHI...................... STEP 2/3
                  PHIOLD = DATAN2(RYOLD(I,J,K,L),RZOLD(I,J,K,L)) 
                  PHINEW = PHIOLD+PHI           
                  COSPHO = DCOS(PHIOLD)
                  SINPHO = DSIN(PHIOLD)
                  COSPHN = DCOS(PHINEW)
                  SINPHN = DSIN(PHINEW)
c                 CALCULATE NEW COORDINATES............................ STEP 3/3 
                  RX(I,J,K,L)=RXOLD(I,J,K,L)+RX(I,J,1,1)  
                  RY(I,J,K,L)=RYOLD(I,J,K,L)*SINPHN/SINPHO+RY(I,J,1,1)
                  RZ(I,J,K,L)=RZOLD(I,J,K,L)*COSPHN/COSPHO+RZ(I,J,1,1)
                  ELSE IF (IAXIS.EQ.2) THEN
c...Rotation around the y-axis        
c                 CALCULATE OLD AND NEW ANGLE PHI...................... STEP 2/3
                  PHIOLD = DATAN2(RZOLD(I,J,K,L),RXOLD(I,J,K,L)) 
                  PHINEW = PHIOLD+PHI           
                  COSPHO = DCOS(PHIOLD)
                  SINPHO = DSIN(PHIOLD)
                  COSPHN = DCOS(PHINEW)
                  SINPHN = DSIN(PHINEW)
c                 CALCULATE NEW COORDINATES............................ STEP 3/3   
                  RX(I,J,K,L)=RXOLD(I,J,K,L)*COSPHN/COSPHO+RX(I,J,1,1)
                  RY(I,J,K,L)=RYOLD(I,J,K,L)+RY(I,J,1,1)
                  RZ(I,J,K,L)=RZOLD(I,J,K,L)*SINPHN/SINPHO+RZ(I,J,1,1)
                  ELSE
c...Rotation around the z-axis
c                 CALCULATE OLD AND NEW ANGLE PHI...................... STEP 2/3
                  PHIOLD = DATAN2(RYOLD(I,J,K,L),RXOLD(I,J,K,L)) 
                  PHINEW = PHIOLD+PHI           
                  COSPHO = DCOS(PHIOLD)
                  SINPHO = DSIN(PHIOLD)
                  COSPHN = DCOS(PHINEW)
                  SINPHN = DSIN(PHINEW)
c                 CALCULATE NEW COORDINATES............................ STEP 3/3   
                  RX(I,J,K,L)=RXOLD(I,J,K,L)*COSPHN/COSPHO+RX(I,J,1,1)
                  RY(I,J,K,L)=RYOLD(I,J,K,L)*SINPHN/SINPHO+RY(I,J,1,1)
                  RZ(I,J,K,L)=RZOLD(I,J,K,L)+RZ(I,J,1,1)
                  ENDIF
   75          CONTINUE
   65       CONTINUE      
   55    CONTINUE
      ENDIF
   45 CONTINUE
c
      RETURN
      END    