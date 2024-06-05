c      
      SUBROUTINE GVERLET
c
c     INTEGRATE Newton EQUATIONS OF MOTION WITH A Verlet "LEAP-FROG" ALGORITHM    
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c
      PARAMETER(NCOMPMAX=5)
      PARAMETER(NMOLMAX=256)
      PARAMETER(NSPECMAX=10)
      PARAMETER(NATSPMAX=100)
      PARAMETER(NLGAUSS=500)
      PARAMETER(NPARTCL=500)
c      
      PARAMETER(BANG=0.5291772083D0)
c
      CHARACTER*150 NLINEG 
      CHARACTER*4   NSYMBL
c
      LOGICAL LG(5),LSCALE            
c                                               
      COMMON/INTGRS/IN(14)
      COMMON/LGCLS /LG
      COMMON/REALS /RL(18)
c           
c...  Positions at step N and N+1..........[after Verlet]     
      COMMON/POSIT /RX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)
c...  Positions at step N stored in a one-dimension array     
      COMMON/POSSD /XX(NPARTCL),YY(NPARTCL),ZZ(NPARTCL)
c...  Positions at step N stored in a three-dimensional array      
      COMMON/POSTD /RXX(NCOMPMAX,NMOLMAX,NPARTCL),
     &              RYY(NCOMPMAX,NMOLMAX,NPARTCL),
     &              RZZ(NCOMPMAX,NMOLMAX,NPARTCL)
c...  Velocities at step N-1/2 and N+1/2...[after Verlet]      
      COMMON/VELOC /VX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              VY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              VZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)
c...  Forces at step N and N+1.............[after Gaussian]           
      COMMON/FORCE /FX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              FY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              FZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX) 
c...  Velocities at step N.................[after Verlet]    
      COMMON/VELOCN/VNX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              VNY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              VNZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX) 
c
      COMMON/PARMOL/NMOL(NMOLMAX),NSPEC(NSPECMAX),
     &NATMSP(NCOMPMAX,NATSPMAX)
      COMMON/ATMASS/ZMASS(NCOMPMAX,NSPECMAX,NATSPMAX)
c...  Atomic numbers      
      COMMON/ATOMZ /NATZ(NCOMPMAX,NSPECMAX,NATSPMAX)
c...  Atomic numbers stored in a one-dimension array      
      COMMON/ATNZSD/ATNZZ(NPARTCL)
c...  Number of atoms p/ molecule for each component      
      COMMON/ATCOMP/NATMOL(NPARTCL)
c...  Chemical symbol of the atomic species suitable for using in a double loop      
      COMMON/ATSYMB/NSYMBL(NCOMPMAX,NPARTCL)
      COMMON/GAUSSL/NLINEG(NLGAUSS)
c
      DIMENSION RTRA(NCOMPMAX,NMOLMAX,NATSPMAX,NATSPMAX)      
c
      EQUIVALENCE (IN(1),NCOMP),(IN(2),NLG),(IN(3),NATOMS)
      EQUIVALENCE (IN(4),NATOMS1),(IN(5),MDSTEP),(IN(6),MAXEQB)
      EQUIVALENCE (IN(7),NSTEP),(IN(8),NMOLCS),(IN(9),LTYPE)
      EQUIVALENCE (IN(11),NGEN),(IN(12),NGA),(IN(13),NGB)
      EQUIVALENCE (LG(3),LSCALE)
      EQUIVALENCE (RL(1),ENRNUC),(RL(2),ENRSCF),(RL(3),STEP) 
      EQUIVALENCE (RL(4),STEPSQ),(RL(5),TKECONV),(RL(6),ACTKE)
      EQUIVALENCE (RL(7),DTEMP),(RL(8),ATEMP),(RL(9),ACTKEAV)
      EQUIVALENCE (RL(10),ENSCFAV),(RL(11),FNATOMS),(RL(12),FEKJM)
      EQUIVALENCE (RL(13),COULV),(RL(15),DIST),(RL(16),DTE)
      EQUIVALENCE (RL(17),DTP),(RL(18),BOX) 
c      
      OPEN(3,FILE='MDPOST.DAT'  ,ACCESS='SEQUENTIAL',STATUS='UNKNOWN') 
      OPEN(4,FILE='MDVELT.DAT'  ,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
      OPEN(5,FILE='gaussr.com'  ,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')      
c
c     WRITE POSITIONS, VELOCITIES AND FORCES WITH SIMILAR FORMAT TO THAT OF ABINIT
c 
c     ....................................
c     Cartesian coordinates (bohr)
c     RX, RY, RZ
c     Cartesian forces (hart/bohr)
c     FX, FY, FZ 
c     Velocities (bohr/(atomic time unit))
c     VX, VY, VZ
c     ....................................
c
c     WRITE POSITIONS AT STEP N TO MDOUTPUT.DAT AND INITIALIZE TO ZERO THE VELOCITIES AT STEP N
c
      IF((NSTEP-1).LE.MAXEQB)THEN
         WRITE(6,*)'Equilibration cartesian coordinates (bohr)'
      ELSE
         WRITE(6,*)'Production cartesian coordinates (bohr)'
      ENDIF
c      
      DO 5 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 15 J=1,JMAX 
            DO 25 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 35 L = 1,LMAX
                  WRITE(6,'(3(1X,E21.14))')RX(I,J,K,L),RY(I,J,K,L),
     &RZ(I,J,K,L)
                  VNX(I,J,K,L) = 0.0D0
                  VNY(I,J,K,L) = 0.0D0
                  VNZ(I,J,K,L) = 0.0D0
   35          CONTINUE
   25       CONTINUE      
   15    CONTINUE
    5 CONTINUE
c
c     CALL GEWALD TO SUBTRACT THE Coulomb NUCLEAR REPULSION ENERGY AND FORCE COMPONENTS AND ADD THE 
c     Coulomb-Ewald ENERGY AND FORCES CALCULATED ONLY BETWEEN ATOMS OF DIFFERENT MOLECULES......Version 2.0 of MDCLG.f
c    
***** CALL GEWALD
c
c     WRITE Ab-Initio FORCES AT STEP N TO MDOUTPUT.DAT
c
      IF((NSTEP-1).LE.MAXEQB)THEN
         WRITE(6,*)'Equilibration cartesian forces (hart/bohr)'
      ELSE
         WRITE(6,*)'Production cartesian forces (hart/bohr)'
      ENDIF
c      
      DO 45 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 55 J=1,JMAX 
            DO 65 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 75 L = 1,LMAX
                  WRITE(6,'(3(1X,E21.14))')FX(I,J,K,L),FY(I,J,K,L),
     &FZ(I,J,K,L)
   75          CONTINUE
   65       CONTINUE      
   55    CONTINUE
   45 CONTINUE
c
c     APPLY VERLET "LEAP-FROG" ALGORITHM
c     CALCULATE VELOCITIES AT STEP (N+1/2) AND (N) AND ACTUAL KINETIC ENERGY FROM THE VELOCITIES AT STEP (N)
c
      ACTKE= 0.0D0
      PX   = 0.0D0
      PY   = 0.0D0
      PZ   = 0.0D0
c     
      IF((NSTEP-1).LE.MAXEQB)THEN
         WRITE(6,*)'Equilibration velocities (bohr/(atomic time unit))'
      ELSE
         WRITE(6,*)'Production velocities (bohr/(atomic time unit))'
      ENDIF 
c      
      DO 85 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 95 J=1,JMAX 
            DO 105 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 115 L = 1,LMAX
c...  Copy velocities at step (N-1/2)
                  VHX = VX(I,J,K,L)
                  VHY = VY(I,J,K,L)
                  VHZ = VZ(I,J,K,L)
c...  Check for conservation of the total linear momentum
                  PX = PX + ZMASS(I,K,L)*VX(I,J,K,L)
                  PY = PY + ZMASS(I,K,L)*VY(I,J,K,L)
                  PZ = PZ + ZMASS(I,K,L)*VZ(I,J,K,L)
c...  Verlet "leap-frog" for velocities  V(N+1/2)=V(N-1/2)+(F(N)/m)*STEP       
                  VX(I,J,K,L)=VX(I,J,K,L)+FX(I,J,K,L)*STEP/ZMASS(I,K,L)
                  VY(I,J,K,L)=VY(I,J,K,L)+FY(I,J,K,L)*STEP/ZMASS(I,K,L)
                  VZ(I,J,K,L)=VZ(I,J,K,L)+FZ(I,J,K,L)*STEP/ZMASS(I,K,L)  
c...  First order velocities at step N   V(N)= 0.5*(V(N-1/2) + V(N+1/2))
                  VNX(I,J,K,L)=0.5D0*(VHX+VX(I,J,K,L))
                  VNY(I,J,K,L)=0.5D0*(VHY+VY(I,J,K,L))
                  VNZ(I,J,K,L)=0.5D0*(VHZ+VZ(I,J,K,L))
c...  Calculate instantaneous actual kinetic energy using velocities at step N
                  ACTKE=ACTKE+0.5D0*ZMASS(I,K,L)*
     &(VNX(I,J,K,L)**2.0D0+VNY(I,J,K,L)**2.0D0+VNZ(I,J,K,L)**2.0D0)
c
c     WRITE ATOMIC VELOCITIES AT STEP N TO MDOUTPUT.DAT
c
                  WRITE(6,'(3(1X,E21.14))')VNX(I,J,K,L),VNY(I,J,K,L),
     &VNZ(I,J,K,L)
c 
  115          CONTINUE
  105       CONTINUE      
   95    CONTINUE
   85 CONTINUE  
c
c     WRITE INSTANTANEOUS AND AVERAGED ENERGIES AND TEMPERATURE IN Hartrees AND Kelvins
c
      DENOM  = DFLOAT(MDSTEP)
c      
*     ATEMP  = 2.0D0*ACTKE/(3.0D0*FNATOMS*BOLTZ)
      ATEMP  = ACTKE*TKECONV
      TDRIFT = ATEMP-DTEMP
      TOTENR = ACTKE + ENRSCF
      TEKJM  = TOTENR*FEKJM
c
c...  CUMULATIVE AVERAGES OF THE KINETIC ENERGY AND SCF ENERGY
c
      ACTKEAV = ACTKEAV + ACTKE
      ENSCFAV = ENSCFAV + ENRSCF
c      
      AVACTKE = ACTKEAV/DENOM
      AVENSCF = ENSCFAV/DENOM
*     AVATEMP = 2.0D0*AVACTKE/(3.0D0*FNATOMS*BOLTZ)
      AVATEMP = AVACTKE*TKECONV
      AVTOTE  = AVACTKE + AVENSCF
c
c...  WRITE NEW ENERGY SUMS TO mdcontrl.dat TO ACCUMULATE FOR THE NEXT TIME-STEP
c     RE-SET TO ZERO THE ENERGY SUMS AT END OF EQUILIBRATION
c 
      IF((NSTEP-1).NE.MAXEQB)THEN
         WRITE(7,*)ACTKEAV,ENSCFAV
      ELSE
         WRITE(7,*)0.0D0,0.0D0
      ENDIF
      
      WRITE(6,39)
      WRITE(6,*)'Components of total energy and temperature 
     &(in Hartree and K):'
      WRITE(6,*)'                                                     '
      WRITE(6,'(A27,E21.12)')'Actual temperature       = ',ATEMP
      WRITE(6,'(A27,E21.12)')'Nuclear kinetic energy   = ',ACTKE
      WRITE(6,'(A27,E21.12)')'Total SCF energy         = ',ENRSCF
      WRITE(6,'(A27,E21.12)')'>>>>>>>>>>>>Total energy = ',TOTENR
      WRITE(6,49)
      WRITE(6,'(A27,E21.12)')'<Actual temperature>     = ',AVATEMP
      WRITE(6,'(A27,E21.12)')'<Nuclear kinetic energy> = ',AVACTKE
      WRITE(6,'(A27,E21.12)')'<Total SCF energy>       = ',AVENSCF
      WRITE(6,'(A27,E21.12)')'>>>>>>>>>><Total energy> = ',AVTOTE
      WRITE(6,*)'                                                     '
      WRITE(6,*)'Other information on the energy and temperature:'
      WRITE(6,*)'                                                     '
      WRITE(6,'(A27,F11.3)')'Desired temperature  [K] = ',DTEMP
      WRITE(6,'(A27,F11.3)')'Actual temperature   [K] = ',ATEMP
      WRITE(6,'(A27,F11.3)')'Temperature drift    [K] = ',TDRIFT
      WRITE(6,'(A27,F11.3)')'<Actual temperature> [K] = ',AVATEMP
      WRITE(6,'(A27,E14.5)')'>>>Total energy [kJ/mol] = ',TEKJM
c
c     ---------------------------------------------------------------------------------HYDROGEN BOND ANALYSIS_START
c      
      WRITE(6,39)
      WRITE(6,*)'Histograms of intra and inter molecular distances and 
     &angles:'
      WRITE(6,*)'                                                     '
c...Calculate time elapsed during simulation in ps [1fs=10¯³ps] for building histograms      
      IF((NSTEP-1).LE.MAXEQB)THEN
         WRITE(6,'(A50)')'Intra-molecular distances (Ang) and angles 
     &(°) Eql'
         TMEPS = 1.0D-3*DTE*DFLOAT(MDSTEP)
         WRITE(6,'(A12,E11.3)')'Time (ps) = ',TMEPS
      ELSE
         WRITE(6,'(A50)')'Intra-molecular distances (Ang) and angles 
     &(°) Prd'
         TMEPS = 1.0D-3*DTP*DFLOAT(MDSTEP)
         WRITE(6,'(A12,E11.3)')'Time (ps) = ',TMEPS
      ENDIF
c
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Calculate internal angles
c      
c...Calculate intra-molecular distances
      DO 125 I = 1,NCOMP     
c...Monatomic molecules     
      IF(NATMOL(I).EQ.1)THEN
         WRITE(6,*)'                                                  '
         WRITE(6,'((A9),(I3))')'Component',I   
         WRITE(6,'(A12)')'------------'           
         WRITE(6,*)'                                                  '
         WRITE(6,'(A20,I3,2X,A25)')'WARNING!!! Component',I,'has no 
     &internal structure' 
      GOTO 125
c...Diatomic molecules - linear_ print out all the distances      
      ELSEIF(NATMOL(I).EQ.2)THEN
         WRITE(6,*)'                                                  '
         WRITE(6,'((A9),(I3))')'Component',I   
         WRITE(6,'(A12)')'------------'       
         WRITE(6,*)'                                                  '    
         WRITE(6,'(2(5X,A5),2(4X,A5))')'  C  ','  M  ','  A  ','  A  '
      ENDIF 
c      
      JMAX  = NMOL(I)
      LMAX = NATMOL(I)
         DO 135 J = 1,JMAX
            DO 145 LA = 1,LMAX-1
            LAB = LA + 1
               DO 155 LB = LAB, LMAX
                  RAXIJ  = RXX(I,J,LA)-RXX(I,J,LB)
                  RAYIJ  = RYY(I,J,LA)-RYY(I,J,LB)
                  RAZIJ  = RZZ(I,J,LA)-RZZ(I,J,LB)
                  RTRASQ = RAXIJ*RAXIJ + RAYIJ*RAYIJ + RAZIJ*RAZIJ     
                  RTRA(I,J,LA,LB) = BANG*DSQRT(RTRASQ)
                  IF(NATMOL(I).EQ.2)THEN                                  
                  WRITE(6,'(4(5X,I4),7X,(A2),(A1),(A2),F15.7)')I,J,LA,
     &LB,NSYMBL(I,LA),'-',NSYMBL(I,LB),RTRA(I,J,LA,LB)  
                  ENDIF
  155          CONTINUE
  145       CONTINUE      
  135    CONTINUE
  125 CONTINUE  
c
c...Calculate intra-molecular angles using the "Law of cosines"
c   Law of cosines - C² = A²+B² - 2ABCos(c) where A,B,C triangle sides and c angle opposite to C
c  
      DO 165 I = 1,NCOMP     
c...Polyatomic molecules - 3 Atoms - linear (e.g., HCN, CO2), angular (H2O, O3)  
c                        - 4 Atoms - trigonal planar (e.g.,BF3, SO3), trigonal pyramidal (e.g., NH3)
c                        - 5 Atoms - tetrahedral (e.g., CH4), square planar (e.g., XeF4)
c                        - 6 Atoms - square pyramidal (e.g., Sb(Ph)5), trigonal bipyramidal (e.g., PCl5)
c                        - 7 Atoms - Octahedral (e.g., SF6)     
         IF(NATMOL(I).GE.3)THEN
            WRITE(6,*)'                                               '
            WRITE(6,'((A9),(I3))')'Component',I   
            WRITE(6,'(A12)')'------------'       
            WRITE(6,*)'                                               '
         ELSE
         GOTO 165
         ENDIF
         IF(NATMOL(I).EQ.3)THEN
         WRITE(6,'(2(5X,A5),5X,3(9X,(A2),(A1),(A2)),9X,4(A2))')
     &'   C ','  M  ',NSYMBL(I,1),'-',NSYMBL(I,2),NSYMBL(I,1),'-',
     &NSYMBL(I,3),NSYMBL(I,2),'-',NSYMBL(I,3),'<)',NSYMBL(I,2),
     &NSYMBL(I,1),NSYMBL(I,3)
         ELSEIF(NATMOL(I).GT.3)THEN
         WRITE(6,'(2(3X,A5),6(4X,(A2),(A1),(A2)),3X,3(3X,4(A2)))')
     &'   C ','  M  ',NSYMBL(I,1),'-',NSYMBL(I,2),NSYMBL(I,1),'-',
     &NSYMBL(I,3),NSYMBL(I,1),'-',NSYMBL(I,4),
     &NSYMBL(I,2),'-',NSYMBL(I,3),NSYMBL(I,3),'-',NSYMBL(I,4),
     &NSYMBL(I,2),'-',NSYMBL(I,4),'<)',NSYMBL(I,2),
     &NSYMBL(I,1),NSYMBL(I,3),'<)',NSYMBL(I,3),NSYMBL(I,1),NSYMBL(I,4),
     &'<)',NSYMBL(I,2),NSYMBL(I,1),NSYMBL(I,4)
         ENDIF
            JMAX  = NMOL(I)
            DO 175 J = 1,JMAX   
                  IF(NATMOL(I).EQ.3)THEN            
                  ASQ   = RTRA(I,J,1,2)**2.0D0
                  BSQ   = RTRA(I,J,1,3)**2.0D0  
                  CSQ   = RTRA(I,J,2,3)**2.0D0
                  ABTWO = 2.0D0*RTRA(I,J,1,2)*RTRA(I,J,1,3)
                  ARGM  = (ASQ+BSQ-CSQ)/ABTWO
                  PHI   = DACOSD(ARGM)
                  WRITE(6,'(2(5X,I4),7X,4(F15.7))')I,J,RTRA(I,J,1,2),
     &RTRA(I,J,1,3),RTRA(I,J,2,3),PHI 
                  ELSEIF(NATMOL(I).GT.3)THEN  
                  ASQ1   = RTRA(I,J,1,2)**2.0D0
                  BSQ1   = RTRA(I,J,1,3)**2.0D0  
                  CSQ1   = RTRA(I,J,2,3)**2.0D0
                  AB1TWO = 2.0D0*RTRA(I,J,1,2)*RTRA(I,J,1,3)
                  ARGM1  = (ASQ1+BSQ1-CSQ1)/AB1TWO
                  PHI1   = DACOSD(ARGM1)
*
                  ASQ2   = RTRA(I,J,1,3)**2.0D0
                  BSQ2   = RTRA(I,J,1,4)**2.0D0  
                  CSQ2   = RTRA(I,J,3,4)**2.0D0
                  AB2TWO = 2.0D0*RTRA(I,J,1,3)*RTRA(I,J,1,4)
                  ARGM2  = (ASQ2+BSQ2-CSQ2)/AB2TWO
                  PHI2   = DACOSD(ARGM2)  
*
                  ASQ3   = RTRA(I,J,1,2)**2.0D0
                  BSQ3   = RTRA(I,J,1,4)**2.0D0  
                  CSQ3   = RTRA(I,J,2,4)**2.0D0
                  AB3TWO = 2.0D0*RTRA(I,J,1,2)*RTRA(I,J,1,4)
                  ARGM3  = (ASQ3+BSQ3-CSQ3)/AB3TWO
                  PHI3   = DACOSD(ARGM3)
                  WRITE(6,'(2(3X,I4),2X,6(1X,F8.3),2X,3(3X,F8.3))')I,J,
     &RTRA(I,J,1,2),RTRA(I,J,1,3),RTRA(I,J,1,4),RTRA(I,J,2,3),
     &RTRA(I,J,3,4),RTRA(I,J,2,4),PHI1,PHI2,PHI3            
                  ENDIF  
  175       CONTINUE
  165 CONTINUE
c
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Count Hydrogen bonds
c               
c...Calculate inter-molecular distances and introduce distance and angle criteria for definition of
c   hydrogen bonding between molecules of the same component for every component         
      WRITE(6,49)
      IF((NSTEP-1).LE.MAXEQB)THEN
         WRITE(6,'(A35)')'Inter-molecular distances Eql (Ang)'
         WRITE(6,'(A12,E11.3)')'Time (ps) = ',TMEPS
      ELSE
         WRITE(6,'(A35)')'Inter-molecular distances Prd (Ang)'
         WRITE(6,'(A12,E11.3)')'Time (ps) = ',TMEPS
      ENDIF
c
*****      WRITE(6,'(2(4X,A5),3(3X,A5))')'  C  ','  M  ','  A  ','  M  ',
*****     &'  A  ' 
c...This loop does not calculate distances for molecules (atoms) of different components      
      DO 185 I = 1,NCOMP
      JMAX  = NMOL(I)
      LMAX = NATMOL(I)
         DO 195 JA = 1,JMAX-1
         JAB = JA + 1
            DO 205 LA = 1,LMAX
               DO 215 JB = JAB,JMAX
                  DO 225 LB = 1, LMAX
                     REXIJ  = RXX(I,JA,LA)-RXX(I,JB,LB)
                     REYIJ  = RYY(I,JA,LA)-RYY(I,JB,LB)
                     REZIJ  = RZZ(I,JA,LA)-RZZ(I,JB,LB)
                     RTERSQ = REXIJ*REXIJ + REYIJ*REYIJ + REZIJ*REZIJ     
                     RTER   = BANG*DSQRT(RTERSQ)
*****                        WRITE(6,'(5(5X,I3),7X,(A2),(A1),(A2),F15.7)')I,
*****     &JA,LA,JB,LB,NSYMBL(I,LA),'-',NSYMBL(I,LB),RTER
  225             CONTINUE
  215          CONTINUE
  205       CONTINUE      
  195    CONTINUE
  185 CONTINUE
c
c               
c...Calculate inter-molecular distances and introduce distance and angle criteria for definition of
c   hydrogen bonding between molecules of different component for the case of a mixture    
c   This loop calculates the inter-molecular distances for molecules (atoms) of different components
c
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Count Hydrogen bonds_not coded
*****      WRITE(6,49)
c     ---------------------------------------------------------------------------------HYDROGEN BOND ANALYSIS_END
c       
      WRITE(6,39)
c
c     WARNINGS!!! WRITE WARNINGS TO MDOUTPUT.DAT IF SOMETHING IS POSSIBLY WRONG
c     .............................................................................  
c     WARN_1. WARN ABOUT POSSIBLE "DIMERIZATION/SUPERPOSITION" BETWEEN ATOMS OF THE CLUSTER
c     WARN_2. WARN ABOUT POSSIBLE "DESAGREGATION/EVAPORATION" OF THE CLUSTER
c     WARN_3. WARN ABOUT POSSIBLE LOST OF CONSERVATION OF THE TOTAL LINEAR MOMENTUM
c
      REPN = 0.0D0
c
c...  Loop over atoms IA
      DO 235 IA = 1, NATOMS1
c...  Loop over atoms JA
         IAP1 = IA + 1
         DO 245 JA = IAP1, NATOMS           
c...  Components of vector distance betten atoms IA and JA                      
            XIJ = XX(IA)-XX(JA)                                                 
            YIJ = YY(IA)-YY(JA)                                                 
            ZIJ = ZZ(IA)-ZZ(JA)                                                     
            RSQ = XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ     
            R   = DSQRT(RSQ)
c...  Calculate the nuclear repulsion energy (nre) to compare with ENRNUC - from Gaussian
c     Move this calculation to GEWALD...........................................................Version 2.0 of MDCLG.f
            REPN = REPN + ATNZZ(IA)*ATNZZ(JA)/R
c     .............................................................................
c...  WARN_1.        
            IF (R.LE.1.5D0) THEN
               WRITE(6,49)
               WRITE(6,'(A48,2(I4))')'WARNING!!! Possible "dimerization
     &" between atoms',IA,JA
               WRITE(6,'(A19,E10.3)')'Distance (bohr) = ',R 
               WRITE(6,49)
            ENDIF
c     .............................................................................
c...  WARN_2.
            IF(LTYPE.LE.1)THEN 
            NLATTP = IDNINT(DFLOAT(NMOLCS-32)/DFLOAT(32)+0.49D0)
            NSPHER = 3+NLATTP
c     RTRUNC is the diameter of the smallest sphere that still contains the system if LPOST TRUE + 5.0   
            RTRUNC = 2.0D0*DIST*DFLOAT(NSPHER-1)+5.0D0
            ELSEIF(LTYPE.GE.2)THEN
            RTRUNC=DSQRT(3.0D0)*BOX+3.0D0    
            ENDIF
            IF (R.GT.RTRUNC) THEN
               WRITE(6,49)
               WRITE(6,'(A55,2(I4))')'WARNING!!! Possible desagregation
     & of the cluster; atoms',IA,JA
               WRITE(6,'(A19,E10.3,A1,E10.3)')'Distance (bohr) = ',R,
     &'>',RTRUNC 
               WRITE(6,49)
            ENDIF    
  245    CONTINUE
  235 CONTINUE           
c     .............................................................................
c     WARN_3.
c
      TLM=PX**2.0D0+PY**2.0D0+PZ**2.0D0
      TLMOMT=DSQRT(TLM)
      IF(TLMOMT.GT.1.0D-5)THEN
      WRITE(6,49)
      WRITE(6,'(A62,E10.3)')'WARNING!!! Total linear momentum
     & (electron masses*bohr/atu) = ',TLMOMT
      WRITE(6,49)
      ENDIF
c
c     ATTENTION_1.
c      
      IF(LSCALE)THEN
         WRITE(6,59)
         WRITE(6,'(A114)')'ATTENTION!!! Scaled velocities for desired 
     &temperature - excpected decrease of the temperature drift 
     &next MD step.'
         WRITE(6,59)
      ENDIF      
c
c     PRINT OUT THE NUCLEAR REPULSION ENERGY OF Gaussian AND OF MDCLG FOR COMPARISON
c     Multiplication of REPN by COULV is redundant COULV=1.0 - multiply only for using the variable COULV 
c     Move this calculation to GEWALD...........................................................Version 2.0 of MDCLG.f 
      REPNUC = COULV*REPN
      WRITE(6,59)
      WRITE(6,'(2(A20,F11.3))')'Gaussian nre [Ha] = ',ENRNUC,
     &'MDCLG nre [Ha] = ',REPNUC
      WRITE(6,59)      
c      
c     .............................................................................      
c
c     CALCULATE POSITIONS AT STEP (N+1)
c
      DO 255 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 265 J=1,JMAX 
            DO 275 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 285 L = 1,LMAX
c...  Velet "leap-frog" for positions R(N+1)= R(N)+V(N+1/2)*STEP
                  RX(I,J,K,L)=RX(I,J,K,L)+VX(I,J,K,L)*STEP
                  RY(I,J,K,L)=RY(I,J,K,L)+VY(I,J,K,L)*STEP
                  RZ(I,J,K,L)=RZ(I,J,K,L)+VZ(I,J,K,L)*STEP
c
c     
c     LTYPE 0 OR LTYPE 1
c     APPLY PERIODIC BOUNDARY CONDITIONS - CENTER OF THE BOX DEFINED AS (0,0,0).................Version 2.0 of MDCLG.f
c     CARTESIAN COORDINATES POSITIVE AND NEGATIVE
c
c     LTYPE 2 OR LTYPE 3
c     APPLY PERIODIC BOUNDARY CONDITIONS - CENTER OF THE BOX DEFINED AS (L/2,L/2,L/2)...........Version 2.0 of MDCLG.f
c     CARTESIAN COORDINATES ALL POSITIVE
c
c     ATTENTION!!! PBC ARE ONLY APPLIED IF A BULK LIQUID IS SIMULATED
c                  NEED TO DEFINE CUBE AND CUBEH FROM DESIRED DENSITY...........................Version 2.0 of MDCLG.f
c                  NEED TO CODE IN THE EWALD SUM FOR POLYATOMIC MOLECULES.......................Version 2.0 of MDCLG.f
c
*      IF (RX(I,J,K,L).LT.CUBEH)RX(I,J,K,L)=RX(I,J,K,L)+CUBE                                              
*      IF (RX(I,J,K,L).GT.CUBEH)RX(I,J,K,L)=RX(I,J,K,L)-CUBE                                                    
c     
*      IF (RY(I,J,K,L).LT.CUBEH)RY(I,J,K,L)=RY(I,J,K,L)+CUBE                                              
*      IF (RY(I,J,K,L).GT.CUBEH)RY(I,J,K,L)=RY(I,J,K,L)-CUBE                                                    
c     
*      IF (RZ(I,J,K,L).LT.CUBEH)RZ(I,J,K,L)=RZ(I,J,K,L)+CUBE                                              
*      IF (RZ(I,J,K,L).GT.CUBEH)RZ(I,J,K,L)=RZ(I,J,K,L)-CUBE 
c   
  285          CONTINUE
  275       CONTINUE      
  265    CONTINUE
  255 CONTINUE
c
c     SCALE VELOCITIES FOR DESIRED TEMPERATURE EVERY KSCALE TIME-STEPS DURING EQUILIBRATION
c
c     ATTENTION!!! VELOCITIES AT STEP N+1/2 ARE SCALED FOR DESIRED TEMPERATURE
c                  THE ACTUAL TEMPERATURE IS CALCULATED USING THE VELOCITIES AT STEP N THROUGH ACTKE
c         
      IF(LSCALE)CALL GSCALE
c 
c     ---------------------------------------------------------------------------------
c  
c     BUILD THE THREE NEW INPUT FILES FOR THE NEXT TIME STEP
c
c     gaussr.com - Positions at step N+1..........[Gausian 03 input]
c     MDVELT.DAT - Velocities at step N+1/2.......[MDCLG input]
c     MDPOST.DAT - Positions at step N+1..........[MDCLG input]
c
c     WRITE GAUSSIAN OPTIONS TO THE NEW GAUSSIAN 03 INPUT FILE gaussr.com 
c 
      DO 295 LL=1,NGA
         WRITE(5,'(A115)')NLINEG(LL)
  295 CONTINUE      
c
c     WRITE HEAD OF MDPOST.DAT
c
      WRITE(3,'((7X,A7),2(17X,A7))')'X(Bohr)','Y(Bohr)','Z(Bohr)'
c
c     WRITE HEAD OF MDVELT.DAT
c
      WRITE(4,'((7X,A12),2(12X,A12))')'VX(Bohr/aut)','VY(Bohr/aut)',
     &'VZ(Bohr/aut)'
c
      DO 305 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 315 J=1,JMAX 
            DO 325 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 335 L = 1,LMAX
                  WRITE(3,*)RX(I,J,K,L),RY(I,J,K,L),RZ(I,J,K,L)
                  WRITE(4,*)VX(I,J,K,L),VY(I,J,K,L),VZ(I,J,K,L)
                  WRITE(5,'(2X,I3,1X,3(1X,F19.14))')NATZ(I,K,L),
     &RX(I,J,K,L),RY(I,J,K,L),RZ(I,J,K,L)   
  335          CONTINUE
  325       CONTINUE      
  315    CONTINUE
  305 CONTINUE
c...  The last line of Gaussian 03 input file must be empty       
      WRITE(5,*)'                                                     '
c
c...  SPECIFY ATOMS'S BASIS SETS IF THE OPTION GEN IS USED IN GAUSSIAN
c
c     ATTENTION!!! NGB CAN BE > NLG IN WHICH CASE GEN WAS NOT USED IN THE GAUSSIAN INPUT FILE
c                  AND NO LINES ARE WRITTEN IN THE FOLLOWING LOOP
c  
      IF(NGB.LT.NLG)THEN
         DO 345 KK=NGB,NLG
            WRITE(5,'(A115)')NLINEG(KK)
  345    CONTINUE
      ENDIF          
c
c     ---------------------------------------------------------------------------------
c
  39  FORMAT(T1,80('-'))
  49  FORMAT(T1,80('.'))
  59  FORMAT(T1,114('.'))
c        
      RETURN
      END   
