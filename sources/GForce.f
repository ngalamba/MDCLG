c      
      SUBROUTINE GFORCE
c
c     READ FORCES, SCF ENERGY TERMS, MULLIKEN ATOMIC CHARGES, ETC... FROM GAUSSIAN OUTPUT FILE 
c     READ POSITIONS AND VELOCITIES FROM PREVIOUS MD STEP 
c
c     ATTENTION!!! SCF ENERGY = KINETIC ENR + POTENTIAL ENR + ELECTRONIC ENR + NUCLEAR REPULSION ENR     
c     ATTENTION!!! THE VELOCITIES READ ARE THE VERLET "LEAP-FROG" HALF MD-STEP
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c
      PARAMETER(NCOMPMAX=5)
      PARAMETER(NMOLMAX=256)
      PARAMETER(NSPECMAX=10)
      PARAMETER(NATSPMAX=100)
      PARAMETER(NPARTCL=500)
c
      LOGICAL LG(5),LNBO
c
      CHARACTER*4  NAMEAT
      CHARACTER*4  NSYMBL
      CHARACTER*15 AN
c      
      CHARACTER*7  PDUMMY
      CHARACTER*7  VDUMMY   
      CHARACTER*7  FDUMMY
      CHARACTER*70 ENRDUMMY
      CHARACTER*70 MLKDUMMY
      CHARACTER*2  ANDUMMY
      CHARACTER*70 NBODUMMY
      CHARACTER*50 MDIPDUMMY
      CHARACTER*7  DIPDUMMY
               
c                                               
      COMMON/INTGRS/IN(14)
      COMMON/LGCLS /LG
      COMMON/REALS /RL(18)     
c      
      COMMON/POSIT /RX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)
      COMMON/POSSD /XX(NPARTCL),YY(NPARTCL),ZZ(NPARTCL) 
      COMMON/POSTD /RXX(NCOMPMAX,NMOLMAX,NPARTCL),
     &              RYY(NCOMPMAX,NMOLMAX,NPARTCL),
     &              RZZ(NCOMPMAX,NMOLMAX,NPARTCL)
      COMMON/VELOC /VX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              VY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              VZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)      
      COMMON/FORCE /FX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              FY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              FZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)         
      COMMON/PARMOL/NMOL(NMOLMAX),NSPEC(NSPECMAX),
     &NATMSP(NCOMPMAX,NATSPMAX)
      COMMON/ATNAME/NAMEAT(NCOMPMAX,NSPECMAX,NATSPMAX)
      COMMON/ATSYMB/NSYMBL(NCOMPMAX,NPARTCL)
      COMMON/ATOMZ /NATZ(NCOMPMAX,NSPECMAX,NATSPMAX)
      COMMON/ATNZSD/ATNZZ(NPARTCL)
      COMMON/MULLIK/CHMULLK(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)
      COMMON/NBOCH /CHNBO(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)
c
      EQUIVALENCE (IN(1),NCOMP),(IN(10),NQUANTUM),(IN(14),NGAUSSV)
      EQUIVALENCE (LG(5),LNBO)
      EQUIVALENCE (RL(1),ENRNUC),(RL(2),ENRSCF)
c
      ENRSCF = 0.0D0
      NSCF = 0
c
c     ATTENTION!!!
c
c     THE FORCES, POSITIONS AND VELOCITIES ARE READ HERE BUT ONLY WRITTEN IN THE ROUTINE GVERLET
c     THE OTHER QUANTITIES ARE READ AND WRITTEN IN THIS ROUTINE EVERY TIME-STEP
c
c     READ FORCE TITTLES OF Gaussian 03 OUTPUT FILE
c
c     CAREFULL!!! IF NEEDED ADJUST FORMATS - DEFAULT FDUMMY, PDUMMY AND VDUMMY SMALL ENOUGH, i.e., CHARACTER*7
c                 TO AVOID FORMAT SPECIFICATION - HENCE VARIABLES NOT USED
c     
      READ(2,*)FDUMMY
      READ(2,*)FDUMMY
      READ(2,*)FDUMMY
c 
c     READ POSITION AND VELOCITY TITTLES OF MD INPUT/OUTPUT FILE
c 
      READ(3,*)PDUMMY 
      READ(4,*)VDUMMY
c
c     READ Ab-Initio FORCES [F(N)], MD POSTITIONS [R(N)] AND MD HALF-STEP VELOCITIES [V(N-1/2)]
c
c     NC AND NA ARE NOT USED - CENTER NUMBER AND ATOMIC NUMBER THAT APPEAR BEFORE THE FORCE IN GAUSSIAN OUTPUT FILE
c
      DO 5 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 15 J=1,JMAX 
            DO 25 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 35 L = 1,LMAX
                  READ(2,*)NC,NA,FX(I,J,K,L),FY(I,J,K,L),FZ(I,J,K,L)
                  READ(3,*)RX(I,J,K,L),RY(I,J,K,L),RZ(I,J,K,L)
                  READ(4,*)VX(I,J,K,L),VY(I,J,K,L),VZ(I,J,K,L)
   35          CONTINUE
   25       CONTINUE      
   15    CONTINUE
    5 CONTINUE   
c    
      CLOSE(3)
      CLOSE(4)
      CLOSE(5) 
c
      IF(NGAUSSV.eq.03)THEN
c     ..................................................................................................... START reading Gaussian 03
c
c     READ GAUSSIAN 03 OUTPUT FILE gaussf.com
c     ---------------------------------------    
c
c     R1. SCF ENERGY TERMS AND CONVERGENCE PARAMETERS
c
c     ENERGY VARIABLES:
c
c     ENRNUC - NUCLEAR REPULSION ENERGY (Ha)
c     ENRSCF - SCF TOTAL ENERGY (Ha)
c     NSCF   - NUMBER OF SCF CYCLES
c     SCFCONV- CONVERGENCE (Ha)
c     ENRKIN - KINETIC ENERGY (Ha)
c     ENRPOT - POTENTIAL ENERGY (Ha)
c     ENRELEC- ELECTRONIC ENERGY (Ha)  
c
         READ(2,'(A31,F21.10)')ENRDUMMY,ENRNUC
         IF(NQUANTUM.EQ.0)THEN
c...Hartree-Fock calculation format      
            READ(2,'(A20,E16.10,A15,I5)')ENRDUMMY,ENRSCF,ENRDUMMY,NSCF
         ELSEIF(NQUANTUM.EQ.1)THEN
c...b3lyp functional calculation format      
c            READ(2,'(A26,E16.10,A15,I5)')ENRDUMMY,ENRSCF,ENRDUMMY,NSCF
c...fluctuating format for different functionals - stopped reading NSCF
            READ(2,'(A26,E16.10)')ENRDUMMY,ENRSCF
         ELSE
            WRITE(*,*)'The type of quantum calulation must be correctly 
     &specified in the input file - see variable NQUANTUM'    
         ENDIF
         READ(2,'(A21,E14.4)')ENRDUMMY,SCFCONV 
         READ(2,'(A30)')ENRDUMMY
         READ(2,'((A4,E19.12),(A4,E19.12),(A4,E19.12))')ENRDUMMY,ENRKIN,
     &ENRDUMMY,ENRPOT,ENRDUMMY,ENRELEC
         READ(2,'(A30)')ENRDUMMY
c
c     R2. MULLIKEN ATOMIC CHARGES 
c
c     NM IS NOT USED         - THE NUMBER BEFORE THE MULLIKEN CHARGES IN GAUSSIAN OUTPUT FILE
c     NC AND AN ARE NOT USED - CENTER NUMBER AND ATOM NAME IN GAUSSIAN OUTPUT FILE
c
         READ(2,'(A25)')MLKDUMMY
         READ(2,*)NM
c        
         DO 45 I = 1,NCOMP
            JMAX = NMOL(I)
            KMAX = NSPEC(I)
            DO 55 J=1,JMAX 
               DO 65 K = 1,KMAX
                  LMAX = NATMSP(I,K)
                  DO 75 L = 1,LMAX
                     READ(2,'(I6,A4,F11.6)')NC,AN,CHMULLK(I,J,K,L)
   75             CONTINUE
   65          CONTINUE      
   55       CONTINUE
   45    CONTINUE
c        
         READ(2,'(A25,F10.5)')MLKDUMMY,SUMMULLK
c
c     R3. NATURAL BOND ORBITAL CHARGES
c
         IF(LNBO)THEN
            READ(2,'(A35)')NBODUMMY
            READ(2,'(A50)')NBODUMMY
            READ(2,'(A55)')NBODUMMY
            READ(2,'(A50)')NBODUMMY
            READ(2,'(A50)')NBODUMMY
            READ(2,'(A50)')NBODUMMY 
            DO 85 I = 1,NCOMP
            JMAX = NMOL(I)
            KMAX = NSPEC(I)
               DO 95 J=1,JMAX 
                  DO 105 K = 1,KMAX
                     LMAX = NATMSP(I,K)
                     DO 115 L = 1,LMAX
                        READ(2,'(A8,I5,F10.5)')AN,NC,CHNBO(I,J,K,L)
  115                CONTINUE
  105             CONTINUE      
   95          CONTINUE
   85       CONTINUE
            READ(2,'(A50)')NBODUMMY
            READ(2,'(A12,F11.5)')NBODUMMY,SUMNBO
         ENDIF
c
c     R4. TOTAL DIPOLE MOMENT
c
         READ(2,'(A45)')MDIPDUMMY
         READ(2,'(4(A7,F10.4))')DIPDUMMY,DIPMX,DIPDUMMY,DIPMY,DIPDUMMY,
     &DIPMZ,DIPDUMMY,DIPMT
c
c
c     R5.
c
c     ..................................................................................................... STOP reading Gaussian 03     
c
      ELSEIF(NGAUSSV.eq.09)THEN
c     ..................................................................................................... START reading Gaussian 09
c
c     READ GAUSSIAN 03 OUTPUT FILE gaussf.com
c     ---------------------------------------    
c
c     R1. SCF ENERGY TERMS AND CONVERGENCE PARAMETERS
c
c     ENERGY VARIABLES:
c
c     ENRNUC - NUCLEAR REPULSION ENERGY (Ha)
c     ENRSCF - SCF TOTAL ENERGY (Ha)
c     NSCF   - NUMBER OF SCF CYCLES
c     SCFCONV- CONVERGENCE (Ha)
c     ENRKIN - KINETIC ENERGY (Ha)
c     ENRPOT - POTENTIAL ENERGY (Ha)
c     ENRELEC- ELECTRONIC ENERGY (Ha)  
c
         READ(2,'(A31,F21.10)')ENRDUMMY,ENRNUC
         
         IF(NQUANTUM.EQ.0)THEN
c...Hartree-Fock calculation format      
            READ(2,'(A20,E16.10,A15,I5)')ENRDUMMY,ENRSCF,ENRDUMMY,NSCF
         ELSEIF(NQUANTUM.EQ.1)THEN
c...b3lyp functional calculation format      
            READ(2,'(A23,E16.10,A15,I5)')ENRDUMMY,ENRSCF,ENRDUMMY,NSCF  
         ELSE
            WRITE(*,*)'The type of quantum calulation must be correctly 
     &specified in the input file - see variable NQUANTUM'    
         ENDIF
         
         READ(2,'(A28,E12.2)')ENRDUMMY,SCFCONV
         
         READ(2,'(A30)')ENRDUMMY
         READ(2,'(A30)')ENRDUMMY
         READ(2,'((A4,E19.12),(A4,E19.12),(A4,E19.12))')ENRDUMMY,ENRKIN,
     &ENRDUMMY,ENRPOT,ENRDUMMY,ENRELEC
c
c     R2. MULLIKEN ATOMIC CHARGES 
c
c     NM IS NOT USED         - THE NUMBER BEFORE THE MULLIKEN CHARGES IN GAUSSIAN OUTPUT FILE
c     NC AND AN ARE NOT USED - CENTER NUMBER AND ATOM NAME IN GAUSSIAN OUTPUT FILE
c
         READ(2,'(A25)')MLKDUMMY
         READ(2,*)NM
c        
         DO 245 I = 1,NCOMP
            JMAX = NMOL(I)
            KMAX = NSPEC(I)
            DO 255 J=1,JMAX 
               DO 265 K = 1,KMAX
                  LMAX = NATMSP(I,K)
                  DO 275 L = 1,LMAX
                     READ(2,'(I6,A4,F11.6)')NC,AN,CHMULLK(I,J,K,L)
  275             CONTINUE
  265          CONTINUE      
  255       CONTINUE
  245    CONTINUE
c        
         READ(2,'(A26,F10.5)')MLKDUMMY,SUMMULLK
c
c     R3. NATURAL BOND ORBITAL CHARGES
c
         IF(LNBO)THEN
            READ(2,'(A35)')NBODUMMY
            READ(2,'(A50)')NBODUMMY
            READ(2,'(A55)')NBODUMMY
            READ(2,'(A50)')NBODUMMY
            READ(2,'(A50)')NBODUMMY
            READ(2,'(A50)')NBODUMMY 
            DO 285 I = 1,NCOMP
            JMAX = NMOL(I)
            KMAX = NSPEC(I)
               DO 295 J=1,JMAX 
                  DO 305 K = 1,KMAX
                     LMAX = NATMSP(I,K)
                     DO 315 L = 1,LMAX
                        READ(2,'(A8,I5,F10.5)')AN,NC,CHNBO(I,J,K,L)
  315                CONTINUE
  305             CONTINUE      
  295          CONTINUE
  285       CONTINUE
            READ(2,'(A50)')NBODUMMY
            READ(2,'(A12,F11.5)')NBODUMMY,SUMNBO
         ENDIF
c
c     R4. TOTAL DIPOLE MOMENT
c
         READ(2,'(A45)')MDIPDUMMY
         READ(2,'(4(A6,F20.4))')DIPDUMMY,DIPMX,DIPDUMMY,DIPMY,DIPDUMMY,
     &DIPMZ,DIPDUMMY,DIPMT
c
c
c     R5.
c
c     ..................................................................................................... STOP reading Gaussian 09     
c
      ENDIF
c     ..................................................................................................... START writting mdoutput       
c
c     WRITE GAUSSIAN OUTPUT QUANTITIES TO MD OUTPUT FILE MDOUTPUT.DAT
c     ---------------------------------------------------------------
c
c     W1. SCF ENERGY TERMS AND CONVERGENCE PARAMETERS
c 
      WRITE(6,*)'                                                     '
      WRITE(6,29)
      WRITE(6,*)'Components of total SCF energy (in Hartree):'
      WRITE(6,*)'                                                     '
      WRITE(6,'(A27,E21.12)')'Kinetic energy           = ',ENRKIN
      WRITE(6,'(A27,E21.12)')'Potential energy         = ',ENRPOT
      WRITE(6,'(A27,E21.12)')'Electronic energy        = ',ENRELEC
      WRITE(6,'(A27,E21.12)')'Nuclear repulsion energy = ',ENRNUC
      WRITE(6,'(A27,E21.12)')'>>>>>>>>Total SCF energy = ',ENRSCF
      WRITE(6,*)'                                                     '
      WRITE(6,*)'Other information on the energy calculation:'
      WRITE(6,*)'                                                     '
      WRITE(6,'(A6,I5,A27)')'At SCF',NSCF,'total energy is converged'
      WRITE(6,'(A14,E14.4)')'Convergence = ',SCFCONV
      WRITE(6,29)
c
c     W2. MULLIKEN CHARGES
c
      IF(.NOT.LNBO)THEN
      WRITE(6,*)'Mulliken atomic charges:'
      WRITE(6,*)'                                                     '
      DO 125 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 135 J=1,JMAX 
            DO 145 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 155 L = 1,LMAX
                  WRITE(6,'(I3,A4,F11.6)')NATZ(I,K,L),NAMEAT(I,K,L),
     &CHMULLK(I,J,K,L) 
  155          CONTINUE
  145       CONTINUE      
  135    CONTINUE
  125 CONTINUE 
      WRITE(6,*)'                                                     '
      WRITE(6,'(A26,F10.5)')'Sum of Mulliken charges = ',SUMMULLK
      WRITE(6,29) 
c      
      ELSEIF(LNBO)THEN   
c
c     W3. MULLIKEN CHARGES AND NATURAL BOND ORBITAL CHARGES
c     
      IF((NSTEP-1).LE.MAXEQB)THEN
      WRITE(6,*)'Equilibration Mulliken and Natural Bond Orbital charge
     &s'
      ELSE
      WRITE(6,*)'Production Mulliken and Natural Bond Orbital charges'
      ENDIF
      WRITE(6,*)'                                                     '
      WRITE(6,'(A23,15X,A18)')'Mulliken atomic charges',
     &'NBO atomic charges'
      WRITE(6,*)'                                                     '
      DO 165 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 175 J=1,JMAX 
            DO 185 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 195 L = 1,LMAX
                  WRITE(6,'(I3,A4,F11.6,21X,F11.6)')NATZ(I,K,L),
     &NAMEAT(I,K,L),CHMULLK(I,J,K,L),CHNBO(I,J,K,L) 
  195          CONTINUE
  185       CONTINUE      
  175    CONTINUE
  165 CONTINUE 
      WRITE(6,*)'                                                     '
      WRITE(6,'(A26,F10.5)')'Sum of Mulliken charges = ',SUMMULLK
      WRITE(6,'(A26,F10.5)')'Sum of NBO charges      = ',SUMNBO
      ENDIF
c
c     W4. TOTAL DIPOLE MOMENT
c
      WRITE(6,*)'                                                     '
      WRITE(6,'(A36)')'Gaussian total dipole moment (Debye)'
      WRITE(6,'(4(A6,F10.4))')' Dx = ',DIPMX,' Dy = ',DIPMY,
     &' Dz = ',DIPMZ,' Dt = ',DIPMT
      WRITE(6,29) 
c
c     W5. 
c      
c     ..................................................................................................... STOP writting mdoutput       
c      
c     TRICKY VARIABLES: PASS THE CARTESIAN COORDINATES FROM A FOUR-DIMENSIONAL ARRAY TO A ONE-DIMENSIONAL ARRAY
c                       ARRAY TO ALLOW SIMPLE CALCULATION OF DISTANCES BETWEEN ATOMS IN THE ROUTINE GVERLET
c     
c                       PASS THE ATOMIC NUMBERS FROM A THREE-DIMENSION ARRAY TO A ONE-DIMENSIONAL ARRAY
c                       TO ALLOW SIMPLE CALCULATION OF NUCLEAR REPULSION ENERGY IN THE ROUTINE GVERLET        
c
c                       PASS THE CARTESIAN COORDINATES FROM A FOUR-DIMENSIONAL ARRAY TO A THREE-DIMENSIONAL ARRAY
c                       THIS ARRAY IS SIMILAR TO THOSE USED IN CLASSICAL MD FOR MIXTURES OF POLYATOMIC MOLECULES
c                       ARRAY USED TO COMPUTE THE INTRA AND INTER MOLECULAR DISTANCES IN THE ROUTINE GVERLET
c
c                       PASS THE NAME OF THE ATOMIC SPECIES OF EACH COMPONENT FROM A THREE-DIMENSIONAL ARRAY TO 
c                       A TWO DIMENSIONAL ARRAY 
c                       USED TO IDENTIFY THE INTRA AND INTER MOLECULAR DISTANCES IN THE ROUTINE GVERLET       
c   
      NN=0
      DO 205 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 215 J=1,JMAX 
            NNN=0
            DO 225 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 235 L = 1,LMAX
                  NN=NN+1
                  XX(NN)    = RX(I,J,K,L)
                  YY(NN)    = RY(I,J,K,L)
                  ZZ(NN)    = RZ(I,J,K,L)
                  ATNZZ(NN) = NATZ(I,K,L)
c                 ..............................
                  NNN=NNN+1
                  RXX(I,J,NNN)  = RX(I,J,K,L)
                  RYY(I,J,NNN)  = RY(I,J,K,L)
                  RZZ(I,J,NNN)  = RZ(I,J,K,L)
                  NSYMBL(I,NNN) = NAMEAT(I,K,L)
  235          CONTINUE
  225       CONTINUE      
  215    CONTINUE
  205 CONTINUE            
c
  29  FORMAT(T1,80('-')) 
c             
      RETURN
      END 