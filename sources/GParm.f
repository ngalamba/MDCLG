      SUBROUTINE GPARM
c
c     READ INPUT PARAMETERS FOR (N,V,E) MD SIMULATION
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c
      PARAMETER(NCOMPMAX=5)
      PARAMETER(NMOLMAX=256)
      PARAMETER(NSPECMAX=10)
      PARAMETER(NATSPMAX=100)
      PARAMETER(NLGAUSS=500)
      PARAMETER(NPARTCL=500)
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Value used by Prof. Benedito = 1836.15 electron masses [proton rest mass] 
c                                         [proton mass]1.672623E-27Kg/[electron mass]9.10939E-31Kg = 1836.152585
c                                         e.g., m(Na) = 22.98977amu*1836.15 = 42212.73                   
*      PARAMETER(PTM=1836.15D0)
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Value used in ABINIT = 1822.88851 electron masses [unified atomic mass constant]
c                                         [amu=]1.6605402E-27Kg/[mass electron]9.10939E-31Kg = 1822.8885 electron masses/amu
c                                         amu = 1.6605402E-27Kg = 1.0E-3Kg/N Avogadro = 1.0E-3Kg/6.0221441990E+23          
c                                         e.g., m(Na) = 22.98977amu*1822.8885 electron masses/amu = 41907.78758 electron masses      
      PARAMETER(PTM=1822.88851D0)
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
c                                        
      PARAMETER(ATU=2.418884D-17)
      PARAMETER(THA=315.773D+3)
      PARAMETER(AVOGN=6.02214199D+23)
      PARAMETER(EHART=4.35974417D-18)
      PARAMETER(BOLTZ=3.16681520D-6)
      PARAMETER(PVACUM=0.0795774715D0)
      PARAMETER(PI=3.141592654D0)
c        
      CHARACTER*4   NAMEAT
      CHARACTER*150 NLINEG
      CHARACTER*25  NAMECH
c 
      LOGICAL LG(5),LPOST,LVEL,LSCALE,LBOMD,LNBO       
      LOGICAL LPRNT  
c
      COMMON/INTGRS/IN(14)
      COMMON/LGCLS /LG 
      COMMON/REALS /RL(18)         
c      
      COMMON/POST0 /RX0(NCOMPMAX,NSPECMAX,NATSPMAX),
     &              RY0(NCOMPMAX,NSPECMAX,NATSPMAX),
     &              RZ0(NCOMPMAX,NSPECMAX,NATSPMAX)   
      COMMON/PARMOL/NMOL(NMOLMAX),NSPEC(NSPECMAX),
     &NATMSP(NCOMPMAX,NATSPMAX)
      COMMON/ATMASS/ZMASS(NCOMPMAX,NSPECMAX,NATSPMAX)
      COMMON/ATNAME/NAMEAT(NCOMPMAX,NSPECMAX,NATSPMAX)
      COMMON/ATOMZ /NATZ(NCOMPMAX,NSPECMAX,NATSPMAX)   
      COMMON/ATNUMB/NUMBAT(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)   
      COMMON/ATCOMP/NATMOL(NPARTCL) 
      COMMON/GAUSSL/NLINEG(NLGAUSS)      
c
      EQUIVALENCE (IN(1),NCOMP),(IN(2),NLG),(IN(3),NATOMS)
      EQUIVALENCE (IN(4),NATOMS1),(IN(5),MDSTEP),(IN(6),MAXEQB)
      EQUIVALENCE (IN(7),NSTEP),(IN(8),NMOLCS),(IN(9),LTYPE)
      EQUIVALENCE (IN(10),NQUANTUM),(IN(11),NGEN)
      EQUIVALENCE (IN(12),NGA),(IN(13),NGB),(IN(14),NGAUSSV)
      EQUIVALENCE (LG(1),LPOST),(LG(2),LVEL),(LG(3),LSCALE)
      EQUIVALENCE (LG(4),LBOMD),(LG(5),LNBO)
      EQUIVALENCE (RL(1),ENRNUC),(RL(2),ENRSCF),(RL(3),STEP)
      EQUIVALENCE (RL(4),STEPSQ),(RL(5),TKECONV),(RL(7),DTEMP)
      EQUIVALENCE (RL(9),ACTKEAV),(RL(10),ENSCFAV),(RL(11),FNATOMS)
      EQUIVALENCE (RL(12),FEKJM),(RL(13),COULV),(RL(14),FNMOLCS)
      EQUIVALENCE (RL(15),DIST),(RL(16),DTE),(RL(17),DTP),(RL(18),BOX)      
c
c...................................................................
c...  READ mdcontrl.dat - MD CONTROL VARIABLES
c...................................................................
c     ATTENTION!!! NSTEP controls the equilibration run and MDSTEP controls the calculation of averages       
      READ(7,*)LPRNT,NSTEP,MDSTEP 
      READ(7,*)LPOST
      READ(7,*)LVEL
      READ(7,*)ACTKEAV,ENSCFAV
      CLOSE(7)           
c
c...................................................................
c...  READ MDINPUT.DAT - LOGICAL VARIABLES
c...................................................................
c
      READ(1,*)NGAUSSV
      READ(1,*)LSCALE
      READ(1,*)LBOMD 
      READ(1,*)LNBO
c
c     READ MDINPUT.DAT - SIMULATION PARAMETERS
c           
      READ(1,*)KSCALE
      READ(1,*)MAXEQB
      READ(1,*)DTE
      READ(1,*)DTP
      READ(1,*)DTEMP
c
           IF (LPOST.OR.LVEL)LBOMD=.FALSE.
           IF (LBOMD)LPOST=.FALSE.
           IF (LBOMD)LVEL =.FALSE.
           IF (MOD(NSTEP,KSCALE).NE.0)LSCALE=.FALSE.
           IF (NSTEP.GT.MAXEQB)LSCALE=.FALSE.     
c
c     READ MDINPUT.DAT - MOLECULAR PARAMETERS
c          CHEMICAL NAME OF THE SYSTEM................................................[e.g., Na+.(H2O)8]
c          NUMBER OF CHEMICAL COMPONENTS..............................................[e.g.,          2]
c          NUMBER OF MOLECULES OF COMPONENT no 1......................................AND...............
c          NUMBER OF DIFFERENT ATOMIC SPECIES OF MOLECULES OF COMPONENT no 1..........[e.g.,        8 2]
c          NUMBER OF ATOMS OF THE ATOMIC SPECIES OF MOLECULES OF COMPONENT no 1.......[e.g.,        1  ]
c                                                                                     [             2  ]
c          NUMBER OF MOLECULES OF COMPONENT no 2......................................AND...............
c          NUMBER OF DIFFERENT ATOMIC SPECIES OF MOLECULES OF COMPONENT no 2..........[e.g.,        1 1]
c          NUMBER OF ATOMS OF THE ATOMIC SPECIES OF MOLECULES OF COMPONENT no 2.......[e.g.,        1  ]
c
c          CALCULATE NUMBER OF ATOMS P/ MOLECULE FOR EACH COMPONENT OF THE CLUSTER - NATMOL
c
      READ(1,*)NAMECH
      READ(1,*)NCOMP
c
      DO 5 LA = 1,NCOMP
         NATMOL(LA)=0.0D0
    5 CONTINUE
       
c
      DO 15 I = 1,NCOMP
         READ(1,*)NMOL(I),NSPEC(I)
         KMAX = NSPEC(I)
         DO 25 K = 1,KMAX
            READ(1,*)NATMSP(I,K)
            NATMOL(I)=NATMOL(I)+NATMSP(I,K)
   25    CONTINUE
   15 CONTINUE 
c
c     READ MDINPUT.DAT - ATOMIC SYMBOLS AND CONSTANTS 
c          CHEMICAL SYMBOL OF THE ATOM...........AND...........ATOMIC NUMBER...........AND..................
c          ATOMIC WEIGHT FOR EACH ATOMIC SPECIES......................................[e.g., O  8  15.9994 ]
c                                                           ..........................[      H  1  1.0079  ]
c                                                           ..........................[      H  1  1.0079  ]
c                                                           ..........................[      Na 11 22.98977]
c
c          CONVERT MASSES TO ATOMIC UNITS OF MASS
c
      DO 35 I = 1,NCOMP
         KMAX = NSPEC(I)
         DO 45 K = 1,KMAX
         LMAX = NATMSP(I,K)
            DO 55 L = 1,LMAX
               READ(1,*)NAMEAT(I,K,L),NATZ(I,K,L),ZMASS(I,K,L)
               ZMASS(I,K,L)=ZMASS(I,K,L)*PTM
   55       CONTINUE      
   45    CONTINUE
   35 CONTINUE
c
c     READ MDINPUT.DAT - GAUSSIAN OPTIONS TO WRITE THE NEW POSITION FILE gaussr.dat 
c     NLG  - TOTAL NUMBER OF GAUSSIAN LINES NEEDED TO BUILD THE GAUSSIAN INPUT FILES
c     NGEN - NUMBER OF LINES USED TO SPECIFY EACG ATOM'S BASIS FUNCTION - THESE LINES MUST
c            COME AFTER THE ATOMS POSITIONS FOLLWED OF AN EMPTY LINE
c 
      READ(1,*)NQUANTUM
      READ(1,*)NLG,NGEN
      DO 65 LL=1,NLG
         READ(1,'(A125)')NLINEG(LL)
   65 CONTINUE  
c
      IF(LPRNT)THEN
         IF(NQUANTUM.EQ.0)THEN
            WRITE(*,*)'                                               '
            WRITE(*,'(A46)')'Hartree-Fock calculation asked for in 
     &this run'
         ELSEIF(NQUANTUM.EQ.1)THEN
            WRITE(*,*)'                                               '
            WRITE(*,'(A59)')'Density Functional Theory calculation 
     &asked for in this run'
         ENDIF
      ENDIF
c
c     READ DISTANCE BETWEEN CONCENTRIC SPHERES, TYPE OF ORDERING AND ATOMIC POSITIONS TO ASSIGN 
c     INITIAL POSITIONS IN THE ROUTINE GPOST IF LPOST .TRUE.
c
      READ(1,*)DIST,LTYPE
      BOX=DIST
c
c    ATTENTION!!! THE FOLLOWING LINES OF MDINPUT.DAT ARE ONLY READ IF LPOST .TRUE.
c                 NOTHING CAN BE READ AFTER THE LAST LINE UNLESS LPOST .TRUE.
c      
      IF(LPOST)THEN
         DO 75 I = 1,NCOMP
         KMAX = NSPEC(I)
            DO 85 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 95 L = 1,LMAX
                  READ(1,*)RX0(I,K,L),RY0(I,K,L),RZ0(I,K,L)   
*check                  write(*,*)RX0(I,K,L),RY0(I,K,L),RZ0(I,K,L) 
   95          CONTINUE
   85       CONTINUE      
   75    CONTINUE
      ENDIF
c
c     READ MDINPUT.DAT - GAUSSIAN # LINE TO WRITE THE FIRST POSITION FILE gaussr.dat IF LPOST .TRUE.
c                        THE GAUSSIAN # LINE READ HERE CANNOT CONTAIN THE OPTION guess=read
c     
      IF(LPOST)THEN
         READ(1,'(A125)')NLINEG(4)
      ENDIF
c
c......................................................................................
c...  CALCULATE SIMULATION PARAMETERS FROM THE INPUT PARAMETERS AND CONSTANT PARAMETERS 
c......................................................................................         
c
      NGA=NLG-NGEN
      NGB=NGA+1
c
      DTEAU = DTE*1.0D-15/ATU
      DTPAU = DTP*1.0D-15/ATU
c      
      IF (NSTEP.LE.MAXEQB)THEN
      STEP  = DTE*1.0D-15/ATU
      ELSE
      STEP  = DTP*1.0D-15/ATU
      ENDIF
c      
      STEPSQ=STEP*STEP
c      
c     ATTENTION!!! This conversion from K to Ha is only used for printing the desired temperature in Ha 
c                  TEMPAU is not used in the calculations     
      TEMPAU=DTEMP/THA
c      
c     Coulomb denominator 1/4*Pi*Vacuum_Permittivity - this is equal to 1.0 in atomic units     
      COULV =1.0D0/(4.0D0*PI*PVACUM)
c      
      NATOMS=0
      NMOLCS=0      
c      
c     COUNT TOTAL NUMBER OF ATOMS OF THE SYSTEM AND TAG THEM WITH NUMBAT  
c   
      DO 105 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 115 J=1,JMAX 
            NMOLCS = NMOLCS+1
            DO 125 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 135 L = 1,LMAX
                  NATOMS=NATOMS+1
                  NUMBAT(I,J,K,L)=NATOMS
  135          CONTINUE
  125       CONTINUE      
  115    CONTINUE
  105 CONTINUE            
c
      NATOMS1=NATOMS-1
      FNATOMS=DFLOAT(NATOMS)
      FNMOLCS=DFLOAT(NMOLCS)
c      
c     Constant factor used to convert the actual kinetic energy [Ha] to the actual temperature [k] 
      TKECONV=2.0D0/(3.0D0*FNATOMS*BOLTZ)
c     Ha to kJ/mol of atoms conversion factor      
      FEKJM  = AVOGN*1.0D-3*EHART/FNATOMS      
c              
c...................................................................
c...  WRITE MDOUTPUT.DAT - MD TITLES AND INPUT PARAMETERS
c...................................................................
c
      IF(LPRNT)THEN
      WRITE(6,'(31X,A55)')'-------------------------------------------
     &----------'
      WRITE(6,'((34X,A35),(2X,A25))')'BORN-OPPENHEIMER MD SIMULATION 
     &OF',NAMECH
      WRITE(6,'(31X,A55)')'-------------------------------------------
     &----------'
      WRITE(6,*)'                                                     '
      WRITE(6,'(A21)')'Simulation parameters'   
      WRITE(6,'(A21)')'---------------------'                                   
      WRITE(6,*)'                                                     '
      WRITE(6,'((A23),(I3))')'Number of components =',NCOMP
      WRITE(6,'((A24),(I4))')'Total number of atoms =',NATOMS
      WRITE(6,'((A32),(I5))')'Number of equilibration steps =',MAXEQB
      WRITE(6,'((A35),(F11.3))')'Time-step for equilibration (fs) =',
     &DTE
      WRITE(6,'((A36),(F11.3))')'Time-step for equilibration (atu) =',
     &DTEAU
      WRITE(6,'((A32),(F11.3))')'Time-step for production (fs) =',DTP
      WRITE(6,'((A33),(F11.3))')'Time-step for production (atu) =',
     &DTPAU
      WRITE(6,'((A18),(F11.3))')'Temperature (K) =',DTEMP
      WRITE(6,'((A19),(F11.7))')'Temperature (Ha) =',TEMPAU
c      
      DO 145 I = 1,NCOMP
         WRITE(6,*)'                                                  '
         WRITE(6,'((A9),(I3))')'Component',I   
         WRITE(6,'(A12)')'------------'                                   
         WRITE(6,*)'                                                  '
         WRITE(6,'((A33),(I3),(A3),(I3))')'Number of molecules of 
     &component',I,' = ',NMOL(I)
         WRITE(6,'((A38),(I3),(A3),(I3))')'Number of atomic species of 
     &component',I,' = ',NSPEC(I)
         KMAX = NSPEC(I)
         DO 155 K = 1,KMAX
            WRITE(6,'((A39),(I3),(A3),(I3))')'Number of atoms p/ 
     &molecule of species',K,' = ',NATMSP(I,K) 
            LMAX = NATMSP(I,K)
            DO 165 L = 1,LMAX
               WRITE(6,'(3(A5),(1X,I3),(1X,A15),(F13.5))')'Atom',
     &NAMEAT(I,K,L),'Z = ',NATZ(I,K,L),'Atomic mass =',ZMASS(I,K,L)                                                                                   
  165       CONTINUE      
  155    CONTINUE
  145 CONTINUE
      WRITE(6,*)'                                                     ' 
      WRITE(6,'(A29)')'Gaussian 03 calculation input'
      WRITE(6,'(A29)')'-----------------------------'
      WRITE(6,*)'                                                     '
      DO 175 LL=1,NLG
         WRITE(6,'(A125)')NLINEG(LL)
  175 CONTINUE 
c
c     WRITE ATOMIC POSITIONS INPUT
c
      IF(LPOST)THEN
         WRITE(6,*)'                                                  ' 
         WRITE(6,'(A29)')'Atomic positions input [Bohr]'
         WRITE(6,'(A29)')'-----------------------------'
         WRITE(6,*)'                                                  '
         DO 185 I = 1,NCOMP
         KMAX = NSPEC(I)
            DO 195 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 205 L = 1,LMAX
                  WRITE(6,'((A5),3(F8.3))')NAMEAT(I,K,L),RX0(I,K,L),
     &RY0(I,K,L),RZ0(I,K,L)  
*check                  write(*,*)NAMEAT(I,K,L),RX0(I,K,L),
*check     &RY0(I,K,L),RZ0(I,K,L)  
  205          CONTINUE
  195       CONTINUE      
  185    CONTINUE 
         IF(LTYPE.LE.1)THEN
         WRITE(6,*)'                                                  ' 
         WRITE(6,'(A57)')'Radii [A] of the spheres used to assign 
     &initial positions'
         WRITE(6,'(A57)')'----------------------------------------
     &-----------------'
         WRITE(6,*)'                                                  '
         ELSEIF(LTYPE.GE.2)THEN
         WRITE(6,*)'                                                  ' 
         WRITE(6,'(A36)')'Lenght [A] of the side of the MD box'
         WRITE(6,'(A36)')'------------------------------------'
         WRITE(6,*)'                                                  '
         ENDIF 
c...The radii of the spheres are written in the subroutine GPOST if LPOST .TRUE.
      ELSE
         WRITE(6,*)'                                                  ' 
         WRITE(6,19) 
      ENDIF 
      ENDIF      
c
c     WARN ABOUT THE END OF EQUILIBRATION AND RESTART STEP COUNTERS
c           
      IF(NSTEP.EQ.(MAXEQB+1))THEN
      MDSTEP = 0
      WRITE(6,*)'                                                     '
      WRITE(6,*)'                                                     '
      WRITE(6,'(47X,A22)')'----------------------'
      WRITE(6,9)
      WRITE(6,'(47X,A22)')'----------------------'
      WRITE(8,9)
      WRITE(*,9)  
      ENDIF                         
c
c...................................................................
c...  WRITE mdcontrl.dat - MD CONTROL VARIABLES
c...................................................................
c
c     UPDATE TIME-STEP AND TURN OFF THE WRITTING OF INPUT PARAMETERS
c

      OPEN(7,FILE='mdcontrl.dat',ACCESS='SEQUENTIAL',STATUS='OLD')
c     
      IF(LBOMD)THEN 
      NSTEP  = NSTEP+1
      MDSTEP = MDSTEP+1 
      ENDIF
      LPRNT  = .FALSE.
      WRITE(7,*)LPRNT,NSTEP,MDSTEP 
c
c     MARK THE CURRENT TIME-STEP IN THE OUTPUT FILE AFTER MDSTEP UPDATED
c 
      IF(LBOMD)THEN
         WRITE(6,*)'                                                  '
         WRITE(6,*)'                                                  '
         IF((NSTEP-1).LE.MAXEQB)THEN
            WRITE(6,*)'MOLECULAR DYNAMICS STEP Eql',MDSTEP
         ELSE
            WRITE(6,*)'MOLECULAR DYNAMICS STEP Prd',MDSTEP
         ENDIF
         WRITE(6,'(A29)')'-----------------------------' 
      ENDIF               
c
c     OUTPUT FILES FORMATS
c                                                          
 9    FORMAT(/47X,'EQUILIBRATION COMPLETE'//)
19    FORMAT(T1,107('_'))      
c      
      RETURN
      END   
